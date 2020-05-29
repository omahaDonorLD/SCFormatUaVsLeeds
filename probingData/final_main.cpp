

//#include "linear_solver.hpp"
//#include <stdlib>
#include <stdio.h>
#include <cstring>
#include <ctime>
#include <string>
#include <cmath>

#include <igraph.h>
#include <assert.h>
#include <float.h>

#include <glpk.h>

using namespace std;

inline
int STREAM_FAIL(const char *FROM_FILE, int AT_LINE, const char *IN_FUNCTION)
{printf("STREAM_FAILURE, line %d, function %s, file %s\n", AT_LINE, IN_FUNCTION, FROM_FILE);return EXIT_FAILURE;};

inline
int MEMO_FAIL(const char *FROM_FILE, int AT_LINE, const char *IN_FUNCTION)
{printf("MEMO_ALLOC_FAILURE, line %d, function %s, file %s\n", AT_LINE, IN_FUNCTION, FROM_FILE);return EXIT_FAILURE;};

typedef struct sln{
	int* g_covers;/* Size : number of ground nodes, gives for each ground the number of uavs in the vicinity.*/
	int** covers;/* Size : number of ground nodes*counters, gives for each ground all the uavs that cover it.*/
	double** uavs;/* matching "contains" list with "uavs_indices". The latter contains the indices 'i'
							of elements in list "uavs" so that one can finds the coordinates of uavs[i] */

	int* counters;// number of elements (size : n_uavs) covered by the uav
	// double** distances;// squared matrix of distances between uavs
	double** dists;// squared matrix of distances between uavs. Type igraph to save space
	//igraph_t* gr;// graph of the solution. Needed for connectivity checking
	int n_uavs;// number of nodes in network : either nodes are uavs (sln : set of uavs), or ground nodes (uav : set of ground nodes)
}sln;


typedef struct block{
	float lowleft[2];// the coordinates of the lower left bound
	float upright[2];// coords of the upper right bound
	char color;// if leaf : 'b' for black (if block contains only 	), 'w' for white (area is empty), 'g' otherwise (node, not leaf, grey node)
	int label;// tells which connected component the block belongs to
	typedef struct block father;
	typedef struct block quadrants[4];// 4 regions of a quadtree's node : NW, NE, SW, SE
}block;

int max_uavs_avail;
int nbr_grnds;
int dim=2;			// dimension of input vectors
double** grnds;		// All ground nodes coordinates
double uavs_range;	// all ranges of available uavs

// limit of map
double bound_1;
double bound_2;

void readData(char** argv);
sln* method1ePasse(double** input_data, int size_input, double threshold);
double k_means(double** data, int n, sln* clusts, double error_tolerance, double range);
double euclDistance(double *node1, double *node2);
void duplicate(double lb, sln* net);
igraph_t* translate(sln* net, double threshold, double* solverSln);
void updateDistMat(sln* net, double range);

void find_covers(sln* net, double range);/* Find every covers for each ground node */

double* solve_linear_model(sln* net, double range, double lb);

void connect_CCs(sln* net, igraph_t* G, double threshold, int** restr_list, int **pairs, int* npairs, bool restrict, int level);
void populate(sln* covers, double threshold, igraph_t* solnG0, int** restr_list, int **pairs, int* npairs);

void freeSln(sln* free_sln);

void readData(char** argv)
{
	FILE* fp;

	// read data (coordinates ground nodes)
	fp=fopen(argv[1],"r");

	// read number of available uavs
	if( fscanf(fp,"%d", &max_uavs_avail) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}
//	uavs_ranges=(double*)malloc((max_uavs_avail+1)*sizeof(double));

	// read range of uavs
	if( fscanf(fp,"%lf", &uavs_range) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}

	int i=0, j=0;
/*
	// read range of each of them
	for(i=1;i<=max_uavs_avail;i++)
		fscanf(fp,"%lf", &uavs_ranges[i]);
*/

	// read number of ground nodes and their coordinates
	if( fscanf(fp,"%d", &nbr_grnds) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}

	// read limits of map
	if( fscanf(fp,"%lf,%lf", &bound_1,&bound_2) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}

	/* allocate memory for ground nodes */
	grnds=(double**)malloc((nbr_grnds+1)*sizeof(double*));
	if(grnds==NULL){ /* memory allocation failure */ MEMO_FAIL(__FILE__, __LINE__, __FUNCTION__); }

	for(i=1;i<=nbr_grnds;i++)
	{
		grnds[i]=(double*)calloc(dim,sizeof(double));// Unique, defines a point : for either a ground node or a uav
		for (j=0;j<dim;j++)
		{
			// skip comma missing for last dim
			if(j==dim-1)	fscanf(fp,"%lf", &grnds[i][j]);
			else fscanf(fp,"%lf,", &grnds[i][j]);
		}
	}

	fclose(fp);

};

double euclDistance(double *node1, double *node2)
{
	int i=0;
	double norm=0;
	for(i=0;i<dim;i++)	norm+=pow(node1[i]-node2[i],2);
    return sqrt(norm);
};

int cmpfuncXaxis(const void *el1, const void *el2)
{
	const double *uav1 = (const double*) el1;
	const double *uav2 = (const double*) el2;
	return uav1[0]-uav2[0];
};

int cmpfuncYaxis(const void *el1, const void *el2)
{
	const double *uav1 = (const double*) el1;
	const double *uav2 = (const double*) el2;
	return uav1[1]-uav2[1];
};

void updateDistMat(sln* net, double range)
{// updates the matrix of distances to avoid constantly having to compute them
	int i,j;
	double dist=0;
	for(i=1;i<=net->n_uavs;i++)
		for(j=i+1;j<=net->n_uavs;j++)
		{
			dist=euclDistance(net->uavs[i],net->uavs[j]);
			net->dists[i][j]=dist;
			net->dists[j][i]=dist;
		}
};

void duplicate(double lb, sln* net)
{/* lb for knowing if ground needs more covers, counts[] the list which gives for each ground node the number of uavs covering it,
	covers[] the list of aforementioned covers */

	int i=0,j=0,ub=0,ind=0;

	for(i=1;i<=nbr_grnds;i++)
	{
		if ( net->g_covers[i] < lb)
		{
printf("Duplicate uavs for gnode %d, which is covered with %d uavs \n", i, net->g_covers[i]);
			ub=net->g_covers[i];
			do
			{
				net->n_uavs++;// increment the number of uavs
				net->counters[net->n_uavs]=1;// Increment number of ground nodes covered by uav (just 1)
				net->covers[i][net->g_covers[i]]=net->n_uavs;// Include new uav as covering ground g
				net->g_covers[i]++;// Increment number of covers for k-th ground node
				/* Duplicate uav coords on new one*/
				ind=net->covers[i][net->g_covers[i]%ub];
				for(j=0;j<dim;j++)
					net->uavs[net->n_uavs][j]=net->uavs[ind][j];
			/* Keep copying until constraint on lb is satisfied */
			}while(net->g_covers[i] < lb);
/*
printf("Duplicated : ");
for(int a=0; a<net->g_covers[i]; a++)
printf(" %d ", net->covers[i][a]);
printf("\n");
*/
		}
	}
printf("\n");
};




void find_covers(sln* net, double range)
{/* Find every covers for each ground node */

	int i=0,j=0;
		
	for(i=1;i<=nbr_grnds;i++)
	{
		for(j=1;j<=net->n_uavs;j++)
		{
			if( euclDistance(net->uavs[j],grnds[i]) <= range && net->covers[i][0] != j)
			{

				net->covers[i][net->g_covers[i]]=j;
				net->g_covers[i]++;
			}
		}
	}

};


double* solve_linear_model(sln* net, double range, double lb)
{

	int i=0,j=0;

	/* find covers for each ground node */
	find_covers(net, range);
	
	/* Find at least one uav that needs to be replicated so that a minimum cover constraint isn't violated : at least lb covers for each target */
	for(i=1;i<=nbr_grnds;i++)
	{
		/* Create duplicate only if needed */
		if( net->g_covers[i] < lb )
		{
			duplicate(lb, net);
			break;
		}
	}

printf("Begin building linear problem with %d ground nodes, and %d uavs\n", nbr_grnds, net->n_uavs);
	/* Build linear problem */
	glp_prob *lp;
	int ia[nbr_grnds*net->n_uavs+1],ja[nbr_grnds*net->n_uavs+1];
	double ar[nbr_grnds*net->n_uavs+1], z=0.;
	
	lp = glp_create_prob();
	glp_set_prob_name(lp, "set cover");
	glp_set_obj_dir(lp, GLP_MIN);

	/* for row names */
	char row_names[nbr_grnds+1][20];
	char col_names[net->n_uavs+1][20];
	char buff[20];

	glp_add_rows(lp, nbr_grnds);
	for(i=1;i<=nbr_grnds;i++)
	{
		sprintf(buff,"r%d",i);
		strcpy(row_names[i],buff);
		glp_set_row_name(lp, i, row_names[i]);
		glp_set_row_bnds(lp, i, GLP_LO, lb, 0.0);
	}

	glp_add_cols(lp, net->n_uavs);
	for(j=1;j<=net->n_uavs;j++)
	{
		sprintf(buff,"x%d",j);
		strcpy(col_names[j],buff);
		glp_set_col_name(lp, j, col_names[j]);
		glp_set_col_bnds(lp, j, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(lp, j, 1.0);
		glp_set_col_kind(lp, j, GLP_IV);
		//glp_set_col_kind(lp, j, GLP_BV);
		//glp_set_col_kind(lp, j, GLP_CV);
	}

printf("RANGE %f, lb %f\n",range,lb);
	int counter=1;
	for(i=1;i<=nbr_grnds;i++)
	{
		for(j=1;j<=net->n_uavs;j++)
		{
			ia[counter] = i;
			ja[counter] = j;
			ar[counter] = ( euclDistance(net->uavs[j],grnds[i]) <= range ? 1.0 : 0.0 );
			counter++;
		}
	}

	glp_load_matrix(lp, nbr_grnds*net->n_uavs, ia, ja, ar);
	

/*
	int ind[net->n_uavs+1];
	double val[net->n_uavs+1];
	for(i=1;i<=nbr_grnds;i++)
	{
		for(j=1;j<net->n_uavs;j++)
		{
			ind[j]=0;
			val[j]=0;
		}
		int len=glp_get_mat_row(lp, i, ind, val);// stores non zero values of row 113 into "ind(indices of cols), val(non zero values)"
//		printf(" LO %f, COL LO %f, COL UP %f %d\n",glp_get_row_lb(lp,113),glp_get_col_lb(lp,88),glp_get_col_ub(lp,88), len);
		printf("gnode %d has %d, covers => ", i, len);
		for(j=1;j<=len;j++)
		printf(" %d, ", ind[j]);
		printf("\n");		
	}
*/

	
//	int result_solver=glp_simplex(lp, NULL);

	glp_iocp parm;
	glp_init_iocp(&parm);
	parm.presolve = GLP_ON;

	int result_solver2=glp_intopt(lp, &parm);
//	z = glp_get_obj_val(lp);
	z = glp_mip_obj_val(lp);
	
	/* Gather results */
	double* soln=(double*)calloc(net->n_uavs+1,sizeof(double));
	double sumrows=0, maxrow=0, minrow=net->n_uavs, tmp=0, lessThanAvrg=0, moreThanvrg=0;// Needed for stats on covers (constraints)
	double* covRes=(double*)calloc(nbr_grnds,sizeof(double));

	for(i=1;i<=nbr_grnds;i++)
	{
//		tmp=glp_get_row_prim(lp, i);
		tmp=glp_mip_row_val(lp, i);
//		printf("cstrnt %d = %f\n", i, tmp);
		sumrows += tmp;
		if(tmp > maxrow)	maxrow=tmp;
		if(tmp < minrow)	minrow=tmp;
		if(tmp >= 3){moreThanvrg++;}else{lessThanAvrg++;}
		covRes[(int)tmp]++;
	}

	for(j=1;j<=net->n_uavs;j++)
	{
//		soln[j]= glp_get_col_prim(lp, j);
		soln[j]= glp_mip_col_val(lp, j);
//		if (soln[j] - ceil(soln[j]) != 0)	printf(" soln %d not an integer = %f\n", j, soln[j]);
	}
	
	glp_print_sol(lp, "glpksols.txt");
	
	glp_delete_prob(lp);

	for(j=1;j<=net->n_uavs;j++)
		if(soln[j]>0)
			printf("[ %d : %f ] ",j, soln[j]);

	
	printf("\nRO == %f, max = %f, min = %f, aver = %f, lessThanAvrg = %f, moreThanvrg = %f\n",sumrows, maxrow, minrow, sumrows/nbr_grnds, lessThanAvrg, moreThanvrg);
	for(i=0;i<=(int)maxrow;i++)	if(covRes[i]>0) printf("cov[%d]=%f\t",i,covRes[i]);
	printf("\n");
	int activeuavs=0;

	FILE* fp;
	fp=fopen("uavs.csv","w");

	for(i=1;i<=net->n_uavs;i++)
	{
		if(soln[i]>0)
		{
			for (j=0;j<dim;j++)
			{
				// skip comma not needed after last dim value
				if(j==dim-1)	fprintf(fp,"%lf\n", net->uavs[i][j]);
				else fprintf(fp,"%lf,", net->uavs[i][j]);
			}
			activeuavs++;
		}
	}
	fclose(fp);

	fp=fopen("uavs_grounds.csv","w");
	int *activecovers=(int*)calloc(nbr_grnds+1,sizeof(int));
	int indice=0;// For ease in reading
	for(i=1;i<=nbr_grnds;i++)
	{
		for(j=0;j<net->g_covers[i];j++)
		{
			indice=net->covers[i][j];// find index of uav
			if(soln[indice]==0.)	continue;// skip if not active
			activecovers[i]++;
			fprintf(fp,"%lf,%lf\n", grnds[i][0], grnds[i][1]);
			fprintf(fp,"%lf,%lf\n\n", net->uavs[indice][0], net->uavs[indice][1]);
		}
	}

	printf(" Obj func : %f and nactive uavs %d\n",z,activeuavs);
	int max=1;
	double sum=activecovers[1];
	printf("Couples (ground, n active uavs covering it) : ");
	for(i=1;i<=nbr_grnds;i++)
	{
		printf(" (%d, %d) ",i,activecovers[i]);
		sum+=activecovers[i];
		if(activecovers[i] > activecovers[max])	max=i;
	}
	printf("\nAverage : %f\n",sum/nbr_grnds);
	printf(" Max degree : %d with %d\n",max,activecovers[max]);
	return soln;
};


sln* blankSolution(int size_input, int nuavs)
{
	int k=0,j=0;
	sln* res=(sln*)malloc(sizeof(sln));
	res->g_covers=(int*)calloc((size_input+1),sizeof(int));
	res->covers=(int**)malloc((size_input+1)*sizeof(int*));
	res->counters=(int*)calloc((size_input+1),sizeof(int));
	res->uavs=(double**)malloc((size_input*3+1)*sizeof(double*));
	res->n_uavs=nuavs;
	res->dists=(double**)malloc((size_input+1)*sizeof(double*));
	/* Initialisation steps */
	for(k=0;k<=size_input;k++)
	{
		res->covers[k]=(int*)calloc((size_input+1),sizeof(int));
		res->uavs[k]=(double*)calloc(dim,sizeof(double));// All uavs to origin
		for(j=0;j<dim;j++)	res->dists[k]=(double*)calloc((size_input+1),sizeof(double));
	}

	return res;
}


sln* method1ePasse(double** input_data, int size_input, double threshold)
{
	int k=2,j=1;
	sln* res=(sln*)malloc(sizeof(sln));
	res->n_uavs=0;
	res->g_covers=(int*)calloc((size_input+1),sizeof(int));
	res->covers=(int**)malloc((size_input+1)*sizeof(int*));
	res->counters=(int*)calloc((size_input+1),sizeof(int));
	//res->gr=(igraph_t*)malloc(sizeof(igraph_t));
	/* Note : allocate large memory for uavs in case duplicates are needed : size_input*3.  */
	res->uavs=(double**)malloc((size_input*3+1)*sizeof(double*));/* Convenient to keep as whole size, makes no need to reallocate each time a uav is removed/added */
	double** cl=(double**)malloc((size_input+1)*sizeof(double*));/* temporary variables, contains sum of elements in cluster */
	res->dists=(double**)malloc((size_input+1)*sizeof(double*));/* matrix of intra-distances between uavs */
// 	if k<=max_uavs_avail and there are less uavs than needed => segmentation fault. Safety:distances** squared matrix  (size_input+1)*(size_input+1)
//	for(k=0;k<=size_input;k++)
//		res->distances[k]=(double*)calloc((size_input+1),sizeof(double));

	/* Initialisation steps */
	for(k=0;k<=size_input;k++)
	{
		res->covers[k]=(int*)calloc((size_input+1),sizeof(int));
		res->uavs[k]=(double*)calloc(dim,sizeof(double));
		for(j=0;j<dim;j++)	res->uavs[k][j]=-1;// if a uav's coords is (-1,-1), then it's not deployed
		cl[k]=(double*)calloc(dim,sizeof(double));
		for(j=0;j<dim;j++)	res->dists[k]=(double*)calloc((size_input+1),sizeof(double));
	}

	/* Init first cluster with first element */
	for(k=0;k<dim;k++)
	{
		res->uavs[1][k]=input_data[1][k];
		cl[1][k]=input_data[1][k];
	}

	res->counters[1]++;
	res->covers[1][0]=1;// label of first element : first cluster
	res->g_covers[1]++;// Increment number of covers for 1st ground node
	res->n_uavs=1;// Number of clusters so far : 1

	double distance_to_closest_centroid=bound_1*bound_2;// Largest distance : the extrem limits of the entire map
	double current_distance=0;

	/* start 1e pass algorithm */
	for(k=2;k<=size_input;k++)
	{
		distance_to_closest_centroid=bound_1*bound_2;// Worst : map's limits
		for(j=1;j<=res->n_uavs;j++)
		{
			current_distance=euclDistance(input_data[k],res->uavs[j]);
			if ( current_distance < distance_to_closest_centroid )
			{
				distance_to_closest_centroid=current_distance;// find the closest cluster
				res->covers[k][0]=j;
			}
		}

		if(distance_to_closest_centroid<threshold)
		{// closest centroid is within admissible range => add element to cluster
			res->counters[res->covers[k][0]]++;/* Increment the number of elements covered by the uav given by the first stored uav in the list of covers of ground node k (covers[k][res->g_covers[k]]) */
			/* update centroids */
			for(j=0;j<dim;j++)
			{
				cl[res->covers[k][0]][j]+=input_data[k][j];
				res->uavs[res->covers[k][0]][j]=cl[res->covers[k][0]][j]/res->counters[res->covers[k][0]];
			}
		}
		else
		{// closest centroid is not within admissible range => Create own new cluster
			res->n_uavs++;
			res->covers[k][0]=res->n_uavs;
			res->g_covers[k]++;// Increment number of covers for k-th ground node
			for(j=0;j<dim;j++)
			{
				cl[res->covers[k][0]][j]=input_data[k][j];
				res->uavs[res->covers[k][0]][j]=input_data[k][j];
			}
			res->counters[res->n_uavs]++;// To [1]
		}
	}

	updateDistMat(res, threshold);

	/* Housekeeping */
	for(k=0;k<=size_input;k++)
	{
		free(cl[k]);
	}
	free(cl);

	return res;
};


double k_means(double** data, int n, sln* clusts, double error_tolerance, double range)
{

    int h, i, j; /* loop counters, of course :) */

    double old_error, error = DBL_MAX; /* sum of squared euclidean distance. DBL_MAX : macro defines the maximum finite floating-point value : 1E+37*/
    double** c1=(double**) malloc((clusts->n_uavs+1)*sizeof(double*)); /* temp centroids */

    assert(data && clusts->n_uavs > 0 && clusts->n_uavs <= n && error_tolerance >= 0); /* for debugging */

    /* initialization : start with given centroids */

     // c1 : temp centroids, clusts->uavs : true centroids
	for (i=0; i<=clusts->n_uavs; i++)
		c1[i]=(double*)calloc(dim, sizeof(double));

	double min_distance = DBL_MAX, distance = 0;

    /** Kmeans loop */
    do {
        /* save error from last step */
		old_error = error;
        error = 0;

        /* clear old counts and temp centroids */
        for (i=1; i <=clusts->n_uavs; clusts->counters[i++] = 0)
        {
			// reinit temp centroids
            for (j=0;j<dim;c1[i][j++] = 0);
        }

        for (h = 1; h <= n; h++)
        {
            /* Find closest cluster */
            min_distance = DBL_MAX;
            for (i=1; i<=clusts->n_uavs; i++)
            {
				distance=0;
				/* considered distance : euclidian squared or squared residuals */
                for(j=0;j<dim;j++)	distance+=pow(data[h][j]-clusts->uavs[i][j],2);
                if (distance < min_distance)
                {
					min_distance=distance;// find the closest cluster
                    clusts->covers[h][0] = i;
                }
            }
            /* update size and temp centroid of the destination cluster */
			for(j=0;j<dim;j++)	c1[clusts->covers[h][0]][j] += data[h][j];
			clusts->counters[clusts->covers[h][0]]++;
            /* update standard error */
            error += min_distance;
        }

		/* update all centroids */
        for (i=1; i<=clusts->n_uavs; i++)
            for (j=0; j<dim; j++)
                clusts->uavs[i][j] = clusts->counters[i] ? c1[i][j]/clusts->counters[i] : c1[i][j];

    } while (fabs(error - old_error) > error_tolerance);/* if for each iteration, the number of changes made are not different from previous */

	updateDistMat(clusts, range);

	/* housekeeping */
	for(h=0;h<=clusts->n_uavs;h++)	free(c1[h]);
	free(c1);

	return error;
}


void connect_CCs(sln* net, igraph_t* Gk, double threshold, int** restr_list, int **pairs, int* npairs, bool restrict, int level)
{// The boolean value "restrict" says wether the restriction strategy should be used (restrict=false since graph G isn't the root graph G0) or not.
	// In the latter case, fill lists : **restr_list, **pairs, and * npairs. Otherwise, use these lists to apply the restrictions

printf("net n_uavs %d threshold %f and restrict %d\n", net->n_uavs, threshold, restrict);
//printf(" vcount %li ecount %li\n", (long int)igraph_vcount(Gk), (long int)igraph_ecount(Gk));

	igraph_integer_t ncomps=-1;// numbers of connected components
	igraph_vector_t labels, compssizes;// labels of connected components for each vertex
	double min_distance=DBL_MAX;// distance of two closest nodes but out of range (belong to two different connected components)
	double current_distance;// used to place a new uav
	int buffn1=-1,buffn2=-1, n1=-1, n2=-1;//two closest but out of range nodes

	// Use of a brute force branching here. Maybe undeterministic could be better.
	// brute force : create one unique connected component by iteratively branching the 2 closest connected components
	// From here, a bit of a mess, does even create a segmentation fault
	int i=0, j=0, k=0, added=0;
	// Note : avoid doing : bool restricted=false; and restrict once only. If done so, it's possible that on the next phases restricted (n1,n2)
	// is included anyway and that we end up with a graph that was already generated. Then : always restrict (n1,n2) in new G.
	
	do
	{
printf("n_uavs %d, k = %d, %d\n", net->n_uavs, k++, ncomps);
/*
if(k++==0)
{// keep topology of initial graph in a file
	FILE* fout;
	fout=fopen("first_graph.csv","w");
	igraph_vector_t ite_edgs_first_graph;
	igraph_vector_init(&ite_edgs_first_graph, 0);
	igraph_get_edgelist(G1, &ite_edgs_first_graph, 0);
	int nfirst=igraph_ecount(G1);
	long int ind1first,ind2first;// needed for the cast
	for (i=0, j=0; j<nfirst; i+=2, j++)
	{
		ind1first=VECTOR(ite_edgs_first_graph)[i];
		ind2first=VECTOR(ite_edgs_first_graph)[i+1];
		fprintf(fout,"%lf,%lf\n", net->uavs[ind1first][0], net->uavs[ind1first][1]);
		fprintf(fout,"%lf,%lf\n", net->uavs[ind2first][0], net->uavs[ind2first][1]);
		fprintf(fout,"\n");
printf("(%ld,%ld)\t", ind1first, ind2first);
//fprintf(bufffp,"%ld-%ld:%lf,%lf:%lf,%lf:W:%lf\n", buff1, buff2, net->uavs[buff1][0], net->uavs[buff1][1], net->uavs[buff2][0], net->uavs[buff2][1],weight);
	}
printf("\n");
	fclose(fout);
	igraph_vector_destroy(&ite_edgs_first_graph);
}
*/
		if(igraph_vector_init(&labels, 0)==IGRAPH_ENOMEM || igraph_vector_init(&compssizes, 0)==IGRAPH_ENOMEM)
			printf(" Memory issues, vector init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);
			
		// compute number of connected components
		igraph_clusters(Gk, &labels, &compssizes, &ncomps, IGRAPH_WEAK);

//printf("CCs Graph %d has %d connected components\n", k++, ncomps);
//printf(" vcount %li ecount %li\n", (long int)igraph_vcount(Gk), (long int)igraph_ecount(Gk));

//printf("start labels\n");
//for(i=0;i<igraph_vcount(&net->gr);i++)	printf("l[%li]:%li\t",i,(long int)VECTOR(labels)[i]);
//printf("end labels\n");

		// free first used memory since initiated again later
		igraph_vector_destroy(&compssizes);

		// if graph not a single connected component, then update matrix of distance with new uavs. Note : exists isolated vertex 0
		if(ncomps>2)
		{//add new node to reach one connected component
// printf("In ncomps %d connected components\n", ncomps);
			min_distance=DBL_MAX;
			// find two closest but out of range nodes
			for (i=1;i<=net->n_uavs;i++)
			{
				for (j=i+1;j<=net->n_uavs;j++)
				{
					// Not interested in nodes in same connected component
					if ( VECTOR(labels)[i] == VECTOR(labels)[j] )	continue;
					/* Find closest clusters but out of range */
					current_distance=net->dists[i][j];
//					if(current_distance<threshold)// skip, in range
//						continue;
					// check only if different connected component
					if(current_distance<min_distance)
					{// keep two nodes
						buffn1=i;
						buffn2=j;
						if( restrict && pairs[level][0] == buffn1 && pairs[level][1] == buffn2 )
						{
printf("i j n1 n2 pairs[level][0] pairs[level][1] level : %d %d %d %d %d %d %d\n", i, j, n1, n2, pairs[level][0], pairs[level][1], level);
printf("level restricted pairs[level][0] == n1, pairs[level][1] == n2 : %d %d %d %d\n", level, restrict, pairs[level][0] == n1, pairs[level][1] == n2);
//printf("in restricted : i j n1 n2 net->n_uavs : %d %d %d %d %d\n", i, j, n1, n2, net->n_uavs);
							continue;// restrict
						}
//printf("not restricted and n1 n2 ncomps : %d %d %d\n", n1, n2, ncomps);
						min_distance=current_distance;
						n1=i;
						n2=j;
					}
				}
			}
			
			if( ! restrict )
			{
				// first iteration, keep track of n1 and n2 for further potential restrictions
				pairs[*npairs][0]=n1;
				pairs[*npairs][1]=n2;
				*npairs=*npairs+1;
// printf("In ! restrict : %d connected components\n", *npairs);
			}
/*
long int zz=igraph_vcount(Gk);
for(long int z=0;z<zz;z++)	printf("%li ",VECTOR(labels)[z]);
printf("\n");
*/
			// free vector of labels
			igraph_vector_destroy(&labels);

			if(n1<0 || n2<0)
			{
				printf(" Something went wrong with graph (Either connectivity or indices. Check), adding new uav will fail %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);
				return;
			}

			// assign coordinates of new created node : middle of segment [n1,n2]
			net->n_uavs++;
			// assign coordinates
			for (j=0;j<dim;j++)
				net->uavs[net->n_uavs][j]=(net->uavs[n1][j]+net->uavs[n2][j])/2;
			// increase number of elements in new cluster : at least itself, since may not cover any ground node, but be used only for connectivity
			// extend size of graph
//			igraph_add_vertices(&net->gr, 1, 0);
			igraph_add_vertices(Gk, 1, 0);			
			// update new distances
			updateDistMat(net, threshold);

			// fill in various lists
			if( ! restrict )
			{
				net->counters[net->n_uavs]++;
				restr_list[added][0]=net->n_uavs;
			}
			int col=1;// gives current column position in restricted list

			// create links with other uavs when needed
			for (j=1;j<net->n_uavs;j++)
			{
				if( net->dists[net->n_uavs][j] < threshold )
				{
					igraph_add_edge(Gk, net->n_uavs, j);
					if( ! restrict )
					{
						restr_list[added][col]=j;
						col++;
					}
				}
			}
			if( ! restrict )	added++;// finished adding links to the row => increment

//			if(net->dists[net->n_uavs][n1]<threshold && solverSln[net->n_uavs] > 0 && solverSln[n1] > 0)
//			if( net->dists[net->n_uavs][n1] < threshold )
//{
//printf("(%d,%d)\n",net->n_uavs, n1);
//				igraph_add_edge(G1, net->n_uavs, n1);
//}
//			if(net->dists[net->n_uavs][n2]<threshold && solverSln[net->n_uavs] > 0 && solverSln[n2] > 0)
//			if( net->dists[net->n_uavs][n2] < threshold )
//{
//printf("(%d,%d)\n",net->n_uavs, n2);
//				igraph_add_edge(G1, net->n_uavs, n2);
//}
//printf("1. dist n_uav-n1 : %lf, thresh : %lf, nuav : %d, n1 : %d\n", net->dists[net->n_uavs][n1], threshold, net->n_uavs, n1);
//printf("2. dist n_uav-n2 : %lf, thresh : %lf, nuav : %d, n2 : %d\n", net->dists[net->n_uavs][n2], threshold, net->n_uavs, n2);
//printf(" %d and %d coords (%lf,%lf) (%lf,%lf)\n", n1, n2, net->uavs[n1][0], net->uavs[n1][1], net->uavs[n2][0], net->uavs[n2][1]);
//printf("dist n1-n2 after update %d - %d\n",(int)net->dists[n2][n1],(int)net->dists[n1][n2]);

//printf("with %li vertices, and new uav %d coords (%lf,%lf)\n", (long int)igraph_vcount(&net->gr), net->n_uavs, net->uavs[net->n_uavs][0], net->uavs[net->n_uavs][1]);
//printf(" vcount %li ecount %li\n", (long int)igraph_vcount(&net->gr), (long int)igraph_ecount(&net->gr));
		}

//	}while(ncomps>2 && k++ < 5);
	}while(ncomps>2);

//printf("final (%f,%f)\n", net->uavs[net->n_uavs][0], net->uavs[net->n_uavs][1]);

/*
	FILE* fout;
	fout=fopen("test_graph.csv","w");
	igraph_vector_t ite_edgs_first_graph;
	igraph_vector_init(&ite_edgs_first_graph, 0);
	igraph_get_edgelist(Gk, &ite_edgs_first_graph, 0);
	int nfirst=igraph_ecount(Gk);
	long int ind1first,ind2first;// needed for the cast
	for (i=0, j=0; j<nfirst; i+=2, j++)
	{
		ind1first=VECTOR(ite_edgs_first_graph)[i];
		ind2first=VECTOR(ite_edgs_first_graph)[i+1];
		fprintf(fout,"%lf,%lf\n", net->uavs[ind1first][0], net->uavs[ind1first][1]);
		fprintf(fout,"%lf,%lf\n", net->uavs[ind2first][0], net->uavs[ind2first][1]);
		fprintf(fout,"\n");
//printf("(%ld,%ld)\t", ind1first, ind2first);
//fprintf(bufffp,"%ld-%ld:%lf,%lf:%lf,%lf:W:%lf\n", buff1, buff2, net->uavs[buff1][0], net->uavs[buff1][1], net->uavs[buff2][0], net->uavs[buff2][1],weight);
	}
//printf("\n");
//printf("Tests vcount %li ecount %li, nuavs %d\n", (long int)igraph_vcount(Gk), (long int)igraph_ecount(Gk), net->n_uavs);
	fclose(fout);
	igraph_vector_destroy(&ite_edgs_first_graph);
*/

};


igraph_t* translate(sln* net, double threshold, double* solverSln)
{// Changes the coordinates of the active UAVs into a graph for the connectivity checking
	// note : computing diameter of graph (O(n*|E|)) cost more than comparing each nodes (O(n^2)) : Complete graph : n*|E| = n*sum_i_in{1...n-1}(i) > O(n^2) for n>3

	long int i=0,j=0;
	igraph_integer_t ncomps=-1;// numbers of connected components
	igraph_vector_t labels, compssizes;// labels of connected components for each vertex
//	double min_distance=DBL_MAX;// distance of two closest nodes but out of range (belong to two different connected components)
//	double current_distance;// used to place a new uav
//	double shared=threshold/4; "shared" is the degree of "merging" of two uavs where the area is shared
//	int n1=-1,n2=-1;//two closest but out of range nodes

	// here, the vectors "labels" are only used to create graph
	int nActivUAVs=0;// used to keep record of the number of uavs
	int indexActivUavs[net->n_uavs];// contains the indices of active uavs
	
	for(i=1;i<=net->n_uavs;i++)
		if(solverSln[i]>0)
		{/* register indices of active uavs (solverSln[i]>0) */
			indexActivUavs[nActivUAVs]=i;
			nActivUAVs++;
		}

	igraph_t* buffgr=(igraph_t*)malloc(sizeof(igraph_t));// buffer graph, temporary, used until remains a unique connected component. 
	
	/* start building graph */
	if(igraph_empty(buffgr, nActivUAVs, IGRAPH_UNDIRECTED)==IGRAPH_EINVAL)
		printf(" Invalid number of vertices, graph init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);
//	if(igraph_empty(&net->gr, net->n_uavs+1, IGRAPH_UNDIRECTED)==IGRAPH_EINVAL)
//		printf(" Invalid number of vertices, graph init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);


printf(" %d active uavs\n", nActivUAVs);
	// add edges to graph : two uavs in each other range
	for(i=0;i<nActivUAVs;i++)
		for(j=i+1;j<nActivUAVs;j++)
			if(net->dists[indexActivUavs[i]][indexActivUavs[j]]<threshold)
{
printf("Enter ---------  (%ld,%ld), (%ld,%ld)\t",i,j,(long int)indexActivUavs[i], (long int)indexActivUavs[j]);
				igraph_add_edge(buffgr, i, j);
}
printf("end init graph with %d nodes\n",nActivUAVs);

		if(igraph_vector_init(&labels, 0)==IGRAPH_ENOMEM || igraph_vector_init(&compssizes, 0)==IGRAPH_ENOMEM)
			printf(" Memory issues, vector init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);			
		// compute number of connected components
		igraph_clusters(buffgr, &labels, &compssizes, &ncomps, IGRAPH_WEAK);


/*
	// add edges to graph, edges : two uavs in range
	for(i=0;i<=net->n_uavs;i++)
		for(j=i+1;j<=net->n_uavs;j++)
//			if(net->dists[i][j]<2*threshold-shared)
			if(net->dists[i][j]<threshold && solverSln[i] > 0 && solverSln[j] > 0)
{
printf("Enter ---------  (%ld,%ld)\t",i,j);
				igraph_add_edge(&net->gr, i, j);
}
printf("end init graph with %d uavs\n",net->n_uavs);

		if(igraph_vector_init(&labels, 0)==IGRAPH_ENOMEM || igraph_vector_init(&compssizes, 0)==IGRAPH_ENOMEM)
			printf(" Memory issues, vector init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);			
		// compute number of connected components
		igraph_clusters(&net->gr, &labels, &compssizes, &ncomps, IGRAPH_WEAK);
*/

printf("Graph 0 has %d connected components\n",ncomps);
		return buffgr;

// Uncomment if have found better result, note : undeterministic better
// concerns creating unique connected components by branching 2 closest connected comp
// From here, a bit of a mess, does even create a segmentation fault
/*
//int k=0;
	do
	{
		
		if(igraph_vector_init(&labels, 0)==IGRAPH_ENOMEM || igraph_vector_init(&compssizes, 0)==IGRAPH_ENOMEM)
			printf(" Memory issues, vector init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);
			
		// compute number of connected components
		igraph_clusters(&net->gr, &labels, &compssizes, &ncomps, IGRAPH_WEAK);
//printf("start labels\n");
//for(i=0;i<igraph_vcount(&net->gr);i++)	printf("l[%li]:%li\t",i,(long int)VECTOR(labels)[i]);
//printf("end labels\n");
		// free first used memory since initiated again later
		igraph_vector_destroy(&compssizes);

		// if graph not one unique connected component then update matrix of distance with new uav (2 as exist isolated vertex 0)
		if(ncomps>2)
		{//add new node to reach one connected component
			min_distance=DBL_MAX;
			// find two closest but out of range nodes
			for (i=1;i<=net->n_uavs;i++)
			{
				for (j=i+1;j<=net->n_uavs && VECTOR(labels)[i] != VECTOR(labels)[j];j++)
				{
*/
					/* Find closest clusters but over range */
/*					current_distance=net->dists[i][j];
//					if(current_distance<threshold)// skip, in range
//						continue;
					// check only if different connected component
//					if(current_distance<min_distance && VECTOR(labels)[i] != VECTOR(labels)[j])
					// Reduncancy!!! If G is well built : all i, j belonging to two different components, then must be out of range
					if(current_distance<min_distance)
					{// keep two nodes
						min_distance=current_distance;
						n1=i;
						n2=j;
					}
				}
			}

			// free vector of labels
			igraph_vector_destroy(&labels);

			if(n1<0 || n2<0)
			{
				printf(" Something went wrong, adding new uav will fail %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);
				return;
			}

			// assign coordinates of new created node
			net->n_uavs++;
			// assign coordinates
			for (j=0;j<dim;j++)	net->uavs[net->n_uavs][j]=(net->uavs[n1][j]+net->uavs[n2][j])/2;
			// increase number of elements in new cluster
			net->counters[net->n_uavs]++;// To 1
			// extend size of graph
			igraph_add_vertices(&net->gr, 1, 0);
			// update new distances
			updateDistMat(net, threshold);
			// create links where needed
//			if(net->dists[net->n_uavs][n1]<2*threshold-shared)
			if(net->dists[net->n_uavs][n1]<threshold && solverSln[net->n_uavs] > 0 && solverSln[n1] > 0)
{
printf("(%d,%d)\n",net->n_uavs, n1);
				igraph_add_edge(&net->gr, net->n_uavs, n1);
}
//			if(net->dists[net->n_uavs][n2]<2*threshold-shared)
			if(net->dists[net->n_uavs][n2]<threshold && solverSln[net->n_uavs] > 0 && solverSln[n2] > 0)
{
printf("(%d,%d)\n",net->n_uavs, n2);
				igraph_add_edge(&net->gr, net->n_uavs, n2);
}
//printf("1. dist n_uav-n1 : %lf, thresh : %lf, nuav : %d, n1 : %d\n", net->dists[net->n_uavs][n1], 2*threshold-shared, net->n_uavs, n1);
printf("1. dist n_uav-n1 : %lf, thresh : %lf, nuav : %d, n1 : %d\n", net->dists[net->n_uavs][n1], threshold, net->n_uavs, n1);
//printf("2. dist n_uav-n2 : %lf, thresh : %lf, nuav : %d, n2 : %d\n", net->dists[net->n_uavs][n2], 2*threshold-shared, net->n_uavs, n2);
printf("2. dist n_uav-n2 : %lf, thresh : %lf, nuav : %d, n2 : %d\n", net->dists[net->n_uavs][n2], threshold, net->n_uavs, n2);
printf(" %d and %d coords (%lf,%lf) (%lf,%lf)\n",n1,n2, net->uavs[n1][0], net->uavs[n1][1], net->uavs[n2][0], net->uavs[n2][1]);
printf("dist n1-n2 after update %d - %d\n",(int)net->dists[n2][n1],(int)net->dists[n1][n2]);

//printf(" %d and %d coords (%lf,%lf) (%lf,%lf)\n",n1,n2, net->uavs[n1][0], net->uavs[n1][1], net->uavs[n2][0], net->uavs[n2][1]);
//printf("with %li vertices, and new uav %d coords (%lf,%lf)\n", (long int)igraph_vcount(&net->gr), net->n_uavs, net->uavs[net->n_uavs][0], net->uavs[net->n_uavs][1]);
//printf(" vcount %li ecount %li\n", (long int)igraph_vcount(&net->gr), (long int)igraph_ecount(&net->gr));
		}
//	}while(ncomps>2 && k++ < 5);
	}while(ncomps>2);

*/

/*
	// write graph into file
	FILE* fp;
	fp=fopen("graphsNewCorrected.csv","w");
	igraph_vector_t ite_edgs;
	igraph_vector_init(&ite_edgs, 0);
	igraph_get_edgelist(&net->gr, &ite_edgs, 0);
	int n=igraph_ecount(&net->gr);
	long int buff1,buff2;// needed for the cast
	for (i=0, j=0; j<n; i+=2, j++)
	{
		buff1=VECTOR(ite_edgs)[i];
		buff2=VECTOR(ite_edgs)[i+1];
		fprintf(fp,"%lf,%lf\n", net->uavs[buff1][0], net->uavs[buff1][1]);
		fprintf(fp,"%lf,%lf\n", net->uavs[buff2][0], net->uavs[buff2][1]);
		fprintf(fp,"\n");
printf("(%ld,%ld)\t",buff1,buff2);
//fprintf(bufffp,"%ld-%ld:%lf,%lf:%lf,%lf:W:%lf\n", buff1, buff2, net->uavs[buff1][0], net->uavs[buff1][1], net->uavs[buff2][0], net->uavs[buff2][1],weight);
	}
printf("\n");
	fclose(fp);
	igraph_vector_destroy(&ite_edgs);
*/
};


void populate(sln* covers, double threshold, igraph_t* solnG0, int** restr_list, int **pairs, int* npairs)
{// note : computing diameter of graph (O(n*|E|)) can cost more than comparing each nodes (O(n^2)) : eg : Complete graph : n*|E| = n*sum_i_in{1...n-1}(i) > O(n^2) for n>3

	long int i=0, j=0, k=0;

	/* the goal of the step is to generate the best NINDIVS by unpiling (Starting by the last added) the list of pairs of
	 * vertices (edges) added to the first graph (solnsGraphs[0]) in order to create the connectivity. The restriction on
	 * the level of pairs unpiled are given by k, and the number of vertices to start with on each iteration is n_in_first_graph-k-1
	 * where n_in_first_graph is the number of vertices in the first connected graph. */
	// long int n_in_first_graph=(long int)igraph_vcount(solnsGraphs[0]);
	// Use of a brute force branching here. Maybe undeterministic could be better.

	// 1st phase : generate 1st connected graph to start with
	solnG0=(igraph_t*)malloc(sizeof(igraph_t));// first connected graph to build.

	// start building graph
	if(igraph_empty(solnG0, covers->n_uavs+1, IGRAPH_UNDIRECTED)==IGRAPH_EINVAL)
		printf(" Invalid number of vertices, graph init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);

printf("G0 with %d active uavs\n", covers->n_uavs);
	// add edges to graph : two uavs in each other range
int count=0;
	for(i=1;i<=covers->n_uavs;i++)
		for(j=i+1;j<=covers->n_uavs;j++)
			if(covers->dists[i][j] < threshold)
{
//printf("Enter ---------  (%ld,%ld) and dist=%f\t",i,j,covers->dists[i][j]);
count++;
				igraph_add_edge(solnG0, i, j);// since j = i+1 to nActivUAVs, then no loop
}
printf("\nend 1st graph and count %d\n",count);
printf("Before connectCCs vcount %li ecount %li and n_uavs %d\n", (long int)igraph_vcount(solnG0), (long int)igraph_ecount(solnG0), covers->n_uavs);

	connect_CCs(covers, solnG0, threshold, restr_list, pairs, npairs, false, 0);
	
printf("Test npairs == %d and connected covers : %d\n", *npairs, covers->n_uavs);

	// Now that more variables are assigned => launch true populating phase
	// solnsGraphs=(igraph_t**)malloc((*npairs)*sizeof(igraph_t*));
	igraph_t* Gk;
	igraph_vector_t edges;
	int nfirst=0;
	long int ind1,ind2;// needed for the cast
	sln* buffnets;
	char path[30];
	char buff[30];
	FILE *fout;

	// Do not forget to record first the data from the first graph (graph of reference) 
		strcpy(path,"./out/G_0");
		fout=fopen(path,"w");
		igraph_vector_init(&edges, 0);
		igraph_get_edgelist(solnG0, &edges, 0);
		nfirst=igraph_ecount(solnG0);
		for (i=0, j=0; j<nfirst; i+=2, j++)
		{
			ind1=VECTOR(edges)[i];
			ind2=VECTOR(edges)[i+1];
			fprintf(fout,"%lf,%lf\n", covers->uavs[ind1][0], covers->uavs[ind1][1]);
			fprintf(fout,"%lf,%lf\n", covers->uavs[ind2][0], covers->uavs[ind2][1]);
			fprintf(fout,"\n");
//	printf("(%ld,%ld)\t", ind1, ind2);
//fprintf(bufffp,"%ld-%ld:%lf,%lf:%lf,%lf:W:%lf\n", buff1, buff2, net->uavs[buff1][0], net->uavs[buff1][1], net->uavs[buff2][0], net->uavs[buff2][1],weight);
		}
//	printf("\n");
		fclose(fout);
		igraph_vector_destroy(&edges);
		
		strcat(path,"_coords");
		fout=fopen(path,"w");
		for(i=1;i<=covers->n_uavs;i++)
			for (j=0;j<dim;j++)
			{
				// skip comma not needed after last dim value
				if(j==dim-1)	fprintf(fout,"%lf\n", covers->uavs[i][j]);
				else fprintf(fout,"%lf,", covers->uavs[i][j]);
			}
		fclose(fout);
		
		find_covers(covers, threshold);
		
int overallredundancy=0;
// Compute redundancy of G_0
for(i=1;i<=nbr_grnds;i++)	overallredundancy+=covers->g_covers[i];
printf("G_0 has overall redundancy of %d\n",overallredundancy);

		if(*npairs%2==0)
		{// we want n even number of generated files, if there are npairs of newly generated files, if added to the starting graph, we are still
			// short of one file => easiest solution : copy starting solution in file named [blabla]_[npairs+1]
			strcpy(path,"./out/G_");
			sprintf(buff, "%d", *npairs+1);
			strcat(path,buff);
			fout=fopen(path,"w");
			igraph_vector_init(&edges, 0);
			igraph_get_edgelist(solnG0, &edges, 0);
			nfirst=igraph_ecount(solnG0);
			for (i=0, j=0; j<nfirst; i+=2, j++)
			{
				ind1=VECTOR(edges)[i];
				ind2=VECTOR(edges)[i+1];
				fprintf(fout,"%lf,%lf\n", covers->uavs[ind1][0], covers->uavs[ind1][1]);
				fprintf(fout,"%lf,%lf\n", covers->uavs[ind2][0], covers->uavs[ind2][1]);
				fprintf(fout,"\n");
	//	printf("(%ld,%ld)\t", ind1, ind2);
	//fprintf(bufffp,"%ld-%ld:%lf,%lf:%lf,%lf:W:%lf\n", buff1, buff2, net->uavs[buff1][0], net->uavs[buff1][1], net->uavs[buff2][0], net->uavs[buff2][1],weight);
			}
	//	printf("\n");
			fclose(fout);
			igraph_vector_destroy(&edges);
			
			strcat(path,"_coords");
			fout=fopen(path,"w");
			for(i=1;i<=covers->n_uavs;i++)
				for (j=0;j<dim;j++)
				{
					// skip comma not needed after last dim value
					if(j==dim-1)	fprintf(fout,"%lf\n", covers->uavs[i][j]);
					else fprintf(fout,"%lf,", covers->uavs[i][j]);
				}
			fclose(fout);
		}

	for (k=0; k < *npairs; k++)
	{
		buffnets=(sln*)malloc(sizeof(sln));
		// some values are not needed for filling (and allocation) since only the position of uavs (and distance between them) are wanted, as well
		// as the connections between them. Other operations (eg:what are the ground nodes covered by a given uav ?) are performed in the MLMPGA
		buffnets->n_uavs=covers->n_uavs-k-1;
//printf("net->n_uavs %d and covers->n_uavs-k-1 %ld and k %ld\n", buffnets->n_uavs, covers->n_uavs-k-1, k);
		buffnets->uavs=(double**)malloc((covers->n_uavs*5)*sizeof(double*));
		buffnets->dists=(double**)malloc((covers->n_uavs*5)*sizeof(double*));
		for(i=0;i<(covers->n_uavs*5);i++)
		{
			buffnets->uavs[i]=(double*)calloc(dim,sizeof(double));// All uavs to origin
			buffnets->dists[i]=(double*)calloc((covers->n_uavs*5),sizeof(double));
			for(j=0;j<dim;j++)
			{
				if(i <= buffnets->n_uavs){	buffnets->uavs[i][j]=covers->uavs[i][j];}
				else{	buffnets->uavs[i][j]=-1;}
			}
		}

		updateDistMat(buffnets, threshold);
		
		Gk=(igraph_t*)malloc(sizeof(igraph_t));

		if(igraph_empty(Gk, buffnets->n_uavs+1, IGRAPH_UNDIRECTED)==IGRAPH_EINVAL)
			printf("Ite %d Invalid number of vertices, graph initialisation failed %s, %d, %s \n", (int)k, __FILE__, __LINE__, __FUNCTION__);

//printf("Phase 1 : Edges found while comparing distances\n");
/*
count=0;
		for(i=1;i<=buffnets->n_uavs;i++)
			for(j=i+1;j<=buffnets->n_uavs;j++)
				if(buffnets->dists[i][j]<threshold)
				{
//printf("(%ld -- %ld)\t",i,j);
count++;
					igraph_add_edge(Gk, i, j);// since j = i+1 to nActivUAVs, then no loop
				}
*/
//printf("\nresult count %d\n", count);


//printf("Phase 2 : Edges found from igraph\n");
		igraph_vector_init(&edges, 0);
		igraph_get_edgelist(solnG0, &edges, 0);
		nfirst=igraph_ecount(solnG0);
count=0;
		// if an edge from graph G0 connects a vertex not in the new graph then skip
		for (i=0, j=0; j<nfirst; i+=2, j++)
		{
			if( VECTOR(edges)[i] <= buffnets->n_uavs && VECTOR(edges)[i+1] <= buffnets->n_uavs )
{
				igraph_add_edge(Gk, VECTOR(edges)[i], VECTOR(edges)[i+1]);
//printf("\n[(%ld -- %ld)\t]", (long int)VECTOR(edges)[i], (long int)VECTOR(edges)[i+1]);
}

//printf("(%ld -- %ld)\t", (long int)VECTOR(edges)[i], (long int)VECTOR(edges)[i+1]);
count++;
//printf("(%ld -- %ld)\t", (long int)VECTOR(edges)[i], (long int)VECTOR(edges)[i+1]);
		}
//printf("\n");
//printf("\ncount new results : %d\n", count);

//printf("tests values vcount %li ecount %li\n", (long int)igraph_vcount(Gk), (long int)igraph_ecount(Gk));

		igraph_vector_destroy(&edges);

printf("K ite %li\n", k);

		connect_CCs(buffnets, Gk, threshold, restr_list, pairs, npairs, true, *npairs-1-k);

		// write into file to save memory
		strcpy(path,"./out/G_");
		sprintf(buff, "%ld", k+1);
		strcat(path,buff);
		fout=fopen(path,"w");
		printf("%s\n",path);

//printf("Ite %li with %d\n", k, buffnets->n_uavs);

		igraph_vector_init(&edges, 0);
		igraph_get_edgelist(Gk, &edges, 0);
		nfirst=igraph_ecount(Gk);
		for (i=0, j=0; j<nfirst; i+=2, j++)
		{
			ind1=VECTOR(edges)[i];
			ind2=VECTOR(edges)[i+1];
			fprintf(fout,"%lf,%lf\n", buffnets->uavs[ind1][0], buffnets->uavs[ind1][1]);
			fprintf(fout,"%lf,%lf\n", buffnets->uavs[ind2][0], buffnets->uavs[ind2][1]);
			fprintf(fout,"\n");
//	printf("(%ld,%ld)\t", ind1, ind2);
//fprintf(bufffp,"%ld-%ld:%lf,%lf:%lf,%lf:W:%lf\n", buff1, buff2, net->uavs[buff1][0], net->uavs[buff1][1], net->uavs[buff2][0], net->uavs[buff2][1],weight);
		}
//	printf("\n");
		fclose(fout);
		igraph_vector_destroy(&edges);
		free(Gk);
		
		strcat(path,"_coords");
		fout=fopen(path,"w");

		for(i=1;i<=buffnets->n_uavs;i++)
			for (j=0;j<dim;j++)
			{
				// skip comma not needed after last dim value
				if(j==dim-1)	fprintf(fout,"%lf\n", buffnets->uavs[i][j]);
				else fprintf(fout,"%lf,", buffnets->uavs[i][j]);
			}
		fclose(fout);

		// Housekeeping
		for(i=0;i<(covers->n_uavs*5);i++)
		{
			//for(j=0;j<dim;j++)	free(buffnets[k]->uavs[i][j]);
			free(buffnets->uavs[i]);
			free(buffnets->dists[i]);
		}
		free(buffnets->uavs);
		free(buffnets->dists);
		free(buffnets);

	}

/*
if(++k==1)
{// keep topology of initial graph
	FILE* fout;
	fout=fopen("first_graph.csv","w");
	igraph_vector_t ite_edgs_first_graph;
	igraph_vector_init(&ite_edgs_first_graph, 0);
	igraph_get_edgelist(&G1, &ite_edgs_first_graph, 0);
	int nfirst=igraph_ecount(&G1);
	long int ind1first,ind2first;// needed for the cast
	for (i=0, j=0; j<nfirst; i+=2, j++)
	{
		ind1first=VECTOR(ite_edgs_first_graph)[i];
		ind2first=VECTOR(ite_edgs_first_graph)[i+1];
		fprintf(fout,"%lf,%lf\n", net->uavs[ind1first][0], net->uavs[ind1first][1]);
		fprintf(fout,"%lf,%lf\n", net->uavs[ind2first][0], net->uavs[ind2first][1]);
		fprintf(fout,"\n");
printf("(%ld,%ld)\t", ind1first, ind2first);
//fprintf(bufffp,"%ld-%ld:%lf,%lf:%lf,%lf:W:%lf\n", buff1, buff2, net->uavs[buff1][0], net->uavs[buff1][1], net->uavs[buff2][0], net->uavs[buff2][1],weight);
	}
printf("\n");
	fclose(fout);
	igraph_vector_destroy(&ite_edgs_first_graph);
}


	// write graph into file
	FILE* fp;
	fp=fopen("genet_clust.csv","w");
	igraph_vector_t ite_edgs;
	igraph_vector_init(&ite_edgs, 0);
	igraph_get_edgelist(&G1, &ite_edgs, 0);
	int n=igraph_ecount(&G1);
	long int ind1,ind2;// needed for the cast
	for (i=0, j=0; j<n; i+=2, j++)
	{
		ind1=VECTOR(ite_edgs)[i];
		ind2=VECTOR(ite_edgs)[i+1];
		fprintf(fp,"%lf,%lf\n", net->uavs[ind1][0], net->uavs[ind1][1]);
		fprintf(fp,"%lf,%lf\n", net->uavs[ind2][0], net->uavs[ind2][1]);
		fprintf(fp,"\n");
printf("(%ld,%ld)\t", ind1, ind2);
//fprintf(bufffp,"%ld-%ld:%lf,%lf:%lf,%lf:W:%lf\n", buff1, buff2, net->uavs[buff1][0], net->uavs[buff1][1], net->uavs[buff2][0], net->uavs[buff2][1],weight);
	}
printf("\n");
	fclose(fp);
	igraph_vector_destroy(&ite_edgs);
*/

};


void freeSln(sln* free_sln)
{
	if(free_sln->g_covers != NULL)
		free(free_sln->g_covers);
	if(free_sln->covers != NULL)
		free(free_sln->covers);
	if(free_sln->counters != NULL)
		free(free_sln->counters);
	//if(free_sln->gr != NULL)
	//	free(free_sln->gr);
	if(free_sln->uavs != NULL)
	{
		int j=0;
		for(;j<=free_sln->n_uavs;j++)
			free(free_sln->uavs[j]);
		free(free_sln->uavs);
	}
	if(free_sln->dists != NULL)
		free(free_sln->dists);
	//igraph_destroy(&free_sln->gr);
	free(free_sln);
};


int main(int argc, char** argv)
{

clock_t begin = clock();


	/* !!! igraph : turn on attribute handling  Ie. even if you don't manipulate attributes explicitly, but create a graph that might have some attributes, eg. read a graph a GraphML file, you need this call before, otherwise the attributes are dropped. */
	igraph_i_set_attribute_table(&igraph_cattribute_table);

	readData(argv);

	//bound_1=1000;
	//bound_2=1000;

//	double radius=(uavs_range/2);
	double radius=uavs_range;

//	sln *res=method1ePasse(grnds, nbr_grnds, radius);
//	sln *trueres=method1ePasse(grnds, nbr_grnds, radius);
	//translate(res, radius);

	sln* res;// true result
	sln* res_plus_1;// true result + 1
	sln* res_plus_2;// true result + 2 and current result
	sln* del;// pointer to result to be freed

	/* Elbow criteria : if 2 next wss are close : deviation(res+1,res+2) > 3/4*(deviation(res,res+1)) */
	int i=1,j=0;

	// to generate the covers
	double prev_deviation=0,wss_minus_1=0,wss=0;/* Caution : wss isn't one of the true result but result+2 */
//	double elbow_ratio=0.7;// if new deviation less than three quarters of previous deviation
	double elbow_ratio=0.7;
// !!! Personal note : if ratio is too low, about 0.001, then there are issues, check if time.!!!!

	bool stop=false;

printf("nbr grnds : %d\n",nbr_grnds);

	do{
		
		/* do until elbow satisfied. Criteria :   */
		res_plus_2=method1ePasse(grnds, nbr_grnds, radius/i);
		wss=k_means(grnds, nbr_grnds, res_plus_2, 0.0001, radius/i);
printf("i %d, wss. %f, n_uavs %d\n", i, wss, res_plus_2->n_uavs);
//		printf("series %d, wss : %f, %d uavs\n", i, wss, res_plus_1->n_uavs);
		if(i<=2)
		{/* the least required steps not yet reached */

			if(i==1)
			{
				i++;
				prev_deviation=wss;
				res=res_plus_2;// just stores result
				continue;// and go to next ite			
			}
			i++;
			prev_deviation=prev_deviation-wss;
			wss_minus_1=wss;
			res_plus_1=res_plus_2;// just stores result
			continue;// and go to next ite			
		}
//printf(", S:%f-%f=%f\n", prev_deviation, wss, prev_wss-wss);
//		buff_next=method1ePasse(grnds, nbr_grnds, radius/(i+1));
//		wss_plus_1=k_means(grnds, nbr_grnds, buff_next, 0.0001, radius/(i+1));
		i++;
		if(wss_minus_1-wss > elbow_ratio*prev_deviation)
		{
printf("Stopping condition holds : wss_minus_1-wss > elbow_ratio*prev_deviation :\n\ti %d, w-1. %f, wss. %f, prevd. %f, ratio. %f, dev %f, f. %d, g. res_plus_1->n_uavs : %d uavs, res->n_uavs : %d uavs\n"
, i, wss_minus_1, wss, prev_deviation, elbow_ratio*prev_deviation, wss_minus_1-wss, wss_minus_1-wss < elbow_ratio*prev_deviation, res_plus_1->n_uavs, res->n_uavs);
			freeSln(res_plus_1);// no longer needed, as only true result "res" is returned
			freeSln(res_plus_2);// same
			stop=true;
		}
		else{/* keep reducing the radius */

			prev_deviation=wss_minus_1-wss;
			wss_minus_1=wss;
			del=res;
			res=res_plus_1;// store result
			res_plus_1=res_plus_2;// store result
			freeSln(del);//Housekeeping
		}
	}while(!stop && wss>0);
//	}while(i<10);


	i--;// out print where it stopped
//	printf("\nFinal series %d, wss : %f, %d uavs\n\n", i, wss, res->n_uavs);
	printf("\nFinal series %d, wss : %f, res uavs %d\n\n", i, wss, res->n_uavs);

/*
	FILE* fp;
	fp=fopen("resTrue.csv","w");
	for(i=1;i<=res->n_uavs;i++)
		for (j=0;j<dim;j++)
		{
			// skip comma not needed after last dim value
			if(j==dim-1)	fprintf(fp,"%lf\n", res->uavs[i][j]);
			else fprintf(fp,"%lf,", res->uavs[i][j]);
		}
	fclose(fp);
*/

	double lb=2.0;

	double* soln=solve_linear_model(res, radius, lb);

/*
clock_t end = clock();
double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("|Time spent %f\n",time_spent);
*/

//	Store the coordinates of the active uavs
//	FILE* fp;
//	fp=fopen("lin_sol1cover.csv","w");
/*
	fp=fopen("lin_solRange250.csv","w");
	for(i=1;i<=res->n_uavs;i++)
	{
		if(soln[i]>0)
			for (j=0;j<dim;j++)
			{
				// skip comma not needed after last dim value
				if(j==dim-1)	fprintf(fp,"%lf\n", res->uavs[i][j]);
				else fprintf(fp,"%lf,", res->uavs[i][j]);
			}
	}
	fclose(fp);
*/
	// keep coordinate of active uavs obtained from the linear model
	int nActivUAVs=0;// used to keep record of the number of uavs
	int indexActivUavs[res->n_uavs+1];// contains the indices of active uavs
	for(i=1;i<=res->n_uavs;i++)
		if(soln[i]>0)
		{/* register indices of active uavs (solverSln[i]>0) */
			nActivUAVs++;
			indexActivUavs[nActivUAVs]=i;
		}

printf("nActivUAVs %d\n", nActivUAVs);
		
	free(soln);

	sln* finalres=(sln*)malloc(sizeof(sln));
	finalres->n_uavs=nActivUAVs;
	finalres->uavs=(double**)malloc((nbr_grnds*2+1)*sizeof(double*));
	finalres->dists=(double**)malloc((nbr_grnds*2+1)*sizeof(double*));
	finalres->counters=(int*)calloc((nbr_grnds*2+1),sizeof(int));
	finalres->covers=(int**)malloc((nbr_grnds*2+1)*sizeof(int*));
	finalres->g_covers=(int*)calloc((nbr_grnds+1),sizeof(int));
	int buffint=0;
	for(i=1;i<=nbr_grnds*2;i++)
	{
		finalres->uavs[i]=(double*)calloc(dim,sizeof(double));// All uavs to origin
		finalres->dists[i]=(double*)calloc((nbr_grnds*2+1),sizeof(double));
		finalres->covers[i]=(int*)malloc((nbr_grnds*2+1)*sizeof(int));
		buffint=indexActivUavs[i];
//printf("buffint %d, i = %d\n",buffint,i);
		for(j=0;j<dim;j++)
		{
			if(i<=nActivUAVs){	finalres->uavs[i][j]=res->uavs[buffint][j];}
			else{finalres->uavs[i][j]=0;}

		}
	}

	FILE* fp;
	fp=fopen("finalres.csv","w");
	for(i=1;i<=finalres->n_uavs;i++)
		for (j=0;j<dim;j++)
		{
			// skip comma not needed after last dim value
			if(j==dim-1)	fprintf(fp,"%lf\n", finalres->uavs[i][j]);
			else fprintf(fp,"%lf,", finalres->uavs[i][j]);
		}
	fclose(fp);

	updateDistMat(finalres, radius);

	freeSln(res);

	// for the graph solutions : 10 individuals for each subpopulation (4) of the ring scheme
	//igraph_t** allGs=(igraph_t**)malloc(10*sizeof(igraph_t*));// useless for now : used when problem modelled as a linear optimisation problem

	//igraph_t* solG0 = translate(res, radius, soln);
	

	// "restr_list" is the restricted list of edges, that depicts the "root" or the list of edges added to the scattered graphs
	// in order to build  the first connected graph. "restr_list" is a matrix where : for i in {1 to the number of ground nodes} :
	// matrix[i][0] == -1, or matrix[i][0] > 0
	// 1. if matrix[i][0] > 0, then matrix[i][0] contains the index of a given uav connected to another existing uav => exists in graph an edge
	// containing value in matrix[i][0] : if exists an edge (u',matrix[i][0]) (or (u',u) where u==matrix[i][0]), then the remaining
	// cells : from matrix[i][1] to matrix[i][m] contain the other m uavs (their indices) to wich u is connected :
	// matrix[i][1,...,m]={u1',u2',...,um'} such that exists (u,u1'),(u,u2'),...,(u,um').
	// to know the value of m for matrix[i], the list matrix[i] is bounded by -1 on the cell following the cell matrix[i][m]
	// <=> the list matrix[i] is such that matrix[i][m+1] == -1
	// eg : for three edges containing uav 5 in the pair (i,j), meaning that there exists (1,5),(3,5),(4,5) in the set of edges
	// then matrix[k] = [5,1,3,4,-1,-1,...,-1]
	// also, for (u',u) with u' in {u1',u2',..., um'} u < u' (to avoid symmetries)
	int** restr_list=(int**)malloc((nbr_grnds+1)*sizeof(int*));
	// initialize the restricted list :
	for (i=0;i<=nbr_grnds;i++)
	{
		restr_list[i]=(int*)malloc((nbr_grnds+1)*sizeof(int));
		for (j=0;j<=nbr_grnds;j++)
			restr_list[i][j]=-1;
	}

	/* keep track of the graph of generated clusters, only concerned about coverage, and less likely to form one connected graph */
	// igraph_t* coverageG=(igraph_t*)malloc(sizeof(igraph_t));
	/* Test if success in creating the graph */
	//if(igraph_empty(coverageG, res->n_uavs+1, IGRAPH_UNDIRECTED)==IGRAPH_EINVAL)
	//	printf(" Invalid number of vertices, graph initialisation failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);

	// keeps track of pairs of used points to create new uavs in order to link separated connected components.
	// usually, at most nbr_grnds-1 edges are needed to build a single connected component
	// but since new positions are constantly created, thus it's not completely possible to assert that only nbr_grnds-1 are needed
	// for the connectivity, thus allocate more space
	int** pairs=(int**)malloc((finalres->n_uavs*3)*sizeof(int*));
	for (i=0;i<(finalres->n_uavs*3);i++)
	{
		pairs[i]=(int*)malloc(2*sizeof(int));
		for (j=0;j<2;j++)	pairs[i][j]=-1;
	}

printf("i=%d and nuavs %d\n",i,finalres->n_uavs);

	igraph_t* solG0=nullptr;
	//igraph_t** solnsGraphs;
	//igraph_t* solnsGraphs;
	//igraph_t** solnsGraphs=(igraph_t**)malloc(NINDIVS*sizeof(igraph_t*));

	int npairs=0;
	//solnsGraphs[0]=build_first_conn_graph(res, radius, restr_list, pairs, &npairs, coverageG);
	//populate(res, radius, solG0, solnsGraphs, restr_list, pairs, &npairs);
	populate(finalres, radius, solG0, restr_list, pairs, &npairs);
// printf("phase 2 : first connected graph, final %d\n", res->n_uavs);


	i=0;
	while(pairs[i][0] > 0 || pairs[i][1] > 0)
	{
		if( (pairs[i][0] > 0 && pairs[i][1] < 0) || (pairs[i][0] > 0 && pairs[i][1] < 0) )
			printf("problem here : not synchronised (%d,%d)\n", pairs[i][0], pairs[i][1]);
//		printf("pair %d : (%f,%f)\n", i, pairs[i][0], pairs[i][1]);
		i++;
	}


/*
	// final connected graph
	fp=fopen("conngraph.csv","w");
	for(i=1;i<=res->n_uavs;i++)
		for (j=0;j<dim;j++)
		{
			// skip comma not needed after last dim value
			if(j==dim-1)	fprintf(fp,"%lf\n", res->uavs[i][j]);
			else fprintf(fp,"%lf,", res->uavs[i][j]);
		}
	fclose(fp);
*/

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("|Time spent %f\n",time_spent);

//	for (i=1;i<=nbr_grnds;i++)
//		printf(" %f ", soln[i]);
//printf("\n");

/*
	FILE* fp;
	fp=fopen("rl.csv","w");
	for(i=1;i<=res->n_uavs;i++)
		for (j=0;j<dim;j++)
		{
			// skip comma not needed after last dim value
			if(j==dim-1)	fprintf(fp,"%lf\n", res->uavs[i][j]);
			else fprintf(fp,"%lf,", res->uavs[i][j]);
		}
	fclose(fp);
//
*/

/*
	fp=fopen("rltrue.csv","w");
	for(i=1;i<=trueres->n_uavs;i++)
		for (j=0;j<dim;j++)
		{
			// skip comma not needed after last dim value
			if(j==dim-1)	fprintf(fp,"%lf\n", trueres->uavs[i][j]);
			else fprintf(fp,"%lf,", trueres->uavs[i][j]);
		}
	fclose(fp);
*/
	

/*
	k_means(grnds, nbr_grnds, res, 0.0001, radius);

	fp=fopen("rlkmeans.csv","w");
	for(i=1;i<=res->n_uavs;i++)
		for (j=0;j<dim;j++)
		{
			// skip comma missing for last dim
			if(j==dim-1)	fprintf(fp,"%lf\n", res->uavs[i][j]);
			else fprintf(fp,"%lf,", res->uavs[i][j]);
		}
	fclose(fp);
*/
//	translate(res, radius, soln);
//	igraph_t graph_sol=translate(res);


//	printf("n uvas : %d and grnds[17][0] : %lf\n", res->n_uavs,grnds[17][0]);
}
