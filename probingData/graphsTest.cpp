

//#include "linear_solver.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <string>
#include <math.h>

#include <igraph.h>
#include <assert.h>
#include <float.h>

#include <glpk.h>


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
	igraph_t gr;// graph of the solution. Needed for connectivity checking
	int n_uavs;// number of nodes in network : either nodes are uavs (sln : set of uavs), or ground nodes (uav : set of ground nodes)
}sln;


typedef struct block{
	float lowleft[2];// the coordinates of the lower left bound
	float upright[2];// coords of the upper right bound
	char color;// if leaf : 'b' for black (if block contains only 	), 'w' for white (area is empty), 'g' otherwise (node, not leaf, grey node)
	int label;// tells which connected component the block belongs to
	struct block father;
	struct block quadrants[4];// 4 regions of a quadtree's node : NW, NE, SW, SE
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
double euclDistance(double *node1, double *node2);
void duplicate(double lb, sln* net);
void translate(sln* net, double threshold);
void updateDistMat(sln* net, double range);

double* solve_linear_model(sln* net, double range, double lb);

void freeSln(sln* free_sln);

void readData(char** argv)
{
	FILE* fp;
	char *tmp[100];// buffer
	int tmp2;// this too


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
			ub=net->g_covers[i];
			do
			{
				net->n_uavs++;// increment the number of uavs
				net->counters[net->n_uavs]=1;// Increment 
				net->covers[i][net->g_covers[i]]=net->n_uavs;// Include new uav as covering ground g
				net->g_covers[i]++;// Increment number of covers for k-th ground node
				/* Duplicate uav coords on new one*/
				ind=net->covers[i][net->g_covers[i]%ub];
				for(j=0;j<dim;j++)
					net->uavs[net->n_uavs][j]=net->uavs[ind][j];
			/* Keep copying until constraint on lb is satisfied */
			}while(net->g_covers[i] < lb);
		}
	}
};

double* solve_linear_model(sln* net, double range, double lb)
{

	int i=0,j=0;

	/* Find all covers for each groun node */
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

	/* Find at least one uav that needs to be replicated so that a given ground node which violates the constraint of at least lb covers */
	for(i=1;i<=nbr_grnds;i++)
	{
		/* Create duplicate if required */
		if( net->g_covers[i] < lb )
		{
			duplicate(lb, net);
			break;
		}
	}

	/* Build problem */
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
for(j=1;j<=net->n_uavs;j++)
{
	ind[j]=0;
	val[j]=0;
}
	int len=glp_get_mat_row(lp, 113, ind, val);
//printf(" LO %f, COL LO %f, COL UP %f\n",glp_get_row_lb(lp,113),glp_get_col_lb(lp,88),glp_get_col_ub(lp,88));
*/
	
	int result_solver=glp_simplex(lp, NULL);
	int result_solver2=glp_intopt(lp, NULL);
	z = glp_get_obj_val(lp);
	
	/* Gather results */
	double* soln=(double*)calloc(net->n_uavs+1,sizeof(double));
	for(j=1;j<=net->n_uavs;j++)
		soln[j]= glp_get_col_prim(lp, j);
	
	glp_delete_prob(lp);

	for(j=1;j<=net->n_uavs;j++)
		if(soln[j]>0)
			printf("[ %d : %f ] ",j, soln[j]);

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
	for(i=2;i<=nbr_grnds;i++)
	{
		printf(" [%d : %d] ",i,activecovers[i]);
		sum+=activecovers[i];
		if(activecovers[i] > max)	max=i;
	}
	printf("\nAverage : %f\n",sum/nbr_grnds);
	printf(" Max degree : %d with %d\n",max,activecovers[max]);
	return soln;
};


sln* method1ePasse(double** input_data, int size_input, double threshold)
{
	int k=2,j=1;
	sln* res=(sln*)malloc(sizeof(sln));
	res->g_covers=(int*)calloc((size_input+1),sizeof(int));
	res->covers=(int**)malloc((size_input+1)*sizeof(int*));
	res->counters=(int*)calloc((size_input+1),sizeof(int));
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


void translate(sln* net, double threshold)
{// note : computing diameter of graph (O(n*|E|)) cost more than comparing each nodes (O(n^2)) : Complete graph : n*|E| = n*sum_i_in{1...n-1}(i) > O(n^2) for n>3

	long int i=0,j=0;
	igraph_integer_t ncomps=-1;// numbers of connected components
	igraph_vector_t labels, compssizes;// labels of connected components for each vertex
	double min_distance=DBL_MAX;// distance of two closest nodes but out of range
	double current_distance;// used to place a new uav, "shared" is the degree of "merging" of two uavs where the area is shared
	double shared=threshold/4;
	int n1=-1,n2=-1;//two closest but out of range nodes

	// here vector "labels" only used to create graph
	if(igraph_empty(&net->gr, net->n_uavs+1, IGRAPH_UNDIRECTED)==IGRAPH_EINVAL)
		printf(" Invalid number of vertices, graph init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);

	// add edges to graph, edges : two uavs in range
	for(i=1;i<=net->n_uavs;i++)
		for(j=i+1;j<=net->n_uavs;j++)
			if(net->dists[i][j]<2*threshold-shared)
{
printf("(%ld,%ld)\t",i,j);
				igraph_add_edge(&net->gr, i, j);
}
printf("end init graph with %d uavs\n",net->n_uavs);

//int k=0;
	do
	{
		
		if(igraph_vector_init(&labels, 0)==IGRAPH_ENOMEM || igraph_vector_init(&compssizes, 0)==IGRAPH_ENOMEM)
			printf(" Memory issues, vector init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);
			
		// compute number of connected components
		igraph_clusters(&net->gr, &labels, &compssizes, &ncomps, IGRAPH_WEAK);
printf("start labels\n");
for(i=0;i<igraph_vcount(&net->gr);i++)	printf("l[%li]:%li\t",i,(long int)VECTOR(labels)[i]);
printf("end labels\n");
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
					/* Find closest clusters but over range */
					current_distance=net->dists[i][j];
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
			if(net->dists[net->n_uavs][n1]<2*threshold-shared)
{
printf("(%d,%d)\n",net->n_uavs, n1);
				igraph_add_edge(&net->gr, net->n_uavs, n1);
}
			if(net->dists[net->n_uavs][n2]<2*threshold-shared)
{
printf("(%d,%d)\n",net->n_uavs, n2);
				igraph_add_edge(&net->gr, net->n_uavs, n2);
}
printf("1. dist n_uav-n1 : %lf, thresh : %lf, nuav : %d, n1 : %d\n", net->dists[net->n_uavs][n1], 2*threshold-shared, net->n_uavs, n1);
printf("2. dist n_uav-n2 : %lf, thresh : %lf, nuav : %d, n2 : %d\n", net->dists[net->n_uavs][n2], 2*threshold-shared, net->n_uavs, n2);
printf(" %d and %d coords (%lf,%lf) (%lf,%lf)\n",n1,n2, net->uavs[n1][0], net->uavs[n1][1], net->uavs[n2][0], net->uavs[n2][1]);
printf("dist n1-n2 after update %d - %d\n",(int)net->dists[n2][n1],(int)net->dists[n1][n2]);

//printf(" %d and %d coords (%lf,%lf) (%lf,%lf)\n",n1,n2, net->uavs[n1][0], net->uavs[n1][1], net->uavs[n2][0], net->uavs[n2][1]);
//printf("with %li vertices, and new uav %d coords (%lf,%lf)\n", (long int)igraph_vcount(&net->gr), net->n_uavs, net->uavs[net->n_uavs][0], net->uavs[net->n_uavs][1]);
//printf(" vcount %li ecount %li\n", (long int)igraph_vcount(&net->gr), (long int)igraph_ecount(&net->gr));
		}
//	}while(ncomps>2 && k++ < 5);
	}while(ncomps>2);


	// write graph into file
	FILE* fp;
	fp=fopen("graphs.csv","w");
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
};

void freeSln(sln* free_sln)
{
	igraph_destroy(&free_sln->gr);
	free(free_sln);
};


int main(int argc, char** argv)
{

	/* !!! igraph : turn on attribute handling  Ie. even if you don't manipulate attributes explicitly, but create a graph that might have some attributes, eg. read a graph a GraphML file, you need this call before, otherwise the attributes are dropped. */
	igraph_i_set_attribute_table(&igraph_cattribute_table);

	readData(argv);

	bound_1=1000;
	bound_2=1000;

//	double radius=(uavs_range/2);
	double radius=(uavs_range/2);

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
	double elbow_ratio=0.7;// if new deviation less than a quarter of previous
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
//printf("i %d, w-1. %f, wss. %f, prevd. %f, ratio. %f, dev %f, f. %d, g. %d uavs\n", i, wss_minus_1, wss, prev_deviation, elbow_ratio*prev_deviation, wss_minus_1-wss, wss_minus_1-wss < elbow_ratio*prev_deviation, res_plus_1->n_uavs);
		if(wss_minus_1-wss > elbow_ratio*prev_deviation)
		{
			free(res_plus_1);// no longer needed, as only true result "res" is returned
			free(res_plus_2);// same
			stop=true;
		}
		else{/* keep reducing the radius */
			prev_deviation=wss_minus_1-wss;
			wss_minus_1=wss;
			del=res;
			res=res_plus_1;// store result
			res_plus_1=res_plus_2;// store result
			free(del);//Housekeeping
		}
	}while(!stop);
//	}while(i<10);

	i--;// out print where it stopped
	printf("\nFinal series %d, wss : %f, %d uavs\n\n", i, wss, res->n_uavs);

	double* soln=solve_linear_model(res, radius, 2.0);

/*
	FILE* fp;
	fp=fopen("lin_sol1cover.csv","w");
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
//	for (i=1;i<=nbr_grnds;i++)
//		printf(" %f ", soln[i]);
//printf("\n");

/*
///*
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
//	translate(res, radius);
//	igraph_t graph_sol=translate(res);


//	printf("n uvas : %d and grnds[17][0] : %lf\n", res->n_uavs,grnds[17][0]);
}
