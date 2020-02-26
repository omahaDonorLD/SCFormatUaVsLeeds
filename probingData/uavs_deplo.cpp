

//#include "linear_solver.hpp"
//#include <stdlib>
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <cmath>

#include <igraph.h>
#include <assert.h>
#include <float.h>
#include <glpk.h>

#include <vector>

using namespace std;

inline
int STREAM_FAIL(const char *FROM_FILE, int AT_LINE, const char *IN_FUNCTION)
{printf("STREAM_FAILURE, line %d, function %s, file %s\n", AT_LINE, IN_FUNCTION, FROM_FILE);return EXIT_FAILURE;};

inline
int MEMO_FAIL(const char *FROM_FILE, int AT_LINE, const char *IN_FUNCTION)
{printf("MEMO_ALLOC_FAILURE, line %d, function %s, file %s\n", AT_LINE, IN_FUNCTION, FROM_FILE);return EXIT_FAILURE;};

typedef struct sln{
	vector<double*> uavs;// vector of pointers to arrays (coordinates of uavs)
	vector<int> *gcovs;// for each ground node, the indices of uavs covering it : an array size 'nbr_grnds' of cells of vectors of indices of uavs
	vector<vector<int>> uavcovs;// for each uav, the ground node it covers
	// adjacency list of distances between uavs. Two structs:
	vector<vector<int>> outdeg;// 1. outdegree relation of uav, note: j < outdeg[j] 
	vector<vector<double>> distuav;// dists[j][outdeg[j]] = distance between uav j and uav outdeg[j]
	//igraph_t* gr;// graph of the solution. Needed for connectivity checking
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
sln* method1ePasse(double** input_data, int nbr_grnds, double threshold);
double k_means(double** data, int n, sln* clusts, double error_tolerance, double range);
double euclDistance(double *node1, double *node2);
void duplicate(double lb, sln* net, int grndi, double range);
igraph_t* translate(sln* net, double threshold, double* solverSln);
void updateDistMat(sln* net);
void addTonetwork(sln* net, double range);

bool uav_in_cover(vector<int> &gcovs, int uavindex);
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
	grnds=(double**)malloc(nbr_grnds*sizeof(double*));
	if(grnds==NULL){ /* memory allocation failure */ MEMO_FAIL(__FILE__, __LINE__, __FUNCTION__); }

	for(i=0;i<nbr_grnds;i++)
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

void updateDistMat(sln* net)
{// updates the matrix of distances to avoid constantly having to compute them
	int i,j;
	double dist=0;
	for(i=0;i<net->uavs.size();i++)
	{
		for(j=0;j<net->outdeg[i].size();j++)
		{
			int adjuav=net->outdeg[i][j];
			dist=euclDistance(net->uavs[i],net->uavs[adjuav]);
			net->distuav[i][adjuav]=dist;
		}
	}
};


void addTonetwork(sln* net, double range)
{// updates the matrix of distances to avoid constantly having to compute them
	int i, newuav=net->uavs.size()-1;
	double dist=0;
	for(i=0;i<newuav;i++)
	{
		dist=euclDistance(net->uavs[i],net->uavs[newuav]);
		if(dist<=range){
			net->outdeg[i].push_back(newuav);
			net->distuav[i].push_back(dist);
		}
	}
};


void duplicate(double lb, sln* net, int grndi, double range)
{/* lb for knowing if ground needs more covers, counts[] the list which gives for each ground node the number of uavs covering it,
	covers[] the list of aforementioned covers */

	int i=0,j=0, n=0, newuavs=0, uavj=0, randindex=0;
	double *buff;

	do
	{
		// duplicate one of the uavs covering ground node i
		randindex=rand()%net->gcovs[grndi].size();
		uavj=net->gcovs[grndi][randindex];

		net->gcovs[grndi].push_back(net->uavs.size());// Include new uav as covering ground g
		buff=new double[dim];
		for(j=0;j<dim;j++)
			buff[j]=net->uavs[uavj][j];
		net->uavs.push_back(buff);
		net->uavcovs[uavj].push_back(grndi);// new uav covers ground node i

		net->outdeg.resize(net->outdeg.size()+1);
		net->distuav.resize(net->distuav.size()+1);
		addTonetwork(net, range);// update distance of network of uavs

	/* Keep duplicating until constraint on lb is satisfied */
	}while(net->gcovs[grndi].size() < lb);

/*
printf("Duplicated : ");
for(int a:net->gcovs[grndi])
printf(" %d ", a);
printf("\n");
*/
};


bool uav_in_cover(vector<int> &gcovs, int uavindex){
	for(int i=0; i<gcovs.size(); i++){
		if(gcovs[i]==uavindex)
			return true;
	}
	return false;
};



void find_covers(sln* net, double range)
{/* Find every covers for each ground node */

	int i=0,j=0;
		
	for(i=0;i<nbr_grnds;i++)
	{
		for(j=0;j<net->uavs.size();j++)
		{
			if( euclDistance(net->uavs[j],grnds[i]) <= range && !uav_in_cover(net->gcovs[i], j))
			{
				net->gcovs[i].push_back(j);
				net->uavcovs[j].push_back(i);
			}
		}
	}
};


double* solve_linear_model(sln* net, double range, double lb)
{

	int i=0,j=0;

	/* find covers for each ground node */
	find_covers(net, range);

	bool changedstate=false;
	/* Find at least one uav that needs to be replicated so that a minimum cover constraint isn't violated : at least lb covers for each target */
	for(i=0;i<nbr_grnds;i++)
	{
//		/* Create duplicate only if needed => if at least one violation of constraint is found then duplicate for all occuring */
		if( net->gcovs[i].size() < lb )
		{
printf("Duplicate uavs for gnode %d, which is covered with %d uavs \n", i, net->gcovs[i].size());
			duplicate(lb, net, i, range);
			if(!changedstate)	changedstate=true;
//			break;
		}
	}
	
	// if uavs added in duplication phase, then search again covers with ground nodes
	if(changedstate)
	{
		find_covers(net, range);
	}

printf("Begin building linear problem with %d ground nodes, and %d uavs\n", nbr_grnds, net->uavs.size());
	/* Build linear problem */
	glp_prob *lp;
	int ia[1+(nbr_grnds*net->uavs.size())],ja[1+(nbr_grnds*net->uavs.size())];
	double ar[1+(nbr_grnds*net->uavs.size())], z=0.;
	
	lp = glp_create_prob();
	glp_set_prob_name(lp, "set cover");
	glp_set_obj_dir(lp, GLP_MIN);

	/* for row names */
	char row_names[nbr_grnds+1][20];
	char col_names[net->uavs.size()+1][20];
	char buff[20];

	glp_add_rows(lp, nbr_grnds);
	for(i=0;i<nbr_grnds;i++)
	{
		sprintf(buff,"r%d",i+1);
		strcpy(row_names[i+1],buff);
		glp_set_row_name(lp, i+1, row_names[i+1]);
		glp_set_row_bnds(lp, i+1, GLP_LO, lb, 0.0);
	}

	glp_add_cols(lp, net->uavs.size());
	for(j=0;j<net->uavs.size();j++)
	{
		sprintf(buff,"x%d",j+1);
		strcpy(col_names[j+1],buff);
		glp_set_col_name(lp, j+1, col_names[j+1]);
		glp_set_col_bnds(lp, j+1, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(lp, j+1, 1.0);
		glp_set_col_kind(lp, j+1, GLP_IV);
		//glp_set_col_kind(lp, j, GLP_BV);
		//glp_set_col_kind(lp, j, GLP_CV);
	}

printf("RANGE %f, lb %f\n",range,lb);
	int counter=1;
	for(i=0;i<nbr_grnds;i++)
	{
		for(j=0;j<net->uavs.size();j++)
		{
			ia[counter] = i+1;
			ja[counter] = j+1;
			ar[counter] = ( euclDistance(net->uavs[j],grnds[i]) <= range ? 1.0 : 0.0 );
			counter++;
		}
	}

	glp_load_matrix(lp, nbr_grnds*net->uavs.size(), ia, ja, ar);

/*
	int ind[net->uavs.size()+1];
	double val[net->uavs.size()+1];
	for(i=0;i<nbr_grnds;i++)
	{
		for(j=0;j<net->uavs.size();j++)
		{
			ind[j+1]=0;
			val[j+1]=0;
		}
		int len=glp_get_mat_row(lp, i+1, ind, val);// stores non zero values of row 113 into "ind(indices of cols), val(non zero values)"
//		printf(" LO %f, COL LO %f, COL UP %f %d\n",glp_get_row_lb(lp,113),glp_get_col_lb(lp,88),glp_get_col_ub(lp,88), len);
		printf("gnode %d has %d, covers => ", i+1, len);
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
	double* soln=(double*)calloc(net->uavs.size(),sizeof(double));
	double sumrows=0, maxrowval=0, minrowval=net->uavs.size(), tmp=0, lessThanAvrg=0, moreThanvrg=0;// Needed for stats on covers (constraints, rows)
	double* covRes=(double*)calloc(nbr_grnds+1,sizeof(double));

	for(i=0;i<nbr_grnds;i++)
	{// for statistics on the values returned by the rows
//		tmp=glp_get_row_prim(lp, i);
		tmp=glp_mip_row_val(lp, i+1);
//		printf("cstrnt %d = %f\n", i, tmp);
		sumrows += tmp;
		if(tmp > maxrowval)	maxrowval=tmp;
		if(tmp < minrowval)	minrowval=tmp;
		if(tmp >= 3){moreThanvrg++;}else{lessThanAvrg++;}
		covRes[(int)tmp]++;
	}

	for(j=0;j<net->uavs.size();j++)
	{
//		soln[j]= glp_get_col_prim(lp, j);
		soln[j]= glp_mip_col_val(lp, j+1);
//		if (soln[j] - ceil(soln[j]) != 0)	printf(" soln %d not an integer = %f\n", j, soln[j]);
	}
	
	glp_print_sol(lp, "glpksols.txt");

	glp_delete_prob(lp);


	for(j=0;j<net->uavs.size();j++)
		if(soln[j]>0)
			printf("[ %d : %f ] ",j, soln[j]);
	
	printf("\nRO == %f, max = %f, min = %f, aver = %f, lessThanAvrg = %f, moreThanvrg = %f\n",sumrows, maxrowval, minrowval, sumrows/nbr_grnds, lessThanAvrg, moreThanvrg);
	for(i=0;i<=(int)maxrowval;i++)	if (covRes[i]>0) printf("cov[%d]=%f\t",i,covRes[i]);
	printf("\n");
	int activeuavs=0;

	FILE* fp;
	fp=fopen("activeduavs.csv","w");
	for(i=0;i<net->uavs.size();i++)
	{
		if(soln[i]>0)
		{
			for (j=0;j<dim;j++)
			{
				// skip comma not needed after last dim value
				fprintf(fp,"%lf", net->uavs[i][j]);
				if(j==dim-1)	fprintf(fp,"\n");
				else	fprintf(fp,",");
			}
			activeuavs++;
		}
	}
	fclose(fp);

	// print in file connections uavs-ground nodes
	fp=fopen("uavs_grounds.csv","w");
	int *activecovers=(int*)calloc(nbr_grnds,sizeof(int));
	for(i=0;i<nbr_grnds;i++)
	{
		for(j=0;j<net->gcovs[i].size();j++)
		{
			int indice=net->gcovs[i][j];// find index of uav
			if(soln[indice]<=0.)	continue;// skip if not active
			activecovers[i]++;
			fprintf(fp,"%lf,%lf\n", grnds[i][0], grnds[i][1]);
			fprintf(fp,"%lf,%lf\n\n", net->uavs[indice][0], net->uavs[indice][1]);
		}
	}

	printf(" Obj func : %f and nactive uavs %d\n",z,activeuavs);
	int max=0;
	double sum=0;
	printf("Couples (ground, n active uavs covering it) : ");
	for(i=0;i<nbr_grnds;i++)
	{
		printf(" (%d, %d) ",i,activecovers[i]);
		sum+=activecovers[i];
		if(activecovers[i] > activecovers[max])	max=i;
	}
	printf("\nAverage : %f\n",sum/nbr_grnds);
	printf(" Max degree : %d with %d\n",max,activecovers[max]);
	
	// house keeping
	delete covRes;
	delete activecovers;

	return soln;
};


sln* blankSolution(int size_input, int nuavs)
{

	int k=0,j=0;
	sln* res=new sln;
/*
	res->g_covers=(int*)calloc((size_input+1),sizeof(int));
	res->covers=(int**)malloc((size_input+1)*sizeof(int*));
	res->counters=(int*)calloc((size_input+1),sizeof(int));
	res->uavs=(double**)malloc((size_input*3+1)*sizeof(double*));
	res->n_uavs=nuavs;
	res->dists=(double**)malloc((size_input+1)*sizeof(double*));
	/* Initialisation steps */
/*
	for(k=0;k<=size_input;k++)
	{
		res->covers[k]=(int*)calloc((size_input+1),sizeof(int));
		res->uavs[k]=(double*)calloc(dim,sizeof(double));// All uavs to origin
		for(j=0;j<dim;j++)	res->dists[k]=(double*)calloc((size_input+1),sizeof(double));
	}
*/
	return res;
}


sln* method1ePasse(double** input_data, int nbr_grnds, double range)
{
	int k=0,j=1;
	// required memory allocation for solution
	sln* res=new sln;
	res->gcovs=new vector<int>[nbr_grnds];
	double *buffdouble;
	//res->gr=(igraph_t*)malloc(sizeof(igraph_t));
	vector<double*> cl;/* temporary variables, will contain sum of elements in cluster, to build clusters */

	/* Init first cluster with first element */
	buffdouble=new double[dim];
	for(k=0;k<dim;k++)// here indices for 'dim' are from 0 to n, sorry...
		buffdouble[k]=input_data[0][k];
	res->uavs.push_back(buffdouble);
	res->outdeg.resize(res->outdeg.size()+1);
	res->distuav.resize(res->distuav.size()+1);
	buffdouble=new double[dim];
	for(k=0;k<dim;k++)// same with temporary value
		buffdouble[k]=input_data[0][k];
	cl.push_back(buffdouble);
	
//	for(k=0;k<nbr_grnds;k++)
//		res->gcovs[k].push_back(-1);

	res->uavcovs.push_back(vector<int>(1,0));// Init first uav to cover first ground node
	res->gcovs[0].push_back(0);// label of first ground node : first cluster

	double distance_to_closest_centroid=bound_1*bound_2;// Largest distance : the extrem limits of the entire map
	double current_distance=0;

	/* start 1e pass algorithm */
	for(k=1; k<nbr_grnds; k++)
	{
		distance_to_closest_centroid=bound_1*bound_2;// Worst : map's limits
		for(j=0;j<res->uavs.size();j++)
		{
			current_distance=euclDistance(input_data[k],res->uavs[j]);
			if ( current_distance < distance_to_closest_centroid )
			{
				distance_to_closest_centroid=current_distance;// find the closest cluster
				if (res->gcovs[k].size()==0) res->gcovs[k].push_back(j);
				else res->gcovs[k][0]=j;
			}
		}

		if(distance_to_closest_centroid<range)
		{// closest centroid is within admissible range => add element to cluster
			res->uavcovs[res->gcovs[k][0]].push_back(k);
			/* update centroids */
			for(j=0;j<dim;j++)
			{
				cl[res->gcovs[k][0]][j]+=input_data[k][j];
				res->uavs[res->gcovs[k][0]][j]=cl[res->gcovs[k][0]][j]/res->uavcovs[res->gcovs[k][0]].size();
			}
		}
		else
		{// closest centroid is not within admissible range => Create own new cluster
			if (res->gcovs[k].size()==0) res->gcovs[k].push_back(res->uavcovs.size());
			else res->gcovs[k][0]=res->uavcovs.size();
			res->uavcovs.push_back(vector<int>(1,k));
			buffdouble=new double[dim];
			for(j=0;j<dim;j++)
				buffdouble[j]=input_data[k][j];
			cl.push_back(buffdouble);
			buffdouble=new double[dim];
			for(j=0;j<dim;j++)
				buffdouble[j]=input_data[k][j];
			res->uavs.push_back(buffdouble);
			res->outdeg.resize(res->outdeg.size()+1);
			res->distuav.resize(res->distuav.size()+1);
			addTonetwork(res, range);// update distance of network of uavs
		}
	}
	
//	updateDistMat(res, threshold); No longer needed, done in "addTonetwork"

	/* Housekeeping */
	for(k=0;k<res->uavs.size();k++)
	{
		delete[] cl[k];
		cl[k]=nullptr;
	}
//	free(cl);

	return res;
};


double k_means(double** data, int ngrounds, sln* clusts, double error_tolerance, double range)
{

    int h, i, j; /* loop counters, of course :) */

    double old_error, error = DBL_MAX; /* sum of squared euclidean distance. DBL_MAX : macro defines the maximum finite floating-point value : 1E+37*/
    double** cl=(double**) malloc((clusts->uavs.size())*sizeof(double*)); /* temp centroids */

//    assert(data && clusts->uavs.size() > 0 && clusts->uavs.size() <= n && error_tolerance >= 0); /* for debugging */
    assert(data && clusts->uavs.size() > 0 && error_tolerance >= 0); /* for debugging */

    /* initialization : start with given centroids */

     // c1 : temp centroids, clusts->uavs : true centroids
	for (i=0; i<clusts->uavs.size(); i++)
		cl[i]=(double*)calloc(dim, sizeof(double));

	double min_distance = DBL_MAX, distance = 0;

    /** Kmeans loop */
    do {
        /* save error from last step */
		old_error = error;
        error = 0;

        /* clear old counts and temp centroids */
        for (i=0; i <clusts->uavs.size(); i++)
        {
			// reinit temp centroids
            for (j=0;j<dim; cl[i][j++] = 0);
			clusts->uavcovs[i].erase(clusts->uavcovs[i].begin(), clusts->uavcovs[i].end());
        }

        for (h = 0; h < ngrounds; h++)
        {
            /* Find closest cluster */
            min_distance = DBL_MAX;
            for (i=0; i<clusts->uavs.size(); i++)
            {
				distance=0;
				/* considered distance : euclidian squared, or squared residuals */
                for(j=0;j<dim;j++)	distance+=pow(data[h][j]-clusts->uavs[i][j],2);
                if (distance < min_distance)
                {
					min_distance=distance;// find the closest cluster
                    clusts->gcovs[h][0] = i;
                }
            }
			/* update size and temp centroid of the destination cluster */
			for(j=0;j<dim;j++)	cl[clusts->gcovs[h][0]][j] += data[h][j];
			clusts->uavcovs[clusts->gcovs[h][0]].push_back(h);
            /* update standard error */
            error += min_distance;
        }

		/* update all centroids */
        for (i=0; i<clusts->uavs.size(); i++)
            for (j=0; j<dim; j++)
                clusts->uavs[i][j] = clusts->uavcovs[i].size() ? cl[i][j]/clusts->uavcovs[i].size() : cl[i][j];

    } while (fabs(error - old_error) > error_tolerance);/* if for each iteration, the number of changes made are not different from previous */

	updateDistMat(clusts);


	/* housekeeping */
	for(h=0;h<clusts->uavs.size();h++)
	{
		delete[] cl[h];
		cl[h]=nullptr;
	}
	delete cl;

	return error;
};


/*
void connect_CCs(sln* net, igraph_t* Gk, double range, vector<int*> &restr_list, int **pairs, int* npairs, bool restrict, int level)
{// The boolean value "restrict" says wether the restriction strategy should be used (then restrict is false, and G isn't root G0), or not.
	// In the latter case, fill lists : **restr_list, **pairs, and * npairs. Otherwise, use these lists to apply the restrictions

printf("net n_uavs %d threshold %f and restrict %d\n", net->uavs.size(), threshold, restrict);
//printf(" vcount %li ecount %li\n", (long int)igraph_vcount(Gk), (long int)igraph_ecount(Gk));

	igraph_integer_t ncomps=-1;// numbers of connected components
	igraph_vector_t labels, compssizes;// labels of connected components for each vertex
	double min_distance=DBL_MAX;// distance of two closest nodes but out of range (belong to two different connected components)
	double current_distance;// used to place a new uav
	int buffn1=-1,buffn2=-1, n1=-1, n2=-1;//two closest but out of range nodes

	// Use of a brute force branching here. (Maybe undeterministic better?)
	// brute force : create one unique connected component by iteratively branching the 2 closest connected components
	// From here, a bit of a mess, check for segmentation faults
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

/*
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
			for (i=0;i<net->uavs.size();i++)
			{
				for (j=0;j<net->outdeg[i].size();j++)
				{
					int uav2=net->outdeg[i][j];
					// Not interested in nodes in same connected component
					if ( VECTOR(labels)[i] == VECTOR(labels)[uav2] )	continue;
					/* Find closest clusters but out of range */
/*
					current_distance=net->distuav[i][uav2];
//					if(current_distance<threshold)// skip, in range
//						continue;
					// check only if different connected component
if(current_distance>=min_distance)
printf("\n\n++++++++++\nProblems: should not verify since skipped uavs in same connected component\n++++++++++\n\n");
					if(current_distance<min_distance)
					{// keep two nodes
						buffn1=i;
						buffn2=uav2;
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
						n2=uav2;
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

			// free vector of labels
			igraph_vector_destroy(&labels);

			if(n1<0 || n2<0)
			{
				printf(" Something went wrong with graph (Either connectivity or indices. Check), adding new uav will fail %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);
				return;
			}

			// assign coordinates of new created node : middle of segment [n1,n2]
			double* buffdouble=new double[dim];
			for (j=0;j<dim;j++)
				buffdouble[j]=(net->uavs[n1][j]+net->uavs[n2][j])/2;
			res->uavs.push_back(buffdouble);
			res->uavcovs.push_back(vector<int>(1,-1));// Init first uav to cover first ground node
			addTonetwork(res, range);

			// increase number of elements in new cluster : at least itself, since may not cover any ground node, but be used only for connectivity
			// extend size of graph
			igraph_add_vertices(Gk, 1, 0);			
			// update new distances
//			updateDistMat(net, threshold);

			// fill in various lists
			if( ! restrict )
			{
				net->counters[net->n_uavs]++;
				restr_list[added][0]=net->uavs.size()-1;
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

/*
};
*/


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
		
		// do until elbow criteria satisfied:
		res_plus_2=method1ePasse(grnds, nbr_grnds, radius/i);
		wss=k_means(grnds, nbr_grnds, res_plus_2, 0.0001, radius/i);
printf("i %d, wss. %f, n_uavs %d\n", i, wss, res_plus_2->uavs.size());
		i++;
		if(i<=3)
		{
			if(i==2)
			{
				prev_deviation=wss;
				res=res_plus_2;// just stores result
				continue;// and go to next ite			
			}
			prev_deviation=prev_deviation-wss;
			wss_minus_1=wss;
			res_plus_1=res_plus_2;// just stores result
			continue;// and go to next ite			
		}
//printf(", S:%f-%f=%f\n", prev_deviation, wss, prev_wss-wss);

		if(wss_minus_1-wss > elbow_ratio*prev_deviation)
		{
printf("Stopping condition holds : wss_minus_1-wss > elbow_ratio*prev_deviation :\n\ti %d, w-1. %f, wss. %f, prevd. %f, ratio. %f, dev %f, f. %d, g. res_plus_1->n_uavs : %d uavs, res->n_uavs : %d uavs\n"
, i, wss_minus_1, wss, prev_deviation, elbow_ratio*prev_deviation, wss_minus_1-wss, wss_minus_1-wss < elbow_ratio*prev_deviation, res_plus_1->uavs.size(), res->uavs.size());
			delete res_plus_1;// no longer needed, as only true result "res" is returned
			res_plus_1=nullptr;
			delete res_plus_2;// same
			res_plus_2=nullptr;
			stop=true;
		}
		else{/* keep reducing the radius */

			prev_deviation=wss_minus_1-wss;
			wss_minus_1=wss;
			del=res;
			res=res_plus_1;// store result
			res_plus_1=res_plus_2;// store result
			delete del;//Housekeeping
			del=nullptr;
		}
	}while(!stop && wss>0);
//	}while(i<10);


	i--;// out print where it stopped
	printf("\nFinal series %d, wss : %f, res uavs %d\n\n", i, wss, res->uavs.size());

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
	return 0;
};