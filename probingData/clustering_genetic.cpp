

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
	int* counters;// size : n_uavs. For each uav, gives the number of targets covered
	int* g_covers;// Size : number of ground nodes. For each ground node gives the number of uavs in its vicinity
	int** covers;// Size : number of ground nodes*counters, gives for each ground all the uavs that cover it
	double** uavs;/* matching "contains" list with "uavs_indices". The latter contains the indices 'i'
							of elements in list "uavs" so that one can finds the coordinates of uavs[i] */

	// double** distances;// square matrix of distances between uavs
	double** dists;// square matrix of distances between uavs
	igraph_t gr;// graph of the solution. When allocated, assumed to have one unique connected components
	int n_uavs;// number of nodes (uavs) in the network
	//int n_uavs;// number of nodes (uavs) in the network : either nodes are uavs (sln : set of uavs), or ground nodes (uav : set of ground nodes)
}sln;


typedef struct block{/* needed for the quadtrees */
	float lowleft[2];// the coordinates of the lower left bound
	float upright[2];// coords of the upper right bound
	char color;/* depicts nature of block :
					- 'b' for black, leaf, block contains only uavs,
					- 'w' for white, leaf, block is empty : no uav
					- 'g' for grey, not leaf, contains other nodes */
	int label;// if block is a leaf, says which connected component the block belongs to
	typedef struct block father;
	typedef struct block quadrants[4];// 4 regions of a quadtree's node : NW, NE, SW, SE
}block;


int max_uavs_avail;
int nbr_grnds;
int dim=2;			// dimension of input vectors
double** grnds;		// All ground nodes coordinates : ngrouns*(1,...,dim)
double uavs_range;	// all ranges of available uavs : can be same for all uavs, or differents

int NINDIVS=60;		// Needed for the genetic version

// limit of map
double bound_1;
double bound_2;

void readData(char** argv);// read data : number of available uavs, range of uavs, number of targets, limit of map, coordinates of ground nodes
sln* method1ePasse(double** input_data, int size_input, double threshold);
double k_means(double** data, int n, sln* clusts, double error_tolerance, double range);
double euclDistance(double *node1, double *node2);
void duplicate(double lb, sln* net);// duplicate uav if needed : number of uavs covering a target not enough
//void translate(sln* net, double threshold, double* solverSln);
/* build the first single connected graph (returned), keep track of the edges used to build the connectivity (restr_list), and the graph
 * of clusters for the coverage requirement (without the additional edges and nodes required for the connectivity */
//igraph_t* build_first_conn_graph(sln* net, double threshold, int** restr_list, double **pairs, int* npairs, igraph_t* coverageG);
void connect_CCs(sln* net, igraph_t* G, double threshold, int** restr_list, int **pairs, int* npairs, bool restrict, int level);
//void populate(sln* covers, double threshold, igraph_t* solnG0, igraph_t** solnsGraphs, int** restr_list, double **pairs, int* npairs);
void populate(sln* covers, double threshold, igraph_t* solnG0, int** restr_list, int **pairs, int* npairs);
void updateDistMat(sln* net, double range);// updates the matrix of distances to avoid constantly having to compute them
void find_covers(sln* net, double range);// Find every covers for each ground node

int cmpfuncXaxis(const void *el1, const void *el2);// compare two uavs on their Xaxis (is u1 left of u2 ?)
int cmpfuncYaxis(const void *el1, const void *el2);// compare two uavs on their Yaxis (is u1 below u2 ?)

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

/**/
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
				net->counters[net->n_uavs]=1;// Increment number of ground nodes covered by uav (just 1)
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

sln* blankSolution(int size_input, int nuavs)
{
	int k=0,j=0;
	sln* res=(sln*)malloc(sizeof(sln));
	res->g_covers=(int*)calloc((size_input+1),sizeof(int));
	res->covers=(int**)malloc((size_input+1)*sizeof(int*));
	res->counters=(int*)calloc((size_input+1),sizeof(int));
	res->uavs=(double**)malloc((size_input+1)*sizeof(double*));
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
	res->g_covers=(int*)calloc((size_input+1),sizeof(int));
	res->covers=(int**)malloc((size_input+1)*sizeof(int*));
	res->counters=(int*)calloc((size_input+1),sizeof(int));
	/* Note : allocate large memory for uavs in case duplicates are needed. However, not useful if number of uavs > number of targets  */
	res->uavs=(double**)malloc((size_input+1)*sizeof(double*));/* Convenient to keep as whole size, makes no need to reallocate each time a uav is removed/added */
	double** cl=(double**)malloc((size_input+1)*sizeof(double*));/* temporary variables, contains sum of elements in cluster */
	res->dists=(double**)malloc((size_input+1)*sizeof(double*));/* matrix of intra-distances between uavs */
// 	if there are more uavs than allocated => segmentation fault. Then:assign larger space:eg: matrix[(size_input+1),(size_input+1)]
//	for(k=0;k<=size_input;k++)
//		res->distances[k]=(double*)calloc((size_input+1),sizeof(double));

	/* Initialisation steps */
	for(k=0;k<=size_input;k++)
	{
		res->covers[k]=(int*)calloc((size_input+1),sizeof(int));
		res->uavs[k]=(double*)calloc(dim,sizeof(double));
		res->dists[k]=(double*)calloc((size_input+1),sizeof(double));
		for(j=0;j<dim;j++)	res->uavs[k][j]=-1;// if a uav's coords is (-1,-1), then it's not deployed
		cl[k]=(double*)calloc(dim,sizeof(double));
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
			res->counters[res->n_uavs]++;// To [1
		}
	}

	/* update distance matrix since there were modifications done on the network */
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

	/* update distance matrix since there were modifications done on the network */
	updateDistMat(clusts, range);

	/* housekeeping */
	for(h=0;h<=clusts->n_uavs;h++)	free(c1[h]);
	free(c1);

	return error;
};


void connect_CCs(sln* net, igraph_t* Gk, double threshold, int** restr_list, int **pairs, int* npairs, bool restrict, int level)
{// The boolean value "restrict" says wether the restriction strategy should be used (restrict=false since graph G isn't the root graph G0) or not.
	// In the latter case, fill lists : **restr_list, **pairs, and * npairs. Otherwise, use these lists to apply the restrictions

//printf("net n_uavs %d threshold %f and restrict %d\n", net->n_uavs, threshold, restrict);
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


/*
//void translate(sln* net, double threshold, double* solverSln)
igraph_t* build_first_conn_graph(sln* net, double threshold, int** restr_list, double **pairs, int* npairs, igraph_t* coverageG)
{// note : computing diameter of graph (O(n*|E|)) can cost more than comparing each nodes (O(n^2)) : eg : Complete graph : n*|E| = n*sum_i_in{1...n-1}(i) > O(n^2) for n>3

	long int i=0,j=0;
	igraph_integer_t ncomps=-1;// numbers of connected components
	igraph_vector_t labels, compssizes;// labels of connected components for each vertex
	double min_distance=DBL_MAX;// distance of two closest nodes but out of range (belong to two different connected components)
	double current_distance;// used to place a new uav
	int n1=-1,n2=-1;//two closest but out of range nodes

//  here, the vectors "labels" are only used to create graph
//	int nActivUAVs=0;// used to keep record of the number of uavs
//	int indexActivUavs[net->n_uavs];// contains the indices of active uavs
	
//	Uncomment following only if linear solver is used
//	for(i=1;i<=net->n_uavs;i++)
//		if(solverSln[i]>0)
//		{// register indices of active uavs ( every uav i with solverSln[i]>0)
//			indexActivUavs[nActivUAVs]=i;// !!! "nActivUAVs" as index is not a mistake : increment and at the current last position, store new index
//			nActivUAVs++;
//		}

	igraph_t* G1=(igraph_t*)malloc(sizeof(igraph_t));// first connected graph to build.
	
	// start building graph
//	if(igraph_empty(&G1, nActivUAVs, IGRAPH_UNDIRECTED)==IGRAPH_EINVAL)
//		printf(" Invalid number of vertices, graph init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);
//	if(igraph_empty(&net->gr, net->n_uavs+1, IGRAPH_UNDIRECTED)==IGRAPH_EINVAL)
//		printf(" Invalid number of vertices, graph init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);


//printf(" %d active uavs\n", nActivUAVs);
	// add edges to graph : two uavs in each other range
//	for(i=0;i<nActivUAVs;i++)
//		for(j=i+1;j<nActivUAVs;j++)
//			if(net->dists[indexActivUavs[i]][indexActivUavs[j]]<threshold)
//{
//printf("Enter ---------  (%ld,%ld), (%ld,%ld)\t",i,j,(long int)indexActivUavs[i], (long int)indexActivUavs[j]);
//				igraph_add_edge(&G1, i, j);// since j = i+1 to nActivUAVs, then no loop
//}
//printf("end init graph with %d nodes\n",nActivUAVs);

		if(igraph_vector_init(&labels, 0)==IGRAPH_ENOMEM || igraph_vector_init(&compssizes, 0)==IGRAPH_ENOMEM)
			printf(" Memory issues, vector init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);			
		// compute number of connected components
		igraph_clusters(&G1, &labels, &compssizes, &ncomps, IGRAPH_WEAK);

printf("Graph has %d connected components\n",ncomps);

	if(igraph_empty(G1, net->n_uavs+1, IGRAPH_UNDIRECTED)==IGRAPH_EINVAL)
		printf(" Invalid number of vertices, graph initialisation failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);


	// Use of a brute force branching here. Maybe undeterministic could be better.
	// brute force : create one unique connected component by iteratively branching the 2 closest connected components
	// From here, a bit of a mess, does even create a segmentation fault
	int k=0, added=0;
	do
	{
		if(igraph_empty(G1, net->n_uavs+1, IGRAPH_UNDIRECTED)==IGRAPH_EINVAL)
			printf(" Invalid number of vertices, graph initialisation failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);

		for(i=1;i<=net->n_uavs;i++)
			for(j=i+1;j<=net->n_uavs;j++)
				if(net->dists[i][j]<threshold)
{
//printf("Enter ---------  (%ld,%ld), (%ld,%ld)\t",i,j,(long int)indexActivUavs[i], (long int)indexActivUavs[j]);
				igraph_add_edge(G1, i, j);// since j = i+1 to nActivUAVs, then no loop
				if (k==0) igraph_add_edge(coverageG, i, j);// store the first graph
}

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

		if(igraph_vector_init(&labels, 0)==IGRAPH_ENOMEM || igraph_vector_init(&compssizes, 0)==IGRAPH_ENOMEM)
			printf(" Memory issues, vector init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);
			
		// compute number of connected components
		igraph_clusters(G1, &labels, &compssizes, &ncomps, IGRAPH_WEAK);
//		igraph_clusters(&net->gr, &labels, &compssizes, &ncomps, IGRAPH_WEAK);

printf("Graph %d has %d connected components\n", k, ncomps);

//printf("start labels\n");
//for(i=0;i<igraph_vcount(&net->gr);i++)	printf("l[%li]:%li\t",i,(long int)VECTOR(labels)[i]);
//printf("end labels\n");

		// free first used memory since initiated again later
		igraph_vector_destroy(&compssizes);

		// if graph not a single connected component then update matrix of distance with new uavs. Note : exists isolated vertex 0
		if(ncomps>2)
		{//add new node to reach one connected component
			min_distance=DBL_MAX;
			// find two closest but out of range nodes
			for (i=1;i<=net->n_uavs;i++)
			{
				for (j=i+1;j<=net->n_uavs;j++)
				{
					// Not interested in nodes in same connected component
					if ( VECTOR(labels)[i] == VECTOR(labels)[j] )	continue;
					// Find closest clusters but out of range
					current_distance=net->dists[i][j];
//					if(current_distance<threshold)// skip, in range
//						continue;
					// check only if different connected component
					if(current_distance<min_distance)
					{// keep two nodes
						min_distance=current_distance;
						n1=i;
						n2=j;
					}
				}
			}
			
			// keep track of n1 and n2
			pairs[*npairs][0]=n1;
			pairs[*npairs][1]=n2;
			*npairs=*npairs+1;

			// free vector of labels
			igraph_vector_destroy(&labels);

			if(n1<0 || n2<0)
			{
				printf(" Something went wrong with graph (Either connectivity or indices. Check), adding new uav will fail %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);
				return NULL;
			}

			// assign coordinates of new created node : middle of segment [n1,n2]
			net->n_uavs++;
			// assign coordinates
			for (j=0;j<dim;j++)
				net->uavs[net->n_uavs][j]=(net->uavs[n1][j]+net->uavs[n2][j])/2;
			// increase number of elements in new cluster : at least itself, since may not cover any ground node, but be used only for connectivity
			net->counters[net->n_uavs]++;// To 1
			// extend size of graph
//			igraph_add_vertices(&net->gr, 1, 0);
			igraph_add_vertices(G1, 1, 0);			
			// update new distances
			updateDistMat(net, threshold);

			// add new uav to restriction list
			restr_list[added][0]=net->n_uavs;
			int col=1;// gives current column position in restricted list
			// create links with other uavs when needed
			for (j=1;j<net->n_uavs;j++)
			{
				if( net->dists[net->n_uavs][j] < threshold )
				{
					igraph_add_edge(G1, net->n_uavs, j);
					restr_list[added][col]=j;
					col++;
				}
			}
			added++;// finished adding links to the row => increment
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

	i=0;
	printf("Start here restricted list\n");
	while(restr_list[i++][0]>=0)
	{
		j=0;
		while(restr_list[i][j]>=0)
		{
			printf("%d\t",restr_list[i][j++]);
		}
		printf("\n");
	}

	// write graph into file
	FILE* fp;
	fp=fopen("genet_clust.csv","w");
	igraph_vector_t ite_edgs;
	igraph_vector_init(&ite_edgs, 0);
	igraph_get_edgelist(G1, &ite_edgs, 0);
	int n=igraph_ecount(G1);
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


	return G1;

};
*/


//void populate(sln* covers, double threshold, igraph_t* solnG0, igraph_t** solnsGraphs, int** restr_list, double **pairs, int* npairs)
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
printf("Before connectCCs vcount %li ecount %li\n", (long int)igraph_vcount(solnG0), (long int)igraph_ecount(solnG0));

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
	igraph_destroy(&free_sln->gr);
	free(free_sln);
};


int main(int argc, char** argv)
{

	clock_t begin = clock();

	/* !!! igraph : turn on attribute handling  Ie. even if you don't manipulate attributes explicitly, but create a graph that might have some attributes, eg. read a graph a GraphML file, you need this call before, otherwise the attributes are dropped. */
	igraph_i_set_attribute_table(&igraph_cattribute_table);

	readData(argv);

	bound_1=1000;
	bound_2=1000;

	double radius=uavs_range;

	sln *res=method1ePasse(grnds, nbr_grnds, radius);
printf("1one pass results with threshold %f and %d uavs\n", radius,res->n_uavs);
	k_means(grnds, nbr_grnds, res, 0.0001, radius);
	printf("phase 1 : clustering, final %d\n", res->n_uavs);


	int i=0,j=0;
	
//print distance matrix
printf("with threshold %f and %d uavs\n", radius,res->n_uavs);
for(i=1;i<=res->n_uavs;i++)
{
	for(j=i+1;j<=res->n_uavs;j++)
		if(res->dists[i][j] < radius)printf("(%d,%d)\t",i,j);
}
printf("\n");


	FILE* fp;
	fp=fopen("resclassiconepass.csv","w");
	for(i=1;i<=res->n_uavs;i++)
		for (j=0;j<dim;j++)
		{
			// skip comma not needed after last dim value
			if(j==dim-1)	fprintf(fp,"%lf\n", res->uavs[i][j]);
			else fprintf(fp,"%lf,", res->uavs[i][j]);
		}
	fclose(fp);

	// for the graph solutions : 60 individuals for each subpopulation (4) of the ring scheme
	igraph_t** allGs=(igraph_t**)malloc(60*sizeof(igraph_t*));// useless for now : used when problem modelled as a linear optimisation problem
	
	// "restr_list" is the restricted list of edges, that depicts the "root" or the list of edges added to the scattered graphs
	// in order to build  the first connected graph. "restr_list" is a matrix where : for i in {1 to the number of ground nodes} :
	// matrix[i][0] == -1, or matrix[i][0] > 0
	// 1. if matrix[i][0] > 0, then matrix[i][0] contains the index of a given uav connected to another existing uav => exists in graph an edge
	// containing value in matrix[i][0] : if an edge (u',matrix[i][0]) (or (u',u) where u==matrix[i][0]), then the remaining
	// cells : matrix[i][1] to matrix[i][m] contain the other m uavs (their indices) with wich u is connected to :
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
	// for the connectivity, thus plan more space
	int** pairs=(int**)malloc((nbr_grnds*10)*sizeof(int*));
	for (i=0;i<(nbr_grnds*10);i++)
	{
		pairs[i]=(int*)malloc(2*sizeof(int));
		for (j=0;j<2;j++)	pairs[i][j]=-1;
	}

	igraph_t* solG0;
	//igraph_t** solnsGraphs;
	//igraph_t* solnsGraphs;
	//igraph_t** solnsGraphs=(igraph_t**)malloc(NINDIVS*sizeof(igraph_t*));

	int npairs=0;
	//solnsGraphs[0]=build_first_conn_graph(res, radius, restr_list, pairs, &npairs, coverageG);
	//populate(res, radius, solG0, solnsGraphs, restr_list, pairs, &npairs);
	populate(res, radius, solG0, restr_list, pairs, &npairs);
// printf("phase 2 : first connected graph, final %d\n", res->n_uavs);


	i=0;
	while(pairs[i][0] > 0 || pairs[i][1] > 0)
	{
		if(pairs[i][0] > 0 && pairs[i][1] < 0 || pairs[i][0] > 0 && pairs[i][1] < 0)
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

};
