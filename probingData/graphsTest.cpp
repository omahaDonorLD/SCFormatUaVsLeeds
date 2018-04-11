
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <string>
#include <math.h>

#include <igraph.h>
#include <assert.h>
#include <float.h>

inline
int STREAM_FAIL(const char *FROM_FILE, int AT_LINE, const char *IN_FUNCTION)
{printf("STREAM_FAILURE, line %d, function %s, file %s\n", AT_LINE, IN_FUNCTION, FROM_FILE);return EXIT_FAILURE;};

inline
int MEMO_FAIL(const char *FROM_FILE, int AT_LINE, const char *IN_FUNCTION)
{printf("MEMO_ALLOC_FAILURE, line %d, function %s, file %s\n", AT_LINE, IN_FUNCTION, FROM_FILE);return EXIT_FAILURE;};

typedef struct sln{
	int* labels;// Size : the total number of elements (ground nodes), stores for each element the cluster it belongs to. Safety as also have info into "covers".
	double** uavs;/* matching "contains" list with "uavs_indices". The latter contains the indices 'i'
							of elements in list "uavs" so that one can finds the coordinates of uavs[i] */
	int* counters;// number of elements (size : n_uavs) covered by the uav
	// double** distances;// squared matrix of distances between uavs
	double** dists;// squared matrix of distances between uavs. Type igraph to save space
	igraph_t gr;// graph of the solution. Needed for connectivity checking
	int n_uavs;// number of nodes in network : either nodes are uavs (sln : set of uavs), or ground nodes (uav : set of ground nodes)
}sln;


int max_uavs_avail;
int nbr_grnds;
int dim=2;			// dimension of input vectors
double** grnds;		// All ground nodes coordinates
double uavs_range;	// all ranges of available uavs

// limit of map
double bound_1;
double bound_2;

void readData(char** argv);
double euclDistance(double *node1, double *node2);
void translate(sln* net, double threshold);
void updateDistMat(sln* net, double range);

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

}


double euclDistance(double *node1, double *node2)
{
	int i=0;
	double norm=0;
	for(i=0;i<dim;i++)	norm+=pow(node1[i]-node2[i],2);
    return sqrt(norm);
}

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
}


bool inRange(double* node1, double* node2, double range)
{
	return ( euclDistance(node1,node2) >= range ? false : true);
};

sln* method1ePasse(double** input_data, int size_input, double threshold)
{
	int k=2,j=1;
	sln* res=(sln*)malloc(sizeof(sln));
	res->labels=(int*)calloc((size_input+1),sizeof(int));
	res->counters=(int*)calloc((size_input+1),sizeof(int));
	res->uavs=(double**)malloc((size_input+1)*sizeof(double*));/* Convenient to keep as whole size, makes no need to reallocate each time a uav is removed/added */
	double** cl=(double**)malloc((size_input+1)*sizeof(double*));/* temporary variables, contains sum of elements in cluster */
	res->dists=(double**)malloc((size_input+1)*sizeof(double*));/* temporary variables, contains sum of elements in cluster */
// 	if k<=max_uavs_avail and there are less uavs than needed then -> segmentation fault. Safety:distances** squared matrix  (size_input+1)*(size_input+1)
//	for(k=0;k<=size_input;k++)
//		res->distances[k]=(double*)calloc((size_input+1),sizeof(double));


	/* Initialisation steps */
	for(k=0;k<=size_input;k++)
	{
		res->uavs[k]=(double*)calloc(dim,sizeof(double));
		for(j=0;j<dim;j++)	res->uavs[k][j]=-1;// if uavs coords is (-1,-1), then not deployed
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
	res->labels[1]=1;// label of first element : first cluster

	res->n_uavs=1;// Number of clusters so far : 1

	double distance_to_closest_centroid=bound_1*bound_2;// Worst : the extrem limits of the entire map
	double current_distance=0;

	/* start 1e pass algorithm */
	for(k=2;k<=size_input;k++)
	{
		distance_to_closest_centroid=bound_1*bound_2;// Worst : extreme considering map
		for(j=1;j<=res->n_uavs;j++)
		{
			current_distance=euclDistance(input_data[k],res->uavs[j]);
			if ( current_distance < distance_to_closest_centroid )
			{
				distance_to_closest_centroid=current_distance;// find the closest cluster
				res->labels[k]=j;
			}
		}


		if(distance_to_closest_centroid<threshold)
		{// closest centroid is within admissible range => add element to cluster
			res->counters[res->labels[k]]++;
			/* update centroids */
			for(j=0;j<dim;j++)
			{
				cl[res->labels[k]][j]+=input_data[k][j];
				res->uavs[res->labels[k]][j]=cl[res->labels[k]][j]/res->counters[res->labels[k]];
			}
		}
		else
		{// closest centroid is not within admissible range => Create own new cluster
			res->n_uavs++;
			res->labels[k]=res->n_uavs;
			for(j=0;j<dim;j++)
			{
				cl[res->labels[k]][j]=input_data[k][j];
				res->uavs[res->labels[k]][j]=input_data[k][j];
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


void k_means(double** data, int n, sln* clusts, double error_tolerance, double range)
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
                    clusts->labels[h] = i;
                }
            }
            /* update size and temp centroid of the destination cluster */
			for(j=0;j<dim;j++)	c1[clusts->labels[h]][j] += data[h][j];
			clusts->counters[clusts->labels[h]]++;
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
				for (j=i+1;j<=net->n_uavs;j++)
				{
					/* Find closest clusters but over range */
					current_distance=net->dists[i][j];
					if(current_distance<threshold)// skip, in range
						continue;
					// check only if different connected component
					if(current_distance<min_distance && VECTOR(labels)[i] != VECTOR(labels)[j])
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

	double radius=uavs_range/4;

	sln *res=method1ePasse(grnds, nbr_grnds, radius);

	int i=0,j=0;
	translate(res, radius);

printf("rl.csv nuavs : %d ",res->n_uavs);
	FILE* fp;
	fp=fopen("rl.csv","w");
	for(i=1;i<=res->n_uavs;i++)
		for (j=0;j<dim;j++)
		{
			// skip comma missing for last dim
			if(j==dim-1)	fprintf(fp,"%lf\n", res->uavs[i][j]);
			else fprintf(fp,"%lf,", res->uavs[i][j]);
		}
	fclose(fp);

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


	printf("n uvas : %d and grnds[17][0] : %lf\n", res->n_uavs,grnds[17][0]);
}
