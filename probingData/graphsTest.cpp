
#include <stdlib.h>
#include <stdio.h>
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
	igraph_matrix_t dists;// squared matrix of distances between uavs. Type igraph to save space
	igraph_t gr;// corresponding graph of the solution. Needed for connectivity checking
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
void translate(sln* net);
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
			dist=( dist > range ? 0 : dist );
			MATRIX(net->dists, i, j)=dist;
			MATRIX(net->dists, j, i)=dist;
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
	igraph_matrix_init(&res->dists,(size_input+1),(size_input+1));
// 	if k<=max_uavs_avail and there are less uavs than needed then -> segmentation fault. Safety:distances** squared matrix  (size_input+1)*(size_input+1)
//	for(k=0;k<=size_input;k++)
//		res->distances[k]=(double*)calloc((size_input+1),sizeof(double));


	/* Initialisation steps */
	for(k=0;k<=size_input;k++)
	{
		res->uavs[k]=(double*)calloc(dim,sizeof(double));
		for(j=0;j<dim;j++)	res->uavs[k][j]=-1;// if uavs coords is (-1,-1), then not deployed
		cl[k]=(double*)calloc(dim,sizeof(double));
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

	updateDistMat(res, threshold*2);

	/* Housekeeping */
	for(k=0;k<=size_input;k++)
	{
		free(cl[k]);
	}
	free(cl);

	return res;
};


void translate(sln* net)
{
	int i=0,j=0;
/*
	int count=0;// number of edges of graph

	typedef struct linklist{
		int node1;
		int node2;
		struct linklist* next;
	}linklist;

	linklist* add=NULL;
	linklist** head=&add;
	i=j=0;

	// Builds a linked list of uavs that can communicate => gather Links for graph)
	for(i=1;i<=net->n_uavs;i++)
		for(j=i+1;j<=net->n_uavs;j++)
			if(inRange(net->uavs[i],net->uavs[j], uavs_range))
			{// build links whenever two uavs can communicate (are in admissible range)
				linklist* new_node=(linklist*) malloc(sizeof(linklist));
				new_node->node1=i;
				new_node->node2=j;
				new_node->next=*head;
				*head=new_node;
				count++;
			}
*/
//	igraph_vector_t edgs;
//	if(igraph_vector_init(&edgs, count*2)==IGRAPH_ENOMEM)
//		printf(" Memory issues, vector init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);

	i=0;
//	linklist* ite=*head;

/*
	while (ite != NULL)
	{
		VECTOR(edgs)[i]=ite->node1;
		VECTOR(edgs)[i+1]=ite->node2;
		ite=ite->next;
		i+=2;
	}
*/
	// igraph : from 0 to n-1 vertices
/*
	int buff=igraph_create(&net->gr, &edgs, net->n_uavs+1, false);
	if ( buff == IGRAPH_EINVEVECTOR || buff == IGRAPH_EINVVID )
			printf(" graph init issues (%s) failed %s, %d, %s \n", buff== IGRAPH_EINVEVECTOR ? "IGRAPH_EINVEVECTOR" : "IGRAPH_EINVVID" , __FILE__, __LINE__, __FUNCTION__);

	j=0;
	for (i=0; i<igraph_vector_size(&edgs); i+=2) {
		printf("%d) [%li-%li] \n", j, (long int) VECTOR(edgs)[i], (long int) VECTOR(edgs)[i+1] );
		j++;
	}
*/


	if ( igraph_weighted_adjacency(&net->gr, &net->dists, IGRAPH_ADJ_MAX, NULL, 1) == IGRAPH_NONSQUARE )
			printf(" Something went wrong for graph init, failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);

	FILE* fp;
	fp=fopen("graphs.csv","w");
	igraph_vector_t ite_edgs;
	igraph_vector_init(&ite_edgs, 0);
	igraph_get_edgelist(&net->gr, &ite_edgs, 0);
	int n=igraph_ecount(&net->gr);
	long int buff1,buff2;// needed for the cast
	double weight;
	for (i=0, j=0; j<n; i+=2, j++)
	{
		buff1=VECTOR(ite_edgs)[i];
		buff2=VECTOR(ite_edgs)[i+1];
		weight=EAN(&net->gr, "weight", j);
//		fprintf(fp,"%lf,%lf,%lf,%lf\n", net->uavs[buff1][0], net->uavs[buff1][0], net->uavs[buff2][0], net->uavs[buff2][0]);
		fprintf(fp,"%lf,%lf\n", net->uavs[buff1][0], net->uavs[buff1][1]);
		fprintf(fp,"%lf,%lf\n", net->uavs[buff2][0], net->uavs[buff2][1]);
		fprintf(fp,"\n");
//		printf("%ld,%ld,%lf\n",buff1,buff2,weight);
	}
	fclose(fp);


/*
	igraph_vector_ptr_t comps_list;
	igraph_vector_ptr_init(&comps_list, 0);// list of components first need to be initiated before call to decompose
*/
	/* 1) -1 : maxcompno, maximum number of components to return. -1 if nolimit
	 * 2)  2 : minelements, minimum number of vertices a component contains to be placed in the components vector. Here 2 skips isolated vertices */
/*
	int someintdontknowyetwhatfor=igraph_decompose(&gr, &comps_list,IGRAPH_WEAK, -1, 2);
	for(i=0;i<igraph_vector_ptr_size(&comps_list);i++)
	{
		igraph_t *buff=VECTOR(comps_list)[i];
		igraph_es_t es;
		igraph_es_all(buff, &es);
		while (!igraph_es_end(buff, &es))
		{
			igraph_real_t *from, *to;
			igraph_get_vertex_attribute(buff, "id", igraph_es_from(buff, &es), (void**) &from, 0);
			igraph_get_vertex_attribute(buff, "id", igraph_es_to(buff, &es), (void**) &to, 0);
			printf("%li %li\n", (long int) *from, (long int) *to);
			igraph_es_next(buff, &es);
		}
	}
*/
//	free_complist(&complist); */

//	igraph_destroy(&gr);
//	igraph_vector_destroy(&edgs);

	/* Housekeeping */
/*
	ite = *head;
	linklist* tmp;
	while (ite != NULL)
	{
		tmp = ite;
		ite=ite->next;
		free(tmp);
	}

	return gr;
*/
};


int main(int argc, char** argv)
{

	/* !!! igraph : turn on attribute handling  Ie. even if you don't manipulate attributes explicitly, but create a graph that might have some attributes, eg. read a graph a GraphML file, you need this call before, otherwise the attributes are dropped. */
	igraph_i_set_attribute_table(&igraph_cattribute_table);

	readData(argv);

	bound_1=1000;
	bound_2=1000;

	double threshold=uavs_range/2;

	sln *res=method1ePasse(grnds, nbr_grnds, threshold);

	int i=0,j=0;

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

	translate(res);
//	igraph_t graph_sol=translate(res);


	printf("%lf\n",grnds[17][0]);
}
