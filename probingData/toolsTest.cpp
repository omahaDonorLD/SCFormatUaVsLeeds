
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

/* Not needed, labels contained in struct solution for each gromax_uav_availund node
typedef struct linklist
{// Choice of LL for ground nodes : more dynamic size and easier to add/remove element at cheap cost both for memory and number of operations
	int id;// index in array of UAVs/Ground nodes
	struct linklist *next;
}linklist;// works for both uav and set of uavs and also for clustering step.
*/

typedef struct sln{
	int* labels;// Size : the total number of elements (ground nodes), stores for each element the cluster it belongs to. Safety as also have info into "covers".
/*	linklist **covers;/* List of n_uavs pointers, each pointing to the head of the linkedlist of indices of ground nodes it covers.
								Linked list since easier to add/remove ground nodes easily and
								keep track of the number of ground nodes within */
	double** uavs;/* matching "contains" list with "uavs_indices". The latter contains the indices 'i'
							of elements in list "uavs" so that one can finds the coordinates of uavs[i] */
	int* counters;// number of elements (size : n_uavs) covered by the uav
	int n_uavs;// number of nodes in network : either nodes are uavs (sln : set of uavs), or ground nodes (uav : set of ground nodes)
}sln;


int max_uavs_avail;
int nbr_grnds;
int dim=2;			// dimension of input vectors
double** grnds;		// All ground nodes coordinates
// Solution above (A unique list of uavs) not relevant as possible to have many solutions at a time
//double** uavs;		// All available uavs : for uavs init to (-1,-1) => not used
//double* uavs_ranges;	// all ranges of available uavs
double uavs_range;	// all ranges of available uavs


// limit of map
double bound_1;
double bound_2;

// Elbow method variables
double elbow_param;	// heuristic for finding "k" needed
//double **sim_mat;	// Euclidian

// One pass method variables


void readData(char** argv);
double euclDistance(double *node1, double *node2);
sln* method1ePasse(double** input_data, int size_input, double threshold);
void k_means(double** data, int n, sln* clusters, double error_tolerance);
igraph_t* translate(sln* net);

//void addToLinkList(linklist** head, int new_index);
//void removeInLinkList(linklist** head, int del_index);

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


/*
void addToLinkList(linklist** head, int new_index)
{// Backward : new element added to the top
	linklist* new_node=(linklist*) malloc(sizeof(linklist));
	new_node->id=new_index;
	new_node->next=*head;
	*head=new_node;
}


void removeInLinkList(linklist** head, int del_index)
{
	linklist* temp = *head;
	linklist* prev;

	// If head node itself holds the key to be deleted
	if (temp != NULL && temp->id == del_index)
	{
		*head=temp->next;// Changed head
		free(temp);// free old head
		return;
	}

	// Search for the key to be deleted, keep track of the previous node since one needs to change 'prev->next'
	while (temp != NULL && temp->id != del_index)
	{
		prev = temp;
		temp = temp->next;
	}

	 // If key not present in linked list
	 if (temp == NULL)
	 {
		 printf("CAREFUL KEY MISSING : %s %d\n", __FUNCTION__, __LINE__);
		 return;
	 }

	 // Unlink the node from linked list
	 prev->next = temp->next;

	 free(temp);  // Free memory
}
*/


sln* method1ePasse(double** input_data, int size_input, double threshold)
{
	int k=2,j=1;
	sln* res=(sln*)malloc(sizeof(sln));
	res->labels=(int*)calloc((size_input+1),sizeof(int));
	res->counters=(int*)calloc((size_input+1),sizeof(int));
	res->uavs=(double**)malloc((size_input+1)*sizeof(double*));/* Convenient to keep as whole size, makes no need to reallocate each time a uav is removed/added */
	double** cl=(double**)malloc((size_input+1)*sizeof(double*));/* temporary variables, contains sum of elements in cluster */

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


//	res->covers=(linklist**)malloc((K+1)*sizeof(linklist*));
//		for(j=1;j<=res->counters[k];j++)
//			addToLinkList(&res->covers[k], elts[k][counters[k]]);

	/* Housekeeping */
	for(k=0;k<=size_input;k++)
	{
		free(cl[k]);
	}
	free(cl);

	return res;
}


void k_means(double** data, int n, sln* clusts, double error_tolerance)
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


	// rebuild the variables 
	// first free previous
/*
	linklist* ite;
	linklist* tmp;
	for(i=0;i<=n;i++)
	{
		ite=clusts->covers[i];
		while (ite != NULL)
		{
			tmp=ite;
			ite=ite->next;
			free(tmp);
		}
		clusts->covers[i]=NULL;
		addToLinkList(&clusts->covers[i], i);

	}
*/

	/* housekeeping */
	for(h=0;h<=clusts->n_uavs;h++)	free(c1[h]);
	free(c1);

}


igraph_t* translate(sln* net)
{
	
};


int main(int argc, char** argv)
{
	readData(argv);

	double threshold=uavs_range;

	sln *res=method1ePasse(grnds, nbr_grnds, threshold);

	FILE* fp;
	fp=fopen("rl.csv","w");
	int i=0,j=0;
	for(i=1;i<=res->n_uavs;i++)
		for (j=0;j<dim;j++)
		{
			// skip comma missing for last dim
			if(j==dim-1)	fprintf(fp,"%lf\n", res->uavs[i][j]);
			else fprintf(fp,"%lf,", res->uavs[i][j]);
		}
	fclose(fp);

	// 2nd tour with larger range :
	threshold=250;
	sln *res2=method1ePasse(res->uavs, res->n_uavs, threshold);

	fp=fopen("rl2.csv","w");
	for(i=1;i<=res2->n_uavs;i++)
		for (j=0;j<dim;j++)
		{
			// skip comma missing for last dim
			if(j==dim-1)	fprintf(fp,"%lf\n", res2->uavs[i][j]);
			else fprintf(fp,"%lf,", res2->uavs[i][j]);
		}
	fclose(fp);

	// Using k means :
	k_means(grnds, nbr_grnds, res, 0.0001);

	fp=fopen("rlkmeans.csv","w");
	for(i=1;i<=res->n_uavs;i++)
		for (j=0;j<dim;j++)
		{
			// skip comma missing for last dim
			if(j==dim-1)	fprintf(fp,"%lf\n", res->uavs[i][j]);
			else fprintf(fp,"%lf,", res->uavs[i][j]);
		}
	fclose(fp);

	printf("%lf\n",grnds[17][0]);
}
