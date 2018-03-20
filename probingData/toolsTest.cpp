
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string>
#include <math.h>


inline
int STREAM_FAIL(const char *FROM_FILE, int AT_LINE, const char *IN_FUNCTION)
{printf("STREAM_FAILURE, line %d, function %s, file %s\n", AT_LINE, IN_FUNCTION, FROM_FILE);return EXIT_FAILURE;};

inline
int MEMO_FAIL(const char *FROM_FILE, int AT_LINE, const char *IN_FUNCTION)
{printf("MEMO_ALLOC_FAILURE, line %d, function %s, file %s\n", AT_LINE, IN_FUNCTION, FROM_FILE);return EXIT_FAILURE;};


/*** Structures
 */
// Each node have coordinates
typedef struct aNode {
	int index;// Unique, defines a point : for either a ground node or a uav
	double x;
	double y;
}aNode;


// A UaV is also a node
// Even though igraph uses specific type, this structure that specifies the coordinates of a UaV in the netwprk
typedef struct aUav{
	aNode identif;
	//long gene=0;
	int covers;	// A UaV contains at least 0 ground nodes
	double range;  	// If "contains" < 0 then is ground node and its range is also < 1
	bool active;  // ground nodes are always unactive. a uav can also be
}aUav;


typedef struct aNode_cluster{
	aNode **C_l;// Worst case N clusters
	aNode *rl;// "centroids"
	int* counters;// number of element in each cluster
	int K;// number of clusters
}aNode_cluster;

int max_uav_avail;
int nbr_grnds;
aNode* GRNDS;		// All ground nodes coordinates
double* UAVs_Range;	// all ranges of available uavs


// limit of map
double bound_1;
double bound_2;

// Elbow method variables
double elbow_param;	// heuristic for finding "k" needed
//double **sim_mat;	// Euclidian

// One pass method variables


void readData(char** argv);
double euclDistance(aNode *uav1, aNode *uav2);
aNode_cluster* method1ePasse(aNode* input_data, int* size_input, double *threshold);

void readData(char** argv)
{
	FILE* fp;
	char *tmp[100];// buffer
	int tmp2;// this too


	// read data (coordinates ground nodes)
	fp=fopen(argv[1],"r");

	// read number of available uavs
	if( fscanf(fp,"%d", &max_uav_avail) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}
	UAVs_Range=(double*)malloc((max_uav_avail+1)*sizeof(double));
	int i=0;
	// read range of each of them
	for(i=1;i<=max_uav_avail;i++)
		fscanf(fp,"%lf", &UAVs_Range[i]);


	// read number of ground nodes and then the coordinates
	if( fscanf(fp,"%d", &nbr_grnds) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}

	/* allocate memory for ground nodes */
	GRNDS=(aNode*)malloc((nbr_grnds+1)*sizeof(aNode));
	if(GRNDS==NULL){ /* memory allocation failure */ MEMO_FAIL(__FILE__, __LINE__, __FUNCTION__); }

	for(i=1;i<=nbr_grnds;i++)
	{
		GRNDS[i].index=i;// Unique, defines a point : for either a ground node or a uav
		fscanf(fp,"%lf,%lf", &GRNDS[i].x, &GRNDS[i].y);
	}

	fclose(fp);

}


double euclDistance(aNode *node1, aNode *node2)
{
    return sqrt(pow(node1->x-node2->x,2)+pow(node1->y-node2->y,2));
}

aNode newCentroid(aNode *Glm, int size)
{
	aNode rlm;
	rlm.x=0;
	rlm.y=0;
	int i=0;
	for(i=1;i<=size;i++)
	{
		rlm.x+=Glm[i].x;
		rlm.y+=Glm[i].y;
	}
	rlm.x/=size;
	rlm.y/=size;


    return rlm;
}



aNode_cluster* method1ePasse(aNode* input_data, int* size_input, double *threshold)
{
	aNode** C_l=(aNode**)malloc((*size_input+1)*sizeof(aNode*));// Worst case N clusters
	aNode* rl=(aNode*)malloc((*size_input+1)*sizeof(aNode));// "centroids"
	int *counters=(int*)malloc((*size_input+1)*sizeof(int));
	int k=2,j=1, index_closest_centroid=0;

	for(k=0;k<=*size_input;k++)
	{
		C_l[k]=(aNode*)malloc((*size_input+1)*sizeof(aNode));// Allocate memory for clusters
		counters[k]=0;
	}

	rl[1]=input_data[1];
	counters[1]++;
	C_l[1][1]=input_data[1];
	int K=1;// Number of clusters

	double distance_to_closest_centroid=bound_1*bound_2;// Worst : the extrem limits of the whole map
	double current_distance=0;

	for(k=2;k<=*size_input;k++)
	{
		distance_to_closest_centroid=bound_1*bound_2;// Worst :
		index_closest_centroid=1;// Default
		for(j=1;j<=K;j++)
		{
			current_distance=euclDistance(&input_data[k],&rl[j]);
			if ( current_distance < distance_to_closest_centroid )
			{
				distance_to_closest_centroid=current_distance;// find the closest cluster
				index_closest_centroid=j;
			}
		}


		if(distance_to_closest_centroid<*threshold)
		{
			counters[index_closest_centroid]++;
			C_l[index_closest_centroid][counters[index_closest_centroid]]=input_data[k];
			rl[index_closest_centroid]=newCentroid(C_l[index_closest_centroid],counters[index_closest_centroid]);// overwrite previous centroid
		}
		else{
			K++;
			rl[K]=input_data[k];// overwrite previous centroid
			counters[K]++;// To [1]
			C_l[K][1]=input_data[k];// C_l[K][counters[K]]
		}
	}

	// Write results
	aNode_cluster *results=(aNode_cluster*)malloc(sizeof(aNode_cluster));
	results->C_l=(aNode**)malloc((K+1)*sizeof(aNode*));// Worst case N clusters
	results->rl=(aNode*)malloc((K+1)*sizeof(aNode));// "centroids"
	results->counters=(int*)malloc((K+1)*sizeof(int));
	results->K=K;
	for(k=1;k<=K;k++)
	{
		results->C_l[k]=(aNode*)malloc((counters[k]+1)*sizeof(aNode));// Worst case N clusters
		results->counters[k]=counters[k];
		results->rl[k]=rl[k];
		for(j=1;j<=counters[k];j++)
				results->C_l[k][j]=C_l[k][j];
	}

	for(k=0;k<=*size_input;k++)	free(C_l[k]);
	free(C_l);
	free(rl);
	free(counters);

	return results;
}

int main(int argc, char** argv)
{
	readData(argv);

	bound_1=1000;
	bound_2=1000;

	double threshold=250/2;

	aNode_cluster *results=method1ePasse(GRNDS, &nbr_grnds, &threshold);

	FILE* fp;
	fp=fopen("rl.csv","w");
	int i=0;
	for(i=1;i<results->K;i++)
		fprintf(fp,"%lf,%lf\n",results->rl[i].x,results->rl[i].y);
	fclose(fp);

	threshold=250;
	aNode_cluster *results2=method1ePasse(results->rl, &results->K, &threshold);

	fp=fopen("rl2.csv","w");
	for(i=1;i<results2->K;i++)
		fprintf(fp,"%lf,%lf\n",results2->rl[i].x,results2->rl[i].y);
	fclose(fp);

	printf("%lf\n",GRNDS[17].x);
}
