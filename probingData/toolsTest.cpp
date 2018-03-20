
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

int max_uav_avail;
int nbr_grnds;
aNode* GRNDS;		// All ground nodes coordinates
double* UAVs_Range;	// all ranges of available uavs

// Elbow method variables
double elbow_param;	// heuristic for finding "k" needed 
//double **sim_mat;	// Euclidian

// One pass method variables
aNode **C_l;// Worst case N clusters
aNode *rl;// "centroids"
int* counters;
int K=1;// Number of clusters


void readData(char** argv);
double euclDistance(aNode *uav1, aNode *uav2);
void method1ePasse();


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



void method1ePasse()
{
	double threshold=250/2;// radius of range uav
	C_l=(aNode**)malloc((nbr_grnds+1)*sizeof(aNode*));// Worst case N clusters
	rl=(aNode*)malloc((nbr_grnds+1)*sizeof(aNode));// "centroids"
	counters=(int*)malloc((nbr_grnds+1)*sizeof(int));
	int k=2,j=1, index_closest_centroid=0;

	for(k=0;k<=nbr_grnds;k++)
	{
		C_l[k]=(aNode*)malloc((nbr_grnds+1)*sizeof(aNode));// Allocate memory for clusters
		counters[k]=0;
	}

	rl[1]=GRNDS[1];
	counters[1]++;
	C_l[1][1]=GRNDS[1];
	K=1;// Number of clusters

	double distance_to_closest_centroid=1000*1000;// Worst : the extrem limits of the whole map
	double current_distance=1000;

	for(k=2;k<=nbr_grnds;k++)
	{
		distance_to_closest_centroid=1000*1000;// Worst :
		index_closest_centroid=1;
		for(j=1;j<=K;j++)
		{
			current_distance=euclDistance(&GRNDS[k],&rl[j]);
			if ( current_distance < distance_to_closest_centroid )
			{
				distance_to_closest_centroid=current_distance;// Worst
				index_closest_centroid=j;
			}
		}
		

		if(distance_to_closest_centroid<threshold)
		{
			counters[index_closest_centroid]++;
			C_l[index_closest_centroid][counters[index_closest_centroid]]=GRNDS[k];
			rl[index_closest_centroid]=newCentroid(C_l[index_closest_centroid],counters[index_closest_centroid]);// overwrite previous centroid
		}
		else{
			K++;
			rl[K]=GRNDS[k];// overwrite previous centroid
			counters[K]++;// To [1]
			C_l[K][1]=GRNDS[k];// C_l[K][counters[K]]
		}
	}

}

int main(int argc, char** argv)
{
	readData(argv);
	method1ePasse();

	FILE* fp;
	fp=fopen("rl.csv","w");
	int i=0;
	for(i=1;i<K;i++)
		fprintf(fp,"%lf,%lf\n",rl[i].x,rl[i].y);
	fclose(fp);

	
	printf("%lf\n",GRNDS[17].x);
}
