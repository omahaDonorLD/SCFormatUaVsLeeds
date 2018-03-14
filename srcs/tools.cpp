
#include "../heads/tools.hpp"

using namespace std;

void readData(char** argv)
{
	FILE* fp;
	char *tmp[100];// Just a buffer
	int tmp2;// this one too

	// parameters
	fp=fopen(argv[1],"r");
	if( !fp ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}

	// Problem parameters
	if( fscanf(fp,"%s", tmp) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// "PRBLM_CRIT", text not needed
	if( fscanf(fp,"%s %d", tmp, &nbr_inds) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// Number of solutions by generations
	if( fscanf(fp,"%s %lf", tmp, &mutaprob) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// Number of generations to perform
	if( fscanf(fp,"%s %lf", tmp, &rank1maxsize) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}
	if( fscanf(fp,"%s %lf", tmp, &crowdSave) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}
	if( fscanf(fp,"%s %lf", tmp, &crowdTotal) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}
	//	fscanf(fp,"%[^:]:%lf", &aha[0], &ahaa, &ahao);// Alternative of reading instream

	// stopping criteria
	if( fscanf(fp,"%s", tmp) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// "STOP_CRIT", text not needed
	if( fscanf(fp,"%s %d", tmp, &max_time) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// Number of generations to perform
	if( fscanf(fp,"%s %d", tmp, &nbr_gen) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// Number of generations to perform
	if( fscanf(fp,"%s %d", tmp, &no_improv_ites) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// Number of iterations allowed to continue with where the improvements of fitness is less than "epsilon"
	if( fscanf(fp,"%s %lf", tmp, &epsilon) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// see previous line
	if( fscanf(fp,"%s %d", tmp, &hyper_vol_CD) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// needed to compute the HyperVolume for NSGA-II
	if( fscanf(fp,"%s %lf", tmp, &evo_part_percent) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}//
	if( fscanf(fp,"%s %lf", tmp, &grasp_part_percent) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}//

	// grasp parameters
	if( fscanf(fp,"%s", tmp) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// "GRASP_CRIT", text not needed
	if( fscanf(fp,"%s %lf", tmp, &grasp_split) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}
	if( fscanf(fp,"%s %d", tmp, &n_delta) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}
	if( fscanf(fp,"%s %lf", tmp, &alpha) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}

	// NSGA-II parameters
	if( fscanf(fp,"%s", tmp) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// "NSGA_CRIT", text not needed
	if( fscanf(fp,"%s %d", tmp, &tmp2) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}
	TESTMODE=(tmp2==1?true : false);
	if( fscanf(fp,"%s %d", tmp, &tmp2) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}
	removeclone=(tmp2==1?true : false);
	if( fscanf(fp,"%s %d", tmp, &tmp2) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}
	parentarena=(tmp2==1?true : false);
	if( fscanf(fp,"%s %d", tmp, &tmp2) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}
	gnuplot=(tmp2==1?true : false);
	if( fscanf(fp,"%s %d", tmp, &tmp2) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}
	results=(tmp2==1?true : false);

	fclose(fp);


	// read data (coordinates ground nodes)
	fp=fopen(argv[2],"r");

	// read number of available uavs
	if( fscanf(fp,"%d", &max_uav_avail) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}	
	UAVs_Range=malloc(max_uav_avail*sizeof(double));
	int i=0;
	// read range of each of them 
	for(i=0;i<max_uav_avail;i++)
		fscanf(fp,"%lf", &UAVs_Range[i]);


	// read number of ground nodes and then the coordinates
	if( fscanf(fp,"%d", &nbr_grnds) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}

	/* allocate memory for ground nodes */
	GRNDS=malloc(nbr_grnds*sizeof(aNode));
	if(GRNDS==NULL){ /* memory allocation failure */ MEMO_FAIL(__LINE__, __FILE__, __FUNCTION__); }

	for(i=0;i<nbr_grnds;i++)
	{
		GRNDS[i].index=i;// Unique, defines a point : for either a ground node or a uav
		fscanf(fp,"%lf,%lf", &GRNDS[i].x, &GRNDS[i].y);
	}

	fclose(fp);
}


void writeData(sln a_sln){
	int i=0;
	FILE *f;
	char path[30];
	char buff[30];
/*
	strcpy(path,"../optimum/optimum250.txt");
	strcpy(buff,argv[4]);
	strtok(buff,"/");
	strcat(path,strtok(NULL,"/"));

	f=fopen(path,"w");
*/	
}


/** 	\brief Check wether the ground node is covered by the uav
 *		\param a pointer on the uav 
 */
bool inRange(aUav* uav, aNode* ground)
{
	return ( sqrt( pow2( uav->idx - ground->x ) + pow2( uav->y - ground->y ) ) > range ? false, true);
};
