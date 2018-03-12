
#include "../heads/tools.hpp"

using namespace std;

void readData(char** input_param, char** input_data)
{
	FILE* fp;
	char *tmp[100];// Just a buffer
	
	// parameters
	// .. add sources for reading parameters
	fscanf(fp,"%[^:]:%lf", &aha[0], &ahaa, &ahao);
	
	// data
	fp=fopen(argv[2],"r");

	if( fscanf(fp,"%d", &nbr_grnds) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}

	/* allocate memory for ground nodes */
	GRNDS=malloc(nbr_grnds*sizeof(aNode));
	if(GRNDS==NULL){ /* memory allocation failure */ MEMO_FAIL(__LINE__, __FILE__, __FUNCTION__); }

	int i=0;
	for(i=0;i<nbr_grnds;i++)
	{
		GRNDS[i].index=i;// Unique, defines a point : for either a ground node or a uav
		fscanf(fp,"%lf,%lf", &GRNDS[i].x, &GRNDS[i].y);
	}

}
