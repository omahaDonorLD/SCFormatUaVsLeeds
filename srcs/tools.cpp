
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
	if( fscanf(fp,"%s %s", tmp, opt_path) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// path of optimal solutions
	if( fscanf(fp,"%s %d", tmp, &nbr_inds) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// Number of solutions by generations
	if( fscanf(fp,"%s %d", tmp, &nbr_subpop) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}// Number of solutions by generations
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


	// read data (number of max avail uavs, number and coordinates of ground nodes)
	fp=fopen(argv[2],"r");

	// read number of available uavs
	if( fscanf(fp,"%d", &max_uavs_avail) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}

	// read range of uavs
	if( fscanf(fp,"%lf", &uavs_range) < 0 ){STREAM_FAIL(__FILE__, __LINE__, __FUNCTION__);}
	int i=0,j=0;

/* For the moment all uavs are supposed to have the same range
//	uavs_ranges=malloc((max_uavs_avail+1)*sizeof(double));
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

	for(i=0;i<=nbr_grnds;i++)
	{
		grnds[i]=(double*)calloc(dim,sizeof(double));// Unique, defines a point : for either a ground node or a uav
		if(i==0)	continue;// only assign values from 1 to nbr_grnds
		for (j=0;j<dim;j++)
		{
			// skip missing comma on last dim
			if(j==dim-1)	fscanf(fp,"%lf", &grnds[i][j]);
			else fscanf(fp,"%lf,", &grnds[i][j]);
		}
	}

	fclose(fp);
};


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
};


double euclDistance(double *node1, double *node2)
{
	int i=0;
	double norm=0;
	for(i=0;i<dim;i++)	norm+=pow(node1[i]-node2[i],2);
    return sqrt(norm);
};

/** 	\brief Check wether the ground node is covered by the uav
 *		\param a pointer on the uav
 */
bool inRange(double* node1, double* node2, double range)
{
	return ( euclDistance(node1,node2) >= range ? false : true);
};


void updateDistMat(sln* net, double range)
{// updates the matrix of distances to avoid constantly having to compute them
	int i,j;
	double dist;
	for(i=1;i<=net->n_uavs;i++)
		for(j=i+1;j<=net->n_uavs;j++)
		{
			dist=euclDistance(net->uavs[i],net->uavs[j]);
			dist=( dist > range ? 0 : dist );
			net->distances[i][j]=net->distances[j][i]=dist;
		}
}



sln* method1ePasse(double** input_data, int size_input, double threshold)
{
	int k=2,j=1;
	sln* res=(sln*)malloc(sizeof(sln));
	res->labels=(int*)calloc((size_input+1),sizeof(int));
	res->counters=(int*)calloc((size_input+1),sizeof(int));
	res->uavs=(double**)malloc((size_input+1)*sizeof(double*));/* Convenient to keep as whole size, makes no need to reallocate each time a uav is removed/added */
	double** cl=(double**)malloc((size_input+1)*sizeof(double*));/* temporary variables, contains sum of elements in cluster */
	// compute and fill the matrix of distances between uavs
	res->distances=(double**)malloc(res->n_uavs*sizeof(double*));
	for(k=1;k<=res->n_uavs;k++)
		res->distances[k]=(double*)calloc(res->n_uavs,sizeof(double));


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

	updateDistMat(res);

	/* Housekeeping */
	for(k=0;k<=size_input;k++)
	{
		free(cl[k]);
	}
	free(cl);

	return res;
};


void k_means(double** data, int n, sln* clusts, double error_tolerance)
{

    int h, i, j; /* loop counters, of course :) */

    double old_error, error = DBL_MAX; /* DBL_MAX : macro, defines the maximum finite floating-point value : 1E+37 */
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

	updateDistMat(clusts);

	/* housekeeping */
	for(h=0;h<=clusts->n_uavs;h++)	free(c1[h]);
	free(c1);

};
