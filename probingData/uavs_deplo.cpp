

//#include "linear_solver.hpp"
//#include <stdlib>
#include <iostream>
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
#include <map>
#include <tuple>
#include <stack>

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
	vector<vector<int>> outdeg;// 1. outdegree relation of uav, note: j < outdeg[j] to avoid symmetries
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
vector<double*>* onepassmethod(double** input_data, int nbr_grnds, double range);
double k_means(double** data, int ngrounds, vector<double*>* clusts, double error_tolerance, double range);
double euclDistance(double *node1, double *node2);
void duplicate(double lb, sln* net, int grndi, double range);
igraph_t* translate(sln* net, double threshold, double* solverSln);
void updateDistMat(sln* net, double range);
void addTonetwork(sln* net, double *toadd, double range);

bool uav_in_cover(vector<int> &gcovs, int uavindex);
void find_covers(sln* net, int uavj, double range);/* Find every covers for each ground node */

map<int,double>* solve_linear_model(sln* net, double range, double lb);

//void connect_CCs(sln* net, igraph_t* Gk, double range, map<int,double>* uavscoverage, vector<int>* uavsconnectivity, stack<tuple<int,int>>* pairs, bool dorestriction);
void connect_CCs(sln* net, igraph_t* Gk, double range, vector<int>* uavsconnectivity, stack<tuple<int,int>>* pairs, bool dorestriction);
//void populate(sln* covers, map<int,double>* uavscoverage, vector<int>* uavsconnectivity, double range, igraph_t* solnG0, stack<tuple<int,int>>* pairs);
void populate(sln* covers, vector<int>* uavsconnectivity, double range, igraph_t* solnG0, stack<tuple<int,int>>* pairs);


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

void updateDistMat(sln* net, double range)
{// updates the matrix of distances to avoid constantly having to compute the euclidian distance
// reminder : vector<vector<int>> outdeg;// 1. outdegree relation of uav, note: j < outdeg[j] to avoid symmetries
//			  vector<vector<double>> distuav;// dists[j][outdeg[j]] = distance between uav j and uav outdeg[j]

	int i,j;
	double dist=0;
	vector<tuple<int,int>> todel;// edges to eventually delete if new distance greater than 'range'

	for(i=0;i<net->uavs.size();i++)
	{
		for(j=0;j<net->outdeg[i].size();j++)
		{
			int adjac_uav=net->outdeg[i][j];
			dist=euclDistance(net->uavs[i],net->uavs[adjac_uav]);
			if (dist>range){
				// remove edge if j and adjac_uav are too far apart, or one of the two uavs is unactive
				todel.push_back(make_tuple(i,j));
			}
			net->distuav[i][j]=dist;// expected to exist, created with calls to 'addTonetwork' 
		}
	}

	if(!todel.empty()){
		for(i=0; i<todel.size(); i++){
			int buffi=get<0>(todel[i]);
			int buffj=get<1>(todel[i]);
			net->outdeg[buffi].erase(net->outdeg[buffi].begin()+buffj);
			net->distuav[buffi].erase(net->distuav[buffi].begin()+buffj);
		}
	}
};


void addTonetwork(sln* net, double *toadd, double range)
{// updates the matrix of distances to avoid constantly having to compute them

	net->uavs.push_back(toadd);
	int i, newuav=net->uavs.size()-1;
	double dist=0;

	// or could also : push_back(std::vector<int>()), then net->outdeg[i].push_back()
	net->outdeg.resize(net->outdeg.size()+1);
	net->distuav.resize(net->distuav.size()+1);
	net->uavcovs.resize(net->uavcovs.size()+1);

	for(i=0;i<newuav;i++)
	{
		dist=euclDistance(net->uavs[i],net->uavs[newuav]);
		if(dist<=range){
			net->outdeg[i].push_back(newuav);
			net->distuav[i].push_back(dist);
		}
	}
	find_covers(net, newuav, range);// find the ground the new uav covers
};


void duplicate(double lb, sln* net, int grndi, double range)
{/* lb for knowing if ground needs more covers, counts[] the list which gives for each ground node the number of uavs covering it,
	covers[] the list of aforementioned covers */

	int j=0, n=0, newuavs=0, uavj=0, randindex=0;
	double *buff;
cout<<"Duplicate uav for node "<<grndi<<endl;
	do
	{
		// duplicate one of the uavs covering ground node i
		randindex=rand()%net->gcovs[grndi].size();
		uavj=net->gcovs[grndi][randindex];
cout<<"Selected "<<net->uavs[uavj][0]<<","<<net->uavs[uavj][1]<<endl;
		buff=new double[dim];
		for(j=0;j<dim;j++)
			buff[j]=net->uavs[uavj][j];
		addTonetwork(net, buff, range);// update distance of network of uavs, does also call to 'find_covers'
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



void find_covers(sln* net, int uavj, double range)
{/* Find every covers for each ground node */

	int i=0,j=0;
		
	for(i=0;i<nbr_grnds;i++)
	{
//		if( euclDistance(net->uavs[uavj],grnds[i]) <= range && !uav_in_cover(net->gcovs[i], uavj))
		if( euclDistance(net->uavs[uavj],grnds[i]) <= range )
		{
if(i==63) cout << " uav: " << uavj << " - " << net->uavs[uavj][0] << "," << net->uavs[uavj][1] << endl;
			net->gcovs[i].push_back(uavj);
			net->uavcovs[uavj].push_back(i);
		}
	}
};


map<int,double>* solve_linear_model(sln* net, double range, double lb)
{

	int i=0,j=0;

//	/* find covers for each ground node */ : done every time a new uav is added
//	find_covers(net, range);

//	bool changedstate=false;
	/* Find at least one uav that needs to be replicated so that a minimum cover constraint isn't violated : at least lb covers for each target */
	for(i=0;i<nbr_grnds;i++)
	{
//		/* Create duplicate if constraint not satisfied => select randomly uavs covering i to duplicate, and update covers */
		if( net->gcovs[i].size() < lb )
		{
printf("Duplicate uavs for gnode %d, which is covered with %d uavs \n", i, net->gcovs[i].size());
for(int z=0;z<net->uavs.size();z++){
//int buffz=net->gcovs[i][z];
cout<<net->uavs[z][0]<<","<<net->uavs[z][1]<<endl;
}
			duplicate(lb, net, i, range);
//			if(!changedstate)	changedstate=true;
//			break;
		}
	}
/*	
	// if uavs added in duplication phase, then search again covers with ground nodes
	if(changedstate)
	{

TO ERASE : Done in addTonetwork		
		find_covers(net, range);
	}
*/

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
	map<int,double> *soln=new map<int,double>;
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
		double buffval=glp_mip_col_val(lp, j+1);
		if(buffval>0)
			(*soln)[j]=buffval;
//		if (soln[j] - ceil(soln[j]) != 0)	printf(" soln %d not an integer = %f\n", j, soln[j]);
	}
	
	glp_print_sol(lp, "glpksols.txt");

	glp_delete_prob(lp);


	for(map<int,double>::iterator it=(*soln).begin(); it!=(*soln).end(); it++)
		printf("[ %d : %f ] ", it->first, it->second);
	
	printf("\nRO == %f, max = %f, min = %f, aver = %f, lessThanAvrg = %f, moreThanvrg = %f\n",sumrows, maxrowval, minrowval, sumrows/nbr_grnds, lessThanAvrg, moreThanvrg);
	for(i=0;i<=(int)maxrowval;i++)	if (covRes[i]>0) printf("cov[%d]=%f\t",i,covRes[i]);
	printf("\n");
	int activeuavs=0;

	FILE* fp;
	fp=fopen("activeuavs.csv","w");
	for(map<int,double>::iterator it=(*soln).begin(); it!=(*soln).end(); it++)
	{
		for (j=0;j<dim;j++)
		{
			// skip comma not needed after last dim value
			fprintf(fp,"%lf", net->uavs[it->first][j]);
			if(j==dim-1)	fprintf(fp,"\n");
			else	fprintf(fp,",");
		}
	}
	fclose(fp);

	// print in file connections uavs-ground nodes
	fp=fopen("uavs_grounds.csv","w");
	int activecovers[nbr_grnds];
	for(i=0;i<nbr_grnds;i++)
	{
		activecovers[i]=0;
		for(map<int,double>::iterator it=(*soln).begin(); it!=(*soln).end(); it++)
		{
			int uavk=it->first;
			if( euclDistance(net->uavs[uavk],grnds[i]) <= range )
			{
				assert(uav_in_cover(net->gcovs[i], uavk));// security safety, otherwise would be a problem
				fprintf(fp,"%lf,%lf\n", grnds[i][0], grnds[i][1]);
				fprintf(fp,"%lf,%lf\n\n", net->uavs[uavk][0], net->uavs[uavk][1]);
				activecovers[i]++;
			}
		}
	}

	printf(" Obj func : %f and nactive uavs %d\n",z, (*soln).size());
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
	
	// housekeeping
	delete covRes;

	return soln;
};


vector<double*>* onepassmethod(double** input_data, int nbr_grnds, double range)
{// first using one pass method, and then fine tuning with k-means

	int k=0,j=1;
	// required memory allocation for solution
	vector<double*>* res=new vector<double*>;
	vector<int> ninres;// gives the n elements in each cluster
	int gcovs[nbr_grnds];
	vector<double*> cl;// temporary centroids, contains sum of elements in cluster, to build clusters

	// Init first cluster with first element
	ninres.push_back(1);
	res->push_back(new double[dim]);
	for(k=0;k<dim;k++)
		(*res)[0][k]=input_data[0][k];

	cl.push_back(new double[dim]);// same with temporary value
	for(k=0;k<dim;k++)
		cl[0][k]=input_data[0][k];

	double largestdist=sqrt(pow(bound_1,2)+pow(bound_2,2));// Largest distance : extrem limits of the map
	double distance_to_closest_centroid=largestdist;
	double current_distance=0;

	// start 1e pass algorithm
	for(k=1; k<nbr_grnds; k++)
	{
		distance_to_closest_centroid=largestdist;
		for(j=0;j<res->size();j++)
		{
			current_distance=euclDistance(input_data[k],(*res)[j]);
			if ( current_distance < distance_to_closest_centroid )
			{
				distance_to_closest_centroid=current_distance;// keep current closest cluster
				gcovs[k]=j;
			}
		}

		if(distance_to_closest_centroid<range)
		{// closest centroid is within admissible range => add element to cluster
			ninres[gcovs[k]]++;
			//update centroids
			for(j=0;j<dim;j++)
			{
				cl[gcovs[k]][j]+=input_data[k][j];
				(*res)[gcovs[k]][j]=cl[gcovs[k]][j]/ninres[gcovs[k]];
			}
		}
		else
		{// closest centroid not within admissible range => Create new cluster
			gcovs[k]=res->size();
			ninres.push_back(1);
			res->push_back(new double[dim]);
			cl.push_back(new double[dim]);
			for(j=0;j<dim;j++){
				(*res)[res->size()-1][j]=input_data[k][j];
				cl[res->size()-1][j]=input_data[k][j];
			}
		}
	}
	
	// Housekeeping
	for(k=0;k<res->size();k++)
	{
		delete[] cl[k];
		cl[k]=nullptr;
	}

	return res;
};


double k_means(double** data, int ngrounds, vector<double*>* clusts, double error_tolerance, double range)
{

    int h, i, j; /* loop counters, of course :) */

    double old_error, error = DBL_MAX; // DBL_MAX : macro defines maximum finite floating-point value : 1E+37
	int ninres[clusts->size()]={};// gives the n elements in each cluster
	int gcovs[nbr_grnds];
	vector<double*> cl(clusts->size());// temporary centroids
	for(h=0;h<clusts->size();h++)
		cl[h]=new double[dim];

    assert(data && clusts->size() > 0 && error_tolerance >= 0); /* for debugging */

    // initialization : can start with any centroid(ordering changes the final result), here first

	double min_distance = DBL_MAX, distance = 0;

    // Kmeans loop
    do {
        // save error from last step
		old_error = error;
        error = 0;

        // clear old counts and temp centroids
        for (i=0; i <clusts->size(); i++)
        {
			// reinit temp centroids
            for (j=0;j<dim; cl[i][j++] = 0){};
			ninres[i]=0;
        }

        for (h = 0; h < ngrounds; h++)
        {
			gcovs[h] = 0;
            // Find closest cluster
            min_distance = DBL_MAX;
            for (i=0; i<clusts->size(); i++)
            {
				distance=0;
				// considered measure : squared residuals
                for(j=0;j<dim;j++)	distance+=pow(data[h][j]-(*clusts)[i][j],2);
                if (distance < min_distance)
                {
					min_distance=distance;// find the closest cluster
                    gcovs[h] = i;
                }
            }
			// update temp centroid of destination cluster
			for(j=0;j<dim;j++)	cl[gcovs[h]][j] += data[h][j];
			ninres[gcovs[h]]++;
            /* update standard error */
            error += min_distance;
        }

		// update centroids
        for (i=0; i<clusts->size(); i++)
            for (j=0; j<dim; j++)
                (*clusts)[i][j] = ninres[i] ? cl[i][j]/ninres[i] : cl[i][j];

    } while (fabs(error - old_error) > error_tolerance);/* if for each iteration, the number of changes made are not different from previous */

	// housekeeping
	for(h=0;h<clusts->size();h++)
	{
		delete[] cl[h];
		cl[h]=nullptr;
	}

	return error;
};


//void connect_CCs(sln* net, igraph_t* Gk, double range, vector<int*> &restr_list, int **pairs, int* npairs, bool restrict, int level)
void connect_CCs(sln* net, igraph_t* Gk, double range, vector<int>* uavsconnectivity, stack<tuple<int,int>>* pairs, bool dorestriction)
{// The boolean value "restrict" says wether the restriction strategy should be used (restrict is false, and G isn't root G0), or not.
	// In restriction not applied, then fill : **restr_list, **pairs, and * npairs. Otherwise, use these lists to apply the restrictions

printf("net n_uavs %d threshold %f and restrict %d\n", net->uavs.size(), range, dorestriction);
//printf(" vcount %li ecount %li\n", (long int)igraph_vcount(Gk), (long int)igraph_ecount(Gk));

	igraph_integer_t ncomps=-1;// numbers of connected components
	igraph_vector_t labels, compssizes;// labels of connected components for each vertex
	double min_distance=DBL_MAX;// distance of two closest nodes but out of range (in differents connected components)
	double current_distance;// used to place a new uav
	int buffn1=-1,buffn2=-1, n1=-1, n2=-1;//two closest but out of range nodes

	// Use of a brute force branching here. (Maybe undeterministic better?)
	// brute force : create one unique connected component by iteratively branching the 2 closest connected components
	// From here, a bit of a mess, pay attention to segmentation faults, clean ASAP
	int i=0, j=0, k=0, added=0;
	// Note : avoid doing : bool restricted=false; and restrict once only. If done so, it's possible that on the next phases restricted (n1,n2)
	// is included anyway and that we end up with a graph that was already generated. Then : always restrict (n1,n2) in new G.

	bool restricted=false;// restriction is performed once only
	do
	{
printf("n_uavs %d, k = %d, %d\n", net->uavs.size(), k++, ncomps);

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
		if(ncomps>1)
		{//add new node to reach one connected component
// printf("In ncomps %d connected components\n", ncomps);
			min_distance=DBL_MAX;
			// find two closest but out of range nodes
			for (i=0;i<net->uavs.size();i++)
			{
				for (j=i+1;j<net->uavs.size();j++)
				{
					// Not interested in nodes in same connected component
					if ( VECTOR(labels)[i] == VECTOR(labels)[j] )	continue;
					/* Find closest clusters but out of range */
					current_distance=euclDistance(net->uavs[i], net->uavs[j]);
					// check only if different connected component
					if(current_distance<min_distance)
					{// keep two nodes
						if( dorestriction && !restricted && get<0>(pairs->top()) == i && get<1>(pairs->top()) == j )
						{// restricted, skip
printf("i uav2 n1 n2 pairs[0] pairs[1] : %d %d %d %d %d %d %d\n", i, j, n1, n2, get<0>(pairs->top()), get<1>(pairs->top()));
printf("tests in \"restricted\", (pairs[0] == n1), (pairs[1] == n2) : %d %d %d\n", get<0>(pairs->top()) == n1, get<1>(pairs->top()) == n2);
							restricted=true;
							continue;
						}
printf("not restricted and n1 n2 ncomps : %d %d %d\n", n1, n2, ncomps);
						min_distance=current_distance;
						n1=i;
						n2=j;
					}
				}
			}
long int zz=igraph_vcount(Gk);
//for(long int z=0;z<zz;z++)	cout<<VECTOR(labels)[z]<<" ";
//cout<<endl;

			// free vector of labels
			igraph_vector_destroy(&labels);

			if(n1<0 || n2<0)
			{
				printf(" Something went wrong with graph, check since adding new uav will fail %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);
				return;
			}
cerr<< "use ("<< net->uavs[n1][0] << "," << net->uavs[n1][1]<< ") and ("<< net->uavs[n2][0] << "," << net->uavs[n2][1]<<")"<<endl;
			// assign coordinates of new created node : middle of segment [n1,n2]
			double* buffdouble=new double[dim];
			for (j=0;j<dim;j++)
				buffdouble[j]=(net->uavs[n1][j]+net->uavs[n2][j])/2;
			addTonetwork(net, buffdouble, range);// update distance of network of uavs

			uavsconnectivity->push_back(net->uavs.size()-1);// keep track of uavs used for connectivity

			// add new vertex and edges to graph
			igraph_add_vertices(Gk, 1, 0);			

cerr<< "Tests |V(G)| = " << (long int)igraph_vcount(Gk) << ", net->uavs.size() = " << net->uavs.size() << " |E(G)|"<< (long int)igraph_ecount(Gk)<<endl;

			int lastuav=net->uavs.size()-1;
			for (i=0;i<net->uavs.size();i++)
			{
				for (j=0;j<net->outdeg[i].size();j++){
					if(net->outdeg[i][j]==lastuav){// add edges with new uav to graph
//cout<< "A LEAST THIS DOES WORK " <<i<<"-"<<j<< endl;
						igraph_add_edge(Gk, i, lastuav);
					}
				}
			}

			if( ! dorestriction )
			{// build first connected sln, keep track of n1 and n2 for further "potential" restrictions
				pairs->push(make_tuple(n1,n2));
			}

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

//	}while(ncomps>2 && k++ < 1);
	}while(ncomps>1);

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


//void populate(sln* covers, map<int,double>* uavscoverage, vector<int>* uavsconnectivity, double range, igraph_t* solnG0, stack<tuple<int,int>>* pairs)
void populate(sln* covers, vector<int>* uavsconnectivity, double range, igraph_t* solnG0, stack<tuple<int,int>>* pairs)
{// note : computing diameter of graph (O(n*|E|)) can cost more than comparing each nodes (O(n^2)) : eg : Complete graph : n*|E| = n*sum_i_in{1...n-1}(i) > O(n^2) for n>3

	long int i=0, j=0;
//	long int k=0;

	/* the goal of the step is to generate the best NINDIVS by unpiling (Start with last) the list of edges added to the first
	 * graph (solnsGraphs[0]) in order to create the unique connected component. The restriction on the current level of edges
	 * to unpile is given by k (if can't use c++ stl's 'stack'). The number of vertices to start with on each iteration is n_in_first_graph-k-1
	 * where n_in_first_graph is the number of vertices in the first connected graph. */
	// long int n_in_first_graph=(long int)igraph_vcount(solnsGraphs[0]);

	// Use of a brute force branching here. Maybe undeterministic could work better.

	// 1st phase : build the 1st connected graph with igraph specs (root). Note : probably not unique connected graph
	solnG0=new igraph_t;// first connected graph to build.


	// start building graph
	if(igraph_empty(solnG0, covers->uavs.size(), IGRAPH_UNDIRECTED)==IGRAPH_EINVAL)
		printf(" Invalid number of vertices, graph init failed %s, %d, %s \n", __FILE__, __LINE__, __FUNCTION__);

printf("G0 with %d active uavs\n", covers->uavs.size());
	// add edges to graph : two uavs in each other range
/*
int count=0;
	for(map<int,double>::iterator it=uavscoverage->begin(); it!=ab.end(); it++)
		cout << it->first << " --- " << it->second << endl;
*/
int count=0;
	for(i=0;i<covers->uavs.size();i++)
	{
		for(j=0;j<covers->outdeg[i].size();j++)
		{// for each uav, if connected to another one, then create link in graph
			count++;
			igraph_add_edge(solnG0, i, covers->outdeg[i][j]);// since j = i+1 to nActivUAVs, then no loop
		}
	}

printf("\nend 1st graph and count %d\n",count);
printf("Before connectCCs vcount %li ecount %li and n_uavs %d\n", (long int)igraph_vcount(solnG0), (long int)igraph_ecount(solnG0), covers->uavs.size());

	connect_CCs(covers, solnG0, range, uavsconnectivity, pairs, false);

printf("Test npairs == %d and should be one connected component (passed connectCCs function)\n", pairs->size(), covers->uavs.size());

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

	// On each iteration, record all data on files to save space
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
		for(i=0;i<covers->uavs.size();i++)
			for (j=0;j<dim;j++)
			{
				// skip comma not needed after last dim value
				if(j==dim-1)	fprintf(fout,"%lf\n", covers->uavs[i][j]);
				else fprintf(fout,"%lf,", covers->uavs[i][j]);
			}
		fclose(fout);

// ________________________________________________________________
/* uncomment here for populating steps
// ________________________________________________________________

TO ERASE : Done in addTonetwork		
		find_covers(covers, threshold);
		
int overallredundancy=0;
// Compute redundancy of G_0
for(i=1;i<=nbr_grnds;i++)	overallredundancy+=covers->g_covers[i];
printf("G_0 has overall redundancy of %d\n",overallredundancy);

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

		updateDistMat(buffnets, range);

		Gk=(igraph_t*)malloc(sizeof(igraph_t));

		if(igraph_empty(Gk, buffnets->n_uavs+1, IGRAPH_UNDIRECTED)==IGRAPH_EINVAL)
			printf("Ite %d Invalid number of vertices, graph initialisation failed %s, %d, %s \n", (int)k, __FILE__, __LINE__, __FUNCTION__);

//printf("Phase 1 : Edges found while comparing distances\n");
//\/\* ------------------------------- Start independent comments
//count=0;
//		for(i=1;i<=buffnets->n_uavs;i++)
//			for(j=i+1;j<=buffnets->n_uavs;j++)
//				if(buffnets->dists[i][j]<threshold)
//				{
//printf("(%ld -- %ld)\t",i,j);
//count++;
//					igraph_add_edge(Gk, i, j);// since j = i+1 to nActivUAVs, then no loop
//				}
//\*\/ ------------------------------- End independent comments
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

// ________________________________________________________________
/* END uncomment here for populating steps
// ________________________________________________________________



// ------------------------------------------------------    note to myself : start Independent comments
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
// ------------------------------------------------------    Note to myself: end independent comments

};



int main(int argc, char** argv)
{

clock_t begin = clock();

	/* !!! igraph : turn on attribute handling  Ie. even if you don't manipulate attributes explicitly, but create a graph that might have some attributes, eg. read a graph a GraphML file, you need this call before, otherwise the attributes are dropped. */
	igraph_i_set_attribute_table(&igraph_cattribute_table);

	readData(argv);

	//bound_1=1000;
	//bound_2=1000;

//	double radius=(uavs_range/2);
	double range=uavs_range;

//	sln *res=onepassmethod(grnds, nbr_grnds, radius);
//	sln *trueres=onepassmethod(grnds, nbr_grnds, radius);
	//translate(res, radius);

	vector<double*>* res;// true result
	vector<double*>* res_plus_1;// true result + 1
	vector<double*>* res_plus_2;// true result + 2 and current result
	vector<double*>* del;// pointer to result to be freed

	/* Elbow criteria : if 2 next wss are close : deviation(res+1,res+2) > 3/4*(deviation(res,res+1)) */
	int i=1,j=0;

	// to generate the covers
	double prev_deviation=0,wss_minus_1=0,wss=0;/* Caution : wss isn't one of the true result but result+2 */
	double elbow_ratio=0.7;// if new deviation less than three quarters of previous deviation
	// !!! Personal note : if ratio is too low, about 0.001, then there are issues, check if time.!!!!

	bool stop=false;

printf("nbr grnds : %d\n",nbr_grnds);

	do{
		
		// do until elbow criteria satisfied:
		res_plus_2=onepassmethod(grnds, nbr_grnds, range/i);
		wss=k_means(grnds, nbr_grnds, res_plus_2, 0.0001, range/i);
printf("i %d, wss. %f, n_uavs %d\n", i, wss, res_plus_2->size());
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
	, i, wss_minus_1, wss, prev_deviation, elbow_ratio*prev_deviation, wss_minus_1-wss, wss_minus_1-wss < elbow_ratio*prev_deviation, res_plus_1->size(), res->size());

			// housekeeping
			for(j=0;j<res_plus_1->size();j++)
			{
				delete[] (*res_plus_1)[j];
				(*res_plus_1)[j]=nullptr;
			}
			delete res_plus_1;// no longer needed, as only true result "res" is returned
			res_plus_1=nullptr;
			for(j=0;j<res_plus_2->size();j++)
			{
				delete[] (*res_plus_2)[j];
				(*res_plus_2)[j]=nullptr;
			}
			delete res_plus_2;// same
			res_plus_2=nullptr;
			stop=true;
		}
		else{/* keep reducing the radius */

			prev_deviation=wss_minus_1-wss;
			wss_minus_1=wss;
			del=res;
			res=res_plus_1;// store result
			res_plus_1=res_plus_2;
			// housekeeping
			for(j=0;j<del->size();j++)
			{
				delete[] (*del)[j];
				(*del)[j]=nullptr;
			}
			delete del;//Housekeeping
			del=nullptr;
		}
	}while(!stop && wss>0);/// !!!!!!!!!!!!!!!!!!!!! CAREFUL!!! CHECK STOPPING CONDITION (wss>0???)
//	}while(i<1);

/*
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

	sln* rawsln=new sln;// create first raw solution
	rawsln->gcovs=new vector<int>[nbr_grnds];
	
	double* buffdouble;
	for(i=0;i<res->size();i++){
		buffdouble=new double[dim];
		for(j=0;j<dim;j++)	buffdouble[j]=(*res)[i][j];
		addTonetwork(rawsln, buffdouble, range);
	}

	// housekeeping
	for(j=0;j<res->size();j++)
	{
		delete[] (*res)[j];
		(*res)[j]=nullptr;
	}
	delete res;//Housekeeping
	res=nullptr;

	double lb=2.0;


	// linear solver : return indices of uavs active for coverage, and the results of their linear program values
	map<int,double>* uavscoverage=solve_linear_model(rawsln, range, lb);

	// copy solution into new and cleaner one, and at the same time store in file the coordinates of active uavs for coverage
	sln* net=new sln;
	net->gcovs=new vector<int>[nbr_grnds];

//	FILE* fp;
//	fp=fopen("finalres.csv","w");
	for(map<int,double>::iterator it=(*uavscoverage).begin(); it!=(*uavscoverage).end(); it++)
	{
		buffdouble=new double[dim];
		for (j=0;j<dim;j++)
		{
			// skip comma not needed after last dim value
//			if(j==dim-1)	fprintf(fp,"%lf\n", rawsln->uavs[it->first][j]);
//			else fprintf(fp,"%lf,", rawsln->uavs[it->first][j]);
			
			buffdouble[j]=rawsln->uavs[it->first][j];
		}
		addTonetwork(net, buffdouble, range);// update distance of network of uavs
	}
//	fclose(fp);
	
	delete rawsln;
	rawsln=nullptr;

	stack<tuple<int,int>>* pairs=new stack<tuple<int,int>>;// used for restrictions
	vector<int> uavsccs;// indices of uavs used to link sparse connected components, empty at start
	igraph_t* solG0=nullptr;

//	populate(net, uavscoverage, &uavsccs, range, solG0, pairs);
	populate(net, &uavsccs, range, solG0, pairs);


	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("|Time spent %f\n",time_spent);

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
