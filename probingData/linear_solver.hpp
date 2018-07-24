#ifndef LIN_SOLV_HPP
#define LIN_SOLV_HPP


//#include "tools.hpp"


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


/* for linear problem */
// Grounds nodes are given as global variables
void solve_linear_model(sln* net, double range, double lb, double* soln);

#endif
