#ifndef TOOLS_HPP
#define TOOLS_HPP


#include <igraph.h>

#include <vector>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>
#include <math.h>


inline
void STREAM_FAIL(char* FROM_FILE, int AT_LINE, char* IN_FUNCTION)
{printf("STREAM_FAILURE, line %d, function %s, file %s\n", AT_LINE, IN_FUNCTION, FROM_FILE);return EXIT_FAILURE;};

inline
void MEMO_FAIL(char* FROM_FILE, int AT_LINE, char* IN_FUNCTION)
{printf("MEMO_ALLOC_FAILURE, line %d, function %s, file %s\n", AT_LINE, IN_FUNCTION, FROM_FILE);return EXIT_FAILURE;};



typedef struct sln{
	int* labels;// Size : the total number of ground nodes, gives for each element the cluster it belongs to.
	double** uavs;// contains the coordinates of uavs[i] in the solution
	double** distances;// squared matrix of distances between uavs
	int* counters;// number of elements (size : n_uavs) covered by the uav
	int n_uavs;// number of nodes in network : either nodes are uavs (sln : set of uavs), or ground nodes (uav : set of ground nodes)
}sln;


/** Variables
 */
// general variables
int nbr_subpop;			// number of subpopulations using different strategies
int nbr_gen;			// number of generations
int nbr_inds;			// number of individuals (slns)
int nbr_grnds;			// number of ground nodes
int max_uavs_avail;		// the total number of uavs available
int p_objctvs;			// 3 : CO, FTO, RO
int m_cnstrnts;			// connectivity cstraint, bounds
double bound_1;			// Limits in 2D space : lower bound for all : 0
double bound_2;
int nbr_clones;			// nombre of clones to remove in one generation

char *opt_path[100];// path of optimal solutions

// stopping criteria and evolutionary parameters
int max_time; // limit on the processing time
int no_improv_ites; // for "no_improv_ites" iterations where the difference in successive fitnesses (improvements) is less than "epsilon", then stop
double epsilon;
double mutaprob;
double evo_part_percent;
double grasp_part_percent;
double grasp_split;// Correspond au pourcentage de population grasp
							// => PopSize*k ici 50% si supérieur ou égale à 1 valeur brute d'individu par direction
int n_delta;// number of grasp directions
double alpha;// RCL restricting coefficient
bool parentarena = 1 ;	// tournament in picking parents (1), or not (0)
double rank1maxsize = 0.5 ;	// Correspond au pourcentage de la population que peut représenter le rang 1, augmentation de la taille de la population si supérieur.
double crowdSave = 0.1 ;	// Correspond au pourcentage de la population sauver par la seconde valeur de crowding, si supérieur à 1 valeur brute d'individu sauvé
double crowdTotal = 0.5 ;	// Correspond au pourcentage de rang que l'on considère pour le calcul du crowding2, ici 50%.
int HyperVolCD = 10 ;
double* graspfit;	// with respect to coefficients grasp, fitness values

// Data visualisation choices
bool TESTMODE=1; // display NSGA-II infos (1), or not (0).
bool gnuplot=1; // plot results
bool results=1; // print results on terminal
bool removeclone=1;	// 1, do not remove generated clones, 0 otherwise

// Set of ground nodes
double** grnds;		// All ground nodes coordinates
//double* uavs_ranges;	// all ranges of available uavs
double uavs_range;	// all ranges of available uavs
/** Functions
 */

// *To develop later if possible*
//long* encodeGene(aUav uav);
//aUav decodeGene(long* gene_indiv);

bool readData(char** argv);
void writeData(sln* a_sln);
double euclDistance(double *node1, double *node2);
bool inRange(double* node1, double* node2);// work both for comparing : uav-uav, uav-ground node
void updateDistMat(sln* net);

// clustering methods
sln* method1ePasse(double** input_data, int size_input, double threshold);
void k_means(double** data, int n, sln* clusts, double error_tolerance);

sln* addUaV(sln to_sln, double* coords);
sln* removeUaV(sln from_sln, int index_uav);
sln* copySln(sln* sln_to_copy);
void freeSln(sln* sln_to_free);
igraph_t* translate(sln* net);// Builds the graph matching the solution
void freeNetwork(igraph_t* net);// Deallocate memory of network

// Linked-list operations main reference : https://www.geeksforgeeks.org/data-structures/linked-list/

void freeIndiv(Individu* indiv_to_copy);
void free2D(double** to_free, int n_rows, int n_cols);




#endif
