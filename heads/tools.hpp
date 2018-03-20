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
*/


/** Variables
 */
// general variables
int nbr_subpop;			// number of subpopulations using different strategies
int nbr_gen;			// number of generations
int nbr_inds;			// number of individuals (slns)
int nbr_grnds;			// number of ground nodes
int nbr_uavs;			// number of UaVs for an individ
int max_uav_avail;		// the total number of uavs available
int nbr_objs;			// 3 : CO, FTO, RO
int nbr_constraints;	// connectivity cstraint, bounds
int bound_1;			// Limits in 2D space : lower bound for all : 0
int bound_2;
int nbr_clones;			// nombre of clones to remove in one generation

// stopping criteria and evolutionary parameters
int max_time; // limit on the processing time
int no_improv_ites; // for "no_improv_ites" iterations where the fitness difference in fitness (improvements) is less than "epsilon", then stop
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
double* GRASPFIT;	// with respect to coefficients grasp, fitness values

// Data visualisation choices
bool TESTMODE=1; // display NSGA-II infos (1), or not (0).
bool gnuplot=1; // plot results
bool results=1; // print results on terminal
bool removeclone=1;	// 1, do not remove generated clones, 0 otherwise

// Set of ground nodes 
aNode* GRNDS;		// All ground nodes coordinates
double* UAVs_Range;	// all ranges of available uavs

/** Functions
 */

// *To develop later if possible*
//long* encodeGene(aUav uav);
//aUav decodeGene(long* gene_indiv);

bool readData(char** argv);
double euclDistance(aNode *uav1, aNode *uav2);
bool inRange(aUav* uav, aNode* ground);
void writeData(sln* a_sln);
void freeIndiv(Individu* indiv_to_copy);
void free2D(double** to_free, int n_rows, int n_cols)


#endif
