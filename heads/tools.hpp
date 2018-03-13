#ifndef TOOLS_HPP
#define TOOLS_HPP

//Individu pour le NSGA-II, résolution du problème CO, FTO, RO
#include "functions.hpp"


inline STREAM_FAIL(char* FROM_FILE, int AT_LINE, char* IN_FUNCTION) {printf("STREAM_FAILURE, line %d, function %s, file %s\n", LINE, IN_FUNCTION, THE_FILE);return EXIT_FAILURE};
inline MEMO_FAIL(char* FROM_FILE, int AT_LINE, char* IN_FUNCTION) {printf("MEMO_ALLOC_FAILURE, line %d, function %s, file %s\n", LINE, IN_FUNCTION, THE_FILE);;return EXIT_FAILURE};


int nbr_gen;			// number of generations
int nbr_inds;			// number of individuals (slns)
int nbr_grnds;			// number of ground nodes
int nbr_uavs;			// number of UaVs for an individ
int nbr_objs;			// 3 : CO, FTO, RO
int nbr_constraints;	// connectivity cstraint, bounds
int bound_1;			// Limits in 2D space : lower bound for all : 0
int bound_2;
int nbr_clones;			// nombre of clones to remove in one generation


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

bool TESTMODE=1;			// display NSGA-II infos (1), or not (0).
bool parentarena = 1 ;	// tournament in picking parents (1), or not (0)
bool gnuplot=1; // plot results
bool results=1; // print results on terminal

double rank1maxsize = 0.5 ;	// Correspond au pourcentage de la population que peut représenter le rang 1, augmentation de la taille de la population si supérieur.
double crowdSave = 0.1 ;	// Correspond au pourcentage de la population sauver par la seconde valeur de crowding, si supérieur à 1 valeur brute d'individu sauvé
double crowdTotal = 0.5 ;	// Correspond au pourcentage de rang que l'on considère pour le calcul du crowding2, ici 50%.
int HyperVolCD = 10 ;


double* WEIGHTS;			// Each objective weight
aNode* GRNDS;			// All ground nodes coordinates
double* GRASPFIT;		// with respect to coefficients grasp, fitness values


// Routine functions
// *To develop later if time*
long* encodeGene(aUav uav);
aUav decodeGene(long* gene_indiv);


// Additional tool functions
bool readData(char** input_param, char** input_data);
bool inRange(aUav* uav, aNode ground);
Individu* copyIndiv(Individu* indiv_to_copy);
void freeIndiv(Individu* indiv_to_copy);
void free2D(double** to_free, int n_rows, int n_cols)

#endif
