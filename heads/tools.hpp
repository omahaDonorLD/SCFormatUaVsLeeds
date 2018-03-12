#ifndef TOOLS_HPP
#define TOOLS_HPP

//Individu pour le NSGA-II, résolution du problème CO, FTO, RO
#include "functions.hpp"
#define STREAM_FAIL(char* FROM_FILE, int AT_LINE, char* IN_FUNCTION) (printf("STREAM_FAILURE, line %d, function %s, file %s\n", LINE, IN_FUNCTION, THE_FILE);)
#define MEMO_FAIL(char* FROM_FILE, int AT_LINE, char* IN_FUNCTION) (printf("MEMO_ALLOC_FAILURE, line %d, function %s, file %s\n", LINE, IN_FUNCTION, THE_FILE);)

int nbr_gen;			// number of generations
int nbr_grnds;			// number of ground nodes
int nbr_uavs;			// number of UaVs for an individ
int nbr_objs;			// 3 : CO, FTO, RO
int nbr_constraints;	// connectivity cstraint, bounds
int bound_1;			// Limits in 2D space : lower bound for all : 0
int bound_2;

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
