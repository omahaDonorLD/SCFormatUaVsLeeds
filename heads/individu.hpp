#ifndef INDIVIDU_HPP
#define INDIVIDU_HPP

//Individu pour le NSGA-II, résolution du problème CO, FTO, RO
#include "functions.hpp"


/** Structures
 */
typedef struct Solution {
		//long** gene;		// *to develop later if time* of size "nbr_uavs"*gene. Gene : sequence of 0 and 1 (int) for each UaV
		int nbr_uavs;		// number of uavs in the solution
		aUav **uavs;		// the uavs in the solution
		igraph_t *network;	// net of uavs, a graph of type "igraph_t" for checking FTO
		int* ranks;			// ranks of each uav
		double* profit;		// Fitness, based on 3 profit value : CO (Coverage), FTO (Fault tole), RO (Redund)
							// last index is total on 3 objectives : profits[nbr_uavs]== Total of the solution
		int* cnstrts_viol; 		// each constraint becomes 1 if violated (minimize). Only valid if sum(cnstrts)==0
		// Necessary to NSGA-II
		double crowding_value;				// crowding1 value
		double crowding_Total;				// crowding2 value
		bool isGRASP;
}sln;


/** Variables
 */
sln** RCList;			// RCList corresponfing to graspfit and alphas


/** Functions
 */
// objectives
void updateFitness(sln* aSolution, int uav_updated);// computes and assign the fitness of the new solution (copy first the unmodified fitness)

// constraints
bool isFeasible(sln* aSolution);			// check if solution is feasible
void makeFeasible(sln* aSolution);			// fix the violated constraints of the solution
void makeConnected(sln* aSolution);		// fix the connectivity constraint of the solution
void coverAllgrnds(sln* aSolution);		// covers all the grounds with the minimum number of uavs

//sln* closestIndiv(aNode* an_indiv, double radius);		// Part of making the covering constraint satisfiable. Finds the closest individuals close to a given one (picked randomly) and with respect to a given radius.

// Init solution. 3 strategies : 1. random, 2. k-means, 3. equally spread
sln* randomize();	// creates randomly an individual contained in search space
sln* clustering();	// creates individuals by clustering (kmeans) the ground nodes in search space ()
sln* uniform();	// creates individuals that are equally spread on the search space : Width/(2*radius-constant). With "cst" so that pairs of UaVs share a space equal to that constant (kind of ;) )

// evolutionary operations
void mutate(sln* an_indiv, int index_uav, double delta);	// move uav in solution wrt delta mutating an individual from 1 to 0
sln* cross(Individu* indiv_1, Individu* indiv_2);	// cross-over only creates 1 child
sln* creategraspfit(double* Coef, double alpha);	// Creates graspfit and RCList from direction Coefficients and alphas
sln* GRASP( double alpha);			// create GRASP individual with tradeoff greedy-random depending on the coefficient "alpha"


// MO tools
bool dominates(Individu* indiv_1, Individu* indiv_2); // Does indiv 1 dominates indiv 2


//void free2D(double** to_free, int n_rows, int n_cols)

#endif
