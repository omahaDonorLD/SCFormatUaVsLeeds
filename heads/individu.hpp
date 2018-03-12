#ifndef INDIVIDU_HPP
#define INDIVIDU_HPP

//Individu pour le NSGA-II, résolution du problème CO, FTO, RO
#include "functions.hpp"


Individu** RCList;			// RCList correspondant à GRASPFIT pour un alpha donné

typedef struct Solution {
		//long** gene;		// *to develop later if time* of size "nbr_uavs"*gene. Gene : sequence of 0 and 1 (int) for each UaV
		int nbr_inds;
		aUav* uavs;			// number of in the solution
		int* rank;
		double* profit;		// Fitness, based on 3 profit value : CO (Coverage), FTO (Fault tole), RO (Redund)
							// last index is total on 3 objectives : profits[nbr_uavs]== Total of the solution
		int* cnstrts_viol; 		// each constraint becomes 1 if violated (minimize). Only valid if sum(cnstrts)==0

		// Necessary to NSGA-II
		double crowding_value;				// crowding1 value
		double crowding_Total;				// crowding2 value
		bool isGRASP;
}sln;


//sln* closestIndiv(aNode* an_indiv, double radius);		// Part of making the covering constraint satisfiable. Finds the closest individuals close to a given one (picked randomly) and with respect to a given radius.


// constraints
bool isFeasible(sln* aSolution);			// check if solution is feasible
void makeFeasible(sln* aSolution);			// fix the violated constraints of the solution
void makeConnected(sln* aSolution);		// fix the connectivity constraint of the solution
void coverAllGrnds(sln* aSolution);		// covers all the grounds with the minimum number of uavs


// objectives
void updateFitness(sln* aSolution, int uav_updated);// computes and assign the fitness of the new solution (copy first the unmodified fitness)

// Init solution. 3 strategies : 1. random, 2. k-means, 3. equally spread
sln* randomize(int bound_1, int bound_2);	// creates randomly an individual
sln* clustering(aNode* GRNDS, int bound_1, int bound_2);	// creates individuals by clustering (kmeans) the ground nodes in search space ()
sln* uniform(int bound_1, int bound_2);	// creates individuals that are equally spread on the search space : Width/(2*radius-constant). With "cst" so that pairs of UaVs share a space equal to that constant (kind of ;) )

// evolutionary operations
void mutate(sln* an_indiv, int index_uav, double delta);	// move uav in solution wrt delta mutating an individual from 1 to 0
sln* cross(Individu* indiv_1, Individu* indiv_2);	// cross-over only creates 1 child
sln* createGRASPFIT(double* Coef, double alpha);	// Creates GRASPFIT and RCList from direction Coefficients and alphas
sln* GRASP( double alpha);			// create GRASP individual with tradeoff greedy-random depending on the coefficient "alpha"


// MO tools
bool dominates(Individu* indiv_1, Individu* indiv_2); // Does indiv 1 dominates indiv 2


// Additional tool functions
bool inRange(aUav* uav, aNode ground);
sln* addUaV(sln sln_to_extend, aUav* uav_to_activate);
sln* removeUaV(sln sln_to_extend, aUav* uav_to_deactivate);
sln* copyIndiv(sln* sln_to_copy);
void freeSln(sln* sln_to_free);
void free2D(double** to_free, int n_rows, int n_cols)

#endif
