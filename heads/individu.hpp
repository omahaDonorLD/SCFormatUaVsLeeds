#ifndef INDIVIDU_HPP
#define INDIVIDU_HPP

//Individu pour le NSGA-II, résolution du problème CO, FTO, RO
#include "functions.hpp"

int nbr_inds;					// Number of individuals composing a solution
int nbr_uavs;					// number of UaVs for an individ
int nbr_objs;					// 3 : CO, FTO, RO
int nbr_constraints;	// connectivity cstraint, bounds
int bound_1;					// Limits in 2D space
int bound_2;

double* WEIGHTS;			// Each objective weight
aNode* GRNDS;			// All ground nodes coordinates
double* GRASPFIT;		// with respect to coefficients grasp, fitness values

typedef struct Individu {
		double* profit;		// Fitness, based on 3 profit value : CO (Coverage), FTO (Fault tole), RO (Redund)
											// profits[nbr_uavs]== Total of the solution

		//vector<int> current_objective_value;	// valeur totale actuelle des objectifs
		//vector<int> current_constraint_value;	// valeur totale actuelle des contraintes
		//int getTotalWeight(int dimension);	// retourne la valeur de la contrainte "dimension"

		//void addObject(int grnd_node);		// rajoute l'item "num_item" au sac à dos si absent
		//void removeObject(int grnd_node);	// supprime l'item "num_item" du sac à dos si présent

		int rank; 							// rang de l'individu

		//void computeObjectives();			// calcule et met à jour les valeurs des objectifs de l'individu
		//void computeConstraints();		// calcule et met à jour les valeurs des contraintes de l'individu



		void restart();		// mutating an individual from 1 to 0
		void restart2();	// mutation by switching to 1 randomly
		void restart3();	// mutation with 1-1 switch


	public: // Everything necessary for NSGA-II


		Individu();
		Individu(Individu* toCopy);			// copy from a different pointer
		~Individu();

		double crowding_value;				// crowding1 value
		double crowding_Total;				// crowding2 value
		bool isGRASP;

		static int getNbrObjectives();
		int getRank();						// returns rank of individual

		void setRank(int in_rank);			// set rank "in_rank" to individual

		static void initProblem(string instance_path);
		static Individu* getChildFrom(Individu* parent1, Individu* parent2);	// Using operator

		int getObjectiveValue(int obj);
		int getContraintValue(int dimension);

		bool dominates(Individu* Indiv2);
		bool isCloneOf(Individu* Indiv2);

		void randomize();					// create randomly an individual
		//static vector<int> RCList;			// RCList correspondant à GRASPFIT pour un alpha donné
		//static void createGRASPFIT(int* Coef, double alpha);	// Permet de créer GRASPFIT et RCList à partir d'un tableau de Coefficient de direction "Coef" et de "alpha"
		void GRASP( double alpha);			// create GRASP individual with tradeoff greedy-ranfom the coefficient "alpha"

		void mutate();		// mutation, restart3


		string toString();	// translate gene of a solution


		void printToScreen();
		static void testData();
};

typedef struct cstrt_info{
	int n_violtd;
	int *index_cnstrts_violtd;
}cstrt_info;



Individu* closestIndiv(Individu* an_indiv, double radius);		// Part of making the covering constraint satisfiable. Finds the closest individuals close to a given one (picked randomly) and with respect to a given radius.

// constraints
bool isFeasible(Individu* aSolution);			// check if a solution is feasible
void makeFeasible(Individu* aSolution);			// fix the violated constraints of the solution
cstrt_info cnstrtVioltd(Individu* an_indiv);	// check the constraints that the individual violates
void makeConnected(Individu* aSolution);		// fix the connectivity constraint of the solution
void coverAllGrnds(Individu* aSolution);		// covers all the grounds with the minimum number of uavs

// objectives
void updateFitness(Individu* aSolution);		// computes the fitness of the solution

#endif
