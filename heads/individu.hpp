#ifndef INDIVIDU_HPP
#define INDIVIDU_HPP

//Individu pour le NSGA-II, résolution du problème sac à dos bi-objectif, multi-dimensionnel/


#include <vector>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>

using namespace std;

// Each node have coordinates
struct aNode {
	bool ground=true;	// Is aNode a ground node (True) or UaV
	float x=0.;
	float y=0.;
};


class Individu
{
	private:// Concerning the problem

		static int nbr_uavs;			// number of UaVs for an individ
		static int nbr_objs;			// 3 : CO, FTO, RO
		static int nbr_constraints;		// connectivity cstraint, bounds
		static int bound_1;				// Limits in 2D space
		static int bound_2;

		static int* range;  			// each UaV's range
		static float* weights;			// Each objective weight
		static aNode* grnds;			// All ground nodes coordinates

		int* numbrsGrnd;		// for each UaV, number of ground nodes covered
		float* profits;		// Fitness, based on 3 profit value : CO (Coverage), FTO (Fault tole), RO (Redund)
										// profits[nbr_objs]== Total on all
		
		bool* active_uavs;				// 0 : uav not needed - 1 : active
		//vector<int> current_objective_value;	// valeur totale actuelle des objectifs
		//vector<int> current_constraint_value;	// valeur totale actuelle des contraintes
		//int getTotalWeight(int dimension);	// retourne la valeur de la contrainte "dimension"
		
		bool constraintRespected(int dimension);// Check if constraint "dimension" is satisfied
		
		bool isFeasible();					// check if solution is feasible
		void makeFeasible();
		
		void addObject(int grnd_node);		// rajoute l'item "num_item" au sac à dos si absent
		void removeObject(int grnd_node);	// supprime l'item "num_item" du sac à dos si présent
		
		int rank; 							// rang de l'individu
		

		
		void computeObjectives();			// calcule et met à jour les valeurs des objectifs de l'individu
		void computeConstraints();			// calcule et met à jour les valeurs des contraintes de l'individu



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
		//static vector<double> GRASPFIT;		// vecteur de valeur des items selon les coefficient grasp
		//static vector<int> RCList;			// RCList correspondant à GRASPFIT pour un alpha donné
		//static void createGRASPFIT(int* Coef, double alpha);	// Permet de créer GRASPFIT et RCList à partir d'un tableau de Coefficient de direction "Coef" et de "alpha"
		void GRASP( double alpha);			// create GRASP individual with tradeoff greedy-ranfom the coefficient "alpha"

		void mutate();		// mutation, restart3

		
		string toString();	// translate gene of a solution


		void printToScreen();	
		static void testData();	
};
#endif
