#ifndef INDIVIDU_HPP
#define INDIVIDU_HPP

//Individu pour le NSGA-II, résolution du problème sac à dos bi-objectif, multi-dimensionnel/


#include <vector>
#include <stdlib.h>
#include <time.h>
#include <string>
#include<iostream>
#include <fstream>
#include <sstream>
#include <climits>

using namespace std;

class Individu
{
	private: //tout ce qui concerne uniquement le problème

		static int nbr_items;				// nombre d'item du sac à doss
		static int nbr_constraints;			// nombre de contraintes
		static int nbr_objectives;			// nombre d'objectif

		static  vector<int> capacities;  	// valeurs max des contraintes
		static  vector<int> weights;		// valeurs actuelles des contraintes par item
		static  vector<int> profits;		// valeurs actuelles des objectifs par item
	
		static int getWeight(int num_item, int dimension);		// retourne le poids d'un item
		static int getProfit(int num_item, int num_objective);	// retourne le profit d'un item
	
		
		vector<bool> picked_objects;			// 0 item non pris - 1 item pris
		vector<int> current_objective_value;	// valeur totale actuelle des objectifs
		vector<int> current_constraint_value;	// valeur totale actuelle des contraintes
	
	
		int getTotalWeight(int dimension);	// retourne la valeur de la contrainte "dimension"
		
		bool constraintRespected(int dimension);				// vérifie si la contrainte "dimension" est respectée
		
		bool isFeasible();					// vérifie si l'individu est viable
		void makeFeasible();				// rend l'individu viable
		
		void addObject(int num_item);		// rajoute l'item "num_item" au sac à dos si absent
		void removeObject(int num_item);	// supprime l'item "num_item" du sac à dos si présent
		
		int rank; 							// rang de l'individu
		

		
		void computeObjectives();			// calcule et met à jour les valeurs des objectifs de l'individu
		void computeConstraints();			// calcule et met à jour les valeurs des contraintes de l'individu



		void restart();		// fait muter un individu en enlevant un 1
		void restart2();	// fait muter un individu en rajoutant tous les 1 possibles de façon aléatoire
		void restart3();	// fait muter un individu en faisant un échange 1-1
		
		

	public: //Tout ce qui est nécessaire à NSGA-II
		
		
		Individu();
		Individu(Individu* toCopy);			// copie avec pointeur différent de l'individu toCopy
		~Individu();
		
		double crowding_value;				// valeur de crowding1
		double crowding_Total;				// valeur de crowding2
		bool isGRASP;
		
		static int getNbrObjectives();		// rend le nombre d'objectif du problème
		int getRank();						// rend le rang de l'individu
		
		void setRank(int in_rank);			// donne le rang "in_rank" à l'individu

		static void initProblem(string instance_path); 							// Initialise le problème
		static Individu* getChildFrom(Individu* parent1, Individu* parent2);	// Permet d'obtenir un enfant à partir de deux parents.

		int getObjectiveValue(int num_objective);	// rend la valeur actuelle de l'objectif "num_objective"
		int getContraintValue(int num_contraint);	// rend la valeur actuelle de la contrainte "num_contraint"

		bool domine(Individu* inIndividu);			// vérifie si l'individu domine "inIndividu"
		bool isCloneOf(Individu* inIndividu);		// vérifie si l'individu est un clone de "inIndividu"

		void randomize();					// crée aléatoirement un individu
		static vector<double> GRASPFIT;		// vecteur de valeur des items selon les coefficient grasp
		static vector<int> RCList;			// RCList correspondant à GRASPFIT pour un alpha donné
		static void createGRASPFIT(int* Coef, double alpha);	// Permet de créer GRASPFIT et RCList à partir d'un tableau de Coefficient de direction "Coef" et de "alpha"
		void GRASP( double alpha);			// crée un individu GRASP avec comme coefficient de compromis glouton aléatoire "alpha"

		void mutate();		// fait muter un individu avec restart3

		
		string toString();	// renvoie le code génétique d'une solution


		void printToScreen();	
		static void testData();	
};
#endif
