#ifndef OUTPUT_HPP
#define OUTPUT_HPP


#include <vector>
#include <stdlib.h>
#include <string>
#include<iostream>
#include <fstream>
#include <sstream>

#include "individu.hpp"

class Output
{
	private:
		string generation_path;
		string pareto_path;
		string plot_path;
		string optimal_solution_path;
		int nbr_obj ;
		
	public:
	
		Output(string in_generation_path, string in_pareto_path, string in_plot_path, string in_optimal_solution_path, int nbrObj);
		~Output();
	
		void writeData(vector<Individu*> population , string output_path);
		void writeRank(int rankmax, int totalrank, double rankAve, double rankAve2, int gener);
		void writeGeneration(vector<Individu*> population, int generation_num);
		void writeRank1(vector<Individu*> rank1, int generation_num);
		void writeSol(vector<Individu*> mySol );
		
		void plotData(int nbr_gen);
		void plotSolutions(FILE *gnuplotpipe, string datafilename, string currparetofile, string optparetofile, string imagefilename, int borders[]);	
		void plotRank(FILE *gnuplotpipe, string rankmax, string totalRank, string rankAve, string imagefilename, int borders[]) ;
		void plotSol(FILE *gnuplotpipe, string mySol, string currparetofile, string imagefilename, int borders[]) ;
};



#endif
