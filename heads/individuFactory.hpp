#ifndef INDIVIDUFACTORY_HPP
#define INDIVIDUFACTORY_HPP

#include <vector>
#include <stdlib.h>
#include <time.h>
#include <string>
#include<iostream>

#include "individu.hpp"

using namespace std;

class IndividuFactory
{
	private:

	public:
		IndividuFactory(string instance_path);
		~IndividuFactory();

		Individu* createRandomSolution();
		Individu* createGraspSolution(double alpha);

};


#endif
