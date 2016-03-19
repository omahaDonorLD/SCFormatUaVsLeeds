#include "individuFactory.hpp"

IndividuFactory::IndividuFactory(string instance_path)
{
	Individu::initProblem(instance_path);
}

IndividuFactory::~IndividuFactory()
{

}

Individu* IndividuFactory::createRandomSolution()
{
	Individu* toReturn = new Individu();
	toReturn->randomize();

	return toReturn;
}

Individu* IndividuFactory::createGraspSolution(double alpha)
{
    Individu* toReturn = new Individu();
    toReturn->GRASP(alpha) ;
    return toReturn;
}

