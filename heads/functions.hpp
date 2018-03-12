#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

// Three objectives to optimize
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>
#include <math.h>

using namespace std;

// Each node have coordinates
typedef struct aNode {
	int index;// Unique, defines a point : for either a ground node or a uav
	double x;
	double y;
}aNode;

typedef struct aUav{
	aNode id;
	//long gene=0;
	int covers;	// A UaV contains at least 0 ground nodes
	double range;  	// If "contains" < 0 then is ground node and its range is also < 1
  bool active;  // ground nodes are always unactive. a uav can also be
}aUav;

int RO(aUav *sln, int n);
int FTO(aUav *sln, int n);
int CO(aUav *sln, int n);

int (*OBJCTV[3])(aUav *sln, int n)={RO,FTO,CO};

#endif
