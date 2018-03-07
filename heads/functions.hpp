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
  int index=-1;// Unique, defines a point : for either
	double x=0.;
	double y=0.;
}aNode;

typedef struct aUav{
  aNode id;
  int covers=-1;	// If value is less than 0 then aNode is a ground node. A UaV contains at least 0 ground nodes
	double range=-1.;  	// If "contains" < 0 then is ground node and its range is also < 1
  bool active=false;  // ground nodes are always unactive. a uav can also be
}aUav;

bool inRange(aUav* uav, aNode ground);
int RO(aUav *sln, int n);
int FTO(aUav *sln, int n);
int CO(aUav *sln, int n);

int (*OBJCTV[3])(aUav *sln, int n)={RO,FTO,CO};

#endif
