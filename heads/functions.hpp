#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "tools.hpp"

int RO(aUav *sln, int n);
int FTO(aUav *sln, int n);
int CO(aUav *sln, int n);

int (*OBJCTV[3])(aUav *sln, int n)={RO,FTO,CO};

double* WEIGHTS;			// Each objective weight

#endif
