#include "../heads/functions.hpp"

bool inRange(aUav* uav, aNode ground)
{
	return ( sqrt( pow2( uav->idx - ground->x ) + pow2( uav->y - ground->y ) ) > range ? false, true);
};

int RO(aNode *sln, int n)
{

};

int FTO(aNode *sln, int n)
{
};

int CO(aNode *sln, int n)
{

};
