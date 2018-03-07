#include "../heads/functions.hpp"

bool inRange(aNode* uav, double range, aNode ground)
{
	return ( sqrt( pow2( uav->x - ground->x ) + pow2( uav->y - ground->y ) ) > range ? false, true);
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
