#ifndef PATTERNFITTING_H
#define PATTERNFITTING_H
#include "Heuristic.h"

using namespace std;

class PatternFitting: public Heuristic
{
public:
	PatternFitting();
	PatternFitting(Instance instance) ;
	~PatternFitting();
	void ImplementAlgorithm(Instance& instance);
};
#endif