#ifndef PATTERNFITTINGBACKLOG_H
#define PATTERNFITTINGBACKLOG_H
#include "Heuristic.h"
#include "DCCP.h"

using namespace std;

class PatternFittingBacklog : public Heuristic
{
public:
	PatternFittingBacklog();
	PatternFittingBacklog(Instance instance);
	~PatternFittingBacklog();
	void ImplementAlgorithm(Instance& instance);
	void CalcObj(Instance& instance);
};
#endif