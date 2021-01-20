#pragma once
#include "Heuristic.h"

using namespace std;

class WW : public Heuristic
{
public:
	WW();
	WW(Instance instance);
	~WW();
	void ImplementAlgorithm(Instance& instance, NumArray2 ySol);
};