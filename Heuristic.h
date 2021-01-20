#ifndef HEURISTIC_H
#define HEURISTIC_H
#include <time.h>
#include "DCCP.h"

class Heuristic
{
protected:
	double algObjValue = -1;
	clock_t duration;

	int numberCPUs;
	int numberProducts;
	int numberPeriods;

	double** x;
	double** s;
	bool** y;

	void CalculateObjectiveFunctionValue(Instance& instance);

public:
	Heuristic();
	Heuristic(Instance& instance);
	~Heuristic();
	void ReportOutput(string fileName);
	virtual void ImplementAlgorithm(Instance& instance) {};
	double GetXvalue(int i, int t);
	double GetSvalue(int j, int t);
	double GetYvalue (int i, int t);
	double GetObj() { return algObjValue; }
	double GetTime() { return (double)duration / CLOCKS_PER_SEC; }
};
#endif