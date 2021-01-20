#include "Heuristic.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>

Heuristic::Heuristic()
{}

Heuristic::Heuristic(Instance& instance)
{
	numberCPUs = instance.U;
	numberPeriods = instance.T;
	numberProducts = instance.P;

	x = new double*[numberCPUs];
	for (int i = 0; i < numberCPUs; i++)
		x[i] = new double[numberPeriods];

	s = new double*[numberProducts];
	for (int j = 0; j < numberProducts; j++)
		s[j] = new double[numberPeriods];

	y = new bool*[numberCPUs];
	for (int i = 0; i < numberCPUs; i++)
		y[i] = new bool[numberPeriods];
}

Heuristic::~Heuristic()
{
	for (int i = 0; i < numberCPUs; i++)
		delete x[i];
	delete x;
 
	for (int j = 0; j < numberProducts; j++)
		delete s[j];
	delete s;

	for (int i = 0; i < numberCPUs; i++)
		delete y[i];
	delete y;
}

void Heuristic::CalculateObjectiveFunctionValue(Instance& instance)
{
	algObjValue = 0;

	for (int t = 0; t < numberPeriods; ++t)
		for (int k = 0; k <= t; ++k)
			for (Product j : instance.Products)
				algObjValue -= j.h[t] * j.d[k];

	for (int i = 0; i < numberCPUs; i++)
	for (int t = 0; t < numberPeriods; t++)
	{
		if (y[i][t]) {
			algObjValue += y[i][t] * instance.CPUs[i].f[t];
			algObjValue += x[i][t] * instance.CPUs[i].c[t];
		}
	}
	//if(algObjValue < 0)
	//{ 
	//	for (int i = 0; i < numberCPUs; i++)
	//		for (int t = 0; t < numberPeriods; t++)
	//		{
	//			if (y[i][t]) {
	//				cout << y[i][t] << '\t' << x[i][t] << endl;
	//			}
	//		}
	//	algObjValue = algObjValue;
	//}

}

void Heuristic::ReportOutput(string fileName)
{
	ofstream file(fileName, ofstream::app);

	file << "Heuristic Objective Value: " << algObjValue << endl;
	file << "Heuristic Duration: " << ((float)duration) / CLOCKS_PER_SEC << endl;

	file << "x(i, t) values: " << endl;
	for (int i = 0; i < numberCPUs; ++i)
	{
		for (int t = 0; t < numberPeriods; ++t)
			file << x[i][t] << " ";
		file << endl;
	}

	file << "y(i, t) values: " << endl;
	for (int i = 0; i < numberCPUs; ++i)
	{
		for (int t = 0; t < numberPeriods; ++t)
			file << y[i][t] << " ";
		file << endl;
	}

	file << "s(j, t) values: " << endl;
	for (int j = 0; j < numberProducts; ++j)
	{
		for (int t = 0; t < numberPeriods; ++t)
			file << s[j][t] << " ";
		file << endl;
	}
}

double Heuristic::GetXvalue(int i, int t)
{
	return x[i][t];
}
double Heuristic::GetSvalue(int j, int t)
{
	return s[j][t];
}
double Heuristic::GetYvalue(int i, int t)
{
	if (y[i][t])
		return 1;
	else
		return 0;
}