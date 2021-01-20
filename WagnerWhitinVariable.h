#pragma once
#include "DCCP.h"
class WagnerWhitinVariable
{
	//vector<int> D;
	//vector<int> C;
	//vector<int> A;
	vector<double> H;
	vector<double> M;
	vector<double> F;
	vector<int> JSTAR;
	vector<int> TEMP_I;

	int N;
	Product* Prod;
	vector<CPU*> Cpus;
public:
	vector<int> I;
	vector<double> Q;
	WagnerWhitinVariable(int PeriodIn, Product* ProductIn, vector<CPU*> &CPUsIn);
	int Compute();
	~WagnerWhitinVariable();
};

