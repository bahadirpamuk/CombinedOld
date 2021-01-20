#pragma once
#include "DCCP.h"
class WagnerWhitin
{
	vector<int> D;
	vector<int> C;
	vector<int> A;
	vector<int> H;

	vector<int> M;
	vector<int> F;
	
	vector<int> S;
	vector<int> JSTAR;

	int N;
	Product* Prod;
	CPU* Cpu;
public:
	vector<int> Q;
	WagnerWhitin(int PeriodIn, Product* ProductIn, CPU* CPUIn);
	int Compute();
	~WagnerWhitin();
};

