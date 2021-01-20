#pragma once
#include "Model.h"

class ELS2 : public Model{
	NumVarArray3* Tijtt;
public:
	ELS2(Instance * Ins_In, Heuristic* h, Parameters* pIn);
	~ELS2();
	void Solve();
	void Output(std::ofstream&);
};