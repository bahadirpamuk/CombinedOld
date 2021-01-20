#pragma once
#include "Model.h"

class ELS1 : public Model {
	NumVarArray3 Tjtt;
public:
	ELS1(Instance * Ins_In, Heuristic* Heur_In, Parameters* pIn);
	~ELS1();
	void Solve();
	void Output(std::ofstream&);
};