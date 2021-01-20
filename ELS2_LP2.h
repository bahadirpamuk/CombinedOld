#pragma once
#include "Model.h"

class ELS2_LP2 : public Model {
	NumVarArray3* Tijtt;
public:

	ELS2_LP2(Instance* Ins_In, Heuristic* Heur_In, Parameters* pIn);
	~ELS2_LP2();
	void Solve();
	void Output(std::ofstream&);
	IloConstraintArray FindUserCuts();
};