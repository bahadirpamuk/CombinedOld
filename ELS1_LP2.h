#pragma once
#include "Model.h"

class ELS1_LP2 : public Model {
	NumVarArray3 Tjtt;
public:

	ELS1_LP2(Instance * Ins, Heuristic* Heur_In, Parameters* pIn);
	~ELS1_LP2();
	void Solve();
	void ExportModel(string Path);
	void Output(std::ofstream&);
	IloConstraintArray FindUserCuts();
};