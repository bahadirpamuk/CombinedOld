#pragma once
#include "Model.h"

class ELS1_LP : public Model {
	NumVarArray3 Tjtt;
	NumVarArray2 Zjt;
public:
	
	ELS1_LP(Instance * Ins, Heuristic* Heur_In, Parameters* pIn);
	~ELS1_LP();
	void Solve();
	void ExportModel(string Path);
	void Output(std::ofstream&);
	IloConstraintArray FindUserCuts();
	IloConstraintArray FindSingleUserCut();
};