#pragma once
#include "BasicModel.h"
class ModelWithW : public Model {
public:
	int userCutCounter = -1;
	int callBackCounter = -1;

	ModelWithW(Instance * Ins_In, Heuristic* Heur_In, Parameters* pIn);
	~ModelWithW();
	void Solve();
	void Output(std::ofstream&);
	IloConstraintArray FindUserCuts();
	//IloConstraintArray FindSingleUserCut();
};