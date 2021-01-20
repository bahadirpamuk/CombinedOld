#pragma once
#include "model.h"

class BasicModel : public Model {
public:
	int userCutCounter = -1;
	int callBackCounter = -1;

	BasicModel(Instance * Ins_In, Heuristic* Heur_In, Parameters* P_In, int T_sub = 0);
	~BasicModel();
	void Solve();
	void Output(std::ofstream&);
	IloConstraintArray FindUserCuts();
	//IloConstraintArray FindSingleUserCut();
};
