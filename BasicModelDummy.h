#pragma once
#include "model.h"

class BasicModelDummy : public Model {
public:
	IloBoolVar Dummy;

	BasicModelDummy(Instance * Ins_In, Heuristic* Heur_In, Parameters* pIn);
	~BasicModelDummy();
	void Solve();
	void Output(std::ofstream&);
	IloConstraintArray FindUserCuts();
	//IloConstraintArray FindSingleUserCut();
};
