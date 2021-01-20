#pragma once
#include "Model.h"
class BasicModelwithZ :
	public Model
{
	NumVarArray2 Zjt;
public:

	BasicModelwithZ(Instance * Ins_In, Heuristic* Heur_In, Parameters* pIn);
	~BasicModelwithZ();
	void Solve();
	void Output(std::ofstream&);
	IloConstraintArray FindUserCuts();
	//IloConstraintArray FindSingleUserCut();
};

