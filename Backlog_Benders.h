#pragma once
#include "DCCP.h"
#include "model.h"
struct VariableInfo
{
char Type;
int i_j;
int t;

VariableInfo(char TypeIn, int i_jIn, int tIn)
{
	Type = TypeIn;
	i_j = i_jIn;
	t = tIn;
}
};
class DualSubProblem : public Model {
	IloNumExpr* ObjectiveExpr;
public:
	NumVarArray2 a_jt;
	NumVarArray2 b_it;
	IloObjective* Objective;
	vector<VariableInfo*> VariableInfos;

	DualSubProblem(Instance * Ins_In);
	~DualSubProblem();
	void Solve();
};

class MasterProblem : public Model{
	IloNumExpr* ObjectiveExpr;
	DualSubProblem* BDSP;
public:
	IloNumVar* q;
	int FeasCuts = 0;
	int OptCuts = 0;

	MasterProblem(Instance * Ins_In);
	~MasterProblem();
	void Solve();
	void Output(std::ofstream&);
};