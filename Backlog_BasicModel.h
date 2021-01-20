#pragma once
#pragma once
#include "model.h"

class Backlog_BasicModel : public Model {
	NumVarArray2 Sjt;
	NumVarArray2 Rjt;

public:
	Backlog_BasicModel(Instance * Ins_In, Heuristic* Heur_In, Parameters *pIn);
	~Backlog_BasicModel();
	void Solve();
	void Output(std::ofstream&);
};