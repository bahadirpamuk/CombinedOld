#pragma once
#include "Model.h"
class BasicModelwithZ_LPCutsOnly : public Model
{
public:
	BasicModelwithZ_LPCutsOnly(Instance * Ins_In, Heuristic* Heur_In);
	~BasicModelwithZ_LPCutsOnly();
	void Solve();
	void Output(std::ofstream&);
};

