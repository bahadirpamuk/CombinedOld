#pragma once
#include "Model.h"

class SFreeModel : public Model{
public:
	SFreeModel(Instance * Ins_In, Heuristic* Heur_In, Parameters* pIn);
	~SFreeModel();
	void Solve();
	void Output(std::ofstream&);
};