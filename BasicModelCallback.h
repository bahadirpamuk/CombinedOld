#pragma once
#include "Model.h"

class BasicModelCallback : public Model {
public:
	BasicModelCallback(Instance * Ins_In, Heuristic* Heur_In, Parameters* pIn);
	~BasicModelCallback();
	void Solve();
	void Output(std::ofstream&);
};
