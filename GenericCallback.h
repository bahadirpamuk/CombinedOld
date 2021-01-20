#pragma once
#include "BasicModelCallback.h"
#include "WW.h"
#include <atomic>
class GenericCallback : public IloCplex::Callback::Function
{
public:
	Model* ref;
	bool CPWW;
	int CurrentNodes[8];
	double CurrentRLXObj[8];
public:
	GenericCallback(Model* _ref, bool CPWW_In); 
	void invoke(const IloCplex::Callback::Context& context);
};

