#include "BasicModel_S-Free.h"
SFreeModel::SFreeModel(Instance* Ins_In, Heuristic* Heur_In, Parameters* pIN) : Model(Ins_In, Heur_In, pIN)
{
	model = new IloModel(env);
	LPRLX = new IloModel(env);

#pragma region Variables	
	IloConstraintArray constraints(env);
	Xit = CreateNumVarArray2(env, Ins->U, Ins->T, "X");
	Yit = CreateBoolVarArray2(env, Ins->U, Ins->T, "Y");
#pragma endregion
#pragma region Constraints
	IloExpr expr(env);
#pragma region Demand Satisfaction
	for (int t = 0; t < Ins->T; ++t)
		for (Product j : Ins->Products)
		{
			expr.clear();

			for (int k = 0; k <= t; ++k)
				for (CPU* i : j.CPUs)
					expr += i->Alphas[j.ID - 1] * Xit[i->ID - 1][k];
			constraints.add(expr >= j.d_t1_t2[0][t]);
		}
#pragma endregion
#pragma region X and Y Link
	for (int t = 0; t < Ins->T; ++t)
		for (CPU i : Ins->CPUs)
		{
			expr.clear();

			//take the maximum dtT / alpha for each J(i)
			double temp = 0;
			int index_j = -1;

			for (Product j : Ins->Products)
				if (i.Alphas[j.ID - 1] > 0 && j.d_t_T[t] / i.Alphas[j.ID - 1] > temp)
				{
					temp = j.d_t_T[t] / i.Alphas[j.ID - 1];
					index_j = j.ID - 1;
				}

			expr += i.Alphas[index_j] * Xit[i.ID - 1][t];
			expr += -Ins->Products[index_j].d_t_T[t] * Yit[i.ID - 1][t];

			constraints.add(expr <= 0);
		}
#pragma endregion
	model->add(constraints); //adding all constraints at once
#pragma endregion
#pragma region Objective
	expr.clear();
	for (int t = 0; t < Ins->T; ++t)
	{
		for (CPU i : Ins->CPUs)
		{
			expr += i.f[t] * Yit[i.ID - 1][t];
			expr += i.c[t] * Xit[i.ID - 1][t];
		}
		for (Product j : Ins->Products)
			for (int k = 0; k <= t; ++k)
				expr -= j.h[t] * j.d[k];
	}
	IloObjective obj = IloMinimize(env, expr);
	model->add(obj); //adding objective function
#pragma endregion
	LPRLX->add(*model);
	for (int i = 0; i < Ins->U; ++i)
		LPRLX->add(IloConversion(env, Yit[i], ILOFLOAT));
	cplex = new IloCplex(env);
#pragma region Export Model
	//cplex->exportModel("BasicModel.lp");
#pragma endregion
};

SFreeModel::~SFreeModel()
{
	delete cplex;
	delete model;
	delete LPRLX;
	env.end();
}

void  SFreeModel::Solve()
{
	LPSolve();
	cplex->extract(*model); // get originalModel
	Model::Solve();
};

void SFreeModel::Output(ofstream& ofc)
{
	Model::Output(ofc);
	ofc << "\n";
}