#include "Backlog_BasicModel.h"
Backlog_BasicModel::Backlog_BasicModel(Instance* Ins_In, Heuristic* Heur_In, Parameters* pIn) : Model(Ins_In, Heur_In, pIn)
{
	model = new IloModel(env);
	LPRLX = new IloModel(env);

#pragma region Variables	
	IloConstraintArray constraints(env);
	Xit = CreateNumVarArray2(env, Ins->U, Ins->T, "X");
	Yit = CreateBoolVarArray2(env, Ins->U, Ins->T, "Y");
	Sjt = CreateNumVarArray2(env, Ins->P, Ins->T, "S");
	if (P->LBonInvVars != 0 || P->UBonInvVars != 0)
	{
		for (Product j : Ins->Products)
			for (int t = 0; t < Ins->T; ++t)
			{
				if (P->LBonInvVars == 1 && P->UBonInvVars == 1)
					Sjt[j.ID - 1][t].setBounds(j.LB[t], j.UB[t]);
				else if (P->LBonInvVars == 1)
					Sjt[j.ID - 1][t].setBounds(j.LB[t], IloInfinity);
				else if (P->UBonInvVars == 1)
					Sjt[j.ID - 1][t].setBounds(0, j.UB[t]);
			}
	}
	if (P->LBonInvVars != 0 || P->UBonInvVars != 0)
	{
		for (Product j : Ins->Products)
			for (int t = 0; t < Ins->T; ++t)
			{
				if (P->LBonInvVars == 1 && P->UBonInvVars == 1)
					Sjt[j.ID - 1][t].setBounds(j.LB[t], j.UB[t]);
				else if (P->LBonInvVars == 1)
					Sjt[j.ID - 1][t].setBounds(j.LB[t], IloInfinity);
				else if (P->UBonInvVars == 1)
					Sjt[j.ID - 1][t].setBounds(0, j.UB[t]);
			}
	}
	Rjt = CreateNumVarArray2(env, Ins->P, Ins->T, "R");
#pragma endregion
#pragma region Constraints
	IloExpr expr(env);
#pragma region Inventory Balance
	for (int t = 0; t < Ins->T; ++t)
		for (Product j : Ins->Products)
		{
			expr.clear();
			if (t != 0) {
				expr += Sjt[j.ID - 1][t - 1];
				expr += -Rjt[j.ID - 1][t - 1];
			}
			expr += -Sjt[j.ID - 1][t];
			expr += Rjt[j.ID - 1][t];
			
			for (CPU i : Ins->CPUs)
				expr += i.Alphas[j.ID - 1] * Xit[i.ID - 1][t];
			constraints.add(expr == j.d[t]);
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
			expr += i.p[t] * Xit[i.ID - 1][t];
		}
		for (Product j : Ins->Products)
		{
			expr += j.h[t] * Sjt[j.ID - 1][t];
			expr += 2 /*hard coded number 2, instances do not have g values*/ * Rjt[j.ID - 1][t];
		}
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

Backlog_BasicModel::~Backlog_BasicModel()
{
	delete cplex;
	delete model;
	delete LPRLX;
	env.end();
}

void  Backlog_BasicModel::Solve()
{
	LPSolve();
	cplex->extract(*model); // get originalModel
	if (P->mode != -1)			//if not default apply auto Cplex
		cplex->setParam(IloCplex::Param::Benders::Strategy, 3);
	Model::Solve();
};

void Backlog_BasicModel::Output(ofstream& ofc)
{
	Model::Output(ofc);
	ofc << "\n";
}
