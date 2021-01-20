#include "BasicModelwithZ.h"
BasicModelwithZ::BasicModelwithZ(Instance* Ins_In, Heuristic* Heur_In, Parameters* pIN) : Model(Ins_In, Heur_In, pIN)
{	model = new IloModel(env);
	LPRLX = new IloModel(env);

	userCutCounter = 0;
	callbackCounter = 0;
	cutTime = 0;

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
	Zjt = CreateNumVarArray2(env, Ins->P, Ins->T, "Z", 0 , 1);
#pragma endregion
#pragma region Constraints
	IloExpr expr(env);
#pragma region Inventory Balance
	for (int t = 0; t < Ins->T; ++t)
		for (Product j : Ins->Products)
		{
			expr.clear();
			if (t != 0)
				expr += Sjt[j.ID - 1][t - 1];
			expr += -Sjt[j.ID - 1][t];
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
#pragma region Z and Y Link 1
	for (int t = 0; t < Ins->T; ++t)
		for (Product j : Ins->Products)
			for (CPU* Ij : j.CPUs)
			{
				expr.clear();
				expr += Yit[Ij->ID - 1][t];
				expr -= Zjt[j.ID - 1][t];
				constraints.add(expr <= 0);
			}
#pragma endregion
#pragma region Z and Y Link 2
	for (int t = 0; t < Ins->T; ++t)
		for (Product j : Ins->Products)
		{
			expr.clear();
			expr -= Zjt[j.ID - 1][t];
			for (CPU* Ij : j.CPUs)
			{
				expr += Yit[Ij->ID - 1][t];
			}
			constraints.add(expr >= 0);
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
			expr += j.h[t] * Sjt[j.ID - 1][t];
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

BasicModelwithZ::~BasicModelwithZ()
{
	delete cplex;
	delete model;
	delete LPRLX;
	env.end();
}

void  BasicModelwithZ::Solve()
{
	LPSolve();
	if (P->mode == 0) {
		AddCuts();
	}
	else {
		cplex->extract(*model); // get originalModel
		Model::Solve();
	}
};

void BasicModelwithZ::Output(ofstream& ofc)
{
	Model::Output(ofc);
	if (P->mode == 0)
	{
		ofc << "\t" << userCutCounter
			<< "\t" << Total_NConst
			<< "\t" << callbackCounter
			<< "\t" << cutTime
			<< "\t" << HeurTime << "\n";
	}
	else
		ofc << "\n";
}

IloConstraintArray BasicModelwithZ::FindUserCuts()
{
	callbackCounter++;//count the number of times this callback is entered

	NumArray2 Val_Xit = CreateNumArray2(env, Ins->U, Ins->T);
	NumArray2 Val_Zjt = CreateNumArray2(env, Ins->P, Ins->T);
	vector<bool> S(Ins->T);
	IloExpr expr(env);
	IloConstraintArray constraints(env);

	for (int i = 0; i < Ins->U; ++i)
		cplex->getValues(Val_Xit[i], Xit[i]);

	for (int j = 0; j < Ins->P; ++j)
		cplex->getValues(Val_Zjt[j], Zjt[j]);

	for (int l = Ins->T - 1; l >= 0; --l)
	{
		for (int j = 0; j < Ins->P; ++j)
		{
			for (int t = 0; t < Ins->T; ++t) //clear the array
				S[t] = false;
			double sum = 0;
			for (int q = 0; q <= l; ++q)
			{
				IloNum sumX = 0;
				for (auto i : Ins->Products[j].CPUs)
					sumX += i->Alphas[j] * Val_Xit[i->ID - 1][q];

				IloNum sumZ = Ins->Products[j].d_t1_t2[q][l] * Val_Zjt[j][q];

				if (S[q] = (sumZ < sumX))
					sum += sumZ;
				else
					sum += sumX;
			}

			if (P->epsilon < Ins->Products[j].d_t1_t2[0/*zero*/][l] - sum)
			{
				expr.clear();
				for (int q = 0; q <= l; ++q)
				{
					if (S[q]) //q in set S
					{
						for (CPU* c : Ins->Products[j].CPUs)
							expr += Ins->Products[j].d_t1_t2[q][l] * Yit[c->ID - 1][q];
					}
					else	//q not in set S
					{
						for (auto i : Ins->Products[j].CPUs)
							expr += i->Alphas[j] * Xit[i->ID - 1][q];
					}
				}
				constraints.add(expr >= Ins->Products[j].d_t1_t2[0][l]); // add the cut
				userCutCounter++;//count the number of constraints added
			}
		}
	}
	Val_Xit.end();
	Val_Zjt.end();
	expr.end();
	return constraints;
}