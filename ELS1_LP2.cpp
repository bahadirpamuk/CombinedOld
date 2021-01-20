#include "ELS1_LP2.h"


void ELS1_LP2::ExportModel(string Path)
{
	cplex->exportModel(Path.c_str());
}

ELS1_LP2::ELS1_LP2(Instance* Ins_In, Heuristic* Heur_In, Parameters* pIN) : Model(Ins_In, Heur_In, pIN)
{
	userCutCounter = 0;
	callbackCounter = 0;
	cutTime = 0;

	model = new IloModel(env);
	LPRLX = new IloModel(env);

#pragma region Variables	
	IloConstraintArray constraints(env);
	Xit = CreateNumVarArray2(env, Ins->U, Ins->T, "X");
	Yit = CreateBoolVarArray2(env, Ins->U, Ins->T, "Y");
	Tjtt = CreateNumVarArray3(env, Ins->P, Ins->T, Ins->T, "Theta");
#pragma endregion
#pragma region Constraints
	IloExpr expr(env);
#pragma region  Demand Satisfaction
	for (int t1 = 0; t1 < Ins->T; ++t1)
		for (Product j : Ins->Products)
		{
			expr.clear();
			for (int t = 0; t <= t1; ++t)
				expr += Tjtt[j.ID - 1][t][t1];
			constraints.add(expr == j.d[t1]);
		}
#pragma endregion
#pragma region X and Theta Link
	for (int t = 0; t < Ins->T; ++t)
		for (Product j : Ins->Products)
		{
			expr.clear();
			for (int t1 = t; t1 < Ins->T; ++t1)
				expr += Tjtt[j.ID - 1][t][t1];
			for (CPU* i : j.CPUs)
				expr -= Xit[i->ID - 1][t] * i->Alphas[j.ID - 1];
			constraints.add(expr <= 0);
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
#pragma endregion
};

ELS1_LP2::~ELS1_LP2()
{
	delete cplex;
	delete model;
	delete LPRLX;
	AddedCuts->end();
	delete AddedCuts;
	env.end();
}

void  ELS1_LP2::Solve()
{
	LPSolve();
	AddCuts();
	if (P->mode != 0 && P->mode != 10 && P->epsilon != -1)
		Model::Solve();
};

IloConstraintArray ELS1_LP2::FindUserCuts()
{
	callbackCounter++;//count the number of times this callback is entered

	NumArray3 Val_Tjtt = CreateNumArray3(env, Ins->P, Ins->T, Ins->T);
	NumArray2 Val_Yit = CreateNumArray2(env, Ins->U, Ins->T);
	vector<bool> S(Ins->T);
	IloExpr expr(env);
	IloConstraintArray constraints(env);

	for (int j = 0; j < Ins->P; ++j)
		for (int t = 0; t < Ins->T; ++t)
			for (int tp = t; tp < Ins->T; ++tp)
				Val_Tjtt[j][t][tp] = cplex->getValue(Tjtt[j][t][tp]);

	for (int i= 0; i < Ins->U; ++i)
		cplex->getValues(Val_Yit[i], Yit[i]);

	for (int l = Ins->T - 1; l >= 0; --l)
	{
		for (int j = 0; j < Ins->P; ++j)
		{
			for (int t = 0; t < Ins->T; ++t) //clear the array
				S[t] = false;
			double sum = 0;
			for (int q = 0; q <= l; ++q)
			{
				IloNum sumT = 0;
				for (int t1 = q; t1 < Ins->T; ++t1)
					sumT += Val_Tjtt[j][q][t1];

				IloNum sumY = 0;
				for(CPU* c : Ins->Products[j].CPUs)
					sumY += Ins->Products[j].d_t1_t2[q][l] * Val_Yit[c->ID - 1][q];

				if (S[q] = (sumY < sumT))
					sum += sumY;
				else
					sum += sumT;
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
						for (int t1 = q; t1 < Ins->T; ++t1)
							expr += Tjtt[j][q][t1];
					}
				}
				constraints.add(expr >= Ins->Products[j].d_t1_t2[0][l]); // add the cut
				userCutCounter++;//count the number of constraints added
			}
		}
	}
	Val_Tjtt.end();
	Val_Yit.end();
	expr.end();
	return constraints;
}

void ELS1_LP2::Output(ofstream& ofc)
{
	Model::Output(ofc);
	ofc << "\t" << userCutCounter
		<< "\t" << Total_NConst
		<< "\t" << callbackCounter
		<< "\t" << cutTime
		<< "\t" << HeurTime << "\n";
}