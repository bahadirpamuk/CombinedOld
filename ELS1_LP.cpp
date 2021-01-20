#include "ELS1_LP.h"

void ELS1_LP::ExportModel(string Path)
{
	cplex->exportModel(Path.c_str());
}

ELS1_LP::ELS1_LP(Instance* Ins_In, Heuristic* Heur_In, Parameters* pIN) : Model(Ins_In, Heur_In,pIN)
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
	Zjt = CreateNumVarArray2(env, Ins->P, Ins->T, "Z" , 0 , 1);
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

ELS1_LP::~ELS1_LP()
{
	delete cplex;
	delete model;
	delete LPRLX;
	AddedCuts->end();
	delete AddedCuts;
	env.end();
}

void  ELS1_LP::Solve()
{
	LPSolve();
	AddCuts();
	if (P->mode != 0 && P->mode != 10 && P->epsilon != -1)
		Model::Solve();
};

IloConstraintArray ELS1_LP::FindUserCuts()
{
	callbackCounter++;//count the number of times this callback is entered

	NumArray3 Val_Tjtt = CreateNumArray3(env, Ins->P, Ins->T, Ins->T);
	NumArray2 Val_Zjt = CreateNumArray2(env, Ins->P, Ins->T);
	vector<bool> S(Ins->T);
	IloExpr expr(env);
	IloConstraintArray constraints(env);
	
	for (int j = 0; j < Ins->P; ++j)
		for (int t = 0; t < Ins->T; ++t)
			for (int tp = t; tp < Ins->T; ++tp)
				Val_Tjtt[j][t][tp] = cplex->getValue(Tjtt[j][t][tp]);

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
				IloNum sumT = 0;
				for (int t1 = q ; t1 < Ins->T ; ++ t1)
						sumT += Val_Tjtt[j][q][t1];
				IloNum sumZ = Ins->Products[j].d_t1_t2[q][l] * Val_Zjt[j][q];

				if (S[q] = (sumZ < sumT))
					sum += sumZ;
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
						expr += Ins->Products[j].d_t1_t2[q][l] * Zjt[j][q];
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
	Val_Zjt.end();
	expr.end();
	return constraints;
}

IloConstraintArray ELS1_LP::FindSingleUserCut()
{
	callbackCounter++;//count the number of times this callback is entered

	NumArray3 Val_Tjtt = CreateNumArray3(env, Ins->P, Ins->T, Ins->T);
	NumArray2 Val_Zjt = CreateNumArray2(env, Ins->P, Ins->T);
	vector<bool> S(Ins->T);
	IloExpr expr(env);
	IloConstraintArray constraints(env);

	IloNum MaxDelta = P->epsilon;
	int jmax = -1;
	int lmax = -1;

	for (int j = 0; j < Ins->P; ++j)
		for (int t = 0; t < Ins->T; ++t)
			for (int tp = t; tp < Ins->T; ++tp)
				Val_Tjtt[j][t][tp] = cplex->getValue(Tjtt[j][t][tp]);

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
				IloNum sumT = 0;
				for (int t1 = q; t1 < Ins->T; ++t1)
					sumT += Val_Tjtt[j][q][t1];
				IloNum sumZ = Ins->Products[j].d_t1_t2[q][l] * Val_Zjt[j][q];

				if (S[q] = (sumZ < sumT))
					sum += sumZ;
				else
					sum += sumT;
			}

			if (MaxDelta < Ins->Products[j].d_t1_t2[0/*zero*/][l] - sum)
			{
				MaxDelta = Ins->Products[j].d_t1_t2[0/*zero*/][l] - sum;
				jmax = j;
				lmax = l;

				expr.clear();
				for (int q = 0; q <= l; ++q)
				{
					if (S[q]) //q in set S
					{
						expr += Ins->Products[j].d_t1_t2[q][l] * Zjt[j][q];
					}
					else	//q not in set S
					{
						for (int t1 = q; t1 < Ins->T; ++t1)
							expr += Tjtt[j][q][t1];
					}
				}
				userCutCounter++;//count the number of constraints added
			}
		}
	}

	if (jmax != -1)
		constraints.add(expr >= Ins->Products[jmax].d_t1_t2[0][lmax]); // add the cut

	Val_Tjtt.end();
	Val_Zjt.end();
	expr.end();
	return constraints;
}

void ELS1_LP::Output(ofstream& ofc)
{
	Model::Output(ofc);
	ofc << "\t" << userCutCounter
		<< "\t" << Total_NConst
		<< "\t" << callbackCounter
		<< "\t" << cutTime
		<< "\t" << HeurTime << "\n";
}