#include "BasicModelwithZ_LPCutsOnly.h"



BasicModelwithZ_LPCutsOnly::BasicModelwithZ_LPCutsOnly(Instance * Ins_In, Heuristic* Heur_In) : Model(Ins_In, Heur_In)
{
	model = new IloModel(env);
#pragma region Variables	
	IloConstraintArray constraints(env);
	Xit = CreateNumVarArray2(env, Ins->U, Ins->T, "X");
	Yit = CreateNumVarArray2(env, Ins->U, Ins->T, "Y");
	Sjt = CreateNumVarArray2(env, Ins->P, Ins->T, "S");
	Zjt = CreateNumVarArray2(env, Ins->P, Ins->T, "Z");
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
			expr += Xit[i.ID - 1][t];
			//take the maximum dtT / alpha for each J(i)
			double temp = 0;

			for (Product j : Ins->Products)
				if (i.Alphas[j.ID - 1] > 0 && j.d_t_T[t] / i.Alphas[j.ID - 1] > temp)
					temp = j.d_t_T[t] / i.Alphas[j.ID - 1];

			expr += -temp * Yit[i.ID - 1][t];

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
}


BasicModelwithZ_LPCutsOnly::~BasicModelwithZ_LPCutsOnly()
{
}
