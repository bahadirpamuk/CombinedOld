#include "BasicModel.h"
#include "GenericCallback.h"
BasicModel::BasicModel(Instance* Ins_In, Heuristic* Heur_In, Parameters* pIN, int T_sub) : Model(Ins_In, Heur_In, pIN)
{	
	if (T_sub == 0)
		T_sub = Ins->T;

	P->T_Sub = T_sub;

	model = new IloModel(env);
	LPRLX = new IloModel(env);

	userCutCounter = 0;
	callBackCounter = 0; 
	cutTime = 0;

#pragma region Variables	
	IloConstraintArray constraints(env);
	Xit = CreateNumVarArray2(env, Ins->U, P->T_Sub, "X");
	Yit = CreateBoolVarArray2(env, Ins->U, P->T_Sub, "Y");
	Sjt = CreateNumVarArray2(env, Ins->P, P->T_Sub, "S");
	if (P->LBonInvVars != 0 || P->UBonInvVars != 0)
	{
		for (Product j : Ins->Products)
			for (int t = 0; t < P->T_Sub; ++t)
			{
				if (P->LBonInvVars == 1 && P->UBonInvVars == 1)
					Sjt[j.ID - 1][t].setBounds(j.LB[t], j.UB[t]);
				else if (P->LBonInvVars == 1)
					Sjt[j.ID - 1][t].setBounds(j.LB[t], IloInfinity);
				else if (P->UBonInvVars == 1)
					Sjt[j.ID - 1][t].setBounds(0, j.UB[t]);
			}
	}
#pragma endregion
#pragma region Constraints
	IloExpr expr(env);
#pragma region Inventory Balance
	for (int t = 0; t < P->T_Sub; ++t)
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
	for (int t = 0; t < P->T_Sub; ++t)
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
	if(P->WWV){
	Wi = CreateBoolVarArray(env, Ins->U, "W");
	#pragma region W Constraints1
	for (int t = 0; t < P->T_Sub; ++t)
		for (CPU i : Ins->CPUs)
		{
			constraints.add(Yit[i.ID - 1][t] <= Wi[i.ID - 1]);
		}
#pragma endregion
#pragma region W Constraints2
		for (CPU i : Ins->CPUs)
		{
			constraints.add(Wi[i.ID - 1] <= IloSum(Yit[i.ID - 1]));
			//expr.clear();
			//expr += Wi[i.ID - 1];
			//for (int t = 0; t < Ins->T; ++t)
			//	expr -= Yit[i.ID - 1][t];
			//constraints.add(expr >= 0);
		}
#pragma endregion
	}
	model->add(constraints); //adding all constraints at once
#pragma endregion
#pragma region Objective
	expr.clear();
	for (int t = 0; t < P->T_Sub; ++t)
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
	{
		if(P->WWV){
			LPRLX->add(IloConversion(env, Wi[i], ILOFLOAT));
		}
		LPRLX->add(IloConversion(env, Yit[i], ILOFLOAT));
	}
	cplex = new IloCplex(env);
#pragma region Export Model
	//cplex->exportModel("BasicModel.lp");
#pragma endregion
};

BasicModel::~BasicModel()
{
	delete cplex;
	delete model;
	delete LPRLX;
	env.end();
}

void  BasicModel::Solve()
{
	vector<double> solns;
	BasicModel* littleModel; 
	high_resolution_clock::time_point startTime;
	if (P->VI1) {
		double tempTimelim = P->TimeLimitMinutes;
		P->TimeLimitMinutes = tempTimelim / (5 * (Ins->T / 6)); //for little models
		startTime = high_resolution_clock::now();
		for (int i = 1; i < (Ins->T / 6) + 1; ++i)
		{
			P->VI1 = false;
			P->littleModel = true;
			littleModel = new BasicModel(Ins, Heur, P, i);
			P->VI1 = true;
			littleModel->cplex->extract(*littleModel->model);
			littleModel->cplex->exportModel("LPExport.lp");
			littleModel->Model::Solve();
			if(cplex->getStatus() == IloCplex::Status::Optimal)
				solns.push_back(littleModel->cplex->getObjValue());
			else
				solns.push_back(littleModel->cplex->getBestObjValue());
			delete littleModel;
		}
		IloExpr expr(env);
		for (int t = 0; t < (Ins->T / 6); ++t)
		{
			//cout << t << "///////////////////////////////////////////////////////////////////////////////////\n";
			for (int t1 = 0; t1 <= t; ++t1)
			{
				for (CPU i : Ins->CPUs)
				{
					expr += i.f[t1] * Yit[i.ID - 1][t1];
					expr += i.p[t1] * Xit[i.ID - 1][t1];
				}
				for (Product j : Ins->Products)
					expr += j.h[t1] * Sjt[j.ID - 1][t1];
			}
			model->add(expr >= solns[t]);
			LPRLX->add(expr >= solns[t]);
		}
		FTime = duration_cast<duration<double>>(high_resolution_clock::now() - startTime).count();
		P->littleModel = false;
		P->TimeLimitMinutes = tempTimelim;
	}
	LPSolve();
	if (P->mode == 0 ) {
		AddCuts();
	}
	else {
		cplex->extract(*model); // get originalModel
		GenericCallback cb((Model*)this, P->CPWW);
		if (P->WWV || P->VI2) {		
			CPXLONG contextmask = IloCplex::Callback::Context::Id::Relaxation;
			cplex->use(&cb, contextmask); //use callback
		}

		if(P->WWV){	
			for (CPU i : Ins->CPUs)
			{
				cplex->setPriority(Wi[i.ID - 1], 1); //default priority is 0, highest priority is branched first
				cplex->setDirection(Wi[i.ID - 1], IloCplex::BranchDirection::BranchDown);
			}
			cplex->setParam(IloCplex::Param::MIP::Strategy::NodeSelect, 0); //depth first search
		}
		Model::Solve();
	}
};

void BasicModel::Output(ofstream& ofc)
{
	Model::Output(ofc);
	if (P->mode == 0)
	{
		ofc << "\t" << userCutCounter
			<< "\t" << Total_NConst
			<< "\t" << callBackCounter
			<< "\t" << cutTime
			<< "\t" << HeurTime << "\n";
	}
	else
		ofc << "\n";
}

IloConstraintArray BasicModel::FindUserCuts()
{
	callBackCounter++;//count the number of times this callback is entered

	NumArray2 Val_Xit = CreateNumArray2(env, Ins->U, P->T_Sub);
	NumArray2 Val_Yit = CreateNumArray2(env, Ins->U, P->T_Sub);
	vector<bool> S(P->T_Sub);
	IloExpr expr(env);
	IloConstraintArray constraints(env);

	for (int i = 0; i < Ins->U; ++i)
		cplex->getValues(Val_Xit[i], Xit[i]);

	for (int i = 0; i < Ins->U; ++i)
		cplex->getValues(Val_Yit[i], Yit[i]);

	for (int l = P->T_Sub - 1; l >= 0; --l)
	{
		for (int j = 0; j < Ins->P; ++j)
		{
			for (int t = 0; t < P->T_Sub; ++t) //clear the array
				S[t] = false;
			double sum = 0;
			for (int q = 0; q <= l; ++q)
			{
				IloNum sumX = 0;
				for (auto i : Ins->Products[j].CPUs)
					sumX += i->Alphas[j] * Val_Xit[i->ID - 1][q];

				IloNum sumY = 0;
				for (CPU* c : Ins->Products[j].CPUs)
					sumY += Ins->Products[j].d_t1_t2[q][l] * Val_Yit[c->ID - 1][q];

				if (S[q] = (sumY < sumX))
					sum += sumY;
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
	Val_Yit.end();
	expr.end();
	return constraints;
}
