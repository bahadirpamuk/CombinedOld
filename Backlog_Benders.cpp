#include "Backlog_Benders.h"
#include <mutex>
mutex theMutex;

MasterProblem::MasterProblem(Instance* Ins_In) : Model(Ins_In)
{
	BDSP = new DualSubProblem(Ins_In);
	model = new IloModel(env);
	ObjectiveExpr = new IloNumExpr(env);
	q = new IloNumVar(env, 0, IloInfinity, "q");	
#pragma region Variables
	Yit = CreateBoolVarArray2(env, Ins->U, Ins->T, "Y");
#pragma endregion
#pragma region Objective
	*ObjectiveExpr += *q;
	for (int t = 0; t < Ins->T; ++t)
		for (CPU i : Ins->CPUs)
			*ObjectiveExpr += i.f[t] * Yit[i.ID - 1][t];

	model->add(IloMinimize(env, *ObjectiveExpr)); //adding objective function
#pragma endregion
	cplex = new IloCplex(env);
#pragma endregion
};

MasterProblem::~MasterProblem()
{
	delete BDSP;
	delete cplex;
	delete model;
	delete ObjectiveExpr;
	delete q;
	//delete LPRLX;
	env.end();
}

ILOLAZYCONSTRAINTCALLBACK3(BendersLazyCallback, Instance, Ins, MasterProblem&, Master, DualSubProblem&, BDSP)
{
	lock_guard<mutex> lock(theMutex);

	for (int t = 0; t < Ins.T; ++t)
		for (int i = 0; i < Ins.U; ++i)
		{
			//take the maximum dtT / alpha for each J(i)
			double temp = 0;
			for (Product j : Ins.Products)
				if (Ins.CPUs[i].Alphas[j.ID - 1] > 0 && j.d_t_T[t] / Ins.CPUs[i].Alphas[j.ID - 1] > temp)
					temp = j.d_t_T[t] / Ins.CPUs[i].Alphas[j.ID - 1];
			//take the maximum dtT / alpha for each J(i)

			IloNum yValue = getValue(Master.Yit[i][t]);
			BDSP.Objective->setLinearCoef(BDSP.b_it[i][t], temp * yValue);
		}
	
	BDSP.Solve();

	double q = getValue(*(Master.q));
	double obj_dual = BDSP.cplex->getObjValue();

	IloCplex::CplexStatus status = BDSP.cplex->getCplexStatus();
	if (BDSP.cplex->getCplexStatus() == IloCplex::Unbounded)
	{
		IloNumArray vals(BDSP.cplex->getEnv());
		IloNumVarArray vars(BDSP.cplex->getEnv());
		BDSP.cplex->getRay(vals, vars);

		IloExpr BendersCut(getEnv());
		for (int currentIndex = 0; currentIndex < vars.getSize(); ++currentIndex)
		{
			IloNumVar currentVar = vars[currentIndex];
			VariableInfo* pVI = static_cast<VariableInfo*>(currentVar.getObject());
			if (pVI->Type == 'a')
				BendersCut += vals[currentIndex] * Ins.Products[pVI->i_j].d[pVI->t];
			else if (pVI->Type == 'b')
			{
				//take the maximum dtT / alpha for each J(i)
				double temp = 0;
				for (Product j : Ins.Products)
					if (Ins.CPUs[pVI->i_j].Alphas[j.ID - 1] > 0 && j.d_t_T[pVI->t] / Ins.CPUs[pVI->i_j].Alphas[j.ID - 1] > temp)
						temp = j.d_t_T[pVI->t] / Ins.CPUs[pVI->i_j].Alphas[j.ID - 1];
				//take the maximum dtT / alpha for each J(i)
				BendersCut += (vals[currentIndex] * temp) * Master.Yit[pVI->i_j][pVI->t];
			}
		}
		//		cout << "Benders cut added!" << endl;
		add(BendersCut <= 0);
		Master.FeasCuts++;
	}
	else if (abs(q - obj_dual) > 0.000001)//optimality cut is needed
	{
		IloExpr BendersCut(getEnv());
		BendersCut += *(Master.q);
		for (int t = 0; t < Ins.T; ++t)
			for (int i = 0; i < Ins.U; ++i)
			{
				//take the maximum dtT / alpha for each J(i)
				double temp = 0;
				for (Product j : Ins.Products)
					if (Ins.CPUs[i].Alphas[j.ID - 1] > 0 && j.d_t_T[t] / Ins.CPUs[i].Alphas[j.ID - 1] > temp)
						temp = j.d_t_T[t] / Ins.CPUs[i].Alphas[j.ID - 1];
				//take the maximum dtT / alpha for each J(i)
				BendersCut -= (temp * BDSP.cplex->getValue(BDSP.b_it[i][t])) * Master.Yit[i][t];		
			}
		double constant = 0;
		for (int t = 0; t < Ins.T; ++t)
			for (int j = 0; j < Ins.P; ++j)
				constant += BDSP.cplex->getValue(BDSP.a_jt[j][t]) * Ins.Products[j].d[t];
		add(BendersCut >= constant);
		Master.OptCuts++;
	}
}

void  MasterProblem::Solve()
{
	cplex->extract(*model);
	cplex->use(BendersLazyCallback(cplex->getEnv(), *Ins, *this, *BDSP));
	cplex->exportModel("MasterProblem.lp");
	cplex->setParam(IloCplex::Param::TimeLimit, 1200);				//set time limit to 20min
	cplex->setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0);	//force to solving to optimality
	cplex->setParam(IloCplex::Threads, 8);							//use all threads
	high_resolution_clock::time_point startTime = high_resolution_clock::now();
	try
	{
		bool success = cplex->solve();
	}
	catch (IloException e)
	{
		cout << e;
	}
	high_resolution_clock::time_point endTime = high_resolution_clock::now();
	MIPTime = duration_cast<duration<double>>(endTime - startTime).count();
	Obj_Val = cplex->getObjValue();
	Rel_Gap_100 = 100 * cplex->getMIPRelativeGap();
	NNodes = cplex->getNnodes();
};

void MasterProblem::Output(ofstream& ofc)
{
	Model::Output(ofc);
	ofc << "\t" << FeasCuts ;
	ofc << "\t" << OptCuts ;
	ofc << "\n";
}



DualSubProblem::DualSubProblem(Instance* Ins_In) : Model (Ins_In)
{
	model = new IloModel(env);
	ObjectiveExpr = new IloNumExpr(env);
	Ins = Ins_In;

#pragma region Variables
	a_jt = CreateNumVarArray2(env, Ins->P, Ins->T, "a", -IloInfinity, IloInfinity);
	b_it = CreateNumVarArray2(env, Ins->U, Ins->T, "b", -IloInfinity, 0);
	for (int t = 0; t < Ins->T; ++t)
	{
		for (int j = 0; j < Ins->P; ++j)
		{
			VariableInfo* pInfo = new VariableInfo('a', j, t);
			VariableInfos.emplace_back(pInfo);
			a_jt[j][t].setObject(pInfo);
		}
		for (int i = 0; i < Ins->U; ++i)
		{
			VariableInfo* pInfo = new VariableInfo('b', i, t);
			VariableInfos.emplace_back(pInfo);
			b_it[i][t].setObject(pInfo);
		}
	}

#pragma endregion
#pragma region Objective
	for (int t = 0; t < Ins->T; ++t) //demand part of the objective
		for (int j = 0; j < Ins->P; ++j)
			*ObjectiveExpr += Ins->Products[j].d[t] * a_jt[j][t];
	//the other part of the objective function will be added at the callback
	//by changing the objective coefficients of b_it...

	Objective = new IloObjective(IloMaximize(env, *ObjectiveExpr));
	model->add(*Objective); //adding objective function
#pragma endregion
#pragma region Constraints
	IloExpr expr(env);
	for (int t = 0; t < Ins->T; ++t)
		for (int i = 0; i < Ins->U; ++i)
		{
			expr.clear();
			expr += b_it[i][t];
			for (int j = 0; j < Ins->P; ++j)
				expr += Ins->CPUs[i].Alphas[j] * a_jt[j][t];
			model->add(expr <= Ins->CPUs[i].c[t]);
		}

	for (int t = 0; t < Ins->T; ++t)
		for (int j = 0; j < Ins->P; ++j)
		{
			if(t == Ins->T - 1)
			{
				model->add( - a_jt[j][t] <= Ins->Products[j].h[t]);
				model->add(   a_jt[j][t] <= 2 /*hard coded number 2, instances do not have g values*/);
			}
			else
			{
				model->add(a_jt[j][t + 1] - a_jt[j][t] <= Ins->Products[j].h[t]);
				model->add(a_jt[j][t] - a_jt[j][t + 1] <= 2 /*hard coded number 2, instances do not have g values*/);

			}
		}
#pragma endregion
	cplex = new IloCplex(env);
	cplex->extract(*model);
#pragma endregion
};

DualSubProblem::~DualSubProblem()
{
	delete cplex;
	delete model;
	delete ObjectiveExpr;
	delete Objective;
	for (int i = 0; i < VariableInfos.size(); ++i)
		delete VariableInfos[i];

	env.end();
}

void  DualSubProblem::Solve()
{
	cplex->extract(*model);
	cplex->exportModel("DualSubProblem.lp");

	cplex->setParam(IloCplex::Param::Preprocessing::Reduce, 0);
	cplex->setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

	try
	{
		bool success = cplex->solve();
	}
	catch (IloException e)
	{
		cout << e;
	}
	//high_resolution_clock::time_point startTime = high_resolution_clock::now();
	//bool success = cplex->solve();
	//high_resolution_clock::time_point endTime = high_resolution_clock::now();
	//MIPTime = duration_cast<duration<double>>(endTime - startTime).count();
	////Obj_Val = cplex->getObjValue();
	////Rel_Gap_100 = 100 * cplex->getMIPRelativeGap();
	////NNodes = cplex->getNnodes();
};


