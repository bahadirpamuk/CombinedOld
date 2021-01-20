#pragma once
#include "DCCP.h"
#include <windows.h>
#include <cmath>
#include <atomic>
#include "Heuristic.h"
#define FracCutAbort 0.002
class Model {
protected:
	
	IloModel* model;
	IloModel* LPRLX;
	IloConstraintArray* AddedCuts;
	Heuristic* Heur;

public:
	IloCplex * cplex;
	NumVarArray2 Xit;
	BoolVarArray2 Yit;
	NumVarArray2 Sjt;
	IloBoolVarArray Wi;

	Instance* Ins;
	IloEnv env;

	Parameters* P;

	double LPTime = 0;
	double MIPTime = 0;
	double Node0Time = 0;
	double HeurTime = 0;
	double FTime = 0;
	double Rlx_Obj_Val = -1;
	double Heur_Obj_Val = -1;
	double Node0_Rlx_Obj_Val = -1;
	double Node0_Feas_Soln = -1;
	double Obj_Val = -1;
	double Rel_Gap_100 = -1;
	int Org_NConst = -1;
	int NNodes = -1;

	double cutTime = 0;
	int Total_NConst = -1;
	atomic_int userCutCounter = -1;
	atomic_int callbackCounter = -1;
	atomic_int heurSolCounter = -1;
	atomic_int heurSolBestObj = -1;

	Model(Instance* Ins_In, Heuristic* Heur_In, Parameters* pIn) {
		Ins = Ins_In;
		Heur = Heur_In;
		P = pIn;
	};
	Model(Instance* Ins_In) {
		Ins = Ins_In;
	};

	virtual void Solve() {
		//we come here as cplex being extracted MIP model
		//HEURISTIC
		if (Heur != NULL && !(P->littleModel))
		{
			Heur->ImplementAlgorithm(*Ins);
			IloNumVarArray startVar(env);
			IloNumArray startVal(env);
			for (int t = 0; t < Ins->T; ++t) {
				for (int i = 0; i < Ins->U; ++i)
				{
					startVar.add(Xit[i][t]);
					startVal.add(Heur->GetXvalue(i, t));
					startVar.add(Yit[i][t]);
					startVal.add(Heur->GetYvalue(i, t));
				}
			}

			cplex->addMIPStart(startVar, startVal);
			Heur_Obj_Val = Heur->GetObj();
			HeurTime = Heur->GetTime();
			startVal.end();
			startVar.end();
		}
		//HEURISTIC

		cplex->setParam(IloCplex::Param::MIP::Limits::Nodes, 0); //solve only node 0
		cplex->setParam(IloCplex::Param::TimeLimit, ((P->TimeLimitMinutes*60) - FTime - LPTime - cutTime - HeurTime) > 0 ? ((P->TimeLimitMinutes * 60) - FTime - LPTime - cutTime - HeurTime) : 0);//set time limit to Timelimit
		cplex->setParam(IloCplex::Param::Threads, P->threadCount);					//thread safe ??
		if (!P->CplexCuts)
		{
			cplex->setParam(IloCplex::Param::MIP::Limits::CutsFactor, 1);   //no cuts, new command replaces the two below..
			//cplex->setParam(IloCplex::Param::MIP::Limits::EachCutLimit, 0); //do not generate cplex cuts except gomory cuts
			//cplex->setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);		//do not generate gomory cuts
		}
		high_resolution_clock::time_point startTime = high_resolution_clock::now();
		try { bool success = cplex->solve(); }
		catch (IloException e) { cout << e; }	
		high_resolution_clock::time_point endTime = high_resolution_clock::now();
		Node0Time = duration_cast<duration<double>>(endTime - startTime).count();//getNode0Time
		Node0_Rlx_Obj_Val = cplex->getBestObjValue(); //getNode0RelaxedObj
		//Node0_Rlx_Obj_Val = cplex->getObjValue(); //BAHADIR bunu KALDIR TRY MODE
		try{
		Node0_Feas_Soln = cplex->getObjValue(); //node0 da eðer bulduysa feas soln.
		}
		catch (IloException e) { cout << e; }
		cplex->setParam(IloCplex::Param::MIP::Limits::Nodes, 9223372036800000000); //solve all nodes
		cplex->setParam(IloCplex::Param::TimeLimit, ((P->TimeLimitMinutes * 60) - FTime - LPTime - cutTime - HeurTime - Node0Time) > 0 ? ((P->TimeLimitMinutes * 60) - FTime - LPTime - cutTime - HeurTime - Node0Time) : 0);				//set time limit to 60min
		cplex->setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0);	//force to solving to optimality
		startTime = high_resolution_clock::now();
		try { bool success = cplex->solve(); }
		catch (IloException e) { cout << e; }
		endTime = high_resolution_clock::now();
		MIPTime = duration_cast<duration<double>>(endTime - startTime).count();
		try { Obj_Val = cplex->getObjValue(); }
		catch (IloException e) { Obj_Val = -1; }
		Rel_Gap_100 = 100 * cplex->getMIPRelativeGap();
		NNodes = cplex->getNnodes();
	};

	void LPSolve() {
		cplex->extract(*LPRLX); // get the relaxed model
		high_resolution_clock::time_point startTime = high_resolution_clock::now();
		bool success = cplex->solve(); //solveRelaxedModel
		high_resolution_clock::time_point endTime = high_resolution_clock::now();
		LPTime = duration_cast<duration<double>>(endTime - startTime).count();  //getLPTime
		cplex->exportModel("LPExport.lp");
		Rlx_Obj_Val = cplex->getObjValue(); //getRelaxedObj
		Org_NConst = cplex->getNrows();		//getoriginal Constraint count	
	};

	virtual void AddCuts() {
		AddedCuts = new IloConstraintArray(env); //collect cuts to add
		//deleted in every sub class destructors

		cplex->setParam(IloCplex::Param::RootAlgorithm, CPX_ALG_DUAL); //use dual simplex in generating cuts

		if (P->epsilon >= 0) // then we need to add cuts...
		{
			high_resolution_clock::time_point startTime = high_resolution_clock::now();
			IloConstraintArray constraints;

			if (P->mode == 0 || P->mode == 1 || P->mode == 2)
				constraints = FindUserCuts();
			else if (P->mode == 10 || P->mode == 11 || P->mode == 12)
				constraints = FindSingleUserCut();

			high_resolution_clock::time_point endTime = high_resolution_clock::now();

			cutTime = duration_cast<duration<double>>(endTime - startTime).count();  //getLPTime

			while (constraints.getSize() != 0)
			{
				LPRLX->add(constraints);
				AddedCuts->add(constraints);
				cplex->exportModel("LPExport.lp");
				startTime = high_resolution_clock::now();
				cplex->solve();
				endTime = high_resolution_clock::now();
				LPTime += duration_cast<duration<double>>(endTime - startTime).count();  //getLPTime
				
				//FindSolution();
				if (cplex->getCplexStatus() == IloCplex::CplexStatus::AbortTimeLim)
					break;
				if (1 - (Rlx_Obj_Val / cplex->getObjValue()) < FracCutAbort)
					break;	//break if specified fraction of improvement is not met
				Rlx_Obj_Val = cplex->getObjValue(); //ADD RESTRICTIONS TO CUT TIMING

				startTime = high_resolution_clock::now();
				if (P->mode == 0 || P->mode == 1 || P->mode == 2)
					constraints = FindUserCuts();
				else if (P->mode == 10 || P->mode == 11 || P->mode == 12)
					constraints = FindSingleUserCut();
				endTime = high_resolution_clock::now();
				cutTime += duration_cast<duration<double>>(endTime - startTime).count();  //getLPTime
			}
		}
		try{ Rlx_Obj_Val = cplex->getObjValue();}

		catch(IloException e)
		{ /*Rlx_Obj_Val = cplex->getBestObjValue();*/ }
		
		cplex->setParam(IloCplex::Param::RootAlgorithm, CPX_ALG_AUTOMATIC); //back to auto after cut generation is done

		if (P->mode != 0 && P->mode != 10 && P->epsilon != -1) {
			if (P->mode == 1 || P->mode == 11)
			{
				(*model).add(*AddedCuts);
				cplex->extract(*model); // get originalModel		
			}
			else if (P->mode == 2 || P->mode == 12)
			{
				cplex->extract(*model); // get originalModel
				cplex->addUserCuts(*AddedCuts);
			}
			Total_NConst = cplex->getNrows();
		}
	};

	virtual IloConstraintArray FindUserCuts() {IloConstraintArray a; return a;};

	virtual IloConstraintArray FindSingleUserCut() { IloConstraintArray a; return a;};

	void FindSolution() {
		////make a copy of LPRLX model
		////loop
		//	//search over Y vars, checking if they are within tolerances
		//	//find the maximum (closest to 1) value among them
		//	//set it to 1
		//	//resolve the model
		////loop
		//high_resolution_clock::time_point startTime = high_resolution_clock::now();

		//IloModel TMPMDL(env);
		//TMPMDL.add(*LPRLX);
		//IloCplex TMPCPX(TMPMDL);
		//TMPCPX.solve(); //it seems possible to get the solution or basis info from existing cplex
		//				//but i dont know how yet
		//		//this is integrality tolerance for original cplex object
		//double toler = cplex->getParam(IloCplex::Param::MIP::Tolerances::Integrality);
		////this is integrality tolerance for original cplex object

		//for (int t = 0; t < Ins->T; ++t)
		//	for (int i = 0; i < Ins->U; ++i)
		//		if (TMPCPX.getValue(Yit[i][t]) > toler)
		//			cout << "Y_" << i << "_" << t << " = " << TMPCPX.getValue(Yit[i][t]) << "\n";

		//bool flag = true;
		//while (flag)
		//{
		//	flag = false;
		//	double max = -1;
		//	IloNumVar MaxVar;
		//	for (int t = 0 ; t < Ins->T ; ++ t)
		//		for (int i = 0 ; i < Ins->U ; ++i)
		//		{
		//			double temp = TMPCPX.getValue(Yit[i][ t]);
		//			if (!(abs(0 - temp) <= toler || abs(1 - temp) <= toler))
		//			{
		//				flag = true;
		//				if (temp > max)
		//				{
		//					max = temp;
		//					MaxVar = Yit[i][t];
		//				}
		//			}
		//		}
		//	if (flag == true) {
		//		TMPMDL.add(1 * MaxVar == 1);
		//		TMPCPX.solve();
		//	}
		//	for (int t = 0; t < Ins->T; ++t)
		//		for (int i = 0; i < Ins->U; ++i)
		//			if (TMPCPX.getValue(Yit[i][t]) > toler)
		//				cout << "Y_" << i << "_" << t << " = " << TMPCPX.getValue(Yit[i][t]) << "\n";
		//}
		//for (int t = 0; t < Ins->T; ++t)
		//	for (int i = 0; i < Ins->U; ++i)
		//		if (TMPCPX.getValue(Yit[i][t]) > toler)
		//			cout << "Y_" << i << "_" << t << " = " << TMPCPX.getValue(Yit[i][t]) << "\n";
		//if (TMPCPX.getObjValue() < heur_best)
		//	heur_best = TMPCPX.getObjValue();
		//TMPMDL.end();
		//TMPCPX.end();
		//HeurTime += duration_cast<duration<double>>(high_resolution_clock::now() - startTime).count();//getNode0Time
	};

	virtual void Output(std::ofstream& ofIn)
	{
			ofIn << FTime << "\t"
			<< LPTime << "\t"
			<< Rlx_Obj_Val << "\t"
			<< HeurTime << "\t"
			<< Heur_Obj_Val << "\t"
			<< Node0Time << "\t"
			<< Node0_Rlx_Obj_Val << "\t"
			<< Node0_Feas_Soln << "\t"
			<< MIPTime << "\t"
			<< Obj_Val << "\t"
			<< Rel_Gap_100 << "\t"
			<< NNodes << "\t"
			<< Org_NConst << "\t"
			<< cutTime << "\t"
			<< Total_NConst << "\t"
			<< userCutCounter << "\t"
			<< callbackCounter << "\t"
			<< heurSolCounter << "\t"
			<< heurSolBestObj;
	};
};