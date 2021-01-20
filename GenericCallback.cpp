#include "GenericCallback.h"
#include "WagnerWhitin.h"
#include "WagnerWhitinVariable.h"

#include "DCCP.h"
#include "date.h"



GenericCallback::GenericCallback(Model* _ref, bool CPWW_In) { CPWW = CPWW_In; ref = _ref; }

void GenericCallback::invoke(const IloCplex::Callback::Context& context)
{
	int CurrentThreadId = context.getIntInfo(IloCplex::Callback::Context::Info::ThreadId);
	if (ref->P->EarlyCutAbortion) {
		if (CurrentNodes[CurrentThreadId] == context.getIntInfo(IloCplex::Callback::Context::Info::NodeCount))
			return;
		else
			CurrentNodes[CurrentThreadId] = context.getIntInfo(IloCplex::Callback::Context::Info::NodeCount);
	}

	NumArray2 xSol = CreateNumArray2(ref->env, ref->Ins->U, ref->Ins->T);
	NumArray2 ySol = CreateNumArray2(ref->env, ref->Ins->U, ref->Ins->T);
	IloNumArray wSol(ref->env, ref->Ins->U);
	vector<CPU*> currentCPUs;
	string removedCPUs; //for debugging
	vector<vector<CPU*>> cpusOut;
	vector<vector<Product*>> productsOut;
	CPU MyCpu;
	double calculatedLB = 0;
	double setLB = 0;

	vector<bool> S(ref->Ins->T);
	IloExpr expr(ref->env);
	WW* WWHeur;

	IloNumVarArray startVar(ref->env);
	IloNumArray startVal(ref->env);

	//string filename = "debug_" + to_string(ref->P->CplexCuts)+"_"+ to_string(ref->P->TimeLimitMinutes) + "_" + ref->P->StartDate + ".txt";
	//std::replace(filename.begin(), filename.end(), ':', '-');

	// Get the current x solution
	std::ofstream fs;
	IloRange bestCut;
	bool bestCutFlag = false;
	double bestCutDiff = 0;
	switch (context.getId()) {
	case IloCplex::Callback::Context::Id::Candidate:
			if (!context.isCandidatePoint()) // The model is always bounded
				throw IloCplex::Exception(-1, "Unbounded solution");
			for (IloInt i = 0; i < ref->Ins->U; ++i) {
				context.getCandidatePoint(ref->Xit[i], xSol[i]);
				context.getCandidatePoint(ref->Yit[i], ySol[i]);
			}
			break;
	case IloCplex::Callback::Context::Id::Relaxation:
		ref->callbackCounter++;//count the number of times this callback is called
		for (IloInt i = 0; i < ref->Ins->U; ++i) {
			context.getRelaxationPoint(ref->Xit[i], xSol[i]);
			context.getRelaxationPoint(ref->Yit[i], ySol[i]);
		}
		if(ref->P->WWV){
			//Basic model callback burada sýkýntý verýcek
			//context.getRelaxationPoint(ref->Wi, wSol);
			for (IloInt i = 0; i < ref->Ins->U; ++i)
				if (!(context.getLocalUB(ref->Wi[i]) < 0.01))
					currentCPUs.push_back(ref->Ins->CPUs_P[i]);
				else
				{
					removedCPUs += std::to_string(i);
					removedCPUs += ",";
				}
			//burada kesin olacak olan CPU larla ilgileniyoruz, aslinda kesin olmayacak olanlari cikarmamiz lazim

			SeperateGraph(cpusOut, productsOut, currentCPUs);
			calculatedLB = 0;
			if (currentCPUs.size() < 20)
				setLB = 0;
			setLB = 0;
			//debugging
			//if (cpusOut.size() > 1)
			//	cout << "DEBUGGING";
			//debugging
			for (int k = 0; k < productsOut.size(); ++k)//for each set of disjoint product set
			{
				setLB = 0;
				for (Product* j : productsOut[k])
				{
					WagnerWhitinVariable WW(ref->Ins->T, j, j->CPUs); //currentCPUs[j] ?
					int result = WW.Compute();
					if (result > setLB)
						setLB = result;
				}
				calculatedLB += setLB;
			}
			//fs.open(filename, std::ios_base::app);
			//fs << context.getIntInfo(IloCplex::Callback::Context::Info::NodeCount) << "\t"
			//	<< removedCPUs << "\t"
			//	<< cpusOut.size() << "\t"
			//	<< productsOut.size() << "\t"
			//	<< calculatedLB << "\t"
			//	<< context.getRelaxationObjective() << "\n";
			//fs.close();
		}

		if (CPWW) {
			WWHeur = new WW(*ref->Ins);

			for (int g = 0; g < ref->Ins->P; ++g) {
				WWHeur->ImplementAlgorithm(*(ref->Ins), ySol);

				for (int t = 0; t < ref->Ins->T; ++t) {
					for (int i = 0; i < ref->Ins->U; ++i)
					{
						startVar.add(ref->Yit[i][t]);
						startVal.add(WWHeur->GetYvalue(i, t));
					}
				}
				context.postHeuristicSolution(startVar, startVal, WWHeur->GetObj(), IloCplex::Callback::Context::SolutionStrategy::Propagate);
				ref->heurSolCounter++;
				double asdasd = WWHeur->GetObj();
				if (WWHeur->GetObj() < ref->heurSolBestObj)
					ref->heurSolBestObj = WWHeur->GetObj();
			}
			delete WWHeur;
		}

		if (ref->P->VI2) {
			for (int l = ref->Ins->T - 1; l >= 0; --l)
			{
				for (int j = 0; j < ref->Ins->P; ++j)
				{
					for (int t = 0; t < ref->Ins->T; ++t) //clear the array
						S[t] = false;
					double sum = 0;
					for (int q = 0; q <= l; ++q)
					{
						IloNum sumX = 0;
						for (auto i : ref->Ins->Products[j].CPUs)
							sumX += i->Alphas[j] * xSol[i->ID - 1][q];

						IloNum sumY = 0;
						for (CPU* c : ref->Ins->Products[j].CPUs)
							sumY += ref->Ins->Products[j].d_t1_t2[q][l] * ySol[c->ID - 1][q];

						if (S[q] = (sumY < sumX))
							sum += sumY;
						else
							sum += sumX;
					}

					if (ref->P->epsilon < ref->Ins->Products[j].d_t1_t2[0/*zero*/][l] - sum)
					{
						expr.clear();
						for (int q = 0; q <= l; ++q)
						{
							if (S[q]) //q in set S
							{
								for (CPU* c : ref->Ins->Products[j].CPUs)
									expr += ref->Ins->Products[j].d_t1_t2[q][l] * ref->Yit[c->ID - 1][q];
							}
							else	//q not in set S
							{
								for (auto i : ref->Ins->Products[j].CPUs)
									expr += i->Alphas[j] * ref->Xit[i->ID - 1][q];
							}
						}

						if (ref->P->singleCut)
						{
							if (bestCutDiff < ref->Ins->Products[j].d_t1_t2[0/*zero*/][l] - sum)
							{
								bestCutDiff = ref->Ins->Products[j].d_t1_t2[0/*zero*/][l] - sum;
								bestCut = expr >= ref->Ins->Products[j].d_t1_t2[0][l];
								bestCutFlag = true;
							}
						}
						else
						{
							context.addUserCut(expr >= ref->Ins->Products[j].d_t1_t2[0][l], IloCplex::UseCutPurge, IloFalse); // add the cut
							ref->userCutCounter++;//count the number of constraints added
						}
					}
				}
			}
			if (ref->P->singleCut && bestCutFlag)
			{
				context.addUserCut(bestCut, IloCplex::UseCutPurge, IloFalse); // add the cut
				ref->userCutCounter++;//count the number of constraints added
			}
		}
		break;
	default:
		throw IloCplex::Exception(-1, "unexpected contextId");
	}
	for (IloInt i = 0; i < ref->Ins->U; ++i)
		xSol[i].end();
	xSol.clear();
	xSol.end();
	for (IloInt i = 0; i < ref->Ins->U; ++i)
		ySol[i].end();
	ySol.clear();
	ySol.end();
	wSol.clear();
	wSol.end();
	expr.end();
	S.clear();
	currentCPUs.clear();
	for (int i = 0; i < cpusOut.size(); i++)
		cpusOut[i].clear();
	cpusOut.clear();
	for (int i = 0; i < productsOut.size(); i++)
		productsOut[i].clear();
	productsOut.clear();
	startVal.end();
	startVar.end();

	currentCPUs.end();
	cpusOut.end();
	productsOut.end();
	bestCut.end();
}
