#include "DCCP.h"
#include <windows.h>
#include <algorithm>

#include "Model.h"
#include "BasicModel.h"
#include "BasicModelwithZ.h"
#include "BasicModel_S-Free.h"
#include "ELS1.h"
#include "ELS2.h"
#include "ELS1_LP.h"
#include "ELS1_LP2.h"
#include "ELS2_LP2.h"
#include "Backlog_BasicModel.h"
#include "Backlog_Benders.h"

#include "PatternFitting.h"
#include "PatternFittingBacklog.h"

#include "BasicModelCallback.h"
#include "BasicModelDummy.h"

#include "WagnerWhitin.h"
#include "WagnerWhitinVariable.h"

//#include "Date.h"

int main(int argc,char *argv[])
{
	int modelID = -1;
	string InstancePath;
	string InstanceName;
	string OutputPath = "";
	string instanceID;

	Parameters p;
	//p.StartDate = date::format("%F %T", std::chrono::system_clock::now());
	bool tryMode = false ;
	if (tryMode)
	{
		modelID = 1;
		InstancePath = "13Test.txt";
		InstanceName = "1_0.5_2_3";
		p.TimeLimitMinutes = 5;
		p.threadCount = 1;
		p.EarlyCutAbortion = false;
		p.epsilon = 0.1;
		p.mode = -1; //zaten -1 di en basta

		p.CPWW = false;
		p.CplexCuts = true;
		p.LBonInvVars = 1;
		p.UBonInvVars = 1;
		p.threadCount = 1;

		p.heur = 1;
		p.VI2 = true;
		p.VI1 = true;
		p.WWV = true;

		//WW TRIALS

		//SMALL EXAMPLE IN EVANS PAPER
		/*
		int N = 4;
		Product prod;
		prod.d.resize(N);
		prod.h.resize(N);
		prod.d = { 60,100,140,200};
		prod.h = {1,1,2,2};

		prod.d_t1_t2.resize(N);
		for (int t1 = 0; t1 < N; ++t1)
		{
			prod.d_t1_t2[t1].resize(N, 0);
			for (int t2 = t1; t2 < N; ++t2)
			{
				double temp = 0;
				for (int i = t1; i <= t2; ++i)
				{
					temp += prod.d[i];
				}
				prod.d_t1_t2[t1][t2] = temp;
			}
		}

		CPU cp_u;
		cp_u.p.resize(N);
		cp_u.f.resize(N);
		cp_u.p = {7,7,8,7};
		cp_u.f = { 150,140,160,160};
		CPU cp_u2;
		cp_u2.p.resize(N);
		cp_u2.f.resize(N);
		cp_u2.p = { 7,4,6,7 };
		cp_u2.f = { 200,190,200,210 };
		vector<CPU*> cpuVec;
		cpuVec.push_back(&cp_u);
		cpuVec.push_back(&cp_u2);
		WagnerWhitinVariable WW(N, &prod, cpuVec);
		WW.Compute();
		*/

		Instance MyInstance;
		//string paths[16] = {	"1_1_0.5_2","1_1_0.5_3" ,"1_1_1_2","1_1_1_3" ,"1_2_0.5_2","1_2_0.5_3" ,"1_2_1_2","1_2_1_3" ,
		//						"2_1_0.5_2","2_1_0.5_3" ,"2_1_1_2","2_1_1_3" , "2_2_0.5_2","2_2_0.5_3" ,"2_2_1_2","2_2_1_3" };
		string paths[10] = { "1_0.5_2_3","2_2_3_3" };
		string seeds[10] = { "13","17", "19", "23", "29", "31", "37", "41", "43","47" };
		string folder = "D://YandexDisk//Thesis//CODES//INSTANCES//";
		//MyInstance.ReadInstance("MILS.txt");

		//bu kisim Callback WW Var'able denemesi icin. LB leri denerken kapattim
		MyInstance.ReadInstance(folder + paths[0] + "//" + seeds[0]+".txt");
		Model* MyDebugModel = new BasicModel(&MyInstance, new PatternFitting(MyInstance), &p);
		//p.mode = 99;
		MyDebugModel->Solve();
		double asdasda = MyDebugModel->cplex->getObjValue();
		//bu kisim Callback WW Var'able denemesi icin. LB leri denerken kapattim

		//MyDebugModel->get;
		//12/092019 'nstance yarat'm'n' deg'st'rd'm 'nfeas'ble c'k'yor
		
		//CPU MyCpu;
		//for (int j = 0; j < 4; ++j) {
		//	//simdi CPU Lari teker teker remove edip ayristirma algoritmasi + WW i kontrol edecegim,
		//	//BB icinde olmasi kesin olan adamlar icin yapiyor, olmayacagi kesin olanlari cikarmamiz lazim.
		//	//bence bu kesin bir problem++++
		//	MyCpu.p.resize(MyInstance.T, INT_MAX);
		//	MyCpu.f.resize(MyInstance.T, INT_MAX);
		//	for (int t = 0; t < MyInstance.T; ++t)
		//		for (CPU* i : MyInstance.Products[j].CPUs)
		//		{
		//			if (i->p[t] / i->Alphas[j] < MyCpu.p[t])
		//				MyCpu.p[t] = i->p[t] / i->Alphas[j];
		//			if (i->f[t] < MyCpu.f[t])
		//				MyCpu.f[t] = i->f[t];
		//		}
		//	cout << "MyCPU : \t f : " << MyCpu.f[0] << "\t p : " << MyCpu.p[0] << "\n";
		//	WagnerWhitin WW(MyInstance.T, &(MyInstance.Products[j]), &MyCpu);
		//	WagnerWhitinVariable WWZ(MyInstance.T, &(MyInstance.Products[j]), MyInstance.Products[j].CPUs);
		//	cout << "Product : " << j << "\t" << WW.Compute() << "\t" << WWZ.Compute() << "\n";
		//	cout << "WW Q :\t";
		//	for (int t = 0; t < MyInstance.T; ++t)
		//		cout << WW.Q[t] << "\t";
		//	cout << endl;
		//	cout << "WWV I :\t";
		//	for (int t = 0; t < MyInstance.T; ++t)
		//		cout << WWZ.I[t] << "\t";
		//	cout << endl;
		//	cout << "WWV Q :\t";
		//	for (int t = 0; t < MyInstance.T; ++t)
		//		cout << WWZ.Q[t] << "\t";
		//	cout << endl;
		//}

		vector<vector<int>> LBs;
		vector<vector<int>> LB2s;

		LBs.resize(10);
		LB2s.resize(10);
		std::ofstream fs;
		fs.open("LBCalc.txt",std::ios_base::app);

		for (int ps = 0; ps < 10; ++ps) {
			LBs[ps].resize(10, 0);
			LB2s[ps].resize(10, 0);
			for (int seed = 0; seed < 10; ++seed)
			{
				cout << "PS:" << ps << "seed:" << seed << "\n";
				string path_ = folder + paths[ps] + "//" + seeds[seed] + ".txt";
				MyInstance.ReadInstance(path_);
				vector<vector<CPU*>> cpusOut;
				vector<vector<Product*>> productsOut;
				SeperateGraph(cpusOut, productsOut, MyInstance.CPUs_P);
				if (cpusOut.size() > 1)
					cout << "DEBUG";

				fs << path_ << "\n";
				for (Product j : MyInstance.Products)
				{
					CPU MyCpu;
					MyCpu.p.resize(MyInstance.T, INT_MAX);
					MyCpu.f.resize(MyInstance.T, INT_MAX);
					for (int t = 0; t < MyInstance.T; ++t)
						for (CPU* i : j.CPUs)
						{
							if (i->p[t] / i->Alphas[j.ID - 1] < MyCpu.p[t])
								MyCpu.p[t] = i->p[t] / i->Alphas[j.ID - 1];
							if (i->f[t] < MyCpu.f[t])
								MyCpu.f[t] = i->f[t];
						}
					WagnerWhitin WW(MyInstance.T, &j, &MyCpu);
					
					int result = WW.Compute();
					fs << "WW: Product " << to_string(j.ID) << "\tLB=\t" << to_string(result) << '\n';
					if (result > LBs[ps][seed])
						LBs[ps][seed] = result;

					WagnerWhitinVariable WWV(MyInstance.T, &j, j.CPUs);
					int result2 = WWV.Compute();
					fs << "WWV: Product " << to_string(j.ID) << "\tLB=\t" << to_string(result2) << '\n';
					if (result2 > LB2s[ps][seed])
						LB2s[ps][seed] = result2;
				}
			}
		}
 		fs << "RESULTS\n";
		for (int ps = 0; ps < 10; ++ps)
		{
			LBs[ps].resize(10, 0);
			fs <<"PS" << ps << "WW\t" << "WWV\t" <<"\n";
			for (int seed = 0; seed < 10; ++seed)
				fs <<"\t"<< LBs[ps][seed] <<"\t"<< LB2s[ps][seed] << '\n';
		}
		fs.close();

		//WW TRIALS

		//LB UB TRIALS
		/*
		Instance MyInstances[10];
		string paths[10] = { "1_0.5_2_3","1_0.5_2_6","2_2_3_3","2_2_3_6","2_2_5_3","2_2_5_6","3_2_3_3","3_2_3_6","3_2_5_3","3_2_5_6" };
		string seeds[10] = {"13","17", "19", "23", "29", "31", "37", "41", "43","47"};
		string folder = "C://Users//bahadir//Desktop//CoPro//INSTANCES//";
		for (int i = 0; i < 10; ++i)
			for (int j = 0; j < 10; ++j)
				MyInstances[i].ReadInstance(folder+paths[i]+"//"+seeds[j]+".txt");
		*/
		//LB UB TRIALS
	}
	else
	{
		modelID = stoi(argv[1]);
		InstancePath = argv[2];
		InstanceName = InstancePath.substr(0, InstancePath.find_last_of('\\'));
		InstanceName = InstanceName.substr(InstanceName.find_last_of('\\') + 1);
		std::replace(InstanceName.begin(), InstanceName.end(), ',', '.');

		//burada ayirdigimiz seyleri kullanmiyoruz su anda
		vector<string> delimStrings = split(InstanceName, "_");
		instanceID = delimStrings[0];
		double TDoub = stod(delimStrings[1]);
		double PDoub = stod(delimStrings[2]);
		double UDoub = stod(delimStrings[3]);
		double densityDoub = stod(delimStrings[4]);
		//burada ayirdigimiz seyleri kullanmiyoruz su anda
		
		
		int temp_index = 0;
		while (argv[3 + temp_index * 2])
		{
			switch (argv[3 + temp_index * 2][1])
			{
			case 'm':
				p.mode = stoi(argv[3 + temp_index * 2 + 1]);
				break;
			case 'e':
				p.epsilon = stod(argv[3 + temp_index * 2 + 1]);
				break;
			case 'h':
				p.heur = stoi(argv[3 + temp_index * 2 + 1]);
				break;
			case 'c':
				if (stoi(argv[3 + temp_index * 2 + 1]) == 1)
					p.CPWW = true;
				else
					p.CPWW = false;
				break;
			case 'x':
				if (stoi(argv[3 + temp_index * 2 + 1]) == 1)
					p.CplexCuts = true;
				else
					p.CplexCuts = false;
				break;
			case 't':
				p.TimeLimitMinutes = stoi(argv[3 + temp_index * 2 + 1]);
				break;
			case 'a':
				if (stoi(argv[3 + temp_index * 2 + 1]) == 1)
					p.EarlyCutAbortion = true;
				else
					p.EarlyCutAbortion = false;
				break;
			case 'l':
				p.LBonInvVars = stoi(argv[3 + temp_index * 2 + 1]);
				break;
			case 'u':
				p.UBonInvVars = stoi(argv[3 + temp_index * 2 + 1]);
				break;
			case 'r':
				p.threadCount = stoi(argv[3 + temp_index * 2 + 1]);
				break;
			case 'w':
				if (stoi(argv[3 + temp_index * 2 + 1]) == 1)
					p.WWV = true;
				else
					p.WWV = false;
				break;
			case '1':
				if (stoi(argv[3 + temp_index * 2 + 1]) == 1)
					p.VI2 = true;
				else
					p.VI2 = false;
				break;
			case '2':
				if (stoi(argv[3 + temp_index * 2 + 1]) == 1)
					p.VI1 = true;
				else
					p.VI1 = false;
				break;
			case 's':
				if (stoi(argv[3 + temp_index * 2 + 1]) == 1)
					p.singleCut = true;
				else
					p.singleCut = false;
				break;
			case 'o':
				OutputPath = argv[3 + temp_index * 2 + 1];
				break;
			}
			temp_index++;
		}
	}
	if(OutputPath == "")
		OutputPath = "D://YandexDisk//Thesis//CODES//OUTPUTS//";
	OutputPath += "MID";
	OutputPath += (modelID < 10 ? "0" : "");
	OutputPath += to_string(modelID) + "_";

	OutputPath += "PS";
	OutputPath += stoi(instanceID) < 10 ? "0" : "";
	OutputPath += instanceID + "_";

	OutputPath += "TL";
	OutputPath += to_string(p.TimeLimitMinutes) + "_";

	OutputPath += "T0" + to_string(p.threadCount) + "_";

	OutputPath += "M" + to_string(p.mode);
	OutputPath += "_";

	OutputPath += (p.CplexCuts ? "CC1_" : "CC0_");

	OutputPath += "E" + to_string(p.epsilon);
	OutputPath += "_";

	OutputPath += "PF" + to_string(p.heur);
	OutputPath += "_";

	OutputPath += (p.CPWW ? "WW1" : "WW0");
	OutputPath += "_";

	OutputPath += (p.EarlyCutAbortion ? "ECA1" : "ECA0");
	OutputPath += "_";

	OutputPath += "LB" + to_string(p.LBonInvVars);
	OutputPath += "_";

	OutputPath += "UB" + to_string(p.UBonInvVars);
	OutputPath += "_";

	OutputPath += "VI1" + to_string(p.VI2);
	OutputPath += "_";

	OutputPath += "SC" + to_string(p.singleCut);
	OutputPath += "_";

	OutputPath += "VI2" + to_string(p.VI1);
	OutputPath += "_";

	OutputPath += "WWV" + to_string(p.WWV);
	OutputPath += ".txt";

	Instance MyInstance;
	MyInstance.ReadInstance(InstancePath);
	Model* MyModel;

	Heuristic* heurMethod;
	switch (p.heur)
	{
	case 0:
		heurMethod = NULL;
		break;
	case 1:
		if (modelID == 30)
			heurMethod = new PatternFittingBacklog(MyInstance);
		else
			heurMethod = new PatternFitting(MyInstance);
		break;
	default:
		heurMethod = NULL;
		break;
	}

	switch (modelID)
	{
	case 99:
		MyModel = new BasicModelCallback(&MyInstance, heurMethod, &p);
		break;
	case 98:
		MyModel = new BasicModelDummy(&MyInstance, heurMethod, &p);
		break;
	case 1:
		MyModel = new BasicModel(&MyInstance, heurMethod, &p);
		break;
	case 2:
		MyModel = new SFreeModel(&MyInstance, heurMethod, &p);
		break;
	case 3:
		MyModel = new ELS1(&MyInstance, heurMethod, &p);
		break;
	case 4:
		MyModel = new ELS2(&MyInstance, heurMethod, &p);
		break;
	case 5:
		MyModel = new BasicModelwithZ(&MyInstance, heurMethod, &p);
		break;
	case 10:
		MyModel = new ELS1_LP(&MyInstance, heurMethod, &p);
		break;
	case 12:
		MyModel = new ELS1_LP2(&MyInstance, heurMethod, &p);
		break;
	case 13:
		MyModel = new ELS2_LP2(&MyInstance, heurMethod, &p);
		break;
	case 30:
		MyModel = new Backlog_BasicModel(&MyInstance, heurMethod, &p);
		break;
	case 31:
		MyModel = new MasterProblem(&MyInstance);
		break;
	}
	MyModel->Solve();
	//std::wstring stemp = std::wstring(_GetDirectoryName(OutputPath).begin(), _GetDirectoryName(OutputPath).end());
	//LPCWSTR sw = stemp.c_str();
	if (tryMode || CreateDirectory(_GetDirectoryName(OutputPath).c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError())
	{
		std::ofstream fs;
		fs.open(OutputPath, std::ios_base::app);
		fs << modelID << "\t"
			<< instanceID << "\t"
			<< p.TimeLimitMinutes << "\t"
			<< p.threadCount << "\t"
			<< p.mode << "\t"
			<< p.CplexCuts << "\t"		
			<< p.CPWW << "\t"
			<< p.EarlyCutAbortion << "\t"
			<< p.LBonInvVars << "\t"
			<< p.UBonInvVars << "\t"
			<< p.epsilon << "\t"
			<< p.VI2 << "\t"
			<< p.singleCut << "\t"
			<< p.VI1 << "\t"
			<< p.heur << "\t"
			<< p.WWV << "\t"
			;
		MyModel->Output(fs);
		fs.close();
	}

	delete MyModel;
	return 0;
}