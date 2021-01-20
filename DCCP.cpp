#include "DCCP.h"
vector<string> split(string str, string token) {
	vector<string>result;
	while (str.size()) {
		int index = str.find(token);
		if (index != string::npos) {
			result.push_back(str.substr(0, index));
			str = str.substr(index + token.size());
			if (str.size() == 0)result.push_back(str);
		}
		else {
			result.push_back(str);
			str = "";
		}
	}
	return result;
}

bool sortfunction(Product* p1, Product* p2)
{ return ((p1-> remainingDemand / p1->currentCPU->Alphas[p1->ID - 1]) > (p2->remainingDemand / p2->currentCPU->Alphas[p2->ID - 1])); }
// An optimized version of Bubble Sort 
//void bubbleSort(vector<Product*> arr, vector<double> alphas)
//{
//	int i, j;
//	bool swapped;
//	for (i = 0; i < arr.size() - 1; i++)
//	{
//		swapped = false;
//		for (j = 0; j < arr.size() - i - 1; j++)
//		{
//			if ((arr[j]->remainingDemand / alphas[(arr[j]->ID)-1]) <  (arr[j + 1]->remainingDemand / alphas[(arr[j + 1]->ID) - 1]))
//			{
//				std::iter_swap(arr[j], arr[j + 1]);
//				swapped = true;
//			}
//		}
//
//		// IF no two elements were swapped by inner loop, then break 
//		if (swapped == false)
//			break;
//	}
//}

ArrString readDelimated(string line)
{
	ArrString result;
	int index = 0;
	if (line[0] == 'p')
		line = line;
	while (index != string::npos)
	{
		result.push_back(line.substr(0, index = line.find_first_of(" \t\n")));
		line = line.substr(index + 1);
	}
	return result;
}

string _GetDirectoryName(string str)
{
	int index = str.find_last_of("/");
	return str.substr(0, index + 1);
}

void Instance::ReadInstance(string PathIn)
{
	Products.clear();
	CPUs.clear();
	Path = PathIn;

	ifstream myfile(Path);
	if (myfile.is_open())
	{
		string line;
		ArrString contents;
		//char delim[3] = { ' ', '\t', '\n' };

		while (getline(myfile, line))
		{
			replace(line.begin(), line.end(), ',', '.');
			contents = readDelimated(line);
			switch (contents[0][0])
			{
			case 'c':
				break;
			case 'p':
				T = stoi(contents[1]);
				U = stoi(contents[2]);
				P = stoi(contents[3]);
				Products.resize(P);
				CPUs.resize(U);
				CPUs_P.resize(U);
				//give them IDs
				for (int i = 0; i < U; ++i) { //fill pointer CPU array
					CPUs[i].ID = i + 1;
					CPUs_P[i] = &CPUs[i];
				}
				for (int i = 0; i < P; ++i)
					Products[i].ID = i + 1;
				break;
			case 'a':
				CPUs[stoi(contents[1]) - 1].Alphas.resize(P);
				for (int i = 0; i < P; ++i) {
					CPUs[stoi(contents[1]) - 1].Alphas[i] = stod(contents[2 + i]);
					if (stod(contents[2 + i]) > 0)
					{
						Products[i].CPUs.push_back(&CPUs[stoi(contents[1]) - 1]);
						CPUs[stoi(contents[1]) - 1].Products.push_back(&Products[i]);
					}
				}
				break;
			case 'f':
				CPUs[stoi(contents[1]) - 1].f.resize(T);
				for (int i = 0; i < T; ++i)
					CPUs[stoi(contents[1]) - 1].f[i] = stod(contents[2 + i]);
				break;
			case 'v':
				CPUs[stoi(contents[1]) - 1].p.resize(T);
				for (int i = 0; i < T; ++i)
					CPUs[stoi(contents[1]) - 1].p[i] = stod(contents[2 + i]);
				break;
			case 'd':
				Products[stoi(contents[1]) - 1].d.resize(T);
				for (int i = 0; i < T; ++i)
					Products[stoi(contents[1]) - 1].d[i] = stod(contents[2 + i]);
				break;
			case 'h':
				Products[stoi(contents[1]) - 1].h.resize(T);
				for (int i = 0; i < T; ++i)
					Products[stoi(contents[1]) - 1].h[i] = stod(contents[2 + i]);
				break;
			}
		}
		myfile.close();
		//calculate d_t_T
		for (int j = 0; j < P; ++j) //for (Product j : Products) create a copy, does not actually change the object
		{
			Products[j].d_t_T.resize(T);
			double temp = 0;
			for (int i = T - 1; i >= 0; --i)
			{
				temp += Products[j].d[i];
				Products[j].d_t_T[i] = temp;
			}
		}
		//calculate d_t_T
		//calculate d_t1_t2
		for (int j = 0; j < P; ++j)
		{
			Products[j].d_t1_t2.resize(T);
			for (int t1 = 0; t1 < T; ++t1)
			{
				Products[j].d_t1_t2[t1].resize(T, 0);
				for (int t2 = t1; t2 < T; ++t2)
				{
					double temp = 0;
					for (int i = t1; i <= t2; ++i)
					{
						temp += Products[j].d[i];
					}
					Products[j].d_t1_t2[t1][t2] = temp;
				}
			}
		}
		//calculate d_t1_t2
		//calculate j.c[t]
		for (int i = 0; i < U; ++i)
		{
			CPUs[i].c.resize(T);
			double h_alpha = 0;
			for (int t = T - 1; t >= 0; --t)
			{
				for (Product* j : CPUs[i].Products)
					h_alpha += j->h[t] * CPUs[i].Alphas[j->ID - 1];
				CPUs[i].c[t] = CPUs[i].p[t] + h_alpha;
			}
		}
		//calculate j.c[t]
		//calculate Lower and Upper bounds on inventory variables
		//calculate LB
		for (int j = 0; j < P; ++j)
		{
			Products[j].LB.resize(T);
			for (int t = 0; t < T; ++t)
			{
				//nonnegative max is omitted, will be dealed with later
				double maxJ = 0;
				for (int j1 = 0; j1 < P; ++j1)
				{
					double minI = INT32_MAX;
					for(CPU* i :Products[j1].CPUs)
						if (minI > Products[j1].d_t1_t2[0][t] * i->Alphas[j] / i->Alphas[j1])
							minI = Products[j1].d_t1_t2[0][t] * i->Alphas[j] / i->Alphas[j1];
					if (maxJ < minI - Products[j].d_t1_t2[0][t])
						maxJ = minI - Products[j].d_t1_t2[0][t];
				}
				if (maxJ > 0)
					Products[j].LB[t] = maxJ;
				else
					Products[j].LB[t] = 0;
			}
		}
		//calculate LB
		//calculate UB
		for (int j = 0; j < P; ++j)
		{
			Products[j].UB.resize(T);
			for (int t = 0; t < T; ++t)
			{
				double maxJI = 0;
				for (int j1 = 0; j1 < P; ++j1)
					for (CPU* i : Products[j1].CPUs)
						if (maxJI < Products[j1].d_t1_t2[0][T - 1] * i->Alphas[j] / i->Alphas[j1])
							maxJI = Products[j1].d_t1_t2[0][T - 1] * i->Alphas[j] / i->Alphas[j1];
				Products[j].UB[t] = maxJI- Products[j].d_t1_t2[0][t];
			}
		}
		//calculate UB
		/*
		std::ofstream fs;
		fs.open("Bounds.txt", std::ios_base::app);
		fs << Path << "\n";
		for (int t = 0; t < T; ++t)
			for (int j = 0; j < P; ++j)
				if(Products[j].LB[t] > 0)
					fs << "L_" << to_string(j + 1) << "_" << to_string(t + 1) << "=\t" << to_string(Products[j].LB[t]) 
					<<"\t"<< "U_" << to_string(j + 1) << "_" << to_string(t + 1) << "=\t" << to_string(Products[j].UB[t]) <<"\n";
		fs.close();
		*/
		//calculate Lower and Upper bounds on inventory variables
	}
}

template < typename T>
std::pair<bool, int > findInVector(const std::vector<T>  & vecOfElements, const T  & element)
{
	std::pair<bool, int > result;

	// Find given element in vector
	auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);

	if (it != vecOfElements.end())
	{
		result.second = distance(vecOfElements.begin(), it);
		result.first = true;
	}
	else
	{
		result.first = false;
		result.second = -1;
	}

	return result;
}

void SeperateGraph(vector<vector<CPU*>> &cpuOutput, vector<vector<Product*>> &productOutput, vector<CPU*> &CPUss)
{
	cpuOutput.resize(1);
	productOutput.resize(1);
	int CurrentSet = 0;
	vector<CPU*> cpuCopyList;
	cpuCopyList.assign((CPUss).begin(),(CPUss).end());
	cpuOutput[CurrentSet].push_back(cpuCopyList[0]);
	for (Product* p : cpuCopyList[0]->Products)
		productOutput[CurrentSet].push_back(p);
	cpuCopyList.erase(cpuCopyList.begin());

	bool cpusUpdated = false;
	bool productsUpdated = true;
	
	while (cpuCopyList.size()>0)
	{
		while (cpusUpdated || productsUpdated)
		{
			if (productsUpdated)
				for (Product* p : productOutput[CurrentSet])
					for (int i = 0 ; i < cpuCopyList.size();++i)
						if(findInVector(cpuCopyList[i]->Products, p).first)
						{
							if(!findInVector(cpuOutput[CurrentSet],cpuCopyList[i]).first)
							{
								cpuOutput[CurrentSet].push_back(cpuCopyList[i]);
								cpusUpdated = true;
								cpuCopyList.erase(cpuCopyList.begin()+i);
								i--;
							}
						}
			productsUpdated = false;
			if (cpusUpdated)
				for (CPU* c : cpuOutput[CurrentSet])
					for (Product* p : c->Products)
					{
						if (!findInVector(productOutput[CurrentSet], p).first)
						{
							productOutput[CurrentSet].push_back(p);
							productsUpdated = true;
						}
					}
			cpusUpdated = false;
		}
		if (cpuCopyList.size() > 0)
		{
			CurrentSet++;
			cpuOutput.resize(CurrentSet + 1);
			productOutput.resize(CurrentSet + 1);		
			cpuOutput[CurrentSet].push_back(cpuCopyList[0]);
			for (Product* p : cpuCopyList[0]->Products)
				productOutput[CurrentSet].push_back(p);
			cpuCopyList.erase(cpuCopyList.begin());
			cpusUpdated = false;
			productsUpdated = true;
		}
	}
	cpuCopyList.end();
}
