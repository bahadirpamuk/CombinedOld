#pragma once
#include "Common.h"
#include <chrono>
using namespace std::chrono;

vector<string> split(string str, string token);

template < typename T>
std::pair<bool, int > findInVector(const std::vector<T>  & vecOfElements, const T  & element);

struct Product;

struct CPU;

void SeperateGraph(vector<vector<CPU*>> &cpuOutput, vector<vector<Product*>> &productOutput, vector<CPU*> &CPUs);

bool sortfunction(Product* p1, Product* p2);

typedef vector<string> ArrString;

ArrString readDelimated(string, char*);

string _GetDirectoryName(string);


struct CPU
{
	int ID;
	vector<double> Alphas;
	vector<double> f;
	vector<double> p;
	vector<double> c;
	vector<Product*> Products;
	int coveredProducts;
	double costToCoverageRatio;
};

struct Product
{
	int ID;
	vector<double> d;
	vector<double> d_t_T;
	vector<vector<double>> d_t1_t2;
	vector<double> h;
	vector<CPU*> CPUs;
	double remainingDemand;
	bool isCovered;
	CPU* currentCPU; //backloglu heuristic icin sort edebilmek icin yaptim
	vector<double> LB; //calculated lower bound on inventory variables wrt. to time
	vector<double> UB; //calculated upper bound on inventory variables wrt. to time
};

struct Instance
{
	int T;
	int U;
	int P;
	vector<Product> Products;
	vector<CPU> CPUs;
	vector<CPU*>CPUs_P;
	string Path;
	void ReadInstance(string);
};

struct Parameters
{
	double epsilon = 0;
	int mode = -1;
	int heur = 0;
	bool CPWW = false;
	bool CplexCuts = true;
	int TimeLimitMinutes = 0;
	bool EarlyCutAbortion = false;
	int LBonInvVars = 0;
	int UBonInvVars = 0;
	int threadCount = 4;
	string StartDate;
	int T_Sub = 0;
	bool WWV = false;
	bool VI2 = false;
	bool VI1 = false;

	bool singleCut = false;

	bool littleModel = false;
};