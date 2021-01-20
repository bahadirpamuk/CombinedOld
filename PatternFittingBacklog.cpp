#include "PATTERNFITTINGBACKLOG.h"

PatternFittingBacklog::PatternFittingBacklog()
{}

PatternFittingBacklog::PatternFittingBacklog(Instance instance) : Heuristic(instance) {}

PatternFittingBacklog::~PatternFittingBacklog()
{

}

void PatternFittingBacklog::ImplementAlgorithm(Instance& instance)
{
	duration = clock();
	for (int t = 0; t < instance.T; t++)
	{
		//double remainingCapacity = instance.capacity[t];
		for (int j = 0; j < instance.P; j++)
		{
			if (t == 0)
				instance.Products[j].remainingDemand = instance.Products[j].d[t];
			else
				instance.Products[j].remainingDemand += instance.Products[j].d[t];

			if (instance.Products[j].remainingDemand <= 0)
				instance.Products[j].isCovered = true;
			else
				instance.Products[j].isCovered = false;

			//cout << "remaining demand " << j << " is " << instance.Products[j]->remainingDemand << " in period " 
			//	<< t << " coverage " << instance.Products[j]->isCovered << endl;
		}

		for (int i = 0; i < instance.U; i++)
		{
			y[i][t] = false;
			x[i][t] = 0;
		}

		bool coverage = false;
		while (!coverage)
		{
			int CPUindex = -1;
			double minRatio = 1000000000000;

			for (int i = 0; i < instance.U; i++)
				//if (!y[i][t]) //a cpu only producing a product is partially produced infinite loop
				{ 
					// calculate instance.CPUs[i]->coveredProducts, we want to cover products with remainingDemand > 0
					instance.CPUs[i].coveredProducts = 0;
					for (Product* jp : instance.CPUs[i].Products)
					{
						if (jp->remainingDemand > 0)
							instance.CPUs[i].coveredProducts++;
						jp->currentCPU = &instance.CPUs[i];
					}
					// calculate x_it
					double xprev = 0;
					if (!y[i][t])
						x[i][t] = 0;
					else
						xprev = x[i][t];
					std::sort(instance.CPUs[i].Products.begin(), instance.CPUs[i].Products.end(), sortfunction);
					double tempx = 0;
					double tempCost = 1000000000;
					for (int j = 0; j < instance.CPUs[i].Products.size(); ++j) //urunler d/a lara gore buyukten kucuge sirali
					{
						Product* tempPro = instance.CPUs[i].Products[j];
						if (tempPro->isCovered)
							break;
						x[i][t] = tempPro->remainingDemand / instance.CPUs[i].Alphas[(tempPro->ID - 1)];

						// determine the overdproduction amounts
						for (Product* j : instance.CPUs[i].Products)
							s[j->ID - 1][t] = x[i][t] * instance.CPUs[i].Alphas[j->ID - 1] - (j->remainingDemand);

						//for (int ss = 0; ss < instance.CPUs[i]->CPUSize; ss++)
						//	cout << "Overproduction of product " << instance.CPUs[i]->alpha[ss][0] << " is " << s[instance.CPUs[i]->alpha[ss][0]][t] << endl;


						if (instance.CPUs[i].coveredProducts > 0 && x[i][t] > 0)
						{
							// calculate instance.CPUs[i]->costToCoverageRatio
							instance.CPUs[i].costToCoverageRatio = 0;

							// add fixed cost
							instance.CPUs[i].costToCoverageRatio += instance.CPUs[i].f[t];

							// add variable cost
							instance.CPUs[i].costToCoverageRatio += instance.CPUs[i].p[t] * x[i][t];

							// add overproduction cost
							for (Product* jp : instance.CPUs[i].Products)
								if (s[jp->ID - 1][t] >  0)
								instance.CPUs[i].costToCoverageRatio += jp->h[t] * s[jp->ID - 1][t];
							// add backlog cost
							if (j > 0) //j 0 dan buyukse bazi urunleri backlog yapiyoruz demektir
								for (int j_back = 0; j_back < j; ++j_back)
								{
									Product* tP_back = instance.CPUs[i].Products[j_back];
									double backlogamount = tP_back->remainingDemand - x[i][t] * instance.CPUs[i].Alphas[tP_back->ID - 1];
									for (int t_back = t; t_back < instance.T; ++t_back)
										instance.CPUs[i].costToCoverageRatio +=
										2 * (instance.T - t_back) * backlogamount;
								}
							if (instance.CPUs[i].costToCoverageRatio > tempCost)
							{
								x[i][t] = tempx;
							}
							else
							{
								tempx = x[i][t];
								tempCost = instance.CPUs[i].costToCoverageRatio;
							}
						}
					}
					x[i][t] += xprev; //eger onceden uretim yaplmssa bundan
					if (instance.CPUs[i].coveredProducts > 0 && x[i][t] > 0)
					{
						// divide coveredProducts
						instance.CPUs[i].costToCoverageRatio /= instance.CPUs[i].coveredProducts;

						if (instance.CPUs[i].costToCoverageRatio < minRatio /* && instance.CPUs[i]->capacityUsage[t] * x[i][t] <= remainingCapacity*/)
						{
							minRatio = instance.CPUs[i].costToCoverageRatio;
							CPUindex = i;
						}
					}
				}

			// pick the CPU with the smallest costToCoverageRatio -> y_it
			if (CPUindex >= 0)
			{
				y[CPUindex][t] = true;
				for (Product* j : instance.CPUs[CPUindex].Products)
				{
					j->remainingDemand -= instance.CPUs[CPUindex].Alphas[j->ID - 1] * x[CPUindex][t];
					if (j->remainingDemand < 0.001)
						j->isCovered = true;
				}

				// check the coverage status of the products 
				coverage = true;
				for (int j = 0; j < instance.P; j++)
					if (instance.Products[j].isCovered == false)
						coverage = false;

				if (!coverage && CPUindex < 0)
					cout << "capacity is not enough to cover the products in period " << t << endl;
			}

			// unselected CPUs should not be produced, i.e., x[i][t] = 0
			for (int i = 0; i < instance.U; i++)
				if (!y[i][t])
					x[i][t] = 0;

			// update inventory levels s_jt 
			for (int j = 0; j < instance.P; j++)
			{
				s[j][t] = 0;
				if (instance.Products[j].remainingDemand < 0)
					s[j][t] = -instance.Products[j].remainingDemand;
			}
		}
	}
	duration = clock() - duration;

	CalcObj(instance);
}

void PatternFittingBacklog::CalcObj(Instance& instance)
{
	algObjValue = 0;
	for (int i = 0; i < numberCPUs; i++)
		for (int t = 0; t < numberPeriods; t++)
		{
			if (y[i][t]) {
				algObjValue += y[i][t] * instance.CPUs[i].f[t];
				algObjValue += x[i][t] * instance.CPUs[i].p[t];
			}
		}
	for (int t = 0; t < numberPeriods; t++)
		for (Product j : instance.Products)
			s[j.ID - 1][t] = 0;

	for (int t = 0; t < numberPeriods; t++)
	{
		for (int i = 0; i < numberCPUs; i++)
			if (y[i][t])
				for (Product* j : instance.CPUs[i].Products)
					s[j->ID - 1][t] += x[i][t] * instance.CPUs[i].Alphas[j->ID - 1];

		for (Product j : instance.Products)
			s[j.ID - 1][t] -= j.d[t];
	}
	for (int j = 0; j < numberProducts; j++)
		for (int t = 0; t < numberPeriods; t++)
			if (s[j][t] > 0)
				algObjValue += s[j][t] * instance.Products[j].h[t];
			else
				algObjValue -= s[j][t] * 2;

};