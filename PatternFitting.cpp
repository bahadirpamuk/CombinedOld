 #include "PatternFitting.h"

PatternFitting::PatternFitting()
{}

PatternFitting::PatternFitting(Instance instance) : Heuristic(instance){}

PatternFitting::~PatternFitting()
{

}

void PatternFitting::ImplementAlgorithm(Instance& instance)
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
			if(!y[i][t])
			{ // calculate x_it				
				x[i][t] = 0;
				for (Product* j : instance.CPUs[i].Products)
				{
					double ratio = j->remainingDemand / instance.CPUs[i].Alphas[j->ID - 1];
					if (ratio > x[i][t])
						x[i][t] = ratio;
				}	
				//cout << "Production amount of CPU " << i << " is " << x[i][t] << endl;

				// determine the overdproduction amounts
				for(Product* j : instance.CPUs[i].Products)
					s[j->ID-1][t] = x[i][t] * instance.CPUs[i].Alphas[j->ID - 1] - (j->remainingDemand);
			
				//for (int ss = 0; ss < instance.CPUs[i]->CPUSize; ss++)
				//	cout << "Overproduction of product " << instance.CPUs[i]->alpha[ss][0] << " is " << s[instance.CPUs[i]->alpha[ss][0]][t] << endl;

				// calculate instance.CPUs[i]->coveredProducts, we want to cover products with remainingDemand > 0
				instance.CPUs[i].coveredProducts = 0;
				for (Product* j : instance.CPUs[i].Products)
					if (j->remainingDemand > 0)
						instance.CPUs[i].coveredProducts++;
				//cout << "CPU " << i << " can cover " << instance.CPUs[i]->coveredProducts << " new products" << endl;

				if (instance.CPUs[i].coveredProducts > 0 && x[i][t] > 0)
				{
					// calculate instance.CPUs[i]->costToCoverageRatio
					instance.CPUs[i].costToCoverageRatio = 0;

					// add fixed cost
					instance.CPUs[i].costToCoverageRatio += instance.CPUs[i].f[t];

					// add variable cost
					instance.CPUs[i].costToCoverageRatio += instance.CPUs[i].p[t] * x[i][t];

					// add overproduction cost
					for (Product* j : instance.CPUs[i].Products)
						instance.CPUs[i].costToCoverageRatio +=
						j->h[t] * s[j->ID-1][t];

					//// multiply capacity usage
					//instance.CPUs[i]->costToCoverageRatio *= (instance.CPUs[i]->capacityUsage[t] * x[i][t] / 10000); // division is for scaling

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
					j->isCovered = true;
					j->remainingDemand -= instance.CPUs[CPUindex].Alphas[j->ID-1] * x[CPUindex][t];
					//cout << "remaining demand " << instance.CPUs[CPUindex]->alpha[ss][0] << " is " 
					//	<< instance.Products[instance.CPUs[CPUindex]->alpha[ss][0]]->remainingDemand << " in period "
					//	<< t << " coverage " << instance.Products[instance.CPUs[CPUindex]->alpha[ss][0]]->isCovered << endl;
				}		
				/*remainingCapacity -= instance.CPUs[CPUindex]->capacityUsage[t] * x[CPUindex][t];*/
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
	duration = clock() - duration;

	CalculateObjectiveFunctionValue(instance);
}