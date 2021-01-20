#include "WW.h"

WW::WW()
{}

WW::WW(Instance instance) : Heuristic(instance){}

WW::~WW()
{

}

void WW::ImplementAlgorithm(Instance& instance, NumArray2 ySol)
{
	duration = clock();
	vector<CPU*> SelectedCPUs(instance.T);

	double * x_partial = new double[numberPeriods];

	for (int t = 0; t < instance.T; t++) {
		x_partial[t] = 0;
		for (int i = 0; i < instance.U; i++)
		{
			y[i][t] = false;
			x[i][t] = 0;
		}
	}

	//Fill remaining demand table
	double ** DTable = new double*[instance.P];
	for (int j = 0; j < instance.P; ++j)
		DTable[j] = new double[instance.T];

	for (int j = 0; j < instance.P; ++j)
		for (int t = 0; t < instance.T; ++t)
			DTable[j][t] = instance.Products[j].d[t];
	//Fill remaining demand table

	//initialize WWTAble
	double ** WWTable = new double*[instance.T];
	for (int t = 0; t < instance.T; ++t)
	{
		WWTable[t] = new double[instance.T];
		for (int t1 = 0; t1 < instance.T; ++t1)
			WWTable[t][t1] = INT_MAX;
	}
	//initialize WWTAble

	//initialize solution array
	int* solution = new int[instance.T];
	//initialize solution array

	//order products randomly
	vector<int> ProductOrder(instance.P);
	for (int j = 0; j < instance.P; ++j)
		ProductOrder[j] = j;
	random_shuffle(ProductOrder.begin(), ProductOrder.end());
	//order products randomly

	//for each randomly shuffled product
	for (int j : ProductOrder)
	{
		for (int t = 0; t < instance.T; t++)
			x_partial[t] = 0;
		//for (int t = 0; t < instance.T; ++t) //clear the array
		//	SelectedCPUs[t] = -1;

		//select 1 CPU per period for procuding selected product
		for (int t = 0; t < instance.T; ++t)
		{
			double temp = -1;
			for (CPU* i : instance.Products[j].CPUs)
				if (ySol[i->ID-1][t] > temp)
				{
					temp = ySol[i->ID - 1][t];
					SelectedCPUs[t] = i;
				}
		}
		//select 1 CPU per period for procuding selected product

		//fill up the WW table


		for (int ti = 0 ; ti < instance.T ; ++ti)
			for (int tj = ti; tj < instance.T; ++tj)
			{
				WWTable[ti][tj] = 0;
				//Previous minimum
				if (ti != 0){
					double temp = INT_MAX;
					for (int ttemp = 0; ttemp < ti; ++ttemp)
						if (WWTable[ttemp][ti - 1] < temp)
							temp = WWTable[ttemp][ti - 1];
					WWTable[ti][tj] += temp;
				}
				//Fixed Cost
				if (!y[SelectedCPUs[ti]->ID - 1 ][ti])
					WWTable[ti][tj] += SelectedCPUs[ti]->f[ti];
				for (int ttemp = ti; ttemp <= tj; ++ttemp)
				{
					if (DTable[j][ttemp] >= 0) //because only the last period could have negative demand, meaning overproduction
						//if this is the case covering last period shouldnt add to the cost
					{
						//Production cost
						WWTable[ti][tj] += ((DTable[j][ttemp] * SelectedCPUs[ttemp]->p[ttemp]) / SelectedCPUs[ttemp]->Alphas[j]);
						//Holding cost
						WWTable[ti][tj] += (ttemp - ti) * instance.Products[j].h[ttemp] * DTable[j][ttemp];
					}
				}
			}
		//fill up the WW table

		////write the WW table
		//for (int ti = 0; ti < instance.T; ++ti, cout<< endl)
		//	for (int tj = 0; tj < instance.T; ++tj)
		//		cout << WWTable[ti][tj] << '\t';
		////write the WW table

		//deduce the solution
			
		double temp = INT_MAX;
		for (int tj = instance.T - 1; tj >= 0; --tj)
			solution[tj] = 0;

		for (int tj = instance.T - 1; tj >=0 ; --tj)
			for (int ti = 0; ti <= tj; ++ti)
				if (WWTable[ti][tj] < temp)
				{
					temp = WWTable[ti][tj];
					solution[tj] = ti;
				}
		//for (int tj = 0; tj < instance.T; ++tj)
		//	cout << solution[tj] << '\t';

		int ProductionPeriod = 0;
		for (int tj = 0; tj < instance.T; ++tj)
		{
			if (DTable[j][tj] <= 0)
			{
				ProductionPeriod = tj + 1;
			}
			else
			{
				if (tj == ProductionPeriod)
				{
					ProductionPeriod = tj;
					y[SelectedCPUs[tj]->ID - 1][tj] = true;

				}
				else if (solution[tj] > solution[tj - 1])
				{
					ProductionPeriod = tj;
					y[SelectedCPUs[tj]->ID - 1][tj] = true;
				}
				x_partial[ProductionPeriod] += DTable[j][tj] / SelectedCPUs[ProductionPeriod]->Alphas[j];
				x[SelectedCPUs[ProductionPeriod]->ID - 1][ProductionPeriod] += DTable[j][tj] / SelectedCPUs[ProductionPeriod]->Alphas[j];
			}
		}

		//deduce the solution
		
		//update remaining demand for all products
		for (int tj = 0; tj < instance.T; ++tj)
			if(x_partial[tj] > 0)
				for (Product* p : SelectedCPUs[tj]->Products)
				{
					DTable[p->ID-1][tj] -= x_partial[tj] * SelectedCPUs[tj]->Alphas[p->ID - 1];
				
						for (int ttemp = tj; ttemp < instance.T; ++ttemp)
						{
							if (DTable[p->ID - 1][ttemp] >= 0 || ttemp + 1 == instance.T)
								break;
							DTable[p->ID - 1][ttemp + 1] += DTable[p->ID - 1][ttemp];
							DTable[p->ID - 1][ttemp] = 0;
						}
				}
		//update remaining demand for all products

		//cout << endl;
		//for (int tj = 0; tj < instance.T; ++tj)
		//	cout << y[0][tj] << '\t';
		//cout << endl;
		//for (int tj = 0; tj < instance.T; ++tj)
		//	cout << x_partial[tj] << '\t';
		//cout << " ////////////////////////////////////////////////////////" << endl;
		//cout << "Selected Product = " << j << endl;
		//cout << "Selected CPUs" << endl;
		//for (int tt = 0; tt < instance.T; ++tt)
		//	cout << SelectedCPUs[tt]->ID - 1 << "\t";
		//cout << endl;
		//cout << "Alphas" << endl;
		//for (int tt = 0; tt < instance.T; ++tt)
		//	cout << SelectedCPUs[tt]->Alphas[j] << "\t";
		//cout << endl;
		//cout << "X_Partial" << endl;
		//for (int tj = 0; tj < instance.T; ++tj)
		//	cout << x_partial[tj] << '\t';
		//cout << endl;
		//cout << " ////////////////////////////////////////////////////////" << endl;
		//for (int jx = 0; jx < instance.P; ++jx, cout << endl)
		//	for (int t = 0; t < instance.T; ++t)
		//		cout << DTable[jx][t] << '\t';
	}
	//for each randomly shuffled product
	//cout << " ////////////////////////////////////////////////////////" << endl;

	//for (int jx = 0; jx < instance.P; ++jx, cout << endl)
	//	for (int t = 0; t < instance.T; ++t)
	//		cout << DTable[jx][t] << '\t';

	//cout << " ////////////////////////////////////////////////////////" << endl;
	//for (int i = 0; i < instance.U; ++i)
	//	for (int t = 0; t < instance.T; ++t)
	//	{
	//		if (y[i][t])
	//			cout << "y" << i << "_" << t << '\t' << y[i][t] << '\t';
	//		if (x[i][t] != 0)
	//			cout << "x" << i << "_" << t << '\t' << x[i][t] << endl;
	//	}

	//cout << endl;
	//for (int tj = 0; tj < instance.T; ++tj)
	//	cout << x_partial[tj] << '\t';
	
	duration = clock() - duration;
	CalculateObjectiveFunctionValue(instance);

	delete[] x_partial;
	delete[] solution;
	
	for (int t = 0; t < instance.T; ++t)
		delete[] WWTable[t];

	for (int j = 0; j < instance.P; ++j)
		delete[] DTable[j];
	
	delete[] DTable;
	delete[] WWTable;
}