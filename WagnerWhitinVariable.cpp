#pragma once
#include "WagnerWhitinVariable.h"

WagnerWhitinVariable::WagnerWhitinVariable(int PeriodIn, Product* ProductIn, vector<CPU*> &CPUsIn)
{
	N = PeriodIn;
	Prod = ProductIn;
	Cpus = CPUsIn;

	for (int t = 0; t < N; ++t)
	{
		H.push_back(Prod->h[t]);
	}
}
WagnerWhitinVariable::~WagnerWhitinVariable()
{
	H.end();
	M.end();
	F.end();
	Q.end();
	JSTAR.end();
	I.end();
	TEMP_I.end();
	Cpus.end();
}

int WagnerWhitinVariable::Compute()
{
	F.resize(N);
	M.resize(N);
	I.resize(N);//cpu selection through periods
	TEMP_I.resize(N);

	JSTAR.resize(N);
	int KM1, JM1;
	double TEMP;
	for (int K = 0; K < N; ++K)
	{
		F[K] = INT_MAX;
		M[K] = Cpus[0]->f[K] + ((Cpus[0]->p[K]/ Cpus[0]->Alphas[Prod->ID-1])* Prod->d[K]);

		for (int i = 0; i < Cpus.size(); ++i)
		{
			TEMP = Cpus[i]->f[K] + ((Cpus[i]->p[K] / Cpus[i]->Alphas[Prod->ID - 1]) * Prod->d[K]);
			if (TEMP < M[K])
			{
				M[K] = TEMP;
				I[K] = i;
			}
		}
	}
	F[0] = M[0];
	JSTAR[0] = 0;
	TEMP = 0;
	double SH = 0;
	for (int K = 1; K < N; ++K)
	{
		KM1 = K - 1;
		for (int J = 0; J <= KM1; ++J)
		{
			SH = 0;
			for (int L = J + 1; L <= K; ++L)
				SH += Prod->d_t1_t2[L][K] * H[L-1];

			JM1 = J - 1;
			M[J] = Cpus[0]->f[J] + ((Cpus[0]->p[J] / Cpus[0]->Alphas[Prod->ID - 1]) * Prod->d_t1_t2[J][K]) + SH;
			TEMP_I[J] = 0;
			
			for (int i = 1; i < Cpus.size(); ++i)
			{
				TEMP = Cpus[i]->f[J] + ((Cpus[i]->p[J] / Cpus[i]->Alphas[Prod->ID - 1]) * Prod->d_t1_t2[J][K]) + SH;
					if (TEMP < M[J])
					{
						M[J] = TEMP;
						TEMP_I[J] = i;
					}
			}
			//M[J] = M[J] + (D[K] * (C[J] + S[J]));

			if (J == 0)
				TEMP = M[J];
			if (J > 0)
				TEMP = F[JM1] + M[J];
			if (TEMP >= F[K])
				continue;
			F[K] = TEMP;
			JSTAR[K] = J;
			I[J] = TEMP_I[J];
		}
		TEMP = F[KM1] + M[K];
		if (TEMP >= F[K])
			continue;
		F[K] = TEMP;
		JSTAR[K] = K;
	}
	//calculate Q
	Q.resize(N, 0);
	for (int K = N - 1; K >= 0; --K)
		for (int J = K; J >= JSTAR[K]; --J)
		{
			Q[JSTAR[K]] += Prod->d[J];
			if (J == JSTAR[K])
				K = JSTAR[K];
		}
	return F[N - 1];
}