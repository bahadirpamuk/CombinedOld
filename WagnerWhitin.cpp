#include "WagnerWhitin.h"

WagnerWhitin::WagnerWhitin(int PeriodIn, Product* ProductIn, CPU* CPUIn)
{
	N = PeriodIn;
	Prod = ProductIn;
	Cpu = CPUIn;

	for (int t = 0; t < N; ++t)
	{
		D.push_back(Prod->d[t]);
		C.push_back(Cpu->p[t]);
		A.push_back(Cpu->f[t]);
		H.push_back(Prod->h[t]);
	}
}
WagnerWhitin::~WagnerWhitin()
{
	D.end();
	C.end();
	A.end();
	H.end();
	M.end();
	F.end();
	Q.end();
	S.end();
	JSTAR.end();
}

int WagnerWhitin::Compute()
{
	F.resize(N);
	S.resize(N);
	M.resize(N);
	JSTAR.resize(N);
	int KM1, JM1, TEMP;
	for (int K = 0; K < N ;++K)
	{
		F[K] = INT_MAX;
		S[K] = 0;
		M[K] = A[K] + (C[K] * D[K]);
	}
	F[0] = M[0];
	JSTAR[0] = 0;
	for (int K = 1; K < N; ++K)
	{
		KM1 = K - 1;
		for (int J = 0; J <= KM1; ++J)
		{
			JM1 = J - 1;
			S[J] = S[J] + H[KM1];
			M[J] = M[J] + (D[K] * (C[J] + S[J]));
			if (J == 0)
				TEMP = M[J];
			if (J > 0)
				TEMP = F[JM1] + M[J];
			if (TEMP >= F[K])
				continue;
			F[K] = TEMP;
			JSTAR[K] = J;
		}
		TEMP = F[KM1] + M[K];
		if (TEMP >= F[K])
			continue;
		F[K] = TEMP;
		JSTAR[K] = K;
	}
	//calculate Q
	Q.resize(N, 0);
	for(int K = N - 1 ; K >=0 ; --K)
		for (int J = K; J >= JSTAR[K]; --J)
		{
			Q[JSTAR[K]] += D[J];
			if(J == JSTAR[K])
				K = JSTAR[K];
		}
	return F[N - 1];
}