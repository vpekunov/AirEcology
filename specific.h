#ifndef __SPECIFIC__
#define __SPECIFIC__

#include "CFGFILE.H"
#include "area3d.h"
#include "KINETIC.H"
#include "solvers.h"

double UwMainPhase = 0;
double _ZeroNu (int i, int x, int y, int z, unsigned char Map);
double _ZeroT (int i, int x, int y, int z, unsigned char Map);
double _ZeroUx (int i, int x, int y, int z, unsigned char Map);
double _ZeroUy (int i, int x, int y, int z, unsigned char Map);
double _ZeroUz (int i, int x, int y, int z, unsigned char Map);
double _ZeroP (int i, int x, int y, int z, unsigned char Map);
double _ZeroKt (int i, int x, int y, int z, unsigned char Map);
double _ZeroEt (int i, int x, int y, int z, unsigned char Map);
double _ZeroFUV (int i, int x, int y, int z, unsigned char Map);
double _ZeroUV (int i, int x, int y, int z, unsigned char Map);
double _ZeroF0 (int i, int x, int y, int z, unsigned char Map);
double _ZeroI (int i, int x, int y, int z, unsigned char Map);
double _ZeroFIR (int i, int x, int y, int z, unsigned char Map);
double _ZeroIR (int i, int x, int y, int z, unsigned char Map);
double _ZeroIRt (int i, int x, int y, int z, unsigned char Map);
double KappaKt = 1.0/0.7194;
double KappaEt = 1.0/0.7194;
double KappaNu = 0.0;
double NuMolNu = 0.0;
double KappaUx = 1;
double KappaUy = 1;
double KappaUV = 0.0;
double NuMolUV = 0.0;
double KappaFUV = 0.0;
double NuMolFUV = 0.0;
double KappaF0 = 0.0;
double NuMolF0 = 0.0;
double KappaI = 0.0;
double NuMolI = 0.0;
double KappaFIR = 0.0;
double NuMolFIR = 0.0;
double KappaIR = 0.0;
double NuMolIR = 0.0;
double KappaIRt = 0.0;
double NuMolIRt = 0.0;
double KappaT = 1.393;
double NuMolT = 1.978E-5;
double KappaUz = 1;
double KappaP = 0.0;
double NuMolP = 0.0;
double _ZeroRok (int i, int x, int y, int z, unsigned char Map);
double _ZeroUkx (int i, int x, int y, int z, unsigned char Map);
double KappaUkx[7] = {0.0};
double NuMolUkx[7] = {0.0};
double _ZeroUky (int i, int x, int y, int z, unsigned char Map);
double KappaUky[7] = {0.0};
double NuMolUky[7] = {0.0};
double _ZeroUkz (int i, int x, int y, int z, unsigned char Map);
double KappaUkz[7] = {0.0};
double NuMolUkz[7] = {0.0};
double _ZeroNk (int i, int x, int y, int z, unsigned char Map);
double _ZeroSO2k (int i, int x, int y, int z, unsigned char Map);
double KappaSO2k[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double NuMolSO2k[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double KappaRok[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double NuMolRok[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double KappaNk[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double NuMolNk[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double RO;
double BT;
double Kx;
double D;
double Ky;
double T0;
double uSolarIntensity;
double SurfaceAlbedo;
double uDiffuseFlux;
double SolarIntensity;
double DiffuseFlux;
double rSolarIntensity;
double rDiffuseFlux;
double IntPlanck[50];
double SurfaceEmissivity;
double tDiffuseFlux;
double ZeroT;
double SurfaceCC;
double SurfaceRO;
double SurfaceLambda;
double UwLight[7];
double InitRok[50];
double InjectStartTime;
double InjectV;
double InjectCellLength;
double InjectN[50];
double Dk[50];
double InitNk[50];
double InitSO2k[50];
double AlphaM;
double AlphaT;
double QeC0[50];
double QeC1[50];
double QeC2[50];
double QeC3[50];
double QeC4[50];
double QeC5[50];
double QaC0[50];
double QaC1[50];
double QaC2[50];
double QaC3[50];
double QaC4[50];
double QaC5[50];
double gC0[50];
double gC1[50];
double gC2[50];
double gC3[50];
double gC4[50];
double gC5[50];
double uQeC0[50];
double uQeC1[50];
double uQeC2[50];
double uQeC3[50];
double uQeC4[50];
double uQeC5[50];
double uQaC0[50];
double uQaC1[50];
double uQaC2[50];
double uQaC3[50];
double uQaC4[50];
double uQaC5[50];
double ugC0[50];
double ugC1[50];
double ugC2[50];
double ugC3[50];
double ugC4[50];
double ugC5[50];
double rQeC0[50];
double rQeC1[50];
double rQeC2[50];
double rQeC3[50];
double rQeC4[50];
double rQeC5[50];
double rQaC0[50];
double rQaC1[50];
double rQaC2[50];
double rQaC3[50];
double rQaC4[50];
double rQaC5[50];
double rgC0[50];
double rgC1[50];
double rgC2[50];
double rgC3[50];
double rgC4[50];
double rgC5[50];
double tQeC0[50];
double tQeC1[50];
double tQeC2[50];
double tQeC3[50];
double tQeC4[50];
double tQeC5[50];
double tQaC0[50];
double tQaC1[50];
double tQaC2[50];
double tQaC3[50];
double tQaC4[50];
double tQaC5[50];
double tgC0[50];
double tgC1[50];
double tgC2[50];
double tgC3[50];
double tgC4[50];
double tgC5[50];
double AirBeta;
double AirAlbedo;
double AirAsymmetry;
double uAirBeta;
double uAirAlbedo;
double uAirAsymmetry;
double rAirBeta;
double rAirAlbedo;
double rAirAsymmetry;
double tAirBeta;
double tAirAlbedo;
double tAirAsymmetry;

#define NumEqs 57

enum {_NumNk0,_NumRok0,_NumSO2k0,_NumUkx0,_NumUky0,_NumUkz0,_NumNk1,_NumRok1,_NumSO2k1,_NumUkx1,_NumUky1,_NumUkz1,_NumNk2,_NumRok2,_NumSO2k2,_NumUkx2,_NumUky2,_NumUkz2,_NumNk3,_NumRok3,_NumSO2k3,_NumUkx3,_NumUky3,_NumUkz3,_NumNk4,_NumRok4,_NumSO2k4,_NumUkx4,_NumUky4,_NumUkz4,_NumNk5,_NumRok5,_NumSO2k5,_NumUkx5,_NumUky5,_NumUkz5,_NumNk6,_NumRok6,_NumSO2k6,_NumUkx6,_NumUky6,_NumUkz6,_NumEt,_NumF0,_NumFIR,_NumFUV,_NumI,_NumIR,_NumIRt,_NumKt,_NumNu,_NumP,_NumT,_NumUV,_NumUx,_NumUy,_NumUz};

#define _NumDn NumEqs /* Вещества */

int PredictControlVar = _NumNu;

#define NumLightPhases 0
#define NumHeavyPhases 8
#define NumPhases (NumLightPhases+NumHeavyPhases)

enum {phMainPhase,phDrops0,phDrops1,phDrops2,phDrops3,phDrops4,phDrops5,phDrops6};

int CarrierPhase = phMainPhase;

#define NumSaves (NumEqs-NumHeavyPhases*(nDims-1))

typedef struct {
 float  * Name[NumEqs];
 float ** Dn; /* Вещества */
 float  * eName[NumEqs];
 float ** eDn; /* Вещества по явному шагу Эйлера */
 int     BoundSize;
 float * Bounds[6];
 float * Parts[3];
 float * PartsV[3];
} StoreStruct;

typedef struct {
  int    IsLight;
  double * _Uw;
  int    _Ux, _Uy, _Uz, _Nu, _T, _Ro;
  const char * Source;
} PhaseDsc;

PhaseDsc PhaseVars[NumPhases] = {
 {0,&UwMainPhase,_NumUx,_NumUy,_NumUz,_NumNu,_NumT,-1,""},
 {0,&UwLight[0],_NumUkx0,_NumUky0,_NumUkz0,-1,-1,_NumRok0,"H2O"},
 {0,&UwLight[1],_NumUkx1,_NumUky1,_NumUkz1,-1,-1,_NumRok1,"H2O"},
 {0,&UwLight[2],_NumUkx2,_NumUky2,_NumUkz2,-1,-1,_NumRok2,"H2O"},
 {0,&UwLight[3],_NumUkx3,_NumUky3,_NumUkz3,-1,-1,_NumRok3,"H2O"},
 {0,&UwLight[4],_NumUkx4,_NumUky4,_NumUkz4,-1,-1,_NumRok4,"H2O"},
 {0,&UwLight[5],_NumUkx5,_NumUky5,_NumUkz5,-1,-1,_NumRok5,"H2O"},
 {0,&UwLight[6],_NumUkx6,_NumUky6,_NumUkz6,-1,-1,_NumRok6,"H2O"}
};

/*!!!!!*/
#define sign(a) ((a) == 0 ? 0 : (a) < 0 ? -1 : +1)

/* Аналитическое решение */
double aGetXY(double D, double RR)
{
 double D2 = D*D;
 double D3 = D*D2;
 double D6 = D3*D3;
 double z = pow(-80.0*D3+108.0*RR+12.0*sqrt(48.0*D6-120.0*D3*RR+81.0*RR*RR),1.0/3.0);
 return z/6.0-4.0/3.0*D2/z-D/3.0;
}

/*!!!!*/
/* Численное решение */
double GetXY(double D, double RR)
{
	double dx;
	double aResult = aGetXY(D, RR);
	double Result = isfinite(aResult) ? aResult : D;
	int    Counter = 0;

	do {
		double z1 = D*D + Result*Result;
		double z2 = D + Result;
		dx = -(z1*z2 - RR) / (2.0*Result*z2 + z1);
		Result += 0.5*dx;
	} while (isfinite(Result) && fabs(dx) >= D*0.0000001 && ++Counter<5);
	return Counter<5 && isfinite(Result) ? Result : isfinite(aResult) ? aResult : D;
}

double CD(double Re)
{
 if (Re<0.001) return CD(0.001);
 else if (Re<0.1) return 24.0/Re;
 else if (Re<1.0) return 24.0/Re*(1.0+3.0/16.0*Re+9.0/160.0*Re*Re*log(2*Re));
 else if (Re<900.0) return 24.0/Re*(1+0.15*pow(Re,0.678));
 else return 0.44;
}

double FFk(double g, double RO, double ROK, double Sigma, double Nu, double D, double dU)
{
 double adU = fabs(dU);
 double Re = D*adU/Nu;
 double We = 0.5*RO*D*dU*dU/Sigma;

 return 3.0/8.0*CD(Re)*(1.0+2.632*We/12.0)*RO/ROK*2.0/D*adU*dU-g;
}

/*!!!!!*/
double GetdU(double g, double RO, double ROK, double Sigma, double Nu, double D, double Init_dU)
{
	double dx;
	double Result = Init_dU;
	int    Counter = 0;

	if (!isfinite(Init_dU) || fabs(Init_dU) > 25) Init_dU = 0.0;

	do {
		double FF = FFk(g, RO, ROK, Sigma, Nu, D, Result);
		double FF1 = FFk(g, RO, ROK, Sigma, Nu, D, Result + 0.00001);
		dx = -FF / (FF1 - FF)*0.00001;
		Result += min(1.5, 0.5*dx);
	} while (fabs(dx) >= 0.0001 && (++Counter<80) && isfinite(Result) && fabs(Result) < 25 + fabs(Init_dU));

	Result = Counter >= 80 || !isfinite(Result) || fabs(Result) > 25 + fabs(Init_dU) ? Init_dU : Result;

	return Result;
}

/* !!! Оставить в specific НАВСЕГДА, а не в функциях эксперимента */
void _GetLU(int NN, int * iRow, double * A, double * LU)
{
	int i, j, k;

	memmove(LU, A, NN*NN*sizeof(double));
	for (i = 0; i<NN; i++)
		iRow[i] = i;
	for (i = 0; i<NN - 1; i++)
	{
		double Big = 0.0;
		int    iBig = -1;

		double Kf;

		for (j = i; j<NN; j++)
		{
			double size = fabs(LU[iRow[j] * NN + i]);

			if (size>Big)
			{
				Big = size;
				iBig = j;
			}
		}
		if (iBig != i)
		{
			int V = iRow[i];
			iRow[i] = iRow[iBig];
			iRow[iBig] = V;
		}
		Kf = 1.0 / LU[iRow[i] * NN + i];

		LU[iRow[i] * NN + i] = Kf;
		for (j = i + 1; j<NN; j++)
		{
			double Fact = Kf*LU[iRow[j] * NN + i];

			LU[iRow[j] * NN + i] = Fact;
			for (k = i + 1; k<NN; k++)
				LU[iRow[j] * NN + k] -= Fact*LU[iRow[i] * NN + k];
		}
	}
	LU[(iRow[NN - 1] + 1)*NN - 1] = 1.0 / LU[(iRow[NN - 1] + 1)*NN - 1];
}

/*!!!  Оставить в specific НАВСЕГДА, а не в функциях эксперимента */
void _SolveLU(int NN, int * iRow, double * LU, double * Y, double * X)
{
	int i, j, k;

	X[0] = Y[iRow[0]];
	for (i = 1; i<NN; i++)
	{
		double V = Y[iRow[i]];

		for (j = 0; j<i; j++)
			V -= LU[iRow[i] * NN + j] * X[j];
		X[i] = V;
	}

	X[NN - 1] *= LU[(iRow[NN - 1] + 1)*NN - 1];

	for (i = 1, j = NN - 2; i<NN; i++, j--)
	{
		double V = X[j];

		for (k = j + 1; k<NN; k++)
			V -= LU[iRow[j] * NN + k] * X[k];
		X[j] = V*LU[iRow[j] * NN + j];
	}
}

/* !!! */
double getLeftSpline(int i, double * Dk, double * XYC, double * XYC2, double * XYC3, double * XYC4, double * XYC5,
	double * XYC6, double * XYC7, double ROK, double * AnC, double * BnC,
	double * Ann, double * Bnn, double * Cnn, double * Dnn) {

	double Pi = 3.1415926535897932;

	double MAT[4][4];
	double LU[4][4]; 
	double FMAT[4]; 
	double XMAT[4];
	int iRow[4];

	MAT[0][0] = XYC3[i - 1];
	MAT[0][1] = XYC2[i - 1];
	MAT[0][2] = XYC[i - 1];
	MAT[0][3] = 1.0;
	FMAT[0] = AnC[i - 1] * XYC[i - 1] + BnC[i - 1];
	MAT[1][0] = XYC3[i];
	MAT[1][1] = XYC2[i];
	MAT[1][2] = XYC[i];
	MAT[1][3] = 1.0;
	FMAT[1] = AnC[i] * XYC[i] + BnC[i];
	MAT[2][0] = 1.0 / 4.0*(XYC4[i] - XYC4[i - 1]);
	MAT[2][1] = 1.0 / 3.0*(XYC3[i] - XYC3[i - 1]);
	MAT[2][2] = 1.0 / 2.0*(XYC2[i] - XYC2[i - 1]);
	MAT[2][3] = XYC[i] - XYC[i - 1];
	FMAT[2] = AnC[i - 1] / 2.0*(Dk[i] * Dk[i] - XYC2[i - 1]) + BnC[i - 1] * (Dk[i] - XYC[i - 1]) +
		AnC[i] / 2.0*(XYC2[i] - Dk[i] * Dk[i]) + BnC[i] * (XYC[i] - Dk[i]);
	MAT[3][0] = Pi / 6.0*ROK*1.0 / 7.0*(XYC7[i] - XYC7[i - 1]);
	MAT[3][1] = Pi / 6.0*ROK*1.0 / 6.0*(XYC6[i] - XYC6[i - 1]);
	MAT[3][2] = Pi / 6.0*ROK*1.0 / 5.0*(XYC5[i] - XYC5[i - 1]);
	MAT[3][3] = Pi / 6.0*ROK*1.0 / 4.0*(XYC4[i] - XYC4[i - 1]);
	FMAT[3] = Pi / 6.0*ROK*(AnC[i - 1] / 5.0*(Dk[i] * Dk[i] * Dk[i] * Dk[i] * Dk[i] - XYC5[i - 1]) + BnC[i - 1] / 4.0*(Dk[i] * Dk[i] * Dk[i] * Dk[i] - XYC4[i - 1]) +
		AnC[i] / 5.0*(XYC5[i] - Dk[i] * Dk[i] * Dk[i] * Dk[i] * Dk[i]) + BnC[i] / 4.0*(XYC4[i] - Dk[i] * Dk[i] * Dk[i] * Dk[i]));

	_GetLU(4, iRow, (double *)MAT, (double *)LU);
	_SolveLU(4, iRow, (double *)LU, FMAT, XMAT);

	Ann[i] = XMAT[0];
	Bnn[i] = XMAT[1];
	Cnn[i] = XMAT[2];
	Dnn[i] = XMAT[3];

	return Ann[i];
}

/* !!! */
double getRightSpline(int i, double * Dk, double * XYC, double * XYC2, double * XYC3, double * XYC4, double * XYC5,
	double * XYC6, double * XYC7, double ROK, double * AnC, double * BnC,
	double * An, double * Bn, double * Cn, double * Dn) {

	double Pi = 3.1415926535897932;

	double MAT[4][4];
	double LU[4][4];
	double FMAT[4];
	double XMAT[4];
	int iRow[4];

	MAT[0][0] = XYC3[i];
	MAT[0][1] = XYC2[i];
	MAT[0][2] = XYC[i];
	MAT[0][3] = 1.0;
	FMAT[0] = AnC[i] * XYC[i] + BnC[i];
	MAT[1][0] = XYC3[i + 1];
	MAT[1][1] = XYC2[i + 1];
	MAT[1][2] = XYC[i + 1];
	MAT[1][3] = 1.0;
	FMAT[1] = AnC[i + 1] * XYC[i + 1] + BnC[i + 1];
	MAT[2][0] = 1.0 / 4.0*(XYC4[i + 1] - XYC4[i]);
	MAT[2][1] = 1.0 / 3.0*(XYC3[i + 1] - XYC3[i]);
	MAT[2][2] = 1.0 / 2.0*(XYC2[i + 1] - XYC2[i]);
	MAT[2][3] = XYC[i + 1] - XYC[i];
	FMAT[2] = AnC[i] / 2.0*(Dk[i + 1] * Dk[i + 1] - XYC2[i]) + BnC[i] * (Dk[i + 1] - XYC[i]) +
		AnC[i + 1] / 2.0*(XYC2[i + 1] - Dk[i + 1] * Dk[i + 1]) + BnC[i + 1] * (XYC[i + 1] - Dk[i + 1]);
	MAT[3][0] = Pi / 6.0*ROK*1.0 / 7.0*(XYC7[i + 1] - XYC7[i]);
	MAT[3][1] = Pi / 6.0*ROK*1.0 / 6.0*(XYC6[i + 1] - XYC6[i]);
	MAT[3][2] = Pi / 6.0*ROK*1.0 / 5.0*(XYC5[i + 1] - XYC5[i]);
	MAT[3][3] = Pi / 6.0*ROK*1.0 / 4.0*(XYC4[i + 1] - XYC4[i]);
	FMAT[3] = Pi / 6.0*ROK*(AnC[i] / 5.0*(Dk[i + 1] * Dk[i + 1] * Dk[i + 1] * Dk[i + 1] * Dk[i + 1] - XYC5[i]) + BnC[i] / 4.0*(Dk[i + 1] * Dk[i + 1] * Dk[i + 1] * Dk[i + 1] - XYC4[i]) +
		AnC[i + 1] / 5.0*(XYC5[i + 1] - Dk[i + 1] * Dk[i + 1] * Dk[i + 1] * Dk[i + 1] * Dk[i + 1]) + BnC[i + 1] / 4.0*(XYC4[i + 1] - Dk[i + 1] * Dk[i + 1] * Dk[i + 1] * Dk[i + 1]));

	_GetLU(4, iRow, (double *)MAT, (double *)LU);
	_SolveLU(4, iRow, (double *)LU, FMAT, XMAT);

	An[i] = XMAT[0];
	Bn[i] = XMAT[1];
	Cn[i] = XMAT[2];
	Dn[i] = XMAT[3];

	return An[i];
}

/* !!!! */ /* NEW */
double calcCoalescense(int Nc, double * Dk, double * X, double * Y, double * XYC,
	double * Ann, double * Bnn, double * Cnn, double * Dnn,
	double * An, double * Bn, double * Cn, double * Dn,
	int Ptr,
	float ** UX, float ** UY, float ** UZ,
	double * UwLight,
	double ROK,
	double * ddN, double * dRo, double * ddRo /* [Nc][Nc+1] */) {
	double Pi = 3.1415926535897932;

	int i, j;

	if (ddRo)
		for (i = 0, Ptr = 0; i < Nc; i++)
			for (j = 0; j <= Nc; j++, Ptr++)
				ddRo[Ptr] = 0.0;

	for (i = 0; i < Nc; i++)
		for (j = i; j < Nc; j++) {
			double h1 = Y[i] - X[i];
			double h2 = Y[j] - X[j];
			double D11 = X[i];
			double D12 = D11 + h1 / 2;
			double D13 = D11 + h1;
			double D21 = X[j];
			double D22 = D21 + h2 / 2;
			double D23 = D21 + h2;

			double N11 = D11 <= XYC[i] ? ((Ann[i] * D11 + Bnn[i])*D11 + Cnn[i])*D11 + Dnn[i] :
				((An[i] * D11 + Bn[i])*D11 + Cn[i])*D11 + Dn[i];
			double N12 = D12 <= XYC[i] ? ((Ann[i] * D12 + Bnn[i])*D12 + Cnn[i])*D12 + Dnn[i] :
				((An[i] * D12 + Bn[i])*D12 + Cn[i])*D12 + Dn[i];
			double N13 = D13 <= XYC[i] ? ((Ann[i] * D13 + Bnn[i])*D13 + Cnn[i])*D13 + Dnn[i] :
				((An[i] * D13 + Bn[i])*D13 + Cn[i])*D13 + Dn[i];
			double N21 = D21 <= XYC[j] ? ((Ann[j] * D21 + Bnn[j])*D21 + Cnn[j])*D21 + Dnn[j] :
				((An[j] * D21 + Bn[j])*D21 + Cn[j])*D21 + Dn[j];
			double N22 = D22 <= XYC[j] ? ((Ann[j] * D22 + Bnn[j])*D22 + Cnn[j])*D22 + Dnn[j] :
				((An[j] * D22 + Bn[j])*D22 + Cn[j])*D22 + Dn[j];
			double N23 = D23 <= XYC[j] ? ((Ann[j] * D23 + Bnn[j])*D23 + Cnn[j])*D23 + Dnn[j] :
				((An[j] * D23 + Bn[j])*D23 + Cn[j])*D23 + Dn[j];

			if (N11 > 1E-5 && N12 > 1E-5 && N13 > 1E-5 && N21 > 1E-5 && N22 > 1E-5 && N23 > 1E-5) {
				double dQ = Pi / 4.0*sqrt((UX[i][Ptr] - UX[j][Ptr])*(UX[i][Ptr] - UX[j][Ptr]) +
										(UY[i][Ptr] - UY[j][Ptr])*(UY[i][Ptr] - UY[j][Ptr]) +
										(UZ[i][Ptr] - UZ[j][Ptr])*(UZ[i][Ptr] - UZ[j][Ptr]))*
					(
					h1 / 6.0*(N11*(h2 / 6.0*(N21*(D11 + D21)*(D11 + D21) + 4 * N22*(D11 + D22)*(D11 + D22) + N23*(D11 + D23)*(D11 + D23))) +
					4 * N12*(h2 / 6.0*(N21*(D12 + D21)*(D12 + D21) + 4 * N22*(D12 + D22)*(D12 + D22) + N23*(D12 + D23)*(D12 + D23))) +
					N13*(h2 / 6.0*(N21*(D13 + D21)*(D13 + D21) + 4 * N22*(D13 + D22)*(D13 + D22) + N23*(D13 + D23)*(D13 + D23))))
					);
				if (dQ > 0.0) {
					double D3;
					int G;
					int found;
					D13 = (N11*D11*D11*D11 + 4 * N12*D12*D12*D12 + N13*D13*D13*D13) / (N11 + 4 * N12 + N13);
					D23 = (N21*D21*D21*D21 + 4 * N22*D22*D22*D22 + N23*D23*D23*D23) / (N21 + 4 * N22 + N23);
					D3 = pow(D13 + D23, 1.0 / 3.0);

					G = 0;
					found = 0;

					while (G < Nc && !found)
						if (Dk[G + 1] > D3)
							found = 1;
						else
							G++;

					ddN[i] -= dQ;
					ddN[j] -= dQ;
					dRo[i] -= dQ*ROK*Pi / 6.0*D13;
					dRo[j] -= dQ*ROK*Pi / 6.0*D23;
					if (G < Nc) {
						ddN[G] += dQ;
						dRo[G] += dQ*ROK*Pi / 6.0*D13;
						dRo[G] += dQ*ROK*Pi / 6.0*D23;
					}
                                        if (ddRo) {
						ddRo[i*(Nc+1)+Nc] -= dQ*ROK*Pi / 6.0*D13;
						ddRo[j*(Nc+1)+Nc] -= dQ*ROK*Pi / 6.0*D23;
						if (G < Nc) {
							ddRo[i*(Nc+1)+G] += dQ*ROK*Pi / 6.0*D13;
							ddRo[j*(Nc+1)+G] += dQ*ROK*Pi / 6.0*D23;
						}
                                        }
				}
			}
		}

	return 0.0;
}

enum { rsAny, rsPositive, rsNegative };

typedef struct {
  const char * FileBase;
  signed char  Projection;
  signed char  SubClass;
  _Solver      Solver;
  char         NeedK;
  char         NeedS;
  _InitVal     _Zero;
  double *     _NuMol;
  double *     _Kappa;
  char         Phase;
  const char * Description;
  char         Restrict;
  unsigned int Flags; // BBBBBBBB xxxxxxxx xxxxxxxx xxxxxxxD
  // BBBBBBBB -- boundary level.
  // D -- fld2RF
} VarDsc;

VarDsc VDefs[NumEqs] = {
 {"nk",-1,0,NULL,1,1,_ZeroNk,&NuMolNk[0],&KappaNk[0],phDrops0,"Концентрация капель [0]", rsPositive,0},
 {"rok",-1,0,NULL,1,1,_ZeroRok,&NuMolRok[0],&KappaRok[0],phDrops0,"Плотность капель [0]", rsPositive,0},
 {"so2k",-1,0,NULL,1,1,_ZeroSO2k,&NuMolSO2k[0],&KappaSO2k[0],phDrops0,"Концентрация SO2 в каплях [0]", rsPositive,0},
 {"uk",0,0,DirectSolver,1,0,_ZeroUkx,&NuMolUkx[0],&KappaUkx[0],phDrops0,"Скорость капель [0]", rsAny,0},
 {"uk",1,0,DirectSolver,1,0,_ZeroUky,&NuMolUky[0],&KappaUky[0],phDrops0,"Скорость капель [0]", rsAny,0},
 {"uk",2,0,DirectSolver,1,0,_ZeroUkz,&NuMolUkz[0],&KappaUkz[0],phDrops0,"Скорость капель [0]", rsAny,0},
 {"nk",-1,1,NULL,1,1,_ZeroNk,&NuMolNk[1],&KappaNk[1],phDrops1,"Концентрация капель [1]", rsPositive,0},
 {"rok",-1,1,NULL,1,1,_ZeroRok,&NuMolRok[1],&KappaRok[1],phDrops1,"Плотность капель [1]", rsPositive,0},
 {"so2k",-1,1,NULL,1,1,_ZeroSO2k,&NuMolSO2k[1],&KappaSO2k[1],phDrops1,"Концентрация SO2 в каплях [1]", rsPositive,0},
 {"uk",0,1,DirectSolver,1,0,_ZeroUkx,&NuMolUkx[1],&KappaUkx[1],phDrops1,"Скорость капель [1]", rsAny,0},
 {"uk",1,1,DirectSolver,1,0,_ZeroUky,&NuMolUky[1],&KappaUky[1],phDrops1,"Скорость капель [1]", rsAny,0},
 {"uk",2,1,DirectSolver,1,0,_ZeroUkz,&NuMolUkz[1],&KappaUkz[1],phDrops1,"Скорость капель [1]", rsAny,0},
 {"nk",-1,2,NULL,1,1,_ZeroNk,&NuMolNk[2],&KappaNk[2],phDrops2,"Концентрация капель [2]", rsPositive,0},
 {"rok",-1,2,NULL,1,1,_ZeroRok,&NuMolRok[2],&KappaRok[2],phDrops2,"Плотность капель [2]", rsPositive,0},
 {"so2k",-1,2,NULL,1,1,_ZeroSO2k,&NuMolSO2k[2],&KappaSO2k[2],phDrops2,"Концентрация SO2 в каплях [2]", rsPositive,0},
 {"uk",0,2,DirectSolver,1,0,_ZeroUkx,&NuMolUkx[2],&KappaUkx[2],phDrops2,"Скорость капель [2]", rsAny,0},
 {"uk",1,2,DirectSolver,1,0,_ZeroUky,&NuMolUky[2],&KappaUky[2],phDrops2,"Скорость капель [2]", rsAny,0},
 {"uk",2,2,DirectSolver,1,0,_ZeroUkz,&NuMolUkz[2],&KappaUkz[2],phDrops2,"Скорость капель [2]", rsAny,0},
 {"nk",-1,3,NULL,1,1,_ZeroNk,&NuMolNk[3],&KappaNk[3],phDrops3,"Концентрация капель [3]", rsPositive,0},
 {"rok",-1,3,NULL,1,1,_ZeroRok,&NuMolRok[3],&KappaRok[3],phDrops3,"Плотность капель [3]", rsPositive,0},
 {"so2k",-1,3,NULL,1,1,_ZeroSO2k,&NuMolSO2k[3],&KappaSO2k[3],phDrops3,"Концентрация SO2 в каплях [3]", rsPositive,0},
 {"uk",0,3,DirectSolver,1,0,_ZeroUkx,&NuMolUkx[3],&KappaUkx[3],phDrops3,"Скорость капель [3]", rsAny,0},
 {"uk",1,3,DirectSolver,1,0,_ZeroUky,&NuMolUky[3],&KappaUky[3],phDrops3,"Скорость капель [3]", rsAny,0},
 {"uk",2,3,DirectSolver,1,0,_ZeroUkz,&NuMolUkz[3],&KappaUkz[3],phDrops3,"Скорость капель [3]", rsAny,0},
 {"nk",-1,4,NULL,1,1,_ZeroNk,&NuMolNk[4],&KappaNk[4],phDrops4,"Концентрация капель [4]", rsPositive,0},
 {"rok",-1,4,NULL,1,1,_ZeroRok,&NuMolRok[4],&KappaRok[4],phDrops4,"Плотность капель [4]", rsPositive,0},
 {"so2k",-1,4,NULL,1,1,_ZeroSO2k,&NuMolSO2k[4],&KappaSO2k[4],phDrops4,"Концентрация SO2 в каплях [4]", rsPositive,0},
 {"uk",0,4,DirectSolver,1,0,_ZeroUkx,&NuMolUkx[4],&KappaUkx[4],phDrops4,"Скорость капель [4]", rsAny,0},
 {"uk",1,4,DirectSolver,1,0,_ZeroUky,&NuMolUky[4],&KappaUky[4],phDrops4,"Скорость капель [4]", rsAny,0},
 {"uk",2,4,DirectSolver,1,0,_ZeroUkz,&NuMolUkz[4],&KappaUkz[4],phDrops4,"Скорость капель [4]", rsAny,0},
 {"nk",-1,5,NULL,1,1,_ZeroNk,&NuMolNk[5],&KappaNk[5],phDrops5,"Концентрация капель [5]", rsPositive,0},
 {"rok",-1,5,NULL,1,1,_ZeroRok,&NuMolRok[5],&KappaRok[5],phDrops5,"Плотность капель [5]", rsPositive,0},
 {"so2k",-1,5,NULL,1,1,_ZeroSO2k,&NuMolSO2k[5],&KappaSO2k[5],phDrops5,"Концентрация SO2 в каплях [5]", rsPositive,0},
 {"uk",0,5,DirectSolver,1,0,_ZeroUkx,&NuMolUkx[5],&KappaUkx[5],phDrops5,"Скорость капель [5]", rsAny,0},
 {"uk",1,5,DirectSolver,1,0,_ZeroUky,&NuMolUky[5],&KappaUky[5],phDrops5,"Скорость капель [5]", rsAny,0},
 {"uk",2,5,DirectSolver,1,0,_ZeroUkz,&NuMolUkz[5],&KappaUkz[5],phDrops5,"Скорость капель [5]", rsAny,0},
 {"nk",-1,6,NULL,1,1,_ZeroNk,&NuMolNk[6],&KappaNk[6],phDrops6,"Концентрация капель [6]", rsPositive,0},
 {"rok",-1,6,NULL,1,1,_ZeroRok,&NuMolRok[6],&KappaRok[6],phDrops6,"Плотность капель [6]", rsPositive,0},
 {"so2k",-1,6,NULL,1,1,_ZeroSO2k,&NuMolSO2k[6],&KappaSO2k[6],phDrops6,"Концентрация SO2 в каплях [6]", rsPositive,0},
 {"uk",0,6,DirectSolver,1,0,_ZeroUkx,&NuMolUkx[6],&KappaUkx[6],phDrops6,"Скорость капель [6]", rsAny,0},
 {"uk",1,6,DirectSolver,1,0,_ZeroUky,&NuMolUky[6],&KappaUky[6],phDrops6,"Скорость капель [6]", rsAny,0},
 {"uk",2,6,DirectSolver,1,0,_ZeroUkz,&NuMolUkz[6],&KappaUkz[6],phDrops6,"Скорость капель [6]", rsAny,0},
 {"e",-1,-1,NULL,1,0,_ZeroEt,&D,&KappaEt,phMainPhase,"Скорость диссипации", rsPositive,0},
 {"f0",-1,-1,SolarSolver,1,1,_ZeroF0,&NuMolF0,&KappaF0,phMainPhase,"Интенсивность VIS-прямого излучения", rsAny,0},
 {"fir",-1,-1,SolarSolver,1,1,_ZeroFIR,&NuMolFIR,&KappaFIR,phMainPhase,"Интенсивность IR-прямого излучения", rsAny,0},
 {"fuv",-1,-1,SolarSolver,1,1,_ZeroFUV,&NuMolFUV,&KappaFUV,phMainPhase,"Интенсивность UV-прямого излучения", rsAny,0},
 /*!!!rsPositive*/{ "i", -1, -1, PoissonSolver, 1, 1, _ZeroI, &NuMolI, &KappaI, phMainPhase, "Интенсивность VIS-диффузного излучения", rsPositive, 0 },
 /*!!!rsPositive*/{"ir",-1,-1,PoissonSolver,1,1,_ZeroIR,&NuMolIR,&KappaIR,phMainPhase,"Интенсивность IR-диффузного излучения", rsPositive,0},
 /*!!!rsPositive*/{ "irt", -1, -1, PoissonSolver, 1, 1, _ZeroIRt, &NuMolIRt, &KappaIRt, phMainPhase, "Интенсивность IR-теплового излучения", rsPositive, 0 },
 {"k",-1,-1,NULL,1,0,_ZeroKt,&D,&KappaKt,phMainPhase,"Энергия пульсаций", rsPositive,0},
 {"nu",-1,-1,DirectSolver,1,0,_ZeroNu,&NuMolNu,&KappaNu,phMainPhase,"Турбулентная вязкость", rsPositive,0},
 {"p",-1,-1,PoissonSolver,1,0,_ZeroP,&NuMolP,&KappaP,phMainPhase,"Давление", rsAny,0},
 {"t",-1,-1,NULL,1,0,_ZeroT,&NuMolT,&KappaT,phMainPhase,"Температура (градусы Цельсия)", rsAny,0},
 /*!!!rsPositive*/{ "uv", -1, -1, PoissonSolver, 1, 1, _ZeroUV, &NuMolUV, &KappaUV, phMainPhase, "Интенсивность UV-диффузного излучения", rsPositive, 0 },
 {"v",0,-1,NULL,1,0,_ZeroUx,&D,&KappaUx,phMainPhase,"Скорость", rsAny,0},
 {"v",1,-1,NULL,1,0,_ZeroUy,&D,&KappaUy,phMainPhase,"Скорость", rsAny,0},
 {"v",2,-1,NULL,1,0,_ZeroUz,&D,&KappaUz,phMainPhase,"Скорость", rsAny,0}
};

int PhaseLinks[NumPhases];

enum {refSO2};

char * SubstRefs[MaxActSubst] = {"SO2"};

char   SLinks[MaxActSubst] = {0};

void AddCustomVars()
{
 AddVar(&tAirAlbedo,"tAirAlbedo",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&tAirBeta,"tAirBeta",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&rAirAsymmetry,"rAirAsymmetry",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&tAirAsymmetry,"tAirAsymmetry",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&SurfaceLambda,"SurfaceLambda",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&SurfaceCC,"SurfaceCC",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&SurfaceRO,"SurfaceRO",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&rAirBeta,"rAirBeta",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&uAirAsymmetry,"uAirAsymmetry",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&AirBeta,"AirBeta",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&RO,"RO",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&AirAlbedo,"AirAlbedo",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&AirAsymmetry,"AirAsymmetry",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&uAirAlbedo,"uAirAlbedo",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&uAirBeta,"uAirBeta",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&ZeroT,"ZeroT",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&rAirAlbedo,"rAirAlbedo",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&Ky,"Ky",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&uSolarIntensity,"uSolarIntensity",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&SurfaceAlbedo,"SurfaceAlbedo",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&D,"D",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&Kx,"Kx",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&BT,"BT",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&tDiffuseFlux,"tDiffuseFlux",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&uDiffuseFlux,"uDiffuseFlux",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&T0,"T0",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&rDiffuseFlux,"rDiffuseFlux",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&SurfaceEmissivity,"SurfaceEmissivity",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&rSolarIntensity,"rSolarIntensity",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(IntPlanck,"IntPlanck",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&DiffuseFlux,"DiffuseFlux",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&SolarIntensity,"SolarIntensity",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddSection("Drops",1);
 AddVar(rgC0,"rgC0",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rQaC5,"rQaC5",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rgC1,"rgC1",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rgC3,"rgC3",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rQaC4,"rQaC4",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rgC5,"rgC5",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rgC4,"rgC4",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rgC2,"rgC2",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rQeC2,"rQeC2",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rQeC4,"rQeC4",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rQeC3,"rQeC3",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tQeC0,"tQeC0",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rQeC5,"rQeC5",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rQaC0,"rQaC0",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rQaC2,"rQaC2",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rQaC1,"rQaC1",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rQaC3,"rQaC3",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tgC4,"tgC4",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tgC1,"tgC1",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tgC0,"tgC0",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tQaC5,"tQaC5",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tgC2,"tgC2",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tgC3,"tgC3",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tgC5,"tgC5",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rQeC1,"rQeC1",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tQaC4,"tQaC4",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tQaC3,"tQaC3",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tQeC4,"tQeC4",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tQeC3,"tQeC3",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tQeC2,"tQeC2",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tQeC5,"tQeC5",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tQaC0,"tQaC0",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tQaC2,"tQaC2",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tQaC1,"tQaC1",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(tQeC1,"tQeC1",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(uQaC2,"uQaC2",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(QeC4,"QeC4",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(QeC3,"QeC3",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(QeC2,"QeC2",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(QeC1,"QeC1",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(QeC5,"QeC5",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(QaC0,"QaC0",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(QaC4,"QaC4",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(QaC3,"QaC3",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(QaC2,"QaC2",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(QaC1,"QaC1",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(QeC0,"QeC0",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&AlphaT,"AlphaT",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&InjectV,"InjectV",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&InjectStartTime,"InjectStartTime",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(InitRok,"InitRok",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(UwLight,"UwLight",fltT,7,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&InjectCellLength,"InjectCellLength",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(InjectN,"InjectN",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&AlphaM,"AlphaM",fltT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(InitSO2k,"InitSO2k",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(InitNk,"InitNk",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(Dk,"Dk",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(QaC5,"QaC5",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(gC0,"gC0",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(uQaC5,"uQaC5",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(uQaC4,"uQaC4",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(uQaC3,"uQaC3",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(uQaC1,"uQaC1",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(ugC0,"ugC0",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(ugC1,"ugC1",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(ugC5,"ugC5",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(ugC4,"ugC4",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(ugC3,"ugC3",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(ugC2,"ugC2",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(uQaC0,"uQaC0",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(uQeC5,"uQeC5",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(gC4,"gC4",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(gC3,"gC3",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(gC2,"gC2",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(gC1,"gC1",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(gC5,"gC5",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(uQeC0,"uQeC0",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(uQeC4,"uQeC4",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(uQeC3,"uQeC3",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(uQeC2,"uQeC2",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(uQeC1,"uQeC1",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(rQeC0,"rQeC0",fltT,50,0,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 nPoissons = 5;
}

#ifdef __PARALLEL__
 #define OffsBuf(S) (((NZ+2)*_Num##S+1)*NX*NY)
 #define DeclareRef(S) float * S = &V[OffsBuf(S)]
 #define DeclareDn(S,i) float * S = &V[((NZ+2)*(_NumDn+i)+1)*NX*NY]
 #define DeclareUpDown(S) \
           float * DH##S = z==0   ? &V[OffsBuf(S)-NY*NX] : &S[(z-1)*NY*NX]; \
           float * UH##S = z==NZs ? &V[OffsBuf(S)+NZ*NY*NX] : &S[(z+1)*NY*NX];
#else
 #define DeclareRef(S) float * S = V->Name[_Num##S]
 #define DeclareDn(S,i) float * S = V->Dn[i]
 #define DeclareUpDown(S) \
           float * DH##S = z==0   ? &S[NZs*NY*NX] : &S[(z-1)*NY*NX]; \
           float * UH##S = z==NZs ? &S[0] : &S[(z+1)*NY*NX];
#endif

#define vector(N) for (i=0; i<N; i++)

#ifdef __MPI2REENT__
namespace UPDOWN {
  float * DH__DIV__[__NPROCS__];
  float * UH__DIV__[__NPROCS__];
};
using namespace UPDOWN;
#endif

#ifdef __PARALLEL__
void CalculateK(WKoeffs * W, float ** Bounds, float * V, float ** Kf, float ** KDn,  float ** Sf, float ** SDn, int Delta, unsigned char * Area, unsigned char * Boundaries, MapItem * Maps)
{
#else
void CalculateK(WKoeffs * W, StoreStruct * V, float ** Kf, float ** KDn, float ** Sf, float ** SDn, unsigned char * Area, unsigned char * Boundaries, MapItem * Maps)
{
 int Delta = 0;
 float ** Bounds = V->Bounds;
#endif
 DeclareRef(Nk0);
 DeclareRef(Rok0);
 DeclareRef(SO2k0);
 DeclareRef(Ukx0);
 DeclareRef(Uky0);
 DeclareRef(Ukz0);
 DeclareRef(Nk1);
 DeclareRef(Rok1);
 DeclareRef(SO2k1);
 DeclareRef(Ukx1);
 DeclareRef(Uky1);
 DeclareRef(Ukz1);
 DeclareRef(Nk2);
 DeclareRef(Rok2);
 DeclareRef(SO2k2);
 DeclareRef(Ukx2);
 DeclareRef(Uky2);
 DeclareRef(Ukz2);
 DeclareRef(Nk3);
 DeclareRef(Rok3);
 DeclareRef(SO2k3);
 DeclareRef(Ukx3);
 DeclareRef(Uky3);
 DeclareRef(Ukz3);
 DeclareRef(Nk4);
 DeclareRef(Rok4);
 DeclareRef(SO2k4);
 DeclareRef(Ukx4);
 DeclareRef(Uky4);
 DeclareRef(Ukz4);
 DeclareRef(Nk5);
 DeclareRef(Rok5);
 DeclareRef(SO2k5);
 DeclareRef(Ukx5);
 DeclareRef(Uky5);
 DeclareRef(Ukz5);
 DeclareRef(Nk6);
 DeclareRef(Rok6);
 DeclareRef(SO2k6);
 DeclareRef(Ukx6);
 DeclareRef(Uky6);
 DeclareRef(Ukz6);
 DeclareRef(Et);
 DeclareRef(F0);
 DeclareRef(FIR);
 DeclareRef(FUV);
 DeclareRef(I);
 DeclareRef(IR);
 DeclareRef(IRt);
 DeclareRef(Kt);
 DeclareRef(Nu);
 DeclareRef(P);
 DeclareRef(T);
 DeclareRef(UV);
 DeclareRef(Ux);
 DeclareRef(Uy);
 DeclareRef(Uz);
 DeclareDn(Source_Drops0,PhaseLinks[phDrops0]);
 DeclareDn(Source_Drops1,PhaseLinks[phDrops1]);
 DeclareDn(Source_Drops2,PhaseLinks[phDrops2]);
 DeclareDn(Source_Drops3,PhaseLinks[phDrops3]);
 DeclareDn(Source_Drops4,PhaseLinks[phDrops4]);
 DeclareDn(Source_Drops5,PhaseLinks[phDrops5]);
 DeclareDn(Source_Drops6,PhaseLinks[phDrops6]);
 float * K_Rok[7] = {Kf[_NumRok0],Kf[_NumRok1],Kf[_NumRok2],Kf[_NumRok3],Kf[_NumRok4],Kf[_NumRok5],Kf[_NumRok6]};
 float * Rok[7] = {Rok0,Rok1,Rok2,Rok3,Rok4,Rok5,Rok6};
 float * S_Rok[7] = {Sf[_NumRok0],Sf[_NumRok1],Sf[_NumRok2],Sf[_NumRok3],Sf[_NumRok4],Sf[_NumRok5],Sf[_NumRok6]};
 float * K_Ukx[7] = {Kf[_NumUkx0],Kf[_NumUkx1],Kf[_NumUkx2],Kf[_NumUkx3],Kf[_NumUkx4],Kf[_NumUkx5],Kf[_NumUkx6]};
 float * Ukx[7] = {Ukx0,Ukx1,Ukx2,Ukx3,Ukx4,Ukx5,Ukx6};
 float * K_Uky[7] = {Kf[_NumUky0],Kf[_NumUky1],Kf[_NumUky2],Kf[_NumUky3],Kf[_NumUky4],Kf[_NumUky5],Kf[_NumUky6]};
 float * Uky[7] = {Uky0,Uky1,Uky2,Uky3,Uky4,Uky5,Uky6};
 float * K_Ukz[7] = {Kf[_NumUkz0],Kf[_NumUkz1],Kf[_NumUkz2],Kf[_NumUkz3],Kf[_NumUkz4],Kf[_NumUkz5],Kf[_NumUkz6]};
 float * Ukz[7] = {Ukz0,Ukz1,Ukz2,Ukz3,Ukz4,Ukz5,Ukz6};
 float * K_Nk[7] = {Kf[_NumNk0],Kf[_NumNk1],Kf[_NumNk2],Kf[_NumNk3],Kf[_NumNk4],Kf[_NumNk5],Kf[_NumNk6]};
 float * Nk[7] = {Nk0,Nk1,Nk2,Nk3,Nk4,Nk5,Nk6};
 float * S_Nk[7] = {Sf[_NumNk0],Sf[_NumNk1],Sf[_NumNk2],Sf[_NumNk3],Sf[_NumNk4],Sf[_NumNk5],Sf[_NumNk6]};
 float * K_SO2k[7] = {Kf[_NumSO2k0],Kf[_NumSO2k1],Kf[_NumSO2k2],Kf[_NumSO2k3],Kf[_NumSO2k4],Kf[_NumSO2k5],Kf[_NumSO2k6]};
 float * SO2k[7] = {SO2k0,SO2k1,SO2k2,SO2k3,SO2k4,SO2k5,SO2k6};
 float * S_SO2k[7] = {Sf[_NumSO2k0],Sf[_NumSO2k1],Sf[_NumSO2k2],Sf[_NumSO2k3],Sf[_NumSO2k4],Sf[_NumSO2k5],Sf[_NumSO2k6]};
 float * K_SO2 = KDn[SLinks[refSO2]];
 float * S_SO2 = SDn[SLinks[refSO2]];
 DeclareDn(SO2,SLinks[refSO2]);
 int i,x,zy,strata,Ptr;

 int BSize = NY*NX*(NZ-2*Delta)*sizeof(float);
 for (i=0; i<NumEqs; i++)
     {
      if (Kf[i]) memset(&(Kf[i][Delta*NY*NX]),0,BSize);
      if (Sf[i]) memset(&(Sf[i][Delta*NY*NX]),0,BSize);
     }
 for (i=0; i<NSubst; i++)
     {
      if (KDn[i]) memset(&(KDn[i][Delta*NY*NX]),0,BSize);
      if (SDn[i]) memset(&(SDn[i][Delta*NY*NX]),0,BSize);
     }
 for (strata=0; strata<=1; strata++)
  #pragma omp parallel if(UseOpenMP)
  #pragma omp for schedule(dynamic,imax(1,NY*NZ/(4*nSMP))) private(zy,x,i,Ptr)
    for (zy=Delta*NY; zy<NY*(NZ-Delta); zy++)
     {
      int z = zy/NY;
      int y = zy%NY;
#ifdef __MPI2REENT__
      DH__DIV__[__id__] = z==0   ? &DIV[NZs*NY*NX] : &DIV[(z-1)*NY*NX];
      UH__DIV__[__id__] = z==NZs ? &DIV[0] : &DIV[(z+1)*NY*NX];
#else
      float * DHDIV = z==0   ? &DIV[NZs*NY*NX] : &DIV[(z-1)*NY*NX];
      float * UHDIV = z==NZs ? &DIV[0] : &DIV[(z+1)*NY*NX];
#endif
      DeclareUpDown(I)
      DeclareUpDown(IR)
      DeclareUpDown(IRt)
      DeclareUpDown(Nu)
      DeclareUpDown(P)
      DeclareUpDown(T)
      DeclareUpDown(UV)
      DeclareUpDown(Ux)
      DeclareUpDown(Uy)
      DeclareUpDown(Uz)
      DeclareUpDown(Ukz0)
      DeclareUpDown(Ukz1)
      DeclareUpDown(Ukz2)
      DeclareUpDown(Ukz3)
      DeclareUpDown(Ukz4)
      DeclareUpDown(Ukz5)
      DeclareUpDown(Ukz6)
      float * DHUkz[7] = {DHUkz0,DHUkz1,DHUkz2,DHUkz3,DHUkz4,DHUkz5,DHUkz6};
      float * UHUkz[7] = {UHUkz0,UHUkz1,UHUkz2,UHUkz3,UHUkz4,UHUkz5,UHUkz6};

      for (x=(z+y+strata)&1, Ptr = zy*NX+((z+y+strata)&1); x<NX; Ptr+=2,x+=2)
          if (Area[Ptr]!=1)
              {
               
		double ddRo[7*(7+1)]; /* NEW */
		double ddSO2kg[7]; /* NEW */
		double ddSO2kf[7]; /* NEW */
               double dUxdx;
               double dUydy;
               double dUzdz;
               double dUxdy;
               double dUydz;
               double dUzdx;
               double dUxdz;
               double dUydx;
               double dUzdy;
               double SSS2;
               double DDD;
               double RNG;
               double RNG3;
               double Pi;
               double ukAlbAsm;
               double utAlpha;
               double utBeta;
               double u_Ft;
               double uFt;
               double kAlbAsm;
               double tAlpha;
               double tBeta;
               double _Ft;
               double Ft;
               double rkAlbAsm;
               double rtAlpha;
               double rtBeta;
               double r_Ft;
               double rFt;
               double tkAlbAsm;
               double ttAlpha;
               double ttBeta;
               double tPlanck;
               double tFt;
               double TK;
               double RO0;
               double CC0;
               double dUkxdx[7];
               double dUkydy[7];
               double dUkzdz[7];
               double dRok[7];
			   double YCell; /* !!!!! REMOVE */
               double dNk[7];
               double H298_SO2;
               double dH_SO2;
               double _DfSO2;
               double ScSO2;
               double tScSO2;
               double Nuss1SO2[7];
               double _HenrySO2[7];
               double LSO2f[7];
               double LSO2g[7];
               double dSO2kf[7];
               double dSO2kg[7];
               double dSO2;
               double KScale;
               double Uref;
               double Z0;
               double ZZ;
               double K0;
               double K1;
               double NuK;
               double NuE;
               double uk1;
               double uk1a;
               double uFluxPi;
               double k1;
               double k1a;
               double FluxPi;
               double rk1;
               double rk1a;
               double rFluxPi;
               double tk1;
               double tk1a;
               double tFluxPi;
               double Mk;
               double Ma;
               double _R;
               double ROK;
               double Dvap0;
               double Lambda0;
               double Dvap1;
               double Lambda1;
               double _Hvap;
               double _Ck;
               double PkF;
               double _Sigma;
               double PkF1;
               double EMPTY[7];
               double NMoles[7];
               double PkF2[7];
               double mCX[7];
               double mCY[7];
               double Mode[7];
               double SPEC1[7];
               double SPEC2[7];
               double State[7];
               double CND1[7];
               double CND2[7];
               double X[7];
               double Y[7];
               double dCX[7];
               double dCY[7];
			   double DC[7]; /*!!!!*/
			   double _XYC[7]; /*!!!*/
			   double XYC[7]; /*!!!*/
			   double XYC2[7]; /*!!!*/
			   double XYC3[7]; /*!!!*/
			   double XYC4[7]; /*!!!*/
			   double XYC5[7]; /*!!!*/
			   double XYC6[7]; /*!!!*/
			   double XYC7[7]; /*!!!*/
			   double _An[7]; /*!!!*/
			   double _Bn[7]; /*!!!*/
			   double _Cn[7]; /*!!!*/
			   double _Dn[7]; /*!!!*/
			   double __An[7]; /*!!!*/
			   double __Bn[7]; /*!!!*/
			   double Cn[7]; /*!!!*/
			   double Dn[7]; /*!!!*/
			   double _Ann[7]; /*!!!*/
			   double _Bnn[7]; /*!!!*/
			   double _Cnn[7]; /*!!!*/
			   double _Dnn[7]; /*!!!*/
			   double Ann[7]; /*!!!*/
			   double Bnn[7]; /*!!!*/
			   double Cnn[7]; /*!!!*/
			   double Dnn[7]; /*!!!*/
			   double CORR_ABCDn[7]; /*!!!*/
			   double CORR_ABCDnn[7]; /*!!!*/
			   double XC1[7]; /*!!!*/
			   double XC2[7]; /*!!!*/
			   double XC3[7]; /*!!!*/
			   double XC4[7]; /*!!!*/
			   double XC5[7]; /*!!!*/
			   double XC6[7]; /*!!!*/
			   double XC7[7]; /*!!!*/
			   double XC8[7]; /*!!!*/
			   double XC9[7]; /*!!!*/
			   double XC10[7]; /*!!!*/
			   double XC11[7]; /*!!!*/
			   double CY1[7]; /*!!!*/
			   double CY2[7]; /*!!!*/
			   double CY3[7]; /*!!!*/
			   double CY4[7]; /*!!!*/
			   double CY5[7]; /*!!!*/
			   double CY6[7]; /*!!!*/
			   double CY7[7]; /*!!!*/
			   double CY8[7]; /*!!!*/
			   double CY9[7]; /*!!!*/
			   double CY10[7]; /*!!!*/
			   double CY11[7]; /*!!!*/
			   double coalesc; /*!!!!*/
			   double ddN[7]; /*!!!!*/
			   double dRo[7]; /*!!!!*/
			   double Nc = 7; /*???????????????????????????????!!!*/
               double X2[7];
               double X3[7];
               double X4[7];
               double X5[7];
               double X6[7];
               double X7[7];
               double X8[7];
               double X9[7];
               double Y2[7];
               double Y3[7];
               double Y4[7];
               double Y5[7];
               double Y6[7];
               double Y7[7];
               double Y8[7];
               double Y9[7];
               double XY1[7];
               double XY2[7];
               double XY3[7];
               double XY4[7];
               double XY5[7];
               double XY6[7];
               double XY7[7];
               double XY8[7];
               double XY9[7];
               double A1[7];
               double B1[7];
               double A2[7];
               double B2[7];
               double nC2[7];
               double Determ[7];
               double AnC[7];
               double BnC[7];
               double CORR[7];
               double An[7];
               double Bn[7];
               double ieC0[7];
               double ieC1[7];
               double ieC2[7];
               double ieC3[7];
               double ieC4[7];
               double ieC5[7];
               double iaC0[7];
               double iaC1[7];
               double iaC2[7];
               double iaC3[7];
               double iaC4[7];
               double iaC5[7];
               double igC0[7];
               double igC1[7];
               double igC2[7];
               double igC3[7];
               double igC4[7];
               double igC5[7];
               double uieC0[7];
               double uieC1[7];
               double uieC2[7];
               double uieC3[7];
               double uieC4[7];
               double uieC5[7];
               double uiaC0[7];
               double uiaC1[7];
               double uiaC2[7];
               double uiaC3[7];
               double uiaC4[7];
               double uiaC5[7];
               double uigC0[7];
               double uigC1[7];
               double uigC2[7];
               double uigC3[7];
               double uigC4[7];
               double uigC5[7];
               double rieC0[7];
               double rieC1[7];
               double rieC2[7];
               double rieC3[7];
               double rieC4[7];
               double rieC5[7];
               double riaC0[7];
               double riaC1[7];
               double riaC2[7];
               double riaC3[7];
               double riaC4[7];
               double riaC5[7];
               double rigC0[7];
               double rigC1[7];
               double rigC2[7];
               double rigC3[7];
               double rigC4[7];
               double rigC5[7];
               double tieC0[7];
               double tieC1[7];
               double tieC2[7];
               double tieC3[7];
               double tieC4[7];
               double tieC5[7];
               double tiaC0[7];
               double tiaC1[7];
               double tiaC2[7];
               double tiaC3[7];
               double tiaC4[7];
               double tiaC5[7];
               double tigC0[7];
               double tigC1[7];
               double tigC2[7];
               double tigC3[7];
               double tigC4[7];
               double tigC5[7];
               double iBeta[7];
               double uiBeta[7];
               double riBeta[7];
               double tiBeta[7];
               double iAbs[7];
               double uiAbs[7];
               double riAbs[7];
               double tiAbs[7];
               double iAsym[7];
               double uiAsym[7];
               double riAsym[7];
               double tiAsym[7];
               double iAlbedo[7];
               double uiAlbedo[7];
               double riAlbedo[7];
               double tiAlbedo[7];
               double eBeta;
               double eScat;
               double kAlbedo;
               double kAsymm;
               double ueBeta;
               double ueScat;
               double ukAlbedo;
               double ukAsymm;
               double reBeta;
               double reScat;
               double rkAlbedo;
               double rkAsymm;
               double teBeta;
               double teScat;
               double tkAlbedo;
               double tkAsymm;
               double kCSubst[7];
			   double _DD[7]; /*!!!*/
			   double DD[7];
			   double dUkx[7];
               double dUky[7];
               double dUkz[7];
               double dUk[7];
               double Rek[7];
               double sRek[7];
               double Prk[7];
               double Sck[7];
               double Nuss1[7];
               double Nuss2[7];
               double PkF3[7];
               double PkF4[7];
               double nX[7];
               double nY[7];
               double LX[7];
               double LY[7];
               double Dz[7];
               double LDz[7];
               double nDz[7];
               double Dzz[7];
               double LDzz[7];
               double nDzz[7];
               double X1[7];
               double Y1[7];
               double LDp[7];
               double Dz2[7];
               double Dz3[7];
               double Dz4[7];
               double Dz5[7];
			   double Dz6[7];/*!!!*/
			   double Dz7[7];/*!!!!*/
			   double Dzz2[7];
               double Dzz3[7];
               double Dzz4[7];
               double Dzz5[7];
			   double Dzz6[7];/*!!!*/
			   double Dzz7[7];/*!!!*/
			   double CTOZERO[7];/*****/
			   double CTOMAX[7]; /*****/
			   double maxSO2k[7]; /*!!!!!*/
			   double RZ[7];
               double PhiM[7];
               double PhiP[7];
			   double _PhiM[7];/*!!!*/
			   double _PhiP[7];/*!!!*/
			   double ppX[7];/*!!!*/
			   double ppY[7];/*!!!!*/
			   double PhiZ[7];
               double RM[7];
               double RP[7];
               double RokM[7];
               double RokP[7];
               double NkM[7];
               double NkP[7];
               double pdNk[7];
               double TOZERO[7];
               double dC;

               dUxdx = ((FastMode ? (0.0) : 0.0)+(SlowMode ? (dFdx(Ux)) : 0.0));
               dUydy = ((FastMode ? (0.0) : 0.0)+(SlowMode ? (dFdy(Uy)) : 0.0));
               dUzdz = ((FastMode ? (0.0) : 0.0)+(SlowMode ? (dFdz(Uz)) : 0.0));
               
               dUxdy = ((FastMode ? (0.0) : 0.0)+(SlowMode ? (dFdy(Ux)) : 0.0));
               dUydz = ((FastMode ? (0.0) : 0.0)+(SlowMode ? (dFdz(Uy)) : 0.0));
               dUzdx = ((FastMode ? (0.0) : 0.0)+(SlowMode ? (dFdx(Uz)) : 0.0));
               
               dUxdz = ((FastMode ? (0.0) : 0.0)+(SlowMode ? (dFdz(Ux)) : 0.0));
               dUydx = ((FastMode ? (0.0) : 0.0)+(SlowMode ? (dFdx(Uy)) : 0.0));
               dUzdy = ((FastMode ? (0.0) : 0.0)+(SlowMode ? (dFdy(Uz)) : 0.0));
               
               Pi = 3.1415926535897932;
               
               TK = 273.15+T[Ptr];
               
               vector(7) dUkxdx[i] = ((FastMode ? (0.0) : 0.0)+(SlowMode ? (dFdx(Ukx[i])) : 0.0));
               vector(7) dUkydy[i] = ((FastMode ? (0.0) : 0.0)+(SlowMode ? (dFdy(Uky[i])) : 0.0));
               vector(7) dUkzdz[i] = ((FastMode ? (0.0) : 0.0)+(SlowMode ? (dFdz(Ukz[i])) : 0.0));
               
			   YCell = floor((ModelTime - InjectStartTime)*InjectV / InjectCellLength); /* !!! REMOVE */
               
               H298_SO2 = 0.012177;
               dH_SO2 = -25942.8;
               
               KScale = 0.208816;
               Z0 = 0.05;
               Mk = 0.018;
               Ma = 0.0288;
               _R = 8.315;
               ROK = 1000.0;
               _Sigma = 0.0761-1.55E-4*T[Ptr];
               
               vector(7) ieC0[i] = QeC0[i];
               vector(7) ieC1[i] = QeC1[i]*1E6;
               vector(7) ieC2[i] = QeC2[i]*1E12;
               vector(7) ieC3[i] = QeC3[i]*1E18;
               vector(7) ieC4[i] = QeC4[i]*1E24;
               vector(7) ieC5[i] = QeC5[i]*1E30;
               vector(7) iaC0[i] = QaC0[i];
               vector(7) iaC1[i] = QaC1[i]*1E6;
               vector(7) iaC2[i] = QaC2[i]*1E12;
               vector(7) iaC3[i] = QaC3[i]*1E18;
               vector(7) iaC4[i] = QaC4[i]*1E24;
               vector(7) iaC5[i] = QaC5[i]*1E30;
               vector(7) igC0[i] = gC0[i];
               vector(7) igC1[i] = gC1[i]*1E6;
               vector(7) igC2[i] = gC2[i]*1E12;
               vector(7) igC3[i] = gC3[i]*1E18;
               vector(7) igC4[i] = gC4[i]*1E24;
               vector(7) igC5[i] = gC5[i]*1E30;
               
               vector(7) uieC0[i] = uQeC0[i];
               vector(7) uieC1[i] = uQeC1[i]*1E6;
               vector(7) uieC2[i] = uQeC2[i]*1E12;
               vector(7) uieC3[i] = uQeC3[i]*1E18;
               vector(7) uieC4[i] = uQeC4[i]*1E24;
               vector(7) uieC5[i] = uQeC5[i]*1E30;
               vector(7) uiaC0[i] = uQaC0[i];
               vector(7) uiaC1[i] = uQaC1[i]*1E6;
               vector(7) uiaC2[i] = uQaC2[i]*1E12;
               vector(7) uiaC3[i] = uQaC3[i]*1E18;
               vector(7) uiaC4[i] = uQaC4[i]*1E24;
               vector(7) uiaC5[i] = uQaC5[i]*1E30;
               vector(7) uigC0[i] = ugC0[i];
               vector(7) uigC1[i] = ugC1[i]*1E6;
               vector(7) uigC2[i] = ugC2[i]*1E12;
               vector(7) uigC3[i] = ugC3[i]*1E18;
               vector(7) uigC4[i] = ugC4[i]*1E24;
               vector(7) uigC5[i] = ugC5[i]*1E30;
               
               vector(7) rieC0[i] = rQeC0[i];
               vector(7) rieC1[i] = rQeC1[i]*1E6;
               vector(7) rieC2[i] = rQeC2[i]*1E12;
               vector(7) rieC3[i] = rQeC3[i]*1E18;
               vector(7) rieC4[i] = rQeC4[i]*1E24;
               vector(7) rieC5[i] = rQeC5[i]*1E30;
               vector(7) riaC0[i] = rQaC0[i];
               vector(7) riaC1[i] = rQaC1[i]*1E6;
               vector(7) riaC2[i] = rQaC2[i]*1E12;
               vector(7) riaC3[i] = rQaC3[i]*1E18;
               vector(7) riaC4[i] = rQaC4[i]*1E24;
               vector(7) riaC5[i] = rQaC5[i]*1E30;
               vector(7) rigC0[i] = rgC0[i];
               vector(7) rigC1[i] = rgC1[i]*1E6;
               vector(7) rigC2[i] = rgC2[i]*1E12;
               vector(7) rigC3[i] = rgC3[i]*1E18;
               vector(7) rigC4[i] = rgC4[i]*1E24;
               vector(7) rigC5[i] = rgC5[i]*1E30;
               
               vector(7) tieC0[i] = tQeC0[i];
               vector(7) tieC1[i] = tQeC1[i]*1E6;
               vector(7) tieC2[i] = tQeC2[i]*1E12;
               vector(7) tieC3[i] = tQeC3[i]*1E18;
               vector(7) tieC4[i] = tQeC4[i]*1E24;
               vector(7) tieC5[i] = tQeC5[i]*1E30;
               vector(7) tiaC0[i] = tQaC0[i];
               vector(7) tiaC1[i] = tQaC1[i]*1E6;
               vector(7) tiaC2[i] = tQaC2[i]*1E12;
               vector(7) tiaC3[i] = tQaC3[i]*1E18;
               vector(7) tiaC4[i] = tQaC4[i]*1E24;
               vector(7) tiaC5[i] = tQaC5[i]*1E30;
               vector(7) tigC0[i] = tgC0[i];
               vector(7) tigC1[i] = tgC1[i]*1E6;
               vector(7) tigC2[i] = tgC2[i]*1E12;
               vector(7) tigC3[i] = tgC3[i]*1E18;
               vector(7) tigC4[i] = tgC4[i]*1E24;
               vector(7) tigC5[i] = tgC5[i]*1E30;
               vector(7) dUkx[i] = Ux[Ptr]-Ukx[i][Ptr];
               vector(7) dUky[i] = Uy[Ptr]-Uky[i][Ptr];
               vector(7) dUkz[i] = Uz[Ptr]-Ukz[i][Ptr];
               
               SSS2 =  2*(dUxdx*dUxdx+dUydy*dUydy+dUzdz*dUzdz)+
                              (dUxdy+dUydx)*(dUxdy+dUydx)+
                              (dUxdz+dUzdx)*(dUxdz+dUzdx)+
                              (dUydz+dUzdy)*(dUydz+dUzdy);
               tPlanck = (IntPlanck[0]+TK*(IntPlanck[1]+TK*(IntPlanck[2]+TK*IntPlanck[3])));
               RO0 = 352.984/TK;
               _DfSO2 = 5.518E-10*pow(TK,1.75);
               Uref = KScale*0.5263;
               ZZ = z==0 && _Z[0]==0.0 ? Z0 : (_Z[z]+_Z[z-1])*0.05+Z0;
               uFluxPi = uDiffuseFlux/Pi;
               FluxPi = DiffuseFlux/Pi;
               rFluxPi = rDiffuseFlux/Pi;
               tFluxPi = tDiffuseFlux/Pi;
               Dvap0 = 0.211E-4*pow(TK/273.0,1.94);
               Lambda0 = 1E-3*(4.39+0.071*TK);
               _Hvap = 2.5E+6*pow(273.15/TK,0.167+3.67E-4*TK);
               _Ck = 4185.5*(0.99618+0.0002874*pow((TK-173)/100.0,5.26)+0.01116*pow(10.0,-0.036*T[Ptr]));
               PkF = (610.7799961+T[Ptr]*(44.36518521+T[Ptr]*(1.428945805+T[Ptr]*(2.650648471E-2+T[Ptr]*(3.031240396E-4+T[Ptr]*(2.034080948E-6+T[Ptr]*6.136820929E-9))))))/_R/TK;
               PkF1 = 4.0*Mk*_Sigma/_R/TK/ROK;
               
			   /*!!!*/vector(7) EMPTY[i] = Nk[i][Ptr]<0.1*TAU || Rok[i][Ptr] / Nk[i][Ptr] <= Pi / 6.0*ROK*(Dk[i] * Dk[i] * Dk[i]);
               /*!!!!*/vector(7) dUk[i] = fabs(UwLight[i]) > 1E-8 ? 0.0 : sqrt(dUkx[i]*dUkx[i]+dUky[i]*dUky[i]+dUkz[i]*dUkz[i]);
               DDD =  Nu[Ptr]*SSS2;
               RNG  = sqrt(SSS2)*Kt[Ptr]/Et[Ptr];
               CC0 = 1.4*(718+0.1167*RO0);
               
               ScSO2 = D/_DfSO2;
               vector(7) _HenrySO2[i] = EMPTY[i] ? 0.0 : H298_SO2*exp(dH_SO2/_R*(1.0/298.0-1/TK))*_R*TK;
               K0 = Uref*ZZ;
               K1 = log(ZZ/Z0);
               NuK = KScale*KScale*(1.09680767+(-0.008122677347+(0.0002099331724+(-0.000002208352551+(0.00000000906293221+(-9.749368416E-12+(-1.187063604E-14)*ZZ)*ZZ)*ZZ)*ZZ)*ZZ)*ZZ);
               NuE = KScale*KScale*KScale*0.3644/ZZ;
               Dvap1 = 2.0*Dvap0/AlphaM*sqrt(2.0*Pi*Mk/_R/TK);
               /*!!!!!!*/vector(7) NMoles[i] = EMPTY[i] ? 0.0 : 1E6*max(0.0, SO2k[i][Ptr])/Nk[i][Ptr];
               RNG3 = RNG*RNG*RNG;
               tScSO2 = pow(ScSO2,1.0/3.0);
               Lambda1 = 2.0*Lambda0/(AlphaT*RO0*CC0)*sqrt(2*Pi*Ma/_R/TK);
               vector(7) PkF2[i] = EMPTY[i] ? 0.0 : 6.0*NMoles[i]*Mk/Pi/ROK;
               vector(7) mCX[i] = EMPTY[i] ? 0.0 : 1E6*Source_Drops0[Ptr]/PkF-exp(PkF1/Dk[i]-PkF2[i]/(Dk[i]*Dk[i]*Dk[i]));
               vector(7) mCY[i] = EMPTY[i] ? 0.0 : 1E6*Source_Drops0[Ptr]/PkF-exp(PkF1/Dk[i+1]-PkF2[i]/(Dk[i+1]*Dk[i+1]*Dk[i+1]));
               vector(7) Mode[i] = mCX[i]*mCY[i]<=0.0 || i==0 && mCX[i]<0.0 && mCY[i]<0.0 ? 0.0 : mCX[i]<0.0 ? -1.0 : +1.0;
               vector(7) SPEC1[i] = i>0 && i<7-1 && Mode[i-1]>0 && Mode[i]>0 && EMPTY[i+1];
               vector(7) SPEC2[i] = i>0 && i<7-1 && EMPTY[i-1] && Mode[i]<0 && Mode[i+1]<0;
               vector(7) State[i] = EMPTY[i] ? 0.0 :
                                  SPEC1[i] || SPEC2[i] ||
                                  i==0 && Mode[i]>0 ||
                                  i==7-1 && Mode[i]<0 ||
                                  i<7-1 && Mode[i]<0 && Mode[i]!=Mode[i+1] ||
                                  i>0 && Mode[i]>0 && Mode[i]!=Mode[i-1] ? -1.0 : 1.0;
               vector(7) CND1[i] = State[i]<0 && (Mode[i]>0 && !SPEC1[i] || SPEC2[i]);
               vector(7) CND2[i] = State[i]<0 && (Mode[i]<0 && !SPEC2[i] || SPEC1[i]);
               vector(7) nC2[i] = (EMPTY[i] || State[i]<0.0) ? 0.0 : 6.0*Rok[i][Ptr]/Pi/ROK;
               /*!!!!*/vector(7) ppX[i] = EMPTY[i] ? Dk[i] : CND1[i] ? amax2(Dk[i],GetXY(Dk[i+1],24.0*Rok[i][Ptr]/(Pi*ROK*Nk[i][Ptr]))) : Dk[i];
               /*!!!!*/vector(7) ppY[i] = EMPTY[i] ? Dk[i+1] : CND2[i] ? amin2(Dk[i+1],GetXY(Dk[i],24.0*Rok[i][Ptr]/(Pi*ROK*Nk[i][Ptr]))) : Dk[i+1];
			   /*!!!!*/vector(7) X[i] = ppX[i] < ppY[i] ? ppX[i] : 0.25*Dk[i] + 0.75*Dk[i + 1];
			   /*!!!!*/vector(7) Y[i] = ppX[i] < ppY[i] ? ppY[i] : 0.75*Dk[i] + 0.25*Dk[i + 1];
			   vector(7) dCX[i] = Mode[i] == 0 ? 0.0 : CND1[i] ? 1E6*Source_Drops0[Ptr] / PkF - exp(PkF1 / X[i] - PkF2[i] / (X[i] * X[i] * X[i])) : mCX[i];
			   vector(7) dCY[i] = Mode[i] == 0 ? 0.0 : CND2[i] ? 1E6*Source_Drops0[Ptr] / PkF - exp(PkF1 / Y[i] - PkF2[i] / (Y[i] * Y[i] * Y[i])) : mCY[i];
			   /*!!!*/ vector(7) _XYC[i] = EMPTY[i] ? 0.5*(Dk[i]+Dk[i+1]) : /*(i == 0 || (State[i] < 0 && Mode[i]>0)) ? X[i] : (i == Nc - 1 || (State[i] < 0 && Mode[i] < 0)) ? Y[i] :*/ 0.5*(X[i] + Y[i]);
			   vector(7) X2[i] = X[i]*X[i];
               vector(7) Y2[i] = Y[i]*Y[i];
               
               vector(7) XY1[i] = (Y[i]-X[i]);
               vector(7) X3[i] = X2[i]*X[i];
               vector(7) Y3[i] = Y2[i]*Y[i];
               vector(7) XY2[i] = (Y2[i]-X2[i])/2.0;
               vector(7) B1[i] = (EMPTY[i] || State[i]<0.0) ? 0.0 : XY1[i];
               vector(7) X4[i] = X3[i]*X[i];
               vector(7) Y4[i] = Y3[i]*Y[i];
               vector(7) XY3[i] = (Y3[i]-X3[i])/3.0;
               
			   vector(7) A1[i] = (EMPTY[i] || State[i]<0.0) ? 0.0 : XY2[i];
               vector(7) X5[i] = X4[i]*X[i];
               vector(7) Y5[i] = Y4[i]*Y[i];
               vector(7) XY4[i] = (Y4[i]-X4[i])/4.0;
               vector(7) X6[i] = X5[i]*X[i];
               vector(7) Y6[i] = Y5[i]*Y[i];
               vector(7) XY5[i] = (Y5[i]-X5[i])/5.0;
               vector(7) B2[i] = (EMPTY[i] || State[i]<0.0) ? 0.0 : XY4[i];
               vector(7) X7[i] = X6[i]*X[i];
               vector(7) Y7[i] = Y6[i]*Y[i];
               vector(7) XY6[i] = (Y6[i]-X6[i])/6.0;
               vector(7) A2[i] = (EMPTY[i] || State[i]<0.0) ? 0.0 : XY5[i];
               vector(7) X8[i] = X7[i]*X[i];
               vector(7) Y8[i] = Y7[i]*Y[i];
               vector(7) XY7[i] = (Y7[i]-X7[i])/7.0;
               vector(7) Determ[i] = (EMPTY[i] || State[i]<0.0) ? 0.0 : A1[i]*B2[i]-B1[i]*A2[i];
               vector(7) X9[i] = X8[i]*X[i];
               vector(7) Y9[i] = Y8[i]*Y[i];
               vector(7) XY8[i] = (Y8[i]-X8[i])/8.0;
               /******/vector(7) AnC[i] = EMPTY[i] ? 0.0 : State[i]<0.0 ? 0.0 : (B2[i]*Nk[i][Ptr]-B1[i]*nC2[i])/Determ[i];
               /******/vector(7) BnC[i] = EMPTY[i] ? 0.0 : State[i]<0.0 ? Nk[i][Ptr]/(Y[i]-X[i]) : (-A2[i]*Nk[i][Ptr]+A1[i]*nC2[i])/Determ[i];
			   vector(7) XY9[i] = (Y9[i] - X9[i]) / 9.0;
               vector(7) CORR[i] = EMPTY[i] ? 0.0 : AnC[i]*X[i]+BnC[i]<0.0 || AnC[i]*Y[i]+BnC[i]<0.0;
               vector(7) __An[i] = EMPTY[i] ? 0.0 : CORR[i] ? 0.0 : AnC[i];
               vector(7) __Bn[i] = EMPTY[i] ? 0.0 : CORR[i] ? Nk[i][Ptr]/(Y[i]-X[i]) : BnC[i];
               
			   /* MOVED !!! */
			   /*!!!!*/vector(7) DC[i] = 0.5*(Dk[i] + Dk[i + 1]);
			   /*****/vector(7) Rek[i] = fabs(UwLight[i]) > 1E-8 ? 0.0 : DC[i] * dUk[i] / D;
			   /*****/vector(7) Prk[i] = D*ROK*_Ck / Lambda0*(1.0 + Lambda1 / DC[i]);
			   /*****/vector(7) Sck[i] = D / Dvap0*(1.0 + Dvap1 / DC[i]);
			   /*****/vector(7) sRek[i] = fabs(UwLight[i]) > 1E-8 ? 0.0 : sqrt(Rek[i]);
			   /*****/vector(7) Nuss1[i] = fabs(UwLight[i]) > 1E-8 ? 2.0 : 2.0 + 0.6*sRek[i] * pow(Sck[i], 1.0 / 3.0);
			   /*****/vector(7) Nuss2[i] = fabs(UwLight[i]) > 1E-8 ? 2.0 : 2.0 + 0.6*sRek[i] * pow(Prk[i], 1.0 / 3.0);
			   vector(7) PkF3[i] = ROK / 2.0 / Nuss1[i] / PkF / Mk;
			   vector(7) PkF4[i] = _Hvap*ROK / 2.0 / Nuss2[i] / TK*(_Hvap*Mk / _R / TK - 1.0);
			   /* MOVED !!!! */
			   vector(7) LX[i] = EMPTY[i] ? 0.0 : 0.5*ROK*dCX[i] / (PkF3[i] / Dvap0*(1.0 + Dvap1 / X[i]) + PkF4[i] / Lambda0*(1.0 + Lambda1 / X[i]));
			   vector(7) LY[i] = EMPTY[i] ? 0.0 : 0.5*ROK*dCY[i] / (PkF3[i] / Dvap0*(1.0 + Dvap1 / Y[i]) + PkF4[i] / Lambda0*(1.0 + Lambda1 / Y[i]));
			   /*!!!*/vector(7) Dz[i] = EMPTY[i] || fabs(Y[i] - Dk[i + 1])>1E-15 ? Y[i] : LX[i]<0 ? Y[i] : max(X[i], min(Y[i], Y[i] - TAU*2.0*LY[i] / ROK / Y[i]));
			   /*!!!*/vector(7) Dzz[i] = EMPTY[i] || fabs(X[i] - Dk[i])>1E-15 ? X[i] : LX[i]<0 ? max(X[i], min(Y[i], X[i] - TAU*2.0*LX[i] / ROK / X[i])) : X[i];
			   /* MOVED !!!! */
			   /*!!!*/ vector(7) XYC[i] = EMPTY[i] ? _XYC[i] : _XYC[i] < Dzz[i] ? Dzz[i] : _XYC[i] > Dz[i] ? Dz[i] : _XYC[i];
			   /*!!!*/ vector(7) XYC2[i] = XYC[i] * XYC[i];
			   /*!!!*/ vector(7) XYC3[i] = XYC2[i] * XYC[i];
			   /*!!!*/ vector(7) XYC4[i] = XYC3[i] * XYC[i];
			   /*!!!*/ vector(7) XYC5[i] = XYC4[i] * XYC[i];
			   /*!!!*/ vector(7) XYC6[i] = XYC5[i] * XYC[i];
			   /*!!!*/ vector(7) XYC7[i] = XYC6[i] * XYC[i];

			   /*!!!*/vector(7) XC1[i] = (XYC[i] - X[i]) / 1.0;
			   /*!!!*/vector(7) XC2[i] = (XYC2[i] - X2[i]) / 2.0;
			   /*!!!*/vector(7) XC3[i] = (XYC3[i] - X3[i]) / 3.0;
			   /*!!!*/vector(7) XC4[i] = (XYC4[i] - X4[i]) / 4.0;
			   /*!!!*/vector(7) XC5[i] = (XYC5[i] - X5[i]) / 5.0;
			   /*!!!*/vector(7) XC6[i] = (XYC6[i] - X6[i]) / 6.0;
			   /*!!!*/vector(7) XC7[i] = (XYC7[i] - X7[i]) / 7.0;
			   /*!!!*/vector(7) XC8[i] = (XYC7[i]*XYC[i] - X8[i]) / 8.0;
			   /*!!!*/vector(7) XC9[i] = (XYC7[i]*XYC2[i] - X9[i]) / 9.0;
			   /*!!!*/vector(7) XC10[i] = (XYC7[i]*XYC3[i] - X9[i]*X[i]) / 10.0;
			   /*!!!*/vector(7) XC11[i] = (XYC7[i]*XYC4[i] - X9[i]*X2[i]) / 11.0;

			   /*!!!*/vector(7) CY1[i] = -(XYC[i] - Y[i]) / 1.0;
			   /*!!!*/vector(7) CY2[i] = -(XYC2[i] - Y2[i]) / 2.0;
			   /*!!!*/vector(7) CY3[i] = -(XYC3[i] - Y3[i]) / 3.0;
			   /*!!!*/vector(7) CY4[i] = -(XYC4[i] - Y4[i]) / 4.0;
			   /*!!!*/vector(7) CY5[i] = -(XYC5[i] - Y5[i]) / 5.0;
			   /*!!!*/vector(7) CY6[i] = -(XYC6[i] - Y6[i]) / 6.0;
			   /*!!!*/vector(7) CY7[i] = -(XYC7[i] - Y7[i]) / 7.0;
			   /*!!!*/vector(7) CY8[i] = -(XYC7[i] * XYC[i] - Y8[i]) / 8.0;
			   /*!!!*/vector(7) CY9[i] = -(XYC7[i] * XYC2[i] - Y9[i]) / 9.0;
			   /*!!!*/vector(7) CY10[i] = -(XYC7[i] * XYC3[i] - Y9[i] * Y[i]) / 10.0;
			   /*!!!*/vector(7) CY11[i] = -(XYC7[i] * XYC4[i] - Y9[i] * Y2[i]) / 11.0;

			   /*!!!*/vector(7) _Dnn[i] = EMPTY[i] ? 0.0 : (i == 0 || EMPTY[i - 1]) ? __Bn[i] : 0.0;
			   /*!!!*/vector(7) _Cnn[i] = EMPTY[i] ? 0.0 : (i == 0 || EMPTY[i-1]) ? __An[i] : 0.0;
			   /*!!!*/vector(7) _Bnn[i] = 0.0;
			   /*!!!*/vector(7) _Ann[i] = EMPTY[i] ? 0.0 : (i == 0 || EMPTY[i-1]) ? 0.0 : getLeftSpline(i, Dk, XYC, XYC2, XYC3, XYC4, XYC5, XYC6, XYC7, ROK, __An, __Bn, _Ann, _Bnn, _Cnn, _Dnn);
			   /*!!!*/vector(7) CORR_ABCDnn[i] = _Ann[i] * X3[i] + _Bnn[i] * X2[i] + _Cnn[i] * X[i] + _Dnn[i] >= 0.0;
			   /*!!!*/vector(7) Dnn[i] = CORR_ABCDnn[i] ? _Dnn[i] : __Bn[i];
			   /*!!!*/vector(7) Cnn[i] = CORR_ABCDnn[i] ? _Cnn[i] : __An[i];
			   /*!!!*/vector(7) Bnn[i] = CORR_ABCDnn[i] ? _Bnn[i] : 0.0;
			   /*!!!*/vector(7) Ann[i] = CORR_ABCDnn[i] ? _Ann[i] : 0.0;

			   /*!!!*/vector(7) _Dn[i] = EMPTY[i] ? 0.0 : (i == Nc-1 || EMPTY[i + 1]) ? __Bn[i] : 0.0;
			   /*!!!*/vector(7) _Cn[i] = EMPTY[i] ? 0.0 : (i == Nc-1 || EMPTY[i + 1]) ? __An[i] : 0.0;
			   /*!!!*/vector(7) _Bn[i] = 0.0;
			   /*!!!*/vector(7) _An[i] = EMPTY[i] ? 0.0 : (i == Nc-1 || EMPTY[i + 1]) ? 0.0 : getRightSpline(i, Dk, XYC, XYC2, XYC3, XYC4, XYC5, XYC6, XYC7, ROK, __An, __Bn, _An, _Bn, _Cn, _Dn);
			   /*!!!*/vector(7) CORR_ABCDn[i] = _An[i] * X3[i] + _Bn[i] * X2[i] + _Cn[i] * X[i] + _Dn[i] >= 0.0;
			   /*!!!*/vector(7) Dn[i] = CORR_ABCDn[i] ? _Dn[i] : __Bn[i];
			   /*!!!*/vector(7) Cn[i] = CORR_ABCDn[i] ? _Cn[i] : __An[i];
			   /*!!!*/vector(7) Bn[i] = CORR_ABCDn[i] ? _Bn[i] : 0.0;
			   /*!!!*/vector(7) An[i] = CORR_ABCDn[i] ? _An[i] : 0.0;

			   /*!!!*/vector(7) kCSubst[i] = EMPTY[i] ? 0.0 : max(0.0, max(ROK/Rok[i][Ptr], 6.0 / (Pi*(
				   An[i] / 7.0*(Y7[i] - XYC7[i]) + Bn[i] / 6.0*(Y6[i] - XYC6[i]) + Cn[i] / 5.0*(Y5[i] - XYC5[i]) + Dn[i] / 4.0*(Y4[i] - XYC4[i]) +
				   Ann[i] / 7.0*(XYC7[i] - X7[i]) + Bnn[i] / 6.0*(XYC6[i] - X6[i]) + Cnn[i] / 5.0*(XYC5[i] - X5[i]) + Dnn[i] / 4.0*(XYC4[i] - X4[i])
				   ))));
			   /*!!!*/vector(7) _DD[i] = EMPTY[i] ? (Dk[i] + Dk[i + 1]) / 2.0 : (An[i] / 5.0*(Y5[i] - XYC5[i]) + Bn[i] / 4.0*(Y4[i] - XYC4[i]) + Cn[i] / 3.0*(Y3[i] - XYC3[i]) + Dn[i] / 2.0*(Y2[i] - XYC2[i]) +
				   Ann[i] / 5.0*(XYC5[i] - X5[i]) + Bnn[i] / 4.0*(XYC4[i] - X4[i]) + Cnn[i] / 3.0*(XYC3[i] - X3[i]) + Dnn[i] / 2.0*(XYC2[i] - X2[i])) / Nk[i][Ptr];
			   /*!!!*/vector(7) DD[i] = _DD[i] < Dk[i] || _DD[i] > Dk[i + 1] ? (Dk[i] + Dk[i + 1]) / 2.0 : _DD[i];
			   /*!!!*/vector(7) nX[i] = EMPTY[i] ? 0.0 : Ann[i] * X3[i] + Bnn[i] * X2[i] + Cnn[i] * X[i] + Dnn[i];
			   /*!!!*/vector(7) nY[i] = EMPTY[i] ? 0.0 : An[i] * Y3[i] + Bn[i] * Y2[i] + Cn[i] * Y[i] + Dn[i];
			   /*!!! MOVED */
			   vector(7) LDp[i] = EMPTY[i] || Mode[i] == 0.0 ? 0.0 : 0.5*ROK*(1E6*Source_Drops0[Ptr] / PkF - exp(PkF1 / DD[i] - PkF2[i] / (DD[i] * DD[i] * DD[i]))) / (PkF3[i] / Dvap0*(1.0 + Dvap1 / DD[i]) + PkF4[i] / Lambda0*(1.0 + Lambda1 / DD[i]));

			   /* !!! MOVED, CORRECTED */
			   vector(7) iBeta[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 :
				   min(1E-2/max(HX[x],max(HY[y],HZ[z])), max(1E-5, Pi / 4.0*
				   (
				   Ann[i] * (ieC0[i] * XC6[i] + ieC1[i] * XC7[i] + ieC2[i] * XC8[i] + ieC3[i] * XC9[i] + ieC4[i] * XC10[i] + ieC5[i] * XC11[i]) +
				   Bnn[i] * (ieC0[i] * XC5[i] + ieC1[i] * XC6[i] + ieC2[i] * XC7[i] + ieC3[i] * XC8[i] + ieC4[i] * XC9[i] + ieC5[i] * XC10[i]) +
				   Cnn[i] * (ieC0[i] * XC4[i] + ieC1[i] * XC5[i] + ieC2[i] * XC6[i] + ieC3[i] * XC7[i] + ieC4[i] * XC8[i] + ieC5[i] * XC9[i]) +
				   Dnn[i] * (ieC0[i] * XC3[i] + ieC1[i] * XC4[i] + ieC2[i] * XC5[i] + ieC3[i] * XC6[i] + ieC4[i] * XC7[i] + ieC5[i] * XC8[i]) +
				   An[i] * (ieC0[i] * CY6[i] + ieC1[i] * CY7[i] + ieC2[i] * CY8[i] + ieC3[i] * CY9[i] + ieC4[i] * CY10[i] + ieC5[i] * CY11[i]) +
				   Bn[i] * (ieC0[i] * CY5[i] + ieC1[i] * CY6[i] + ieC2[i] * CY7[i] + ieC3[i] * CY8[i] + ieC4[i] * CY9[i] + ieC5[i] * CY10[i]) +
				   Cn[i] * (ieC0[i] * CY4[i] + ieC1[i] * CY5[i] + ieC2[i] * CY6[i] + ieC3[i] * CY7[i] + ieC4[i] * CY8[i] + ieC5[i] * CY9[i]) +
				   Dn[i] * (ieC0[i] * CY3[i] + ieC1[i] * CY4[i] + ieC2[i] * CY5[i] + ieC3[i] * CY6[i] + ieC4[i] * CY7[i] + ieC5[i] * CY8[i])
				   )))
				   ) : 0.0));
			   vector(7) uiBeta[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 :
				   min(1E-2/max(HX[x],max(HY[y], HZ[z])), max(1E-5, Pi / 4.0*
				   (
				   Ann[i] * (uieC0[i] * XC6[i] + uieC1[i] * XC7[i] + uieC2[i] * XC8[i] + uieC3[i] * XC9[i] + uieC4[i] * XC10[i] + uieC5[i] * XC11[i]) +
				   Bnn[i] * (uieC0[i] * XC5[i] + uieC1[i] * XC6[i] + uieC2[i] * XC7[i] + uieC3[i] * XC8[i] + uieC4[i] * XC9[i] + uieC5[i] * XC10[i]) +
				   Cnn[i] * (uieC0[i] * XC4[i] + uieC1[i] * XC5[i] + uieC2[i] * XC6[i] + uieC3[i] * XC7[i] + uieC4[i] * XC8[i] + uieC5[i] * XC9[i]) +
				   Dnn[i] * (uieC0[i] * XC3[i] + uieC1[i] * XC4[i] + uieC2[i] * XC5[i] + uieC3[i] * XC6[i] + uieC4[i] * XC7[i] + uieC5[i] * XC8[i]) +
				   An[i] * (uieC0[i] * CY6[i] + uieC1[i] * CY7[i] + uieC2[i] * CY8[i] + uieC3[i] * CY9[i] + uieC4[i] * CY10[i] + uieC5[i] * CY11[i]) +
				   Bn[i] * (uieC0[i] * CY5[i] + uieC1[i] * CY6[i] + uieC2[i] * CY7[i] + uieC3[i] * CY8[i] + uieC4[i] * CY9[i] + uieC5[i] * CY10[i]) +
				   Cn[i] * (uieC0[i] * CY4[i] + uieC1[i] * CY5[i] + uieC2[i] * CY6[i] + uieC3[i] * CY7[i] + uieC4[i] * CY8[i] + uieC5[i] * CY9[i]) +
				   Dn[i] * (uieC0[i] * CY3[i] + uieC1[i] * CY4[i] + uieC2[i] * CY5[i] + uieC3[i] * CY6[i] + uieC4[i] * CY7[i] + uieC5[i] * CY8[i])
				   )))
				   ) : 0.0));
			   vector(7) riBeta[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 :
				   min(1E-2/max(HX[x],max(HY[y], HZ[z])), max(1E-5, Pi / 4.0*
				   (
				   Ann[i] * (rieC0[i] * XC6[i] + rieC1[i] * XC7[i] + rieC2[i] * XC8[i] + rieC3[i] * XC9[i] + rieC4[i] * XC10[i] + rieC5[i] * XC11[i]) +
				   Bnn[i] * (rieC0[i] * XC5[i] + rieC1[i] * XC6[i] + rieC2[i] * XC7[i] + rieC3[i] * XC8[i] + rieC4[i] * XC9[i] + rieC5[i] * XC10[i]) +
				   Cnn[i] * (rieC0[i] * XC4[i] + rieC1[i] * XC5[i] + rieC2[i] * XC6[i] + rieC3[i] * XC7[i] + rieC4[i] * XC8[i] + rieC5[i] * XC9[i]) +
				   Dnn[i] * (rieC0[i] * XC3[i] + rieC1[i] * XC4[i] + rieC2[i] * XC5[i] + rieC3[i] * XC6[i] + rieC4[i] * XC7[i] + rieC5[i] * XC8[i]) +
				   An[i] * (rieC0[i] * CY6[i] + rieC1[i] * CY7[i] + rieC2[i] * CY8[i] + rieC3[i] * CY9[i] + rieC4[i] * CY10[i] + rieC5[i] * CY11[i]) +
				   Bn[i] * (rieC0[i] * CY5[i] + rieC1[i] * CY6[i] + rieC2[i] * CY7[i] + rieC3[i] * CY8[i] + rieC4[i] * CY9[i] + rieC5[i] * CY10[i]) +
				   Cn[i] * (rieC0[i] * CY4[i] + rieC1[i] * CY5[i] + rieC2[i] * CY6[i] + rieC3[i] * CY7[i] + rieC4[i] * CY8[i] + rieC5[i] * CY9[i]) +
				   Dn[i] * (rieC0[i] * CY3[i] + rieC1[i] * CY4[i] + rieC2[i] * CY5[i] + rieC3[i] * CY6[i] + rieC4[i] * CY7[i] + rieC5[i] * CY8[i])
				   )))
				   ) : 0.0));
			   vector(7) tiBeta[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 :
				   min(1E-2/max(HX[x],max(HY[y], HZ[z])), max(1E-5, Pi / 4.0*
				   (
				   Ann[i] * (tieC0[i] * XC6[i] + tieC1[i] * XC7[i] + tieC2[i] * XC8[i] + tieC3[i] * XC9[i] + tieC4[i] * XC10[i] + tieC5[i] * XC11[i]) +
				   Bnn[i] * (tieC0[i] * XC5[i] + tieC1[i] * XC6[i] + tieC2[i] * XC7[i] + tieC3[i] * XC8[i] + tieC4[i] * XC9[i] + tieC5[i] * XC10[i]) +
				   Cnn[i] * (tieC0[i] * XC4[i] + tieC1[i] * XC5[i] + tieC2[i] * XC6[i] + tieC3[i] * XC7[i] + tieC4[i] * XC8[i] + tieC5[i] * XC9[i]) +
				   Dnn[i] * (tieC0[i] * XC3[i] + tieC1[i] * XC4[i] + tieC2[i] * XC5[i] + tieC3[i] * XC6[i] + tieC4[i] * XC7[i] + tieC5[i] * XC8[i]) +
				   An[i] * (tieC0[i] * CY6[i] + tieC1[i] * CY7[i] + tieC2[i] * CY8[i] + tieC3[i] * CY9[i] + tieC4[i] * CY10[i] + tieC5[i] * CY11[i]) +
				   Bn[i] * (tieC0[i] * CY5[i] + tieC1[i] * CY6[i] + tieC2[i] * CY7[i] + tieC3[i] * CY8[i] + tieC4[i] * CY9[i] + tieC5[i] * CY10[i]) +
				   Cn[i] * (tieC0[i] * CY4[i] + tieC1[i] * CY5[i] + tieC2[i] * CY6[i] + tieC3[i] * CY7[i] + tieC4[i] * CY8[i] + tieC5[i] * CY9[i]) +
				   Dn[i] * (tieC0[i] * CY3[i] + tieC1[i] * CY4[i] + tieC2[i] * CY5[i] + tieC3[i] * CY6[i] + tieC4[i] * CY7[i] + tieC5[i] * CY8[i])
				   )))
				   ) : 0.0));
			   vector(7) iAbs[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 :
				   min(iBeta[i], max(0.0, Pi / 4.0*
				   (
				   Ann[i] * (iaC0[i] * XC6[i] + iaC1[i] * XC7[i] + iaC2[i] * XC8[i] + iaC3[i] * XC9[i] + iaC4[i] * XC10[i] + iaC5[i] * XC11[i]) +
				   Bnn[i] * (iaC0[i] * XC5[i] + iaC1[i] * XC6[i] + iaC2[i] * XC7[i] + iaC3[i] * XC8[i] + iaC4[i] * XC9[i] + iaC5[i] * XC10[i]) +
				   Cnn[i] * (iaC0[i] * XC4[i] + iaC1[i] * XC5[i] + iaC2[i] * XC6[i] + iaC3[i] * XC7[i] + iaC4[i] * XC8[i] + iaC5[i] * XC9[i]) +
				   Dnn[i] * (iaC0[i] * XC3[i] + iaC1[i] * XC4[i] + iaC2[i] * XC5[i] + iaC3[i] * XC6[i] + iaC4[i] * XC7[i] + iaC5[i] * XC8[i]) +
				   An[i] * (iaC0[i] * CY6[i] + iaC1[i] * CY7[i] + iaC2[i] * CY8[i] + iaC3[i] * CY9[i] + iaC4[i] * CY10[i] + iaC5[i] * CY11[i]) +
				   Bn[i] * (iaC0[i] * CY5[i] + iaC1[i] * CY6[i] + iaC2[i] * CY7[i] + iaC3[i] * CY8[i] + iaC4[i] * CY9[i] + iaC5[i] * CY10[i]) +
				   Cn[i] * (iaC0[i] * CY4[i] + iaC1[i] * CY5[i] + iaC2[i] * CY6[i] + iaC3[i] * CY7[i] + iaC4[i] * CY8[i] + iaC5[i] * CY9[i]) +
				   Dn[i] * (iaC0[i] * CY3[i] + iaC1[i] * CY4[i] + iaC2[i] * CY5[i] + iaC3[i] * CY6[i] + iaC4[i] * CY7[i] + iaC5[i] * CY8[i])
				   )))
				   ) : 0.0));
			   vector(7) uiAbs[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 :
				   min(uiBeta[i], max(0.0, Pi / 4.0*
				   (
				   Ann[i] * (uiaC0[i] * XC6[i] + uiaC1[i] * XC7[i] + uiaC2[i] * XC8[i] + uiaC3[i] * XC9[i] + uiaC4[i] * XC10[i] + uiaC5[i] * XC11[i]) +
				   Bnn[i] * (uiaC0[i] * XC5[i] + uiaC1[i] * XC6[i] + uiaC2[i] * XC7[i] + uiaC3[i] * XC8[i] + uiaC4[i] * XC9[i] + uiaC5[i] * XC10[i]) +
				   Cnn[i] * (uiaC0[i] * XC4[i] + uiaC1[i] * XC5[i] + uiaC2[i] * XC6[i] + uiaC3[i] * XC7[i] + uiaC4[i] * XC8[i] + uiaC5[i] * XC9[i]) +
				   Dnn[i] * (uiaC0[i] * XC3[i] + uiaC1[i] * XC4[i] + uiaC2[i] * XC5[i] + uiaC3[i] * XC6[i] + uiaC4[i] * XC7[i] + uiaC5[i] * XC8[i]) +
				   An[i] * (uiaC0[i] * CY6[i] + uiaC1[i] * CY7[i] + uiaC2[i] * CY8[i] + uiaC3[i] * CY9[i] + uiaC4[i] * CY10[i] + uiaC5[i] * CY11[i]) +
				   Bn[i] * (uiaC0[i] * CY5[i] + uiaC1[i] * CY6[i] + uiaC2[i] * CY7[i] + uiaC3[i] * CY8[i] + uiaC4[i] * CY9[i] + uiaC5[i] * CY10[i]) +
				   Cn[i] * (uiaC0[i] * CY4[i] + uiaC1[i] * CY5[i] + uiaC2[i] * CY6[i] + uiaC3[i] * CY7[i] + uiaC4[i] * CY8[i] + uiaC5[i] * CY9[i]) +
				   Dn[i] * (uiaC0[i] * CY3[i] + uiaC1[i] * CY4[i] + uiaC2[i] * CY5[i] + uiaC3[i] * CY6[i] + uiaC4[i] * CY7[i] + uiaC5[i] * CY8[i])
				   )))
				   ) : 0.0));
			   vector(7) riAbs[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 :
				   min(riBeta[i], max(0.0, Pi / 4.0*
				   (
				   Ann[i] * (riaC0[i] * XC6[i] + riaC1[i] * XC7[i] + riaC2[i] * XC8[i] + riaC3[i] * XC9[i] + riaC4[i] * XC10[i] + riaC5[i] * XC11[i]) +
				   Bnn[i] * (riaC0[i] * XC5[i] + riaC1[i] * XC6[i] + riaC2[i] * XC7[i] + riaC3[i] * XC8[i] + riaC4[i] * XC9[i] + riaC5[i] * XC10[i]) +
				   Cnn[i] * (riaC0[i] * XC4[i] + riaC1[i] * XC5[i] + riaC2[i] * XC6[i] + riaC3[i] * XC7[i] + riaC4[i] * XC8[i] + riaC5[i] * XC9[i]) +
				   Dnn[i] * (riaC0[i] * XC3[i] + riaC1[i] * XC4[i] + riaC2[i] * XC5[i] + riaC3[i] * XC6[i] + riaC4[i] * XC7[i] + riaC5[i] * XC8[i]) +
				   An[i] * (riaC0[i] * CY6[i] + riaC1[i] * CY7[i] + riaC2[i] * CY8[i] + riaC3[i] * CY9[i] + riaC4[i] * CY10[i] + riaC5[i] * CY11[i]) +
				   Bn[i] * (riaC0[i] * CY5[i] + riaC1[i] * CY6[i] + riaC2[i] * CY7[i] + riaC3[i] * CY8[i] + riaC4[i] * CY9[i] + riaC5[i] * CY10[i]) +
				   Cn[i] * (riaC0[i] * CY4[i] + riaC1[i] * CY5[i] + riaC2[i] * CY6[i] + riaC3[i] * CY7[i] + riaC4[i] * CY8[i] + riaC5[i] * CY9[i]) +
				   Dn[i] * (riaC0[i] * CY3[i] + riaC1[i] * CY4[i] + riaC2[i] * CY5[i] + riaC3[i] * CY6[i] + riaC4[i] * CY7[i] + riaC5[i] * CY8[i])
				   )))
				   ) : 0.0));
			   vector(7) tiAbs[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 :
				   min(tiBeta[i], max(0.0, Pi / 4.0*
				   (
				   Ann[i] * (tiaC0[i] * XC6[i] + tiaC1[i] * XC7[i] + tiaC2[i] * XC8[i] + tiaC3[i] * XC9[i] + tiaC4[i] * XC10[i] + tiaC5[i] * XC11[i]) +
				   Bnn[i] * (tiaC0[i] * XC5[i] + tiaC1[i] * XC6[i] + tiaC2[i] * XC7[i] + tiaC3[i] * XC8[i] + tiaC4[i] * XC9[i] + tiaC5[i] * XC10[i]) +
				   Cnn[i] * (tiaC0[i] * XC4[i] + tiaC1[i] * XC5[i] + tiaC2[i] * XC6[i] + tiaC3[i] * XC7[i] + tiaC4[i] * XC8[i] + tiaC5[i] * XC9[i]) +
				   Dnn[i] * (tiaC0[i] * XC3[i] + tiaC1[i] * XC4[i] + tiaC2[i] * XC5[i] + tiaC3[i] * XC6[i] + tiaC4[i] * XC7[i] + tiaC5[i] * XC8[i]) +
				   An[i] * (tiaC0[i] * CY6[i] + tiaC1[i] * CY7[i] + tiaC2[i] * CY8[i] + tiaC3[i] * CY9[i] + tiaC4[i] * CY10[i] + tiaC5[i] * CY11[i]) +
				   Bn[i] * (tiaC0[i] * CY5[i] + tiaC1[i] * CY6[i] + tiaC2[i] * CY7[i] + tiaC3[i] * CY8[i] + tiaC4[i] * CY9[i] + tiaC5[i] * CY10[i]) +
				   Cn[i] * (tiaC0[i] * CY4[i] + tiaC1[i] * CY5[i] + tiaC2[i] * CY6[i] + tiaC3[i] * CY7[i] + tiaC4[i] * CY8[i] + tiaC5[i] * CY9[i]) +
				   Dn[i] * (tiaC0[i] * CY3[i] + tiaC1[i] * CY4[i] + tiaC2[i] * CY5[i] + tiaC3[i] * CY6[i] + tiaC4[i] * CY7[i] + tiaC5[i] * CY8[i])
				   )))
				   ) : 0.0));
			   vector(7) iAsym[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 :
				   min(1.0, max(0.0, (
				   Ann[i] * (igC0[i] * XC4[i] + igC1[i] * XC5[i] + igC2[i] * XC6[i] + igC3[i] * XC7[i] + igC4[i] * XC8[i] + igC5[i] * XC9[i]) +
				   Bnn[i] * (igC0[i] * XC3[i] + igC1[i] * XC4[i] + igC2[i] * XC5[i] + igC3[i] * XC6[i] + igC4[i] * XC7[i] + igC5[i] * XC8[i]) +
				   Cnn[i] * (igC0[i] * XC2[i] + igC1[i] * XC3[i] + igC2[i] * XC4[i] + igC3[i] * XC5[i] + igC4[i] * XC6[i] + igC5[i] * XC7[i]) +
				   Dnn[i] * (igC0[i] * XC1[i] + igC1[i] * XC2[i] + igC2[i] * XC3[i] + igC3[i] * XC4[i] + igC4[i] * XC5[i] + igC5[i] * XC6[i]) +
				   An[i] * (igC0[i] * CY4[i] + igC1[i] * CY5[i] + igC2[i] * CY6[i] + igC3[i] * CY7[i] + igC4[i] * CY8[i] + igC5[i] * CY9[i]) +
				   Bn[i] * (igC0[i] * CY3[i] + igC1[i] * CY4[i] + igC2[i] * CY5[i] + igC3[i] * CY6[i] + igC4[i] * CY7[i] + igC5[i] * CY8[i]) +
				   Cn[i] * (igC0[i] * CY2[i] + igC1[i] * CY3[i] + igC2[i] * CY4[i] + igC3[i] * CY5[i] + igC4[i] * CY6[i] + igC5[i] * CY7[i]) +
				   Dn[i] * (igC0[i] * CY1[i] + igC1[i] * CY2[i] + igC2[i] * CY3[i] + igC3[i] * CY4[i] + igC4[i] * CY5[i] + igC5[i] * CY6[i])
				   ) / Nk[i][Ptr]))
				   ) : 0.0));
			   vector(7) uiAsym[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 :
				   min(1.0, max(0.0, (
				   Ann[i] * (uigC0[i] * XC4[i] + uigC1[i] * XC5[i] + uigC2[i] * XC6[i] + uigC3[i] * XC7[i] + uigC4[i] * XC8[i] + uigC5[i] * XC9[i]) +
				   Bnn[i] * (uigC0[i] * XC3[i] + uigC1[i] * XC4[i] + uigC2[i] * XC5[i] + uigC3[i] * XC6[i] + uigC4[i] * XC7[i] + uigC5[i] * XC8[i]) +
				   Cnn[i] * (uigC0[i] * XC2[i] + uigC1[i] * XC3[i] + uigC2[i] * XC4[i] + uigC3[i] * XC5[i] + uigC4[i] * XC6[i] + uigC5[i] * XC7[i]) +
				   Dnn[i] * (uigC0[i] * XC1[i] + uigC1[i] * XC2[i] + uigC2[i] * XC3[i] + uigC3[i] * XC4[i] + uigC4[i] * XC5[i] + uigC5[i] * XC6[i]) +
				   An[i] * (uigC0[i] * CY4[i] + uigC1[i] * CY5[i] + uigC2[i] * CY6[i] + uigC3[i] * CY7[i] + uigC4[i] * CY8[i] + uigC5[i] * CY9[i]) +
				   Bn[i] * (uigC0[i] * CY3[i] + uigC1[i] * CY4[i] + uigC2[i] * CY5[i] + uigC3[i] * CY6[i] + uigC4[i] * CY7[i] + uigC5[i] * CY8[i]) +
				   Cn[i] * (uigC0[i] * CY2[i] + uigC1[i] * CY3[i] + uigC2[i] * CY4[i] + uigC3[i] * CY5[i] + uigC4[i] * CY6[i] + uigC5[i] * CY7[i]) +
				   Dn[i] * (uigC0[i] * CY1[i] + uigC1[i] * CY2[i] + uigC2[i] * CY3[i] + uigC3[i] * CY4[i] + uigC4[i] * CY5[i] + uigC5[i] * CY6[i])
				   ) / Nk[i][Ptr]))
				   ) : 0.0));
			   vector(7) riAsym[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 :
				   min(1.0, max(0.0, (
				   Ann[i] * (rigC0[i] * XC4[i] + rigC1[i] * XC5[i] + rigC2[i] * XC6[i] + rigC3[i] * XC7[i] + rigC4[i] * XC8[i] + rigC5[i] * XC9[i]) +
				   Bnn[i] * (rigC0[i] * XC3[i] + rigC1[i] * XC4[i] + rigC2[i] * XC5[i] + rigC3[i] * XC6[i] + rigC4[i] * XC7[i] + rigC5[i] * XC8[i]) +
				   Cnn[i] * (rigC0[i] * XC2[i] + rigC1[i] * XC3[i] + rigC2[i] * XC4[i] + rigC3[i] * XC5[i] + rigC4[i] * XC6[i] + rigC5[i] * XC7[i]) +
				   Dnn[i] * (rigC0[i] * XC1[i] + rigC1[i] * XC2[i] + rigC2[i] * XC3[i] + rigC3[i] * XC4[i] + rigC4[i] * XC5[i] + rigC5[i] * XC6[i]) +
				   An[i] * (rigC0[i] * CY4[i] + rigC1[i] * CY5[i] + rigC2[i] * CY6[i] + rigC3[i] * CY7[i] + rigC4[i] * CY8[i] + rigC5[i] * CY9[i]) +
				   Bn[i] * (rigC0[i] * CY3[i] + rigC1[i] * CY4[i] + rigC2[i] * CY5[i] + rigC3[i] * CY6[i] + rigC4[i] * CY7[i] + rigC5[i] * CY8[i]) +
				   Cn[i] * (rigC0[i] * CY2[i] + rigC1[i] * CY3[i] + rigC2[i] * CY4[i] + rigC3[i] * CY5[i] + rigC4[i] * CY6[i] + rigC5[i] * CY7[i]) +
				   Dn[i] * (rigC0[i] * CY1[i] + rigC1[i] * CY2[i] + rigC2[i] * CY3[i] + rigC3[i] * CY4[i] + rigC4[i] * CY5[i] + rigC5[i] * CY6[i])
				   ) / Nk[i][Ptr]))
				   ) : 0.0));
			   vector(7) tiAsym[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 :
				   min(1.0, max( 0.0, (
				   Ann[i] * (tigC0[i] * XC4[i] + tigC1[i] * XC5[i] + tigC2[i] * XC6[i] + tigC3[i] * XC7[i] + tigC4[i] * XC8[i] + tigC5[i] * XC9[i]) +
				   Bnn[i] * (tigC0[i] * XC3[i] + tigC1[i] * XC4[i] + tigC2[i] * XC5[i] + tigC3[i] * XC6[i] + tigC4[i] * XC7[i] + tigC5[i] * XC8[i]) +
				   Cnn[i] * (tigC0[i] * XC2[i] + tigC1[i] * XC3[i] + tigC2[i] * XC4[i] + tigC3[i] * XC5[i] + tigC4[i] * XC6[i] + tigC5[i] * XC7[i]) +
				   Dnn[i] * (tigC0[i] * XC1[i] + tigC1[i] * XC2[i] + tigC2[i] * XC3[i] + tigC3[i] * XC4[i] + tigC4[i] * XC5[i] + tigC5[i] * XC6[i]) +
				   An[i] * (tigC0[i] * CY4[i] + tigC1[i] * CY5[i] + tigC2[i] * CY6[i] + tigC3[i] * CY7[i] + tigC4[i] * CY8[i] + tigC5[i] * CY9[i]) +
				   Bn[i] * (tigC0[i] * CY3[i] + tigC1[i] * CY4[i] + tigC2[i] * CY5[i] + tigC3[i] * CY6[i] + tigC4[i] * CY7[i] + tigC5[i] * CY8[i]) +
				   Cn[i] * (tigC0[i] * CY2[i] + tigC1[i] * CY3[i] + tigC2[i] * CY4[i] + tigC3[i] * CY5[i] + tigC4[i] * CY6[i] + tigC5[i] * CY7[i]) +
				   Dn[i] * (tigC0[i] * CY1[i] + tigC1[i] * CY2[i] + tigC2[i] * CY3[i] + tigC3[i] * CY4[i] + tigC4[i] * CY5[i] + tigC5[i] * CY6[i])
				   ) / Nk[i][Ptr]))
				   ) : 0.0));
			   /* !!! END MOVED */

			   /*!!!*/vector(7) iAlbedo[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 : min(1.0, max(0.0, (iBeta[i] - iAbs[i]) / iBeta[i]))) : 0.0));
			   /*!!!*/vector(7) uiAlbedo[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 : min(1.0, max(0.0, (uiBeta[i] - uiAbs[i]) / uiBeta[i]))) : 0.0));
			   /*!!!*/vector(7) riAlbedo[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 : min(1.0, max(0.0, (riBeta[i] - riAbs[i]) / riBeta[i]))) : 0.0));
			   /*!!!*/vector(7) tiAlbedo[i] = ((FastMode ? (0.0) : 0.0) + (SlowMode ? (EMPTY[i] ? 0.0 : min(1.0, max(0.0, (tiBeta[i] - tiAbs[i]) / tiBeta[i]))) : 0.0));
               
			   /*!!!*/eBeta = min(1E-2/max(HX[x],max(HY[y], HZ[z])), AirBeta + ((iBeta[0]) + (iBeta[1]) + (iBeta[2]) + (iBeta[3]) + (iBeta[4]) + (iBeta[5]) + (iBeta[6])));
               
			   /*!!!*/ueBeta = min(1E-2/max(HX[x],max(HY[y], HZ[z])), uAirBeta + ((uiBeta[0]) + (uiBeta[1]) + (uiBeta[2]) + (uiBeta[3]) + (uiBeta[4]) + (uiBeta[5]) + (uiBeta[6])));
               
			   /*!!!*/reBeta = min(1E-2/max(HX[x],max(HY[y], HZ[z])), rAirBeta + ((riBeta[0]) + (riBeta[1]) + (riBeta[2]) + (riBeta[3]) + (riBeta[4]) + (riBeta[5]) + (riBeta[6])));
               
			   /*!!!*/teBeta = min(1E-2/max(HX[x],max(HY[y], HZ[z])), tAirBeta + ((tiBeta[0]) + (tiBeta[1]) + (tiBeta[2]) + (tiBeta[3]) + (tiBeta[4]) + (tiBeta[5]) + (tiBeta[6])));

			   /*!!!*/eScat = min(eBeta, AirAlbedo*AirBeta + ((iAlbedo[0] * iBeta[0]) + (iAlbedo[1] * iBeta[1]) + (iAlbedo[2] * iBeta[2]) + (iAlbedo[3] * iBeta[3]) + (iAlbedo[4] * iBeta[4]) + (iAlbedo[5] * iBeta[5]) + (iAlbedo[6] * iBeta[6])));
			   /*!!!*/ueScat = min(ueBeta, uAirAlbedo*uAirBeta + ((uiAlbedo[0] * uiBeta[0]) + (uiAlbedo[1] * uiBeta[1]) + (uiAlbedo[2] * uiBeta[2]) + (uiAlbedo[3] * uiBeta[3]) + (uiAlbedo[4] * uiBeta[4]) + (uiAlbedo[5] * uiBeta[5]) + (uiAlbedo[6] * uiBeta[6])));
			   /*!!!*/reScat = min(reBeta, rAirAlbedo*rAirBeta + ((riAlbedo[0] * riBeta[0]) + (riAlbedo[1] * riBeta[1]) + (riAlbedo[2] * riBeta[2]) + (riAlbedo[3] * riBeta[3]) + (riAlbedo[4] * riBeta[4]) + (riAlbedo[5] * riBeta[5]) + (riAlbedo[6] * riBeta[6])));
			   /*!!!*/teScat = min(teBeta, tAirAlbedo*tAirBeta + ((tiAlbedo[0] * tiBeta[0]) + (tiAlbedo[1] * tiBeta[1]) + (tiAlbedo[2] * tiBeta[2]) + (tiAlbedo[3] * tiBeta[3]) + (tiAlbedo[4] * tiBeta[4]) + (tiAlbedo[5] * tiBeta[5]) + (tiAlbedo[6] * tiBeta[6])));

			   /*!!!!!*/kAlbedo = min(1.0, eScat / eBeta);
			   /*!!!!!*/kAsymm = min(1.0, (AirAlbedo*AirBeta*AirAsymmetry + ((iAlbedo[0] * iBeta[0] * iAsym[0]) + (iAlbedo[1] * iBeta[1] * iAsym[1]) + (iAlbedo[2] * iBeta[2] * iAsym[2]) + (iAlbedo[3] * iBeta[3] * iAsym[3]) + (iAlbedo[4] * iBeta[4] * iAsym[4]) + (iAlbedo[5] * iBeta[5] * iAsym[5]) + (iAlbedo[6] * iBeta[6] * iAsym[6]))) / eScat);
               /*!!!!!*/ukAlbedo = min(1.0, ueScat/ueBeta);
			   /*!!!!!*/ukAsymm = min(1.0, (uAirAlbedo*uAirBeta*uAirAsymmetry + ((uiAlbedo[0] * uiBeta[0] * uiAsym[0]) + (uiAlbedo[1] * uiBeta[1] * uiAsym[1]) + (uiAlbedo[2] * uiBeta[2] * uiAsym[2]) + (uiAlbedo[3] * uiBeta[3] * uiAsym[3]) + (uiAlbedo[4] * uiBeta[4] * uiAsym[4]) + (uiAlbedo[5] * uiBeta[5] * uiAsym[5]) + (uiAlbedo[6] * uiBeta[6] * uiAsym[6]))) / ueScat);
			   /*!!!!!*/rkAlbedo = min(1.0, reScat / reBeta);
			   /*!!!!!*/rkAsymm = min(1.0, (rAirAlbedo*rAirBeta*rAirAsymmetry + ((riAlbedo[0] * riBeta[0] * riAsym[0]) + (riAlbedo[1] * riBeta[1] * riAsym[1]) + (riAlbedo[2] * riBeta[2] * riAsym[2]) + (riAlbedo[3] * riBeta[3] * riAsym[3]) + (riAlbedo[4] * riBeta[4] * riAsym[4]) + (riAlbedo[5] * riBeta[5] * riAsym[5]) + (riAlbedo[6] * riBeta[6] * riAsym[6]))) / reScat);
			   /*!!!!!*/tkAlbedo = min(1.0, teScat / teBeta);
			   /*!!!!!*/tkAsymm = min(1.0, (tAirAlbedo*tAirBeta*tAirAsymmetry + ((tiAlbedo[0] * tiBeta[0] * tiAsym[0]) + (tiAlbedo[1] * tiBeta[1] * tiAsym[1]) + (tiAlbedo[2] * tiBeta[2] * tiAsym[2]) + (tiAlbedo[3] * tiBeta[3] * tiAsym[3]) + (tiAlbedo[4] * tiBeta[4] * tiAsym[4]) + (tiAlbedo[5] * tiBeta[5] * tiAsym[5]) + (tiAlbedo[6] * tiBeta[6] * tiAsym[6]))) / teScat);
			   ukAlbAsm = ukAlbedo*ukAsymm;
               utAlpha = ueBeta*(1-ukAlbedo);
               u_Ft = 3.0*ueBeta*ukAlbedo/4.0/Pi;
               kAlbAsm = kAlbedo*kAsymm;
               tAlpha = eBeta*(1-kAlbedo);
               _Ft = 3.0*eBeta*kAlbedo/4.0/Pi;
               rkAlbAsm = rkAlbedo*rkAsymm;
               rtAlpha = reBeta*(1-rkAlbedo);
               r_Ft = 3.0*reBeta*rkAlbedo/4.0/Pi;
               
               tkAlbAsm = tkAlbedo*tkAsymm;
               ttAlpha = teBeta*(1-tkAlbedo);

			   utBeta = ueBeta*(1 - ukAlbAsm);
               uFt = u_Ft*FUV[Ptr];
               tBeta = eBeta*(1-kAlbAsm);
               Ft = _Ft*F0[Ptr];
               rtBeta = reBeta*(1-rkAlbAsm);
               rFt = r_Ft*FIR[Ptr];
               ttBeta = teBeta*(1-tkAlbAsm);
               tFt = 3.0*ttAlpha*tPlanck;
               
               uk1 = 2.0/3.0/utBeta;
               
               k1 = 2.0/3.0/tBeta;
               
               rk1 = 2.0/3.0/rtBeta;
               
               tk1 = 2.0/3.0/ttBeta;

			   vector(7) X1[i] = EMPTY[i] ? X[i] : LX[i]<0 ? X[i] : X[i] + TAU*2.0*LX[i] / ROK / X[i];
               vector(7) Y1[i] = EMPTY[i] ? Y[i] : LX[i]<0 ? Y[i]+TAU*2.0*LY[i]/ROK/Y[i] : Y[i];

			   /*?????*/vector(7) Nuss1SO2[i] = fabs(UwLight[i]) > 1E-8 ? 2.0 : 2.0 + 0.6*sRek[i] * tScSO2;
			   /*?????*/vector(7) LSO2f[i] = EMPTY[i] ? 0.0 : -Nuss1SO2[i] * _DfSO2*1E6*kCSubst[i] / _HenrySO2[i];
			   /*?????*/vector(7) LSO2g[i] = EMPTY[i] ? 0.0 : Nuss1SO2[i] * _DfSO2*1E6*max(0.0, SO2[Ptr]);
			   
			   /*!!!!!*/vector(7) maxSO2k[i] = EMPTY[i] ? 0.0 : max(0.0, SO2[Ptr])*_HenrySO2[i] / kCSubst[i];

			   uk1a = uk1*(1 + SurfaceAlbedo);
               k1a = k1*(1+SurfaceAlbedo);
               rk1a = rk1*(1+SurfaceAlbedo);
               tk1a = tk1*(2.0-SurfaceEmissivity);
			   /*!!! MOVED*/vector(7) Dz2[i] = Dz[i] * Dz[i];
			   /*!!! MOVED*/vector(7) Dzz2[i] = Dzz[i] * Dzz[i];
			   /*!!! MOVED*/vector(7) Dz3[i] = Dz2[i] * Dz[i];
			   /*!!! MOVED*/vector(7) Dzz3[i] = Dzz2[i] * Dzz[i];
			   vector(7) LDz[i] = EMPTY[i] || Mode[i] == 0.0 ? 0.0 : LX[i]<0 ? LY[i] : 0.5*ROK*(1E6*Source_Drops0[Ptr] / PkF - exp(PkF1 / Dz[i] - PkF2[i] / (Dz[i] * Dz[i] * Dz[i]))) / (PkF3[i] / Dvap0*(1.0 + Dvap1 / Dz[i]) + PkF4[i] / Lambda0*(1.0 + Lambda1 / Dz[i]));
			   /*!!!*/vector(7) nDz[i] = EMPTY[i] ? 0.0 : LX[i]<0 ? nY[i] : An[i] * Dz3[i] + Bn[i] * Dz2[i] + Cn[i] * Dz[i] + Dn[i];
               vector(7) LDzz[i] = EMPTY[i] || Mode[i]==0.0 ? 0.0 : LX[i]<0 ? 0.5*ROK*(1E6*Source_Drops0[Ptr]/PkF-exp(PkF1/Dzz[i]-PkF2[i]/(Dzz[i]*Dzz[i]*Dzz[i])))/(PkF3[i]/Dvap0*(1.0+Dvap1/Dzz[i])+PkF4[i]/Lambda0*(1.0+Lambda1/Dzz[i])) : LX[i];
			   /*!!!*/vector(7) nDzz[i] = EMPTY[i] ? 0.0 : LX[i]<0 ? Ann[i] * Dzz3[i] + Bnn[i] * Dzz2[i] + Cnn[i] * Dzz[i] + Dnn[i] : nX[i];
			   /*!!! MOVED */vector(7) Dz4[i] = Dz3[i] * Dz[i];
			   /*!!! MOVED */vector(7) Dzz4[i] = Dzz3[i] * Dzz[i];
			   /*!!! MOVED */vector(7) Dz5[i] = Dz4[i] * Dz[i];
			   /*!!! MOVED */vector(7) Dzz5[i] = Dzz4[i] * Dzz[i];
			   /*!!! MOVED */vector(7) Dz6[i] = Dz5[i] * Dz[i];
			   /*!!! MOVED */vector(7) Dzz6[i] = Dzz5[i] * Dzz[i];
			   /*!!! MOVED */vector(7) Dz7[i] = Dz6[i] * Dz[i];
			   /*!!! MOVED */vector(7) Dzz7[i] = Dzz6[i] * Dzz[i];
			   /*!!!*/vector(7) _PhiM[i] = EMPTY[i] ? 0.0 : Pi*0.5*(nDzz[i] * LDzz[i] * Dzz[i] + nX[i] * LX[i] * X[i])*(Dzz[i] - X[i]);
			   /*!!!*/vector(7) PhiM[i] = EMPTY[i] ? 0.0 : fabs(Rok[i][Ptr]*(Dzz[i] - X[i])/(Y[i]-X[i]) / TAU) < fabs(_PhiM[i]) ? sign(0.5*(LX[i]+LDzz[i]))*max(0.0, Rok[i][Ptr]*(Dzz[i] - X[i])/(Y[i]-X[i]) / TAU) : _PhiM[i];
			   /*!!!*/vector(7) _PhiP[i] = EMPTY[i] ? 0.0 : Pi*0.5*(nY[i] * LY[i] * Y[i] + nDz[i] * LDz[i] * Dz[i])*(Y[i] - Dz[i]);
			   /*!!!*/vector(7) PhiP[i] = EMPTY[i] ? 0.0 : fabs(Rok[i][Ptr]*(Y[i] - Dz[i])/(Y[i]-X[i]) / TAU) < fabs(_PhiP[i]) ? sign(0.5*(LY[i]+LDz[i]))*max(0.0, Rok[i][Ptr]*(Y[i] - Dz[i])/(Y[i]-X[i]) / TAU) : _PhiP[i];
			   /* WOW SCHEME!!!!!!!!!!!!!!!! ERROR CORRECTED!!!!!!!!! */
			   /*!!!*/vector(7) NkM[i] = EMPTY[i] ? 0.0 : max(0.0, min(Nk[i][Ptr]*(Dzz[i]-X[i])/(Y[i]-X[i]), Ann[i] / 4.0*(Dzz4[i] - X4[i]) + Bnn[i] / 3.0*(Dzz3[i] - X3[i]) + Cnn[i] / 2.0*(Dzz2[i] - X2[i]) + Dnn[i] * (Dzz[i] - X[i])));
			   /* WOW SCHEME!!!!!!!!!!!!!!!! ERROR CORRECTED!!!!!!!!! */
			   /*!!!!!*/vector(7) NkP[i] = EMPTY[i] ? 0.0 : max(0.0, min(Nk[i][Ptr]*(Y[i]-Dz[i])/(Y[i]-X[i]), An[i] / 4.0*(Y4[i] - Dz4[i]) + Bn[i] / 3.0*(Y3[i] - Dz3[i]) + Cn[i] / 2.0*(Y2[i] - Dz2[i]) + Dn[i] * (Y[i] - Dz[i])));
			   /*!!!  WOW SCHEME ERROR CORRECTED!*/vector(7) RZ[i] = EMPTY[i] ? 0.0 : max(0.0, min(Pi*Nk[i][Ptr]*(Dz[i] - Dzz[i]), Pi*(
				   An[i] / 5.0*(Dz5[i] - XYC5[i]) + Bn[i] / 4.0*(Dz4[i] - XYC4[i]) + Cn[i] / 3.0*(Dz3[i] - XYC3[i]) + Dn[i] / 2.0*(Dz2[i] - XYC2[i]) +
				   Ann[i] / 5.0*(XYC5[i] - Dzz5[i]) + Bnn[i] / 4.0*(XYC4[i] - Dzz4[i]) + Cnn[i] / 3.0*(XYC3[i] - Dzz3[i]) + Dnn[i] / 2.0*(XYC2[i] - Dzz2[i])
				   )));
			   /*!!! WOW SCHEME ERROR CORRECTED!*/vector(7) RM[i] = EMPTY[i] ? 0.0 : max(0.0, min(Pi*Nk[i][Ptr]*(Dzz[i]-X[i]), Pi*(Ann[i] / 5.0*(Dzz5[i] - X5[i]) + Bnn[i] / 4.0*(Dzz4[i] - X4[i]) + Cnn[i] / 3.0*(Dzz3[i] - X3[i]) + Dnn[i] / 2.0*(Dzz2[i] - X2[i]))));
			   /*!!! WOW SCHEME ERROR CORRECTED!*/vector(7) RP[i] = EMPTY[i] ? 0.0 : max(0.0, min(Pi*Nk[i][Ptr]*(Y[i]-Dz[i]), Pi*(An[i] / 5.0*(Y5[i] - Dz5[i]) + Bn[i] / 4.0*(Y4[i] - Dz4[i]) + Cn[i] / 3.0*(Y3[i] - Dz3[i]) + Dn[i] / 2.0*(Y2[i] - Dz2[i]))));
               vector(7) pdNk[i] = 1.0/TAU*(-NkP[i]-NkM[i])+
                                (i>0.0 ? 1.0/TAU*NkP[i-1] : 0.0)+
                                (i<7-1 ? 1.0/TAU*NkM[i+1] : 0.0);

			   vector(7) PhiZ[i] = RZ[i] * LDp[i];
               vector(7) TOZERO[i] = X1[i]>=Dk[i+1] || Y1[i]<=Dk[i] || (Nk[i][Ptr]<1.0 && pdNk[i]<=0.0);
               vector(7) dNk[i] = TOZERO[i] ? -Nk[i][Ptr]/TAU : pdNk[i];
               vector(7) dSO2kg[i] = TOZERO[i] ? -1E6*SO2k[i][Ptr]/TAU :
                                    RZ[i]*LSO2g[i]+
                                    (i>0 ? RP[i-1]*LSO2g[i-1] : 0.0)+
                                    (i<7-1 ? RM[i+1]*LSO2g[i+1] : 0.0);
			   /* WOW SCHEME!!!!!!!!!!!!!!!! ERROR CORRECTED!!!!!!!!! */
			   /*!!!*/vector(7) RokM[i] = EMPTY[i] ? 0.0 : max(0.0, min(Rok[i][Ptr], Pi / 6.0*ROK*(Ann[i] / 7.0*(Dzz7[i] - X7[i]) + Bnn[i] / 6.0*(Dzz6[i] - X6[i]) + Cnn[i] / 5.0*(Dzz5[i] - X5[i]) + Dnn[i] / 4.0*(Dzz4[i] - X4[i]))));
			   /* WOW SCHEME!!!!!!!!!!!!!!!! ERROR CORRECTED!!!!!!!!! */
			   /*!!!*/vector(7) RokP[i] = EMPTY[i] ? 0.0 : max(0.0, min(Rok[i][Ptr], Pi / 6.0*ROK*(An[i] / 7.0*(Y7[i] - Dz7[i]) + Bn[i] / 6.0*(Y6[i] - Dz6[i]) + Cn[i] / 5.0*(Y5[i] - Dz5[i]) + Dn[i] / 4.0*(Y4[i] - Dz4[i]))));
               dC = ((PhiZ[0])+(PhiZ[1])+(PhiZ[2])+(PhiZ[3])+(PhiZ[4])+(PhiZ[5])+(PhiZ[6]))+((PhiP[1-1])+(PhiP[2-1])+(PhiP[3-1])+(PhiP[4-1])+(PhiP[5-1])+(PhiP[6-1]))+((PhiM[0+1])+(PhiM[1+1])+(PhiM[2+1])+(PhiM[3+1])+(PhiM[4+1])+(PhiM[5+1]));
               vector(7) dRok[i] = TOZERO[i] ? -Rok[i][Ptr]/TAU :
                                  PhiZ[i]+1.0/TAU*(-RokP[i]-RokM[i])+
                                  (i>0.0 ? PhiP[i-1]+1.0/TAU*RokP[i-1] : 0.0)+
                                  (i<7-1 ? PhiM[i+1]+1.0/TAU*RokM[i+1] : 0.0);
               vector(7) dSO2kf[i] = TOZERO[i] ? 0.0 :
                                    RZ[i]*LSO2f[i]+1E6/(TAU*ROK)*(-RokP[i]-RokM[i])+
                                    (i>0 ? RP[i-1]*LSO2f[i-1]+1E6/TAU*RokP[i-1]/ROK : 0.0)+
                                    (i<7-1 ? RM[i+1]*LSO2f[i+1]+1E6/TAU*RokM[i+1]/ROK : 0.0);

			   /*!!!!*/ vector(7) ddN[i] = 0.0;
			   /*!!!!*/ vector(7) dRo[i] = 0.0;
			   /*!!!!*/ coalesc = calcCoalescense(Nc, Dk, X, Y, XYC, Ann, Bnn, Cnn, Dnn, An, Bn, Cn, Dn, Ptr, Ukx, Uky, Ukz, UwLight, ROK, ddN, dRo, ddRo); /* NEW CORR */

/* NEW */		vector(7) ddSO2kg[i] = (
                              (i != 0 ? SO2k[0][Ptr]*ddRo[0*(7+1)+i] : 0.0) +
                              (i != 1 ? SO2k[1][Ptr]*ddRo[1*(7+1)+i] : 0.0) +
                              (i != 2 ? SO2k[2][Ptr]*ddRo[2*(7+1)+i] : 0.0) +
                              (i != 3 ? SO2k[3][Ptr]*ddRo[3*(7+1)+i] : 0.0) +
                              (i != 4 ? SO2k[4][Ptr]*ddRo[4*(7+1)+i] : 0.0) +
                              (i != 5 ? SO2k[5][Ptr]*ddRo[5*(7+1)+i] : 0.0) +
                              (i != 6 ? SO2k[6][Ptr]*ddRo[6*(7+1)+i] : 0.0)
                        )/ROK;
/* NEW */		vector(7) ddSO2kf[i] = (ddRo[i*(7+1)+i] + ddRo[i*(7+1)+7])/ROK;

/* NEW */		   /*!!!!!*/ vector(7) CTOZERO[i] = SO2k[i][Ptr] + TAU*(1E-6*dSO2kg[i] + ddSO2kg[i] + SO2k[i][Ptr] * (1E-6*dSO2kf[i] + ddSO2kf[i])) < 0;
/* NEW */		   /*!!!!!*/ vector(7) CTOMAX[i] = SO2k[i][Ptr] + TAU*(1E-6*dSO2kg[i] + ddSO2kg[i] + SO2k[i][Ptr] * (1E-6*dSO2kf[i] + ddSO2kf[i])) > maxSO2k[i];

			   /*!!! CHANGED && MOVED !!!*/dSO2 = (
				   (CTOZERO[0] ? 0.0 : (RZ[0] * (LSO2f[0] * max(0.0, SO2k[0][Ptr]) + LSO2g[0]))) +
				   (CTOZERO[1] ? 0.0 : (RZ[1] * (LSO2f[1] * max(0.0, SO2k[1][Ptr]) + LSO2g[1]))) +
				   (CTOZERO[2] ? 0.0 : (RZ[2] * (LSO2f[2] * max(0.0, SO2k[2][Ptr]) + LSO2g[2]))) +
				   (CTOZERO[3] ? 0.0 : (RZ[3] * (LSO2f[3] * max(0.0, SO2k[3][Ptr]) + LSO2g[3]))) +
				   (CTOZERO[4] ? 0.0 : (RZ[4] * (LSO2f[4] * max(0.0, SO2k[4][Ptr]) + LSO2g[4]))) +
				   (CTOZERO[5] ? 0.0 : (RZ[5] * (LSO2f[5] * max(0.0, SO2k[5][Ptr]) + LSO2g[5]))) +
				   (CTOZERO[6] ? 0.0 : (RZ[6] * (LSO2f[6] * max(0.0, SO2k[6][Ptr]) + LSO2g[6])))
				   ) +
			       (
				   (CTOZERO[1-1] ? 0.0 : (RP[1 - 1] * (LSO2f[1 - 1] * max(0.0, SO2k[1 - 1][Ptr]) + LSO2g[1 - 1]))) +
				   (CTOZERO[2-1] ? 0.0 : (RP[2 - 1] * (LSO2f[2 - 1] * max(0.0, SO2k[2 - 1][Ptr]) + LSO2g[2 - 1]))) +
				   (CTOZERO[3-1] ? 0.0 : (RP[3 - 1] * (LSO2f[3 - 1] * max(0.0, SO2k[3 - 1][Ptr]) + LSO2g[3 - 1]))) +
				   (CTOZERO[4-1] ? 0.0 : (RP[4 - 1] * (LSO2f[4 - 1] * max(0.0, SO2k[4 - 1][Ptr]) + LSO2g[4 - 1]))) +
				   (CTOZERO[5-1] ? 0.0 : (RP[5 - 1] * (LSO2f[5 - 1] * max(0.0, SO2k[5 - 1][Ptr]) + LSO2g[5 - 1]))) +
				   (CTOZERO[6-1] ? 0.0 : (RP[6 - 1] * (LSO2f[6 - 1] * max(0.0, SO2k[6 - 1][Ptr]) + LSO2g[6 - 1])))
				   ) +
			       (
				   (CTOZERO[0+1] ? 0.0 : (RM[0 + 1] * (LSO2f[0 + 1] * max(0.0, SO2k[0 + 1][Ptr]) + LSO2g[0 + 1]))) +
				   (CTOZERO[1+1] ? 0.0 : (RM[1 + 1] * (LSO2f[1 + 1] * max(0.0, SO2k[1 + 1][Ptr]) + LSO2g[1 + 1]))) +
				   (CTOZERO[2+1] ? 0.0 : (RM[2 + 1] * (LSO2f[2 + 1] * max(0.0, SO2k[2 + 1][Ptr]) + LSO2g[2 + 1]))) +
				   (CTOZERO[3+1] ? 0.0 : (RM[3 + 1] * (LSO2f[3 + 1] * max(0.0, SO2k[3 + 1][Ptr]) + LSO2g[3 + 1]))) +
				   (CTOZERO[4+1] ? 0.0 : (RM[4 + 1] * (LSO2f[4 + 1] * max(0.0, SO2k[4 + 1][Ptr]) + LSO2g[4 + 1]))) +
				   (CTOZERO[5+1] ? 0.0 : (RM[5 + 1] * (LSO2f[5 + 1] * max(0.0, SO2k[5 + 1][Ptr]) + LSO2g[5 + 1])))
				   );

			   Kf[_NumP][Ptr] += ((FastMode ? (0.0) : 0.0) + (SlowMode ? (-RO*(dUxdx*dUxdx + dUydy*dUydy + dUzdz*dUzdz +
                                           2.0*(dUydx*dUxdy+dUydz*dUzdy+dUxdz*dUzdx)
                                           -0.1*DIV[Ptr]/TAU
                                           -d2NFdx2(Nu,DIV)
                                           -d2NFdy2(Nu,DIV)
                                           -d2NFdz2(Nu,DIV))) : 0.0));
               
               Kf[_NumP][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (RO*BT*9.81*dFdz(T)) : 0.0));
               
               Kf[_NumUx][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (-1.0/RO*dFdx(P)) : 0.0));
               
               Kf[_NumUy][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (-1.0/RO*dFdy(P)) : 0.0));
               
               Kf[_NumUz][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (-1.0/RO*dFdz(P)) : 0.0));
               
               Kf[_NumUz][Ptr] +=  BT*9.81*(T[Ptr]-T0);
               
               Kf[_NumKt][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (DDD-Et[Ptr]-BT*9.81*Nu[Ptr]/0.7194*dFdz(T)) : 0.0));
               
               Kf[_NumEt][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? ( (1.42*DDD-(1.68+0.0845*RNG3*(1-RNG/4.38)/(1+0.012*RNG3))*Et[Ptr])*Et[Ptr]/Kt[Ptr]) : 0.0));
               
               Kf[_NumNu][Ptr] +=  0.0845 * Kt[Ptr] * Kt[Ptr]/Et[Ptr];
               
               Kf[_NumFUV][Ptr] +=  ueBeta;
               
               Kf[_NumUV][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (-uFt*utBeta) : 0.0));
               
               Kf[_NumF0][Ptr] +=  eBeta;
               
               Kf[_NumI][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (-Ft*tBeta) : 0.0));
               
               Kf[_NumFIR][Ptr] +=  reBeta;
               
               Kf[_NumIR][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (-rFt*rtBeta) : 0.0));
               
               Kf[_NumIRt][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (-tFt*ttBeta) : 0.0));
               
               Kf[_NumT][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? ((tAlpha*(4.0*Pi*I[Ptr]+F0[Ptr])+rtAlpha*(4.0*Pi*IR[Ptr]+FIR[Ptr])+4.0*Pi*(ttAlpha*IRt[Ptr]-tFt/3.0))/RO0/CC0) : 0.0));

			   vector(7) { /*!!!*/
                K_Rok[i][Ptr] +=  ((FastMode ? (dRok[i]+dRo[i]) : 0.0)+(SlowMode ? (0.0) : 0.0));
               
               }
               
			   vector(7) { /*!!!!*/
				   K_Rok[i][Ptr] += ((FastMode ? (0.0) : 0.0) + (SlowMode ? (_Z[0] == 0.0 && z == 1 && InjectStartTime <= ModelTime ?
					   (
					   EMPTY[i] ? InjectN[i] / 24.0 / (Dk[i + 1] - Dk[i])*Pi*ROK*(Dk[i + 1] * Dk[i + 1] * Dk[i + 1] * Dk[i + 1] - Dk[i] * Dk[i] * Dk[i] * Dk[i]) :
					   InjectN[i] / 6.0 / XY1[i] * Pi*ROK*XY4[i]
					   ) : 0.0) : 0.0));
               
               }
               
               vector(7) {
                K_Ukx[i][Ptr] +=  Ux[Ptr];
               
               }
               
               vector(7) {
                K_Uky[i][Ptr] +=  Uy[Ptr];
               
               }
               
			   vector(7) { /*!!!!!!!*/
				   K_Ukz[i][Ptr] += ((FastMode ? (0.0) : 0.0) + (SlowMode ? (fabs(UwLight[i]) > 1E-8 ? Uz[Ptr] - UwLight[i] : Uz[Ptr] - GetdU(9.81, RO0, ROK, _Sigma, D, DD[i], Uz[Ptr] - Ukz[i][Ptr])) : 0.0));
               
               }

			   vector(7) { /*!!!!!*/
                K_Nk[i][Ptr] +=  ((FastMode ? (dNk[i]+ddN[i]) : 0.0)+(SlowMode ? (0.0) : 0.0));
               
               }
               
			   vector(7) { /*!!!!!*/
                K_Nk[i][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (_Z[0]==0.0 && z==1 && InjectStartTime<=ModelTime ? InjectN[i] : 0.0) : 0.0));
               
               }
               
/* NEW */		   vector(7) { /*!!!!!*/
				   K_SO2k[i][Ptr] += ((FastMode ?
					   (CTOZERO[i] ? -SO2k[i][Ptr] / TAU : CTOMAX[i] ? maxSO2k[i]/TAU : (1E-6*dSO2kg[i] + ddSO2kg[i])) :
						   0.0)
					   + (SlowMode ? (0.0) : 0.0));
               }
               
               K_SO2[Ptr] += (_Z[0]==0.0 && z==1 ? 1.114E-14 : 0.0); /* !!!!! */
               
			   /*!!!!!*/K_SO2[Ptr] += ((FastMode ? (-1E-6*dSO2) : 0.0) + (SlowMode ? (0.0) : 0.0));

               Sf[_NumFUV][Ptr] +=  uSolarIntensity;
               
               Sf[_NumUV][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (3.0*utAlpha*utBeta) : 0.0));
               
               Sf[_NumF0][Ptr] +=  SolarIntensity;
               
               Sf[_NumI][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (3.0*tAlpha*tBeta) : 0.0));
               
               Sf[_NumFIR][Ptr] +=  rSolarIntensity;
               
               Sf[_NumIR][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (3.0*rtAlpha*rtBeta) : 0.0));
               
               Sf[_NumIRt][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (3.0*ttAlpha*ttBeta) : 0.0));
               
			   vector(7) { /*!!!!!*/
                S_Rok[i][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (0.0/*-dUkxdx[i]-dUkydy[i]-dUkzdz[i]*/) : 0.0));
               
               }
               
			   vector(7) { /*!!!!!*/
                S_Nk[i][Ptr] +=  ((FastMode ? (0.0) : 0.0)+(SlowMode ? (0.0/*-dUkxdx[i]-dUkydy[i]-dUkzdz[i]*/) : 0.0));
               
               }
               
/* NEW */		   vector(7) { /*!!!!!!!*/
				   S_SO2k[i][Ptr] += ((FastMode ? ((CTOZERO[i] ? 0.0 : CTOMAX[i] ? -1.0/TAU : (1E-6*dSO2kf[i] + ddSO2kf[i]))) : 0.0)
					+ (SlowMode ? (0.0) : 0.0));
               
               }
               
			   if (IsBound && !IsExchng) Ux[Ptr] = IsCond(Ux, 'B') ?
                                      2.5*Kx*Uref*K1 :
                                      Ux[Ptr];
               
               if (IsBound && !IsExchng) Uy[Ptr] =  IsCond(Uy,'B') ?
                                       2.5*Ky*Uref*K1 :
                                       Uy[Ptr];
               
               if (IsBound && !IsExchng) Kt[Ptr] =  IsCond(Nu,'B') ? NuK : Kt[Ptr];
               
               if (IsBound && !IsExchng) if (IsCond(Et,'A'))
                                   {
                                    if (IsForw)   _fw(Et[Ptr]) = pow(0.0845,0.75)*pow(_fw(Kt[Ptr]),1.5)/0.41/HY[YM];
                                    if (IsBack)   _bw(Et[Ptr]) = pow(0.0845,0.75)*pow(_bw(Kt[Ptr]),1.5)/0.41/HY[y];
                                    if (IsBottom) _bt(Et[Ptr]) = pow(0.0845,0.75)*pow(_bt(Kt[Ptr]),1.5)/0.41/HZ[z];
                                    if (IsTop)    _tp(Et[Ptr]) = pow(0.0845,0.75)*pow(_tp(Kt[Ptr]),1.5)/0.41/HZ_ZM;
                                    if (IsLeft)   _lf(Et[Ptr]) = pow(0.0845,0.75)*pow(_lf(Kt[Ptr]),1.5)/0.41/HX[x];
                                    if (IsRight)  _rg(Et[Ptr]) = pow(0.0845,0.75)*pow(_rg(Kt[Ptr]),1.5)/0.41/HX[XM];
                                   }
               
               if (IsBound && !IsExchng) Et[Ptr] =  IsCond(Nu,'B') ? NuE : Et[Ptr];
               
               if (IsBound && !IsExchng) if (IsCond(Nu,'A'))
                                   {
                                    if (IsForw)   _fw(Nu[Ptr]) = 0.0845*pow(_fw(Kt[Ptr]),2)/_fw(Et[Ptr]);
                                    if (IsBack)   _bw(Nu[Ptr]) = 0.0845*pow(_bw(Kt[Ptr]),2)/_bw(Et[Ptr]);
                                    if (IsBottom) _bt(Nu[Ptr]) = 0.0845*pow(_bt(Kt[Ptr]),2)/_bt(Et[Ptr]);
                                    if (IsTop)    _tp(Nu[Ptr]) = 0.0845*pow(_tp(Kt[Ptr]),2)/_tp(Et[Ptr]);
                                    if (IsLeft)   _lf(Nu[Ptr]) = 0.0845*pow(_lf(Kt[Ptr]),2)/_lf(Et[Ptr]);
                                    if (IsRight)  _rg(Nu[Ptr]) = 0.0845*pow(_rg(Kt[Ptr]),2)/_rg(Et[Ptr]);
                                   }
               
               if (IsBound && !IsExchng) Nu[Ptr] =  IsCond(Nu,'B') ? 0.0845*NuK*NuK/NuE : Nu[Ptr];
               
               if (IsBound && !IsExchng)       if (IsClosed(UV))
                       {
                        if (IsForw)   _fw(UV[Ptr]) = (uk1a*UV[ZYMX]+HY[y-1]*SolarY*_fw(FUV[Ptr])*SurfaceAlbedo/Pi)/(uk1a+HY[y-1]*(1-SurfaceAlbedo));
                        if (IsBack)   _bw(UV[Ptr]) = (uk1a*UV[ZYPX]-HY[y]*SolarY*  _bw(FUV[Ptr])*SurfaceAlbedo/Pi)/(uk1a+HY[y]*(1-SurfaceAlbedo));
                        if (IsBottom) _bt(UV[Ptr]) = (uk1a*UV[ZPYX]-HZ[z]*SolarZ*  _bt(FUV[Ptr])*SurfaceAlbedo/Pi)/(uk1a+HZ[z]*(1-SurfaceAlbedo));
                        if (IsTop)    _tp(UV[Ptr]) = (uk1a*UV[ZMYX]+HZ_ZM*SolarZ*  _tp(FUV[Ptr])*SurfaceAlbedo/Pi)/(uk1a+HZ_ZM*(1-SurfaceAlbedo));
                        if (IsLeft)   _lf(UV[Ptr]) = (uk1a*UV[ZYXP]-HX[x]*SolarX*  _lf(FUV[Ptr])*SurfaceAlbedo/Pi)/(uk1a+HX[x]*(1-SurfaceAlbedo));
                        if (IsRight)  _rg(UV[Ptr]) = (uk1a*UV[ZYXM]+HX[x-1]*SolarX*_rg(FUV[Ptr])*SurfaceAlbedo/Pi)/(uk1a+HX[x-1]*(1-SurfaceAlbedo));
                       }
                     else
                       {
                        if (IsForw)   _fw(UV[Ptr]) = (uk1*UV[ZYMX]+HY[y-1]*uFluxPi)/(uk1+HY[y-1]);
                        if (IsBack)   _bw(UV[Ptr]) = (uk1*UV[ZYPX]+HY[y]*uFluxPi)/(uk1+HY[y]);
                        if (IsBottom) _bt(UV[Ptr]) = (uk1*UV[ZPYX]+HZ[z]*uFluxPi)/(uk1+HZ[z]);
                        if (IsTop)    _tp(UV[Ptr]) = (uk1*UV[ZMYX]+HZ_ZM*uFluxPi)/(uk1+HZ_ZM);
                        if (IsLeft)   _lf(UV[Ptr]) = (uk1*UV[ZYXP]+HX[x]*uFluxPi)/(uk1+HX[x]);
                        if (IsRight)  _rg(UV[Ptr]) = (uk1*UV[ZYXM]+HX[x-1]*uFluxPi)/(uk1+HX[x-1]);
                       }
			   if (IsBound && !IsExchng)       if (IsClosed(I))
                       {
                        if (IsForw)   _fw(I[Ptr]) = (k1a*I[ZYMX]+HY[y-1]*SolarY*_fw(F0[Ptr])*SurfaceAlbedo/Pi)/(k1a+HY[y-1]*(1-SurfaceAlbedo));
                        if (IsBack)   _bw(I[Ptr]) = (k1a*I[ZYPX]-HY[y]*SolarY*  _bw(F0[Ptr])*SurfaceAlbedo/Pi)/(k1a+HY[y]*(1-SurfaceAlbedo));
                        if (IsBottom) _bt(I[Ptr]) = (k1a*I[ZPYX]-HZ[z]*SolarZ*    _bt(F0[Ptr])*SurfaceAlbedo/Pi)/(k1a+HZ[z]*(1-SurfaceAlbedo));
                        if (IsTop)    _tp(I[Ptr]) = (k1a*I[ZMYX]+HZ_ZM*SolarZ*  _tp(F0[Ptr])*SurfaceAlbedo/Pi)/(k1a+HZ_ZM*(1-SurfaceAlbedo));
                        if (IsLeft)   _lf(I[Ptr]) = (k1a*I[ZYXP]-HX[x]*SolarX* _lf(F0[Ptr])*SurfaceAlbedo/Pi)/(k1a+HX[x]*(1-SurfaceAlbedo));
                        if (IsRight)  _rg(I[Ptr]) = (k1a*I[ZYXM]+HX[x-1]*SolarX*_rg(F0[Ptr])*SurfaceAlbedo/Pi)/(k1a+HX[x-1]*(1-SurfaceAlbedo));
			           }
                     else
                       {
                        if (IsForw)   _fw(I[Ptr]) = (k1*I[ZYMX]+HY[y-1]*FluxPi)/(k1+HY[y-1]);
                        if (IsBack)   _bw(I[Ptr]) = (k1*I[ZYPX]+HY[y]*FluxPi)/(k1+HY[y]);
                        if (IsBottom) _bt(I[Ptr]) = (k1*I[ZPYX]+HZ[z]*FluxPi)/(k1+HZ[z]);
                        if (IsTop)    _tp(I[Ptr]) = (k1*I[ZMYX]+HZ_ZM*FluxPi)/(k1+HZ_ZM);
                        if (IsLeft)   _lf(I[Ptr]) = (k1*I[ZYXP]+HX[x]*FluxPi)/(k1+HX[x]);
                        if (IsRight)  _rg(I[Ptr]) = (k1*I[ZYXM]+HX[x-1]*FluxPi)/(k1+HX[x-1]);
                       }
               
               if (IsBound && !IsExchng)       if (IsClosed(IR))
                       {
                        if (IsForw)   _fw(IR[Ptr]) = (rk1a*IR[ZYMX]+HY[y-1]*SolarY*_fw(FIR[Ptr])*SurfaceAlbedo/Pi)/(rk1a+HY[y-1]*(1-SurfaceAlbedo));
                        if (IsBack)   _bw(IR[Ptr]) = (rk1a*IR[ZYPX]-HY[y]*SolarY*  _bw(FIR[Ptr])*SurfaceAlbedo/Pi)/(rk1a+HY[y]*(1-SurfaceAlbedo));
                        if (IsBottom) _bt(IR[Ptr]) = (rk1a*IR[ZPYX]-HZ[z]*SolarZ*    _bt(FIR[Ptr])*SurfaceAlbedo/Pi)/(rk1a+HZ[z]*(1-SurfaceAlbedo));
                        if (IsTop)    _tp(IR[Ptr]) = (rk1a*IR[ZMYX]+HZ_ZM*SolarZ*  _tp(FIR[Ptr])*SurfaceAlbedo/Pi)/(rk1a+HZ_ZM*(1-SurfaceAlbedo));
                        if (IsLeft)   _lf(IR[Ptr]) = (rk1a*IR[ZYXP]-HX[x]*SolarX* _lf(FIR[Ptr])*SurfaceAlbedo/Pi)/(rk1a+HX[x]*(1-SurfaceAlbedo));
                        if (IsRight)  _rg(IR[Ptr]) = (rk1a*IR[ZYXM]+HX[x-1]*SolarX*_rg(FIR[Ptr])*SurfaceAlbedo/Pi)/(rk1a+HX[x-1]*(1-SurfaceAlbedo));
                       }
                     else
                       {
                        if (IsForw)   _fw(IR[Ptr]) = (rk1*IR[ZYMX]+HY[y-1]*rFluxPi)/(rk1+HY[y-1]);
                        if (IsBack)   _bw(IR[Ptr]) = (rk1*IR[ZYPX]+HY[y]*rFluxPi)/(rk1+HY[y]);
                        if (IsBottom) _bt(IR[Ptr]) = (rk1*IR[ZPYX]+HZ[z]*rFluxPi)/(rk1+HZ[z]);
                        if (IsTop)    _tp(IR[Ptr]) = (rk1*IR[ZMYX]+HZ_ZM*rFluxPi)/(rk1+HZ_ZM);
                        if (IsLeft)   _lf(IR[Ptr]) = (rk1*IR[ZYXP]+HX[x]*rFluxPi)/(rk1+HX[x]);
                        if (IsRight)  _rg(IR[Ptr]) = (rk1*IR[ZYXM]+HX[x-1]*rFluxPi)/(rk1+HX[x-1]);
                       }
               
               if (IsBound && !IsExchng)       if (IsClosed(IRt))
                       {
                        if (IsForw)   _fw(IRt[Ptr]) = (tk1a*IRt[ZYMX]+HY[y-1]*SurfaceEmissivity*(IntPlanck[0]+(273.15+_fw(T[Ptr]))*(IntPlanck[1]+(273.15+_fw(T[Ptr]))*(IntPlanck[2]+(273.15+_fw(T[Ptr]))*IntPlanck[3])))/Pi)/(tk1a+HY[y-1]*SurfaceEmissivity);
                        if (IsBack)   _bw(IRt[Ptr]) = (tk1a*IRt[ZYPX]+HY[y]*  SurfaceEmissivity*(IntPlanck[0]+(273.15+_bw(T[Ptr]))*(IntPlanck[1]+(273.15+_bw(T[Ptr]))*(IntPlanck[2]+(273.15+_bw(T[Ptr]))*IntPlanck[3])))/Pi)/(tk1a+HY[y]*SurfaceEmissivity);
                        if (IsBottom) _bt(IRt[Ptr]) = (tk1a*IRt[ZPYX]+HZ[z]*    SurfaceEmissivity*(IntPlanck[0]+(273.15+_bt(T[Ptr]))*(IntPlanck[1]+(273.15+_bt(T[Ptr]))*(IntPlanck[2]+(273.15+_bt(T[Ptr]))*IntPlanck[3])))/Pi)/(tk1a+HZ[z]*SurfaceEmissivity);
                        if (IsTop)    _tp(IRt[Ptr]) = (tk1a*IRt[ZMYX]+HZ_ZM*  SurfaceEmissivity*(IntPlanck[0]+(273.15+_tp(T[Ptr]))*(IntPlanck[1]+(273.15+_tp(T[Ptr]))*(IntPlanck[2]+(273.15+_tp(T[Ptr]))*IntPlanck[3])))/Pi)/(tk1a+HZ_ZM*SurfaceEmissivity);
                        if (IsLeft)   _lf(IRt[Ptr]) = (tk1a*IRt[ZYXP]+HX[x]* SurfaceEmissivity*(IntPlanck[0]+(273.15+_lf(T[Ptr]))*(IntPlanck[1]+(273.15+_lf(T[Ptr]))*(IntPlanck[2]+(273.15+_lf(T[Ptr]))*IntPlanck[3])))/Pi)/(tk1a+HX[x]*SurfaceEmissivity);
                        if (IsRight)  _rg(IRt[Ptr]) = (tk1a*IRt[ZYXM]+HX[x-1]*SurfaceEmissivity*(IntPlanck[0]+(273.15+_rg(T[Ptr]))*(IntPlanck[1]+(273.15+_rg(T[Ptr]))*(IntPlanck[2]+(273.15+_rg(T[Ptr]))*IntPlanck[3])))/Pi)/(tk1a+HX[x-1]*SurfaceEmissivity);
                       }
                     else
                       {
                        if (IsForw)   _fw(IRt[Ptr]) = (tk1*IRt[ZYMX]+HY[y-1]*tFluxPi)/(tk1+HY[y-1]);
                        if (IsBack)   _bw(IRt[Ptr]) = (tk1*IRt[ZYPX]+HY[y]*tFluxPi)/(tk1+HY[y]);
                        if (IsBottom) _bt(IRt[Ptr]) = (tk1*IRt[ZPYX]+HZ[z]*tFluxPi)/(tk1+HZ[z]);
                        if (IsTop)    _tp(IRt[Ptr]) = (tk1*IRt[ZMYX]+HZ_ZM*tFluxPi)/(tk1+HZ_ZM);
                        if (IsLeft)   _lf(IRt[Ptr]) = (tk1*IRt[ZYXP]+HX[x]*tFluxPi)/(tk1+HX[x]);
                        if (IsRight)  _rg(IRt[Ptr]) = (tk1*IRt[ZYXM]+HX[x-1]*tFluxPi)/(tk1+HX[x-1]);
                       }

			   if (IsBound && !IsExchng)    if (IsClosed(T) && SlowMode)
                     {
                      if (IsForw)
                         {
                          double TKK = 273.15+_fw(T[Ptr]);
                          double LambdaT = 1E-3*(4.39+0.071*TKK);
               
                          _fw(T[Ptr]) += TAU/SurfaceCC/SurfaceRO*(SurfaceLambda*(_d2Fdx2(_forw,T)  +_d2Fdz2(_forw,T))+((1-SurfaceAlbedo)*  (Pi*(_fw(I[Ptr])+_fw(IR[Ptr])+_fw(UV[Ptr])-k1*dFdyf(I)-rk1*dFdyf(IR)-uk1*dFdyf(UV))+(_fw(F0[Ptr])+_fw(FIR[Ptr])+_fw(FUV[Ptr]))*SolarY)+SurfaceEmissivity*(Pi*(_fw(IRt[Ptr])-tk1*dFdyf(IRt))-(IntPlanck[0]+TKK*(IntPlanck[1]+TKK*(IntPlanck[2]+TKK*IntPlanck[3]))))+LambdaT/HY[YM]*(T[ZYMX]-_fw(T[Ptr])))/HY[YM]);
					  }
                      if (IsBack)
                         {
                          double TKK = 273.15+_bw(T[Ptr]);
                          double LambdaT = 1E-3*(4.39+0.071*TKK);
               
                          _bw(T[Ptr]) += TAU/SurfaceCC/SurfaceRO*(SurfaceLambda*(_d2Fdx2(_back,T)  +_d2Fdz2(_back,T))+((1-SurfaceAlbedo)*  (Pi*(_bw(I[Ptr])+_bw(IR[Ptr])+_bw(UV[Ptr])+k1*dFdyb(I)+rk1*dFdyb(IR)+uk1*dFdyb(UV))-(_bw(F0[Ptr])+_bw(FIR[Ptr])+_bw(FUV[Ptr]))*SolarY)+SurfaceEmissivity*(Pi*(_bw(IRt[Ptr])+tk1*dFdyb(IRt))-(IntPlanck[0]+TKK*(IntPlanck[1]+TKK*(IntPlanck[2]+TKK*IntPlanck[3]))))+LambdaT/HY[y]*(T[ZYPX]- _bw(T[Ptr])))/HY[y]);
					  }
                      if (IsBottom)
                         {
                          double TKK = 273.15+_bt(T[Ptr]);
                          double LambdaT = 1E-3*(4.39+0.071*TKK);
               
                          _bt(T[Ptr]) += TAU/SurfaceCC/SurfaceRO*(SurfaceLambda*(_d2Fdx2(_bottom,T)+_d2Fdy2(_bottom,T))+((1-SurfaceAlbedo)*(Pi*(_bt(I[Ptr])+_bt(IR[Ptr])+_bt(UV[Ptr])+k1*dFdzb(I)+rk1*dFdzb(IR)+uk1*dFdzb(UV))-(_bt(F0[Ptr])+_bt(FIR[Ptr])+_bt(FUV[Ptr]))*SolarZ)+SurfaceEmissivity*(Pi*(_bt(IRt[Ptr])+tk1*dFdzb(IRt))-(IntPlanck[0]+TKK*(IntPlanck[1]+TKK*(IntPlanck[2]+TKK*IntPlanck[3]))))+LambdaT/HZ[z]*(T[ZPYX]-   _bt(T[Ptr])))/HZ[z]);
					  }
                      if (IsTop)
                         {
                          double TKK = 273.15+_tp(T[Ptr]);
                          double LambdaT = 1E-3*(4.39+0.071*TKK);
               
                          _tp(T[Ptr]) += TAU/SurfaceCC/SurfaceRO*(SurfaceLambda*(_d2Fdx2(_top,T)   +_d2Fdy2(_top,T))+((1-SurfaceAlbedo)*   (Pi*(_tp(I[Ptr])+_tp(IR[Ptr])+_tp(UV[Ptr])-k1*dFdzt(I)-rk1*dFdzt(IR)-uk1*dFdzt(UV))+(_tp(F0[Ptr])+_tp(FIR[Ptr])+_tp(FUV[Ptr]))*SolarZ)+SurfaceEmissivity*(Pi*(_tp(IRt[Ptr])-tk1*dFdzt(IRt))-(IntPlanck[0]+TKK*(IntPlanck[1]+TKK*(IntPlanck[2]+TKK*IntPlanck[3]))))+LambdaT/HZ_ZM*(T[ZMYX]- _tp(T[Ptr])))/HZ_ZM);
					  }
                      if (IsLeft)
                         {
                          double TKK = 273.15+_lf(T[Ptr]);
                          double LambdaT = 1E-3*(4.39+0.071*TKK);
               
                          _lf(T[Ptr]) += TAU/SurfaceCC/SurfaceRO*(SurfaceLambda*(_d2Fdy2(_left,T)  +_d2Fdz2(_left,T))+((1-SurfaceAlbedo)*  (Pi*(_lf(I[Ptr])+_lf(IR[Ptr])+_lf(UV[Ptr])+k1*dFdxl(I)+rk1*dFdxl(IR)+uk1*dFdxl(UV))-(_lf(F0[Ptr])+_lf(FIR[Ptr])+_lf(FUV[Ptr]))*SolarX)+SurfaceEmissivity*(Pi*(_lf(IRt[Ptr])+tk1*dFdxl(IRt))-(IntPlanck[0]+TKK*(IntPlanck[1]+TKK*(IntPlanck[2]+TKK*IntPlanck[3]))))+LambdaT/HX[x]*(T[ZYXP]-_lf(T[Ptr])))/HX[x]);
					  }
                      if (IsRight)
                         {
                          double TKK = 273.15+_rg(T[Ptr]);
                          double LambdaT = 1E-3*(4.39+0.071*TKK);
               
                          _rg(T[Ptr]) += TAU/SurfaceCC/SurfaceRO*(SurfaceLambda*(_d2Fdy2(_right,T) +_d2Fdz2(_right,T))+((1-SurfaceAlbedo)* (Pi*(_rg(I[Ptr])+_rg(IR[Ptr])+_rg(UV[Ptr])-k1*dFdxr(I)-rk1*dFdxr(IR)-uk1*dFdxr(UV))+(_rg(F0[Ptr])+_rg(FIR[Ptr])+_rg(FUV[Ptr]))*SolarX)+SurfaceEmissivity*(Pi*(_rg(IRt[Ptr])-tk1*dFdxr(IRt))-(IntPlanck[0]+TKK*(IntPlanck[1]+TKK*(IntPlanck[2]+TKK*IntPlanck[3]))))+LambdaT/HX[XM]*(T[ZYXM]-_rg(T[Ptr])))/HX[XM]);
					  }
                     }
               KDn[PhaseLinks[phDrops0]][Ptr] +=  ((FastMode ? (-1E-6/Mk*dC) : 0.0)+(SlowMode ? (0.0) : 0.0));
			  }
     }
 #if defined(__IMITATE__) && !defined(WIN32)
   Check();
 #endif
}

double _ZeroNu (int i, int x, int y, int z, unsigned char Map)
{
 double Result;
 double KT = _ZeroKt(i,x,y,z,Map);
 double ET = _ZeroEt(i,x,y,z,Map);
 
 Result = 0.0845*KT*KT/ET;

 return Result;
}
double _ZeroT (int i, int x, int y, int z, unsigned char Map)
{
 return ZeroT - (_Z[z]-_Z[0])*0.01; /*!!!!!*/
}
double _ZeroUx (int i, int x, int y, int z, unsigned char Map)
{
 double Result;
 double KScale = 0.208816;
 double Uref = KScale*0.5263;
 double Z0 = 0.05;
 double ZZ = z==0 && _Z[0]==0.0 ? Z0 : (_Z[z]+_Z[z-1])*0.05+Z0;
 double K0 = Uref*ZZ;
 double K1 = log(ZZ/Z0);
 
 Result = 2.5*Kx*Uref*K1;

 return Result;
}
double _ZeroUy (int i, int x, int y, int z, unsigned char Map)
{
 double Result;
 double KScale = 0.208816;
 double Uref = KScale*0.5263;
 double Z0 = 0.05;
 double ZZ = z==0 && _Z[0]==0.0 ? Z0 : (_Z[z]+_Z[z-1])*0.05+Z0;
 double K0 = Uref*ZZ;
 double K1 = log(ZZ/Z0);
 
 Result = 2.5*Ky*Uref*K1;

 return Result;
}
double _ZeroUz (int i, int x, int y, int z, unsigned char Map)
{
 return 0;
}
double _ZeroP (int i, int x, int y, int z, unsigned char Map)
{
 return 0;
}
double _ZeroKt (int i, int x, int y, int z, unsigned char Map)
{
 double Result;
 double KScale = 0.208816;
 double Z0 = 0.05;
 double ZZ = z==0 && _Z[0]==0.0 ? Z0 : (_Z[z]+_Z[z-1])*0.05+Z0;
 
 Result = KScale*KScale*(1.09680767+(-0.008122677347+(0.0002099331724+(-0.000002208352551+(0.00000000906293221+(-9.749368416E-12+(-1.187063604E-14)*ZZ)*ZZ)*ZZ)*ZZ)*ZZ)*ZZ);

 return Result;
}
double _ZeroEt (int i, int x, int y, int z, unsigned char Map)
{
 double Result;
 double KScale = 0.208816;
 double Z0 = 0.05;
 double ZZ = z==0 && _Z[0]==0.0 ? Z0 : (_Z[z]+_Z[z-1])*0.05+Z0;
 
 Result = KScale*KScale*KScale*0.3644/ZZ;

 return Result;
}
double _ZeroFUV (int i, int x, int y, int z, unsigned char Map)
{
 return 0;
}
double _ZeroUV (int i, int x, int y, int z, unsigned char Map)
{
 return 0;
}
double _ZeroF0 (int i, int x, int y, int z, unsigned char Map)
{
 return 0;
}
double _ZeroI (int i, int x, int y, int z, unsigned char Map)
{
 return 0;
}
double _ZeroFIR (int i, int x, int y, int z, unsigned char Map)
{
 return 0;
}
double _ZeroIR (int i, int x, int y, int z, unsigned char Map)
{
 return 0;
}
double _ZeroIRt (int i, int x, int y, int z, unsigned char Map)
{
 return 0;
}
double _ZeroRok (int i, int x, int y, int z, unsigned char Map)
{
 return InitRok[i];
}
double _ZeroUkx (int i, int x, int y, int z, unsigned char Map)
{
 return 0;
}
double _ZeroUky (int i, int x, int y, int z, unsigned char Map)
{
 return 0;
}
double _ZeroUkz (int i, int x, int y, int z, unsigned char Map)
{
 return 0;
}
double _ZeroNk (int i, int x, int y, int z, unsigned char Map)
{
 return InitNk[i];
}
double _ZeroSO2k (int i, int x, int y, int z, unsigned char Map)
{
 return InitSO2k[i];
}
#endif
