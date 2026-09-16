#ifndef __SOLVERS__
#define __SOLVERS__

#include "area3d.h"
#include "KINETIC.H"

extern int UseOpenMP;

typedef int (*_Solver) (
                        float ** Bounds,
                        int Iteration, float * H, float * K, float * S, float * H1,
                        unsigned char * Map, unsigned char * Area, unsigned char * Boundaries,
                        DescLine * DescH, DescLine * DescP, DescLine * DescV,
                        float * _HX,  float * _HY,  float * _HZ,
                        float * _HXX, float * _HYY, float * _HZZ,
                        float * _HXS, float * _HYS, float * _HZS,
                        double * L, double * M, double * G,
                        unsigned int Flags,
                        char Projection
);

void CorrectBounds (float ** Bounds,
                    float * H, int NumZ, int NumY, int NumX, int ZStep, int YStep, int XStep,
                    DescLine * Desc, unsigned char * Map, unsigned char * Area,
                    unsigned char Projection, unsigned char Direction,
                    int UseOpenMP);

#define _XP ((x+1)      % NumX)
#define _XM ((x+NumX-1) % NumX)

#define _ZYX  (z*ZStep+y*YStep+x*XStep)

#define _ZYXP (z*ZStep+y*YStep+_XP*XStep)
#define _ZYXM (z*ZStep+y*YStep+_XM*XStep)

#ifdef __PARALLEL__
#define HX_XM (x==0 ? HXbf : HX[x-1])
#else
#define HX_XM HX[(x+NumX-1) % NumX]
#endif

#define QRAN_SHIFT 15
#define QRAN_MASK ((1 << QRAN_SHIFT) - 1)
#define QRAN_MAX QRAN_MASK
#define QRAN_A 1664525
#define QRAN_C 1013904223

#define omp_RAND_MAX (omp_in_parallel() ? QRAN_MAX : RAND_MAX)

unsigned int srands[1024] = { 0 };

void omp_srand(unsigned int seed) {
	int i;
	for (i = 0; i < sizeof(srands)/sizeof(srands[0]); i++)
		srands[i] = (seed + i) % (QRAN_MAX+1);
	srand(seed);
}

unsigned int omp_rand() {
	if (omp_in_parallel()) {
		int id = omp_get_thread_num();

		srands[id] = ((long long) QRAN_A * srands[id] + QRAN_C) & 0xFFFFFFFF;

		return (srands[id] >> 16) & QRAN_MAX;
	} else
		return rand();
}

int DirectSolver (float ** Bounds,
                  int Iteration, float * H, float * K, float * S, float * H1,
                  unsigned char * Map, unsigned char * Area, unsigned char * Boundaries,
                  DescLine * DescH, DescLine * DescP, DescLine * DescV,
                  float * HX,  float * HY,  float * HZ,
                  float * HXX, float * HYY, float * HZZ,
                  float * HXS, float * HYS, float * HZS,
                  double * L, double * M, double * G,
                  unsigned int Flags,
                  char Projection)
{
 DeclareStandard
 int zy;

 #pragma omp parallel if(UseOpenMP)
 #pragma omp for schedule(dynamic,imax(1,NY*NZ/(2*nSMP))) private(zy,z,y,x,Ptr)
 for (zy=0; zy<NZ*NY; zy++)
    {
     z = zy/NY;
     y = zy%NY;
     for (x=0, Ptr=NX*zy; x<NX; Ptr++,x++)
         if (Map[Ptr]==Zero)
            H[Ptr] = K==NULL ? 0.0f : K[Ptr];
    }

 return 1;
}

int EmptySolver (float ** Bounds,
                  int Iteration, float * H, float * K, float * S, float * H1,
                  unsigned char * Map, unsigned char * Area, unsigned char * Boundaries,
                  DescLine * DescH, DescLine * DescP, DescLine * DescV,
                  float * HX,  float * HY,  float * HZ,
                  float * HXX, float * HYY, float * HZZ,
                  float * HXS, float * HYS, float * HZS,
                  double * L, double * M, double * G,
                  unsigned int Flags,
                  char Projection)
{
 return 1;
}

int RozhkovSolver (float ** Bounds,
                   int Iteration, float * H, float * K, float * S, float * H1,
                   unsigned char * Map, unsigned char * Area, unsigned char * Boundaries,
                   DescLine * DescH, DescLine * DescP, DescLine * DescV,
                   float * HX,  float * HY,  float * HZ,
                   float * HXX, float * HYY, float * HZZ,
                   float * HXS, float * HYS, float * HZS,
                   double * L, double * M, double * G,
                   unsigned int Flags,
                   char Projection)
{
 DeclareStandard
 int zy;

 #pragma omp parallel if(UseOpenMP)
 #pragma omp for schedule(dynamic,imax(1,NY*NZ/(2*nSMP))) private(zy,z,y,x,Ptr)
 for (zy=0; zy<NZ*NY; zy++)
    {
     z = zy/NY;
     y = zy%NY;
     for (x=0, Ptr=NX*zy; x<NX; Ptr++,x++)
         if (Map[Ptr]==Zero)
            H[Ptr] = (float)((H[Ptr]+TAU*(K==NULL ? 0.0f : K[Ptr]))/(1.0-TAU*(S==NULL ? 0.0f : S[Ptr])));
    }

 return 1;
}

#define MaxNPoissons 256

int PoissonSolver (float ** Bounds,
                  int Iteration, float * H, float * K, float * S, float * H1,
                  unsigned char * Map, unsigned char * Area, unsigned char * Boundaries,
                  DescLine * DescH, DescLine * DescP, DescLine * DescV,
                  float * HX,  float * HY,  float * HZ,
                  float * HXX, float * HYY, float * HZZ,
                  float * HXS, float * HYS, float * HZS,
                  double * L, double * M, double * G,
                  unsigned int Flags,
                  char Projection)
{
 /* Make global. Create HMax, Delta for each proc & so on. Fill by zeros... */
 static double Theta[MaxNPoissons];
 static double r2[MaxNPoissons];
 static int Pcnt = -1;
 char BND = Flags>>shBoundary;

 int x,y,z,zy;
 int Ptr = 0;
 double Pi = 3.1415926535897932;
 double lMax[256]   = {0.0,0.0,0.0,0.0};
 double lDelta[256] = {0.0,0.0,0.0,0.0};
 double Delta = 0.0;
 double HMax = 0.0;

 for (x=0; x<nSMP; x++)
     {
      lMax[x] = 0.0;
      lDelta[x] = 0.0;
     }

 if (Pcnt<0)
    {
     for (x=0; x<nPoissons; x++)
         {
          Theta[x] = -1.0;
          r2[x] = -1.0;
         }
     Pcnt = nPoissons-1;
    }

 if (Iteration==1) Pcnt = (Pcnt+1) % nPoissons;

 if (r2[Pcnt]>=0.0) /* Chebyshev acceleration */
    Theta[Pcnt] = Iteration==1 ? 1.0 : Iteration==2 ? 2.0*(r2[Pcnt]-2.0)/(3.0*r2[Pcnt]-4.0) : (r2[Pcnt]*Theta[Pcnt]-4.0)/(r2[Pcnt]+r2[Pcnt]*Theta[Pcnt]-4.0);
 else
    {
     if (Theta[Pcnt]<0.0)
        { /* Search Optimal Theta (with recursive calls) => r2 */
         int BestIters = 100000;
         double BestTheta = 1.0;
         double Step = 0.1;
         int CanBack = 0;
         Theta[Pcnt] = 1.0;
         do {
             int GO;

             Iteration = 2; /* To prevent Pcnt incrementation */
             memset(H1,0,BoardSize);
             do {
                 GO = PoissonSolver(Bounds,Iteration++,H1,NULL,(S!=NULL ? H1 : NULL),NULL,Map,Area,Boundaries,DescH,DescP,DescV,HX,HY,HZ,HXX,HYY,HZZ,HXS,HYS,HZS,L,M,G,Flags,Projection);
                 CorrectBounds(Bounds,H1, NZ,NX,NY, NY*NX,1,NX, DescP, Map, Area, Projection, 1, UseOpenMP);
                 CorrectBounds(Bounds,H1, NY,NX,NZ, NX,1,NY*NX, DescV, Map, Area, Projection, 2, UseOpenMP);
                 CorrectBounds(Bounds,H1, NZ,NY,NX, NY*NX,NX,1, DescH, Map, Area, Projection, 0, UseOpenMP);
             } while (!GO);
             if (GO>=0 && Iteration<=BestIters)
                {
                 BestIters = Iteration;
                 BestTheta = Theta[Pcnt];
                 CanBack = 1;
                }
             else
                {
                 if (CanBack) Step *= -0.1;
                 BestIters = 100000;
                }
             Theta[Pcnt] += Step;
         } while (Theta[Pcnt]>=1.0 && Theta[Pcnt]<2.0 && fabs(Step)>=0.01);
         if (BestTheta<1.0) BestTheta = 1.0;
         if (BestTheta>1.99) BestTheta = 1.99;

         r2[Pcnt] = 1.0-(2.0/BestTheta-1.0)*(2.0/BestTheta-1.0);

         Theta[Pcnt] = 1.0;
        }
    }

 #pragma omp parallel if(UseOpenMP) private(zy,z,y,x,Ptr) shared(Theta,Pcnt,lMax,lDelta)
 {
  #ifdef __OPENMP__
    int      ThreadID = omp_get_thread_num();
  #else
    int      ThreadID = 0;
  #endif
  #pragma omp for schedule(dynamic,imax(1,NY*NZ/imax(20,4*nSMP)))
  for (zy=0; zy<NY*NZ; zy++) /* Odd loop */
    {
     z = zy/NY;
     y = zy%NY;
     for (x=(z+y)&1, Ptr = zy*NX+((z+y)&1); x<NX; Ptr+=2,x+=2)
       if (Map[Ptr]==Zero && !(Boundaries[Ptr] && Boundaries[Ptr]<=BND))
         { /* Non-equal steps */
          double _HXP = *VAL(Bounds,_right,&H[ZYXP]);
          double _HXM = *VAL(Bounds,_left,&H[ZYXM]);
          double _HYP = *VAL(Bounds,_forw,&H[ZYPX]);
          double _HYM = *VAL(Bounds,_back,&H[ZYMX]);
          double _HZP = *VAL(Bounds,_top,&H[ZPYX]);
          double _HZM = *VAL(Bounds,_bottom,&H[ZMYX]);
          double HHH = H[Ptr]==FictZero ? 0.0 : H[Ptr];
          double NewH = H[Ptr]==FictZero ? 0.0 :
                        (1-Theta[Pcnt])*HHH+
                        Theta[Pcnt]/
                             (2.0/HXX[x]+2.0/HYY[y]+2.0/HZ[z]/HZ_ZM+
                              (H==S ? 0.01*sin(12.0*Pi*zy/(NY*NZ))*cos(12.0*Pi*x/NX) : S==NULL ? 0.0 : S[Ptr]))*
                          (
                           2.0/HXS[x]*(_HXP/HX[x]+_HXM/HX[(x+NXs)%NX])+
                           2.0/HYS[y]*(_HYP/HY[y]+_HYM/HY[(y+NYs)%NY])+
                           2.0/HZS[z]*(_HZP/HZ[z]+_HZM/HZ_ZM)-
                           (H1==NULL ? -sin(12.0*Pi*zy/(NY*NZ))*cos(12.0*Pi*x/NX) : K==NULL ? 0.0 : K[Ptr])
                          );
          double delta = fabs(NewH-HHH);
          double fH = fabs(HHH);

          if (fH>lMax[ThreadID]) lMax[ThreadID] = fH;
          if (delta>lDelta[ThreadID]) lDelta[ThreadID] = delta;
          H[Ptr] = NewH;
         }
    }
 }
 #pragma omp parallel if(UseOpenMP) private(zy,z,y,x,Ptr) shared(Theta,Pcnt,lMax,lDelta)
 {
  #ifdef __OPENMP__
    int      ThreadID = omp_get_thread_num();
  #else
    int      ThreadID = 0;
  #endif
  #pragma omp for schedule(dynamic,imax(1,NY*NZ/imax(20,4*nSMP)))
  for (zy=0; zy<NY*NZ; zy++) /* Odd loop */
    {
     z = zy/NY;
     y = zy%NY;
     for (x=(z+y+1)&1, Ptr = zy*NX+((z+y+1)&1); x<NX; Ptr+=2,x+=2)
        if (Map[Ptr]==Zero && !(Boundaries[Ptr] && Boundaries[Ptr]<=BND))
          { /* Non-equal steps */
          double _HXP = *VAL(Bounds,_right,&H[ZYXP]);
          double _HXM = *VAL(Bounds,_left,&H[ZYXM]);
          double _HYP = *VAL(Bounds,_forw,&H[ZYPX]);
          double _HYM = *VAL(Bounds,_back,&H[ZYMX]);
          double _HZP = *VAL(Bounds,_top,&H[ZPYX]);
          double _HZM = *VAL(Bounds,_bottom,&H[ZMYX]);
          double HHH = H[Ptr]==FictZero ? 0.0 : H[Ptr];
          double NewH = H[Ptr]==FictZero ? 0.0 :
                        (1-Theta[Pcnt])*HHH+
                        Theta[Pcnt]/
                             (2.0/HXX[x]+2.0/HYY[y]+2.0/HZ[z]/HZ_ZM+
                              (H==S ? 0.01*sin(12.0*Pi*zy/(NY*NZ))*cos(12.0*Pi*x/NX) : S==NULL ? 0.0 : S[Ptr]))*
                          (
                           2.0/HXS[x]*(_HXP/HX[x]+_HXM/HX[(x+NXs)%NX])+
                           2.0/HYS[y]*(_HYP/HY[y]+_HYM/HY[(y+NYs)%NY])+
                           2.0/HZS[z]*(_HZP/HZ[z]+_HZM/HZ_ZM)-
                           (H1==NULL ? -sin(12.0*Pi*zy/(NY*NZ))*cos(12.0*Pi*x/NX) : K==NULL ? 0.0 : K[Ptr])
                          );
           double delta = fabs(NewH-HHH);
           double fH = fabs(HHH);

           if (fH>lMax[ThreadID]) lMax[ThreadID] = fH;
           if (delta>lDelta[ThreadID]) lDelta[ThreadID] = delta;
           H[Ptr] = NewH;
          }
    }
 }
 HMax  = lMax[0];
 Delta = lDelta[0];
 #ifdef __OPENMP__
 for (x=1; x<nSMP; x++)
     {
      if (lMax[x]>HMax) HMax = lMax[x];
      if (lDelta[x]>Delta) Delta = lDelta[x];
     }
 #endif
 if (HMax<1E-8) HMax = 1.0;
// printf("%lf:%lf=%lf ",Delta,HMax,Delta/HMax);
 Delta /= HMax;
 return HMax > 10000 ? -1 : (H1==NULL ? Delta<=1E-2 : Delta<=1E-3);
}

 /* Решаем уравнение Пуассона для поля критичности. lapl(GRAV) = ERRS */
void SolveCheckPoisson(float * GRAV, float * ERRS) {
 int  x, y, z;
 long Ptr;
 char GO;
 do {
    const double Theta = 1.83;
    const double H = 1.0;
    const double H2 = H*H;
    double lMax[256]   = {0.0,0.0,0.0,0.0};
    double lDelta[256] = {0.0,0.0,0.0,0.0};
    double Delta = 0.0;
    double HMax = 0.0;
    #pragma omp parallel if (UseOpenMP) private(z,y,x,Ptr) shared(lMax, lDelta) 
    {
     #ifdef __OPENMP__
     int ThreadID = omp_get_thread_num();
     #else
     int ThreadID = 0;
     #endif
     int zy;
     #pragma omp for schedule(guided)
     for (zy=0; zy<NY*NZ; zy++) /* Odd loop */
       {
        z = zy/NY;
        y = zy%NY;
        for (x=(z+y)&1, Ptr = zy*NX+((z+y)&1); x<NX; Ptr+=2,x+=2)
          if (x > 0 && x < NXs && y > 0 && y < NYs && z > 0 && z < NZs)
            { /* Non-equal steps */
             double _HXP = GRAV[ZYXP];
             double _HXM = GRAV[ZYXM];
             double _HYP = GRAV[ZYPX];
             double _HYM = GRAV[ZYMX];
             double _HZP = GRAV[ZPYX];
             double _HZM = GRAV[ZMYX];
             double HHH = GRAV[Ptr];
             double NewH = (1-Theta)*HHH+
                           Theta/(2.0+2.0+2.0)*
                             (
                              _HXP+_HXM+
                              _HYP+_HYM+
                              _HZP+_HZM+
                              H2*ERRS[Ptr]
                             );
             double delta = fabs(NewH-HHH);
             double fH = fabs(HHH);

             if (fH>lMax[ThreadID]) lMax[ThreadID] = fH;
             if (delta>lDelta[ThreadID]) lDelta[ThreadID] = delta;
             GRAV[Ptr] = NewH;
            }
       }
    }
    #pragma omp parallel if(UseOpenMP) private(z,y,x,Ptr) shared(lMax,lDelta)
    {
     #ifdef __OPENMP__
     int ThreadID = omp_get_thread_num();
     #else
     int ThreadID = 0;
     #endif
     int zy;
     #pragma omp for schedule(guided)
     for (zy=0; zy<NY*NZ; zy++) /* Odd loop */
       {
        z = zy/NY;
        y = zy%NY;
        for (x=(z+y+1)&1, Ptr = zy*NX+((z+y+1)&1); x<NX; Ptr+=2,x+=2)
           if (x > 0 && x < NXs && y > 0 && y < NYs && z > 0 && z < NZs)
             { /* Non-equal steps */
              double _HXP = GRAV[ZYXP];
              double _HXM = GRAV[ZYXM];
              double _HYP = GRAV[ZYPX];
              double _HYM = GRAV[ZYMX];
              double _HZP = GRAV[ZPYX];
              double _HZM = GRAV[ZMYX];
              double HHH = GRAV[Ptr];
              double NewH = (1-Theta)*HHH+
                           Theta/(2.0+2.0+2.0)*
                             (
                              _HXP+_HXM+
                              _HYP+_HYM+
                              _HZP+_HZM+
                              H2*ERRS[Ptr]
                             );
              double delta = fabs(NewH-HHH);
              double fH = fabs(HHH);

              if (fH>lMax[ThreadID]) lMax[ThreadID] = fH;
              if (delta>lDelta[ThreadID]) lDelta[ThreadID] = delta;
              GRAV[Ptr] = NewH;
             }
       }
    }
    HMax  = lMax[0];
    Delta = lDelta[0];
#ifdef __OPENMP__
    for (x=1; x<nSMP; x++)
        {
         if (lMax[x]>HMax) HMax = lMax[x];
         if (lDelta[x]>Delta) Delta = lDelta[x];
        }
#endif
    if (HMax<1E-8) HMax = 1.0;
    // fprintf(StdOutput,"%lf:%lf=%lf ",Delta,HMax,Delta/HMax); fflush(StdOutput);
    Delta /= HMax;
    GO = Delta>1E-4;
 } while (GO);
}

#define nRays 5

int SolarSolver (float ** Bounds,
                 int Iteration, float * H, float * K, float * S, float * H1,
                 unsigned char * Map, unsigned char * Area, unsigned char * Boundaries,
                 DescLine * DescH, DescLine * DescP, DescLine * DescV,
                 float * HX,  float * HY,  float * HZ,
                 float * HXX, float * HYY, float * HZZ,
                 float * HXS, float * HYS, float * HZS,
                 double * L, double * M, double * G,
                 unsigned int Flags,
                 char Projection)
{
 DeclareStandard
 int zy;
 double x01,x02,y01,y02,z01,z02;
 double dx = -SolarX;
 double dy = -SolarY;
 double dz = -SolarZ;
 char zrx = fabs(dx)<1E-3;
 char zry = fabs(dy)<1E-3;
 char zrz = fabs(dz)<1E-3;
 char px = !zrx && dx>0.0;
 char py = !zry && dy>0.0;
 char pz = !zrz && dz>0.0;

 omp_srand(184415/*(unsigned int) &SolarSolver*/);
 
 #pragma omp for schedule(dynamic,imax(1,NY*NZ/imax(20,4*nSMP))) private(zy,z,y,x,Ptr,x01,x02,y01,y02,z01,z02)
 for (zy=0; zy<NY*NZ; zy++)
    {
     int j;

     z = zy/NY;
     z01=0.0;
     z02=HZ[0];
     for (j=0; j<z; j++, z01=z02, z02 += j<NZ ? HZ[j] : 0.0);
     
     y = zy%NY;
     y01=0.0;
     y02=HY[0];
     for (j=0; j<y; j++, y01=y02, y02 += j<NY ? HY[j] : 0.0);
     Ptr=zy*NX;

     for (x=0, x01=0.0, x02=HX[0]; x<NX; Ptr++, x++, x01=x02, x02 += x<NX ? HX[x] : 0.0)
         if (Map[Ptr]==Zero ||
             Map[Ptr]==BoundClosed &&
               (IsTop && SolarZ>0.0 || IsBottom && SolarZ<0.0 ||
                IsLeft && SolarX<0.0 || IsRight && SolarX>0.0 ||
                IsForw && SolarY>0.0 || IsBack && SolarY<0.0
               )
            )
            {
             double kx = IsLeft ? 0.99 : IsRight ? 0.01 : 0.0;
             double ky = IsBack ? 0.99 : IsForw ? 0.01 : 0.0;
             double kz = IsBottom ? 0.99 : IsTop ? 0.01 : 0.0;
             int Rays;

             H[Ptr] = 0.0;
             for (Rays=0; Rays<nRays; Rays++)
               {
                int XX = x, YY = y, ZZ = z;
                double X1 = x01, X2 = x02, Y1 = y01, Y2 = y02, Z1 = z01, Z2 = z02;
                double xx = X1+kx*(X2-X1)+(1-kx)*((double) omp_rand()/omp_RAND_MAX);
                double yy = Y1+ky*(Y2-Y1)+(1-ky)*((double) omp_rand()/omp_RAND_MAX);
                double zz = Z1+kz*(Z2-Z1)+(1-kz)*((double) omp_rand()/omp_RAND_MAX);
                double F0 = 0.0;
                int    PTR = Ptr;
                int    Ptrs[512] = {PTR};
                double Ts[512];
                int    NP = 1;
                char   STOP = 0;

                do {
                  double tx = zrx ? 1E8 : px ? (X2-xx)/dx : (X1-xx)/dx;
                  double ty = zry ? 1E8 : py ? (Y2-yy)/dy : (Y1-yy)/dy;
                  double tz = zrz ? 1E8 : pz ? (Z2-zz)/dz : (Z1-zz)/dz;
                  double tr = amin3(tx,ty,tz);
                  int    Dirs[3];
                  char   cntDirs = 0;

                  if (fabs(tr-tx)<1E-15)
                     if (px)
                        {
                         PTR++; X1 = X2;
                         if (++XX<NX) X2 += HX[XX];
                         else STOP = 1;
                         Dirs[cntDirs++] = _right;
                        }
                     else
                        {
                         PTR--; X2 = X1;
                         if (--XX>=0) X1 -= HX[XX];
                         else STOP = 1;
                         Dirs[cntDirs++] = _left;
                        }
                  if (!STOP && fabs(tr-ty)<1E-15)
                     if (py)
                        {
                         PTR+=NX; Y1 = Y2;
                         if (++YY<NY) Y2 += HY[YY];
                         else STOP = 1;
                         Dirs[cntDirs++] = _forw;
                        }
                     else
                        {
                         PTR-=NX; Y2 = Y1;
                         if (--YY>=0) Y1 -= HY[YY];
                         else STOP = 1;
                         Dirs[cntDirs++] = _back;
                        }
                  if (!STOP && fabs(tr-tz)<1E-15)
                     if (pz)
                        {
                         PTR+=NX*NY; Z1 = Z2;
                         if (++ZZ<NZ) Z2 += HZ[ZZ];
                         else STOP = 1;
                         Dirs[cntDirs++] = _top;
                        }
                     else
                        {
                         PTR-=NX*NY; Z2 = Z1;
                         if (--ZZ>=0) Z1 -= HZ[ZZ];
                         else STOP = 1;
                         Dirs[cntDirs++] = _bottom;
                        }
                  Ts[NP-1] = tr;
                  if (!STOP)
                     if (Map[PTR]!=Zero)
                        {
                         if (Map[PTR]==BoundClosed || cntDirs==0)
                            F0 = 0.0;
                         else
                            {
                             int k;

                             for (k=0, F0=0.0; k<cntDirs; k++)
                                 F0 += *VAL(Bounds,Dirs[k],&H[PTR]);
                             F0 /= cntDirs;
                            }
                         STOP = 1;
                        }
                     else
                        {
                         xx += tr*dx; yy += tr*dy; zz += tr*dz;
                         Ptrs[NP++] = PTR;
                        }
                  else
                     F0 = S==NULL ? 0.0 : S[Ptr];
                } while(!STOP);
                if (F0>1E-5)
                   {
                    double OptPath = 0.0;

                    if (K)
                       for (; NP>0; NP--)
                         OptPath += K[Ptrs[NP-1]]*Ts[NP-1];
                    H[Ptr] += F0*exp(-OptPath);
                   }
               }
             H[Ptr] /= nRays;
            }
    }

 return 1;
}

#endif
