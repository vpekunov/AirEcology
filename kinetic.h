#ifndef __KINETIC__
#define __KINETIC__

#include <stdlib.h>
#include <math.h>

/* Максимальное число веществ, участвующих в реакции */
#define MaxActSubst 40
/* Максимальная длина имени вещества */
#define SubstNameLength 21

#define Adams_Rozhkov_Method 0
#define Rosenbrock_Method 1

extern int nSMP;

typedef struct {
  int    NumLeft;
  int    NumRight;
  int    Indexes[MaxActSubst];
  double k[MaxActSubst];
} Reaction;

typedef struct {
 int      LossH;    /* Признак малости H */
 double   LastH;    /* Последнее H */
 double   ReachTau; /* Достигнутое значение TAU */
 float    Tmp;      /* Температура проведения реакций */
} ErrorKinetic;

typedef double Derivs[7];

typedef char   SubstName[SubstNameLength];

#ifndef M_LOG10E
#define M_LOG10E 0.434294481903251827651
#endif

/* Constants */
/* Константы STIFF */
double dfltzr = 1.0E-18;
int    maxite = 3;
int    maxfai = 3;
double rmxini = 1.0E+4;
double rmxnor = 10.0;
double rmxfai = 2.0;
int    idelay = 10;
double rhcorr = 0.25;
double rherr3 = 0.1;
double thrshl = 0.0;
double rctest = 0.3;
double bias1  = 1.3;
double bias2  = 1.2;
double bias3  = 1.4;

double ZeroK = 273.0;
double R     = 1.987E-3;

double Koeffs[5][3] =
  {
   {1.0,2.0,3.0},
   {1.0,4.5,6.0},
   {0.5,7.3333333,9.1666667},
   {0.1666667,10.4166667,12.5},
   {0.04166667,13.7,1.0}
  };

int MaxIterations = 30000; /* Максимальное число итераций */

double InitTime = 0.0; /* Начальное время интегрирования */
double CalcEps = 1E-2; /* Относительная точность интегрирования ST1/eps */
double MinH  = 1.0E-18;/* Минимальный шаг интегрирования ST1/hmin */
/* End of Constants */

double DefaultT = 27.0;

/* Shared variables */
double InitH = 1E-7;  /* Начальный шаг интегрирования h */
/* Имена веществ [NSubst] */
SubstName Names[MaxActSubst];
/* Скорости витания веществ [NSubst] */
double    Uw[MaxActSubst];
/* Карта участия веществ в реакциях */
char      Activity[MaxActSubst];
/* Индексы веществ, участвующих в химических реакциях */
int       TranMap[MaxActSubst];
/* Ei : Энергии активации KINET/ea [NReact] */
double * EA;
/* Ai : Lg( предэкспоненты констант скорости ) [NReact] */
double * _A;
/* Ni : Показатели степеней в температурных множителях KINET/tn [NReact] */
double * Tn;
/* Реакции KINET/lr [NReact] */
Reaction * LR;

int NReact;  /* Число реакций KINET/m */
int NSubst;  /* Число веществ KINET/n = ST1/n */
int NASubst; /* Число веществ, участвующих в реакции */
/* End of shared variables */

typedef struct {
  ErrorKinetic KinErrorInfo; /* Блок информации по ошибке */

  int LossPrecision;        /* Флаг потери точности при вычислениях */
  int Iters; /* Действительное число итераций */

  float  Tk;   /* Температура в Кельвинах KINET/tk */
  double _EndTime; /* Конечное время интегрирования */

  /* Начальная концентрация веществ [NASubst] */
  float * Conc0;
  /* Конечная концентрация [NASubst] */
  float * Conc1;

  /* Рабочие переменные [NASubst] */
  Derivs * Y;
  /* KINET */
  int    LgAScale; /* ml = Lg(aml) */
  double AScale; /* Масштабный множитель aml */
  double * C;    /* Вычисляемые константы скоростей реакций c [NReact] */
  double Tem;    /* Температура в ккал/моль tem */
  /* ST1 */
  double Time;
  double H;
  int    MaxDer;
  int    jStart;
  /* ST2 */
  double * Ymax; /* [NASubst] */
  /* ST3 */
  double * Err; /* [NASubst] */
  /* ST4 */
  double * Df; /* [NASubst][NASubst] */
  /* ST5 */
  double * FSave1; /* [NASubst] */
  double * FSave2; /* [NASubst] */
  /* STIFF */
  int    * iRow; /* [NASubst] */

  double El[6];
  double Tq[4];
  double edn;
  double e;
  double eup;
  double bnd;

  int   Nq;
  int   l;
  int   lmax;
  int   idoub;
  double rmax;
  double rc;

  double rh;

  double trend;
  double oldl0;
  double hold;

  int   evalja;
  int   conver;
} KineticContext;

KineticContext * Contexts;

void SetTk(KineticContext * C, float Tmp)
{
 C->Tk = (float) (ZeroK+Tmp);
}

double amax2(double First, double Second)
{
 return First>Second ? First : Second;
}

double amax3(double First, double Second, double Third)
{
 return First>Second ? amax2(First,Third) : amax2(Second,Third);
}

double amin2(double First, double Second)
{
 return First<Second ? First : Second;
}

double amin3(double First, double Second, double Third)
{
 return First<Second ? amin2(First,Third) : amin2(Second,Third);
}

int AllocateKineticByReactions()
{
 int doubleNR = NReact*sizeof(double);
 int C = 1;
 int i;

 memset(Activity,0,sizeof(Activity));
 /* Ei : Энергии активации KINET/ea [NReact] */
 EA = (double *) malloc(doubleNR);
 /* Ai : Lg( предэкспоненты констант скорости ) [NReact] */
 _A  = (double *) malloc(doubleNR);
 /* Ni : Показатели степеней в температурных множителях KINET/tn [NReact] */
 Tn = (double *) malloc(doubleNR);
 /* Реакции [NReact] */
 LR = (Reaction *) malloc(NReact*sizeof(Reaction));

 for (i=0; i<nSMP; i++)
     C = C && (Contexts[i].C = (double *) malloc(doubleNR));

 return (EA && _A && Tn && LR && C);
}

/* Выводим из рассмотрения вещества, фактически не участвующие в реакциях */
void ReTranslate()
{
 char Map[MaxActSubst];
 char Delta = 0;
 int  i,j;

 NASubst = 0;
 for (i=0; i<NSubst; i++)
     if (Activity[i])
        {
         Map[i] = Delta;
         TranMap[NASubst++] = i;
        }
     else
        Delta++;
 for (i=0; i<NReact; i++)
     for (j=0; j<LR[i].NumLeft+LR[i].NumRight; j++)
         LR[i].Indexes[j] -= Map[LR[i].Indexes[j]];
}

int AllocateKineticBySubsts()
{
 int floatNS  = NASubst*sizeof(float);
 int doubleNS = NASubst*sizeof(double);
 int C = 1;
 int i;

 for (i=0; i<nSMP; i++)
    {
     C = C && (Contexts[i].Conc0 = (float *) malloc(floatNS));
     C = C && ( Contexts[i].Conc1 = (float *) malloc(floatNS));
     C = C && ( Contexts[i].Y     = (Derivs *) malloc(NASubst*sizeof(Derivs)));
     C = C && ( Contexts[i].Ymax  = (double *) malloc(doubleNS));
     C = C && ( Contexts[i].Err   = (double *) malloc(doubleNS));
     C = C && ( Contexts[i].Df    = (double *) malloc(NASubst*doubleNS));
     C = C && ( Contexts[i].FSave1 = (double *) malloc(doubleNS));
     C = C && ( Contexts[i].FSave2 = (double *) malloc(doubleNS));
     C = C && ( Contexts[i].iRow   = (int *) malloc(NASubst*sizeof(int)));
    }

 return C;
}

void FreeKinetic()
{
 int i;

 free(EA);
 free(_A);
 free(Tn);
 free(LR);

 for (i=0; i<nSMP; i++)
     {
      free(Contexts[i].Conc0);
      free(Contexts[i].Conc1);
      free(Contexts[i].Y);
      free(Contexts[i].C);
      free(Contexts[i].Ymax);
      free(Contexts[i].Err);
      free(Contexts[i].Df);
      free(Contexts[i].FSave1);
      free(Contexts[i].FSave2);
      free(Contexts[i].iRow);
     }
}

int Translate(int NumEq, char * Equation)
{
 int Count = 0;

 int Long  = strlen(Equation);

 int Left  = 0;
 int All   = 0;

 if (!Long || Equation[Long-1]=='+' || Equation[Long-1]=='-')
    return 0;

 while (Count<Long)
   {
    int    BegCount = Count;
    double How  = 0.0;
    int    iHow = 1;

    int CurLong;
    int End;

    int Found = 0;
    int Index = 0;

    if (Equation[Count]=='-') Count++;
    while (Count<Long && Equation[Count]>='0' && Equation[Count]<='9') Count++;
    if (Count<Long && Equation[Count]=='.') Count++;
    while (Count<Long && Equation[Count]>='0' && Equation[Count]<='9') Count++;

    if (Count>BegCount)
       {
        char Buf = Equation[Count];

        Equation[Count] = 0;
        How = atof(&Equation[BegCount]);
        Equation[Count] = Buf;
        if (How==0.0) return 0;
        if (floor(How)==How && How>0.0)
           {
            iHow = (int)How;
            How = 0.0;
           }
       }

    if (All+iHow>MaxActSubst) return 0;

    End = Count;

    while (Equation[End]!='+' && Equation[End]!='-' && End<Long) End++;

    CurLong = End-Count;

    while (Index<NSubst && !Found)
       {
        if (strlen(Names[Index])==CurLong)
           {
            int CurN = 0;
            int CurE = Count;

            while (CurN<CurLong && Names[Index][CurN]==Equation[CurE])
              {
               CurN++;
               CurE++;
              }
            if (CurN>=CurLong) Found = 1;
           }

        if (!Found) Index++;
       }

    if (!Found) return 0;

    while (iHow--)
      {
       LR[NumEq].k[All] = How;
       LR[NumEq].Indexes[All++] = Index;
      }
    Activity[Index] = 1;

    if (End<Long-1)
       if (Equation[End]=='-')
          if (Equation[++End]=='>' && !Left)
             Left = All;
          else return 0;
    Count = End+1;
   }
 LR[NumEq].NumLeft  = Left;
 LR[NumEq].NumRight = All-Left;

 return (Left);
}

/* ВЫЧИСЛЕНИЕ ЗНАЧЕНИЙ ЯКОБИАНА ИНТЕГРИРУЕМОЙ СИСТЕМЫ ОДУ ПО */
/* ЗНАЧЕНИЯМ ВРЕМЕНИ t И РЕШЕНИЙ y, СООТВЕТСТВУЮЩИХ ЭТОМУ МОМЕНТУ */
/* ВРЕМЕНИ. ЭЛЕМЕНТЫ ЯКОБИАНА ЗАСЫЛАЮТСЯ В МАССИВ Df[50,50) */
/* PEDERV */
void Jacobian(KineticContext * C)
{
 int i,j,k;

 int PtrIJ = 0;

 for (i=0;i<NASubst;i++)
     for (j=0;j<NASubst;j++,PtrIJ++)
         C->Df[PtrIJ]=0.0;
/* ЦИКЛ ПО ЧИСЛУ УРАВНЕНИЙ */
 for (i=0;i<NReact;i++)
     {
      int Left  = LR[i].NumLeft;
      int Right = LR[i].NumRight;

      for (j=0;j<Left;j++)
          {
           int Cur = LR[i].Indexes[j];

           double W  = C->C[i];

           if (LR[i].k[j]!=0.0) W*=LR[i].k[j]*pow(C->Y[Cur][0],LR[i].k[j]-1.0);
           for (k=0;k<Left;k++)
               if (k!=j)
                  if (LR[i].k[k]==0.0)
                     W*=C->Y[LR[i].Indexes[k]][0];
                  else
                     W*=pow(C->Y[LR[i].Indexes[k]][0],LR[i].k[k]);
           for (k=0;k<Left;k++)
               C->Df[NASubst*LR[i].Indexes[k]+Cur] -= W;
           for (k=0;k<Right;k++)
               C->Df[NASubst*LR[i].Indexes[Left+k]+Cur] += W;
          }
     }
#if defined(__IMITATE__) && !defined(WIN32)
 Check();
#endif
}

/* ВЫЧИСЛЯЕТ ЗНАЧЕНИЯ ПРАВЫХ ЧАСТЕЙ СИСТЕМЫ ОДУ ПО ЗНАЧЕНИЯМ */
/* ВРЕМЕНИ t И РЕШЕНИЙ y0, СООТВЕТСТВУЮЩИХ ЭТОМУ МОМЕНТУ ВРЕМЕНИ, */
/* И ЗАСЫЛАЕТ ИХ В МАССИВ f */
/* DIFFUN */
void DifFun(KineticContext * CC,double * Y0,int Ystep,double * F)
{
 int i,j;

 for (i=0;i<NASubst;i++)
     {
      if (Y0[Ystep*i]<0.0) Y0[Ystep*i] = 0.0;
      F[i] = 0.0;
     }

 for (i=0;i<NReact;i++)
     {
      int Left  = LR[i].NumLeft;
      int Right = LR[i].NumRight;

      double W = CC->C[i];

      for (j=0;j<Left;j++)
          if (LR[i].k[j]==0.0)
             W*=Y0[Ystep*LR[i].Indexes[j]];
          else
             W*=pow(Y0[Ystep*LR[i].Indexes[j]],LR[i].k[j]);
      /* Здесь W содержит значение скорости соответствующей реакции.
         Делить на AScale */
      for (j=0;j<Left;j++)
          F[LR[i].Indexes[j]] -= W;
      for (j=0;j<Right;j++)
          F[LR[i].Indexes[Left+j]] += W;
     }
#if defined(__IMITATE__) && !defined(WIN32)
 Check();
#endif
}

void JacobianD(KineticContext * CC)
{
 int i,j;
 double F0[MaxActSubst];
 double F1[MaxActSubst];
 const double DX = 0.0000001;
 int PtrIJ;

 DifFun(CC,(double *) CC->Y,7,F0);
 for (i=0;i<NASubst;i++)
     {
      CC->Y[i][0] += DX;
      DifFun(CC,(double *) CC->Y,7,F1);
      CC->Y[i][0] -= DX;
      for (j=0,PtrIJ=i; j<NASubst; j++,PtrIJ+=NASubst)
          CC->Df[PtrIJ] = (F1[j]-F0[j])/DX;
     }
#if defined(__IMITATE__) && !defined(WIN32)
 Check();
#endif
}

/* ПЕРЕМАСШТАБИРОВАНИЕ РЕШЕНИЙ СИСТЕМЫ И ИХ ПРОИЗВОДНЫХ */
/* ПРИ ИЗМЕНЕНИИ ШАГА ИНТЕГРИРОВАНИЯ */
/* RESCAL */
void Rescale(KineticContext * CC)
{
 int i,j;

 double R1=1.0;

 CC->rh = amax2(CC->rh,MinH/fabs(CC->H));
 CC->rh = amin3(CC->rh,CC->_EndTime/fabs(CC->H),CC->rmax);
 for (j=1;j<CC->l;j++)
     {
      R1 *= CC->rh;
      for (i=0;i<NASubst;i++)
          CC->Y[i][j] *= R1;
     }
 CC->H *= CC->rh;
 CC->rc *= CC->rh;
 CC->idoub = CC->l+1;
#if defined(__IMITATE__) && !defined(WIN32)
 Check();
#endif
}

/* ВЫЧИСЛЕНИЕ ЗНАЧЕНИЯ ПРЕДИКТОРА */
void Predictor(KineticContext * CC, int Nq)
{
 int i,j,j1,j2;

 for (j1=0;j1<Nq;j1++)
     for (j2=j1;j2<Nq;j2++)
         {
          j = Nq-j2+j1-1;

          for (i=0;i<NASubst;i++)
              CC->Y[i][j] += CC->Y[i][j+1];
         }
#if defined(__IMITATE__) && !defined(WIN32)
 Check();
#endif
}

void RePredictor(KineticContext * CC, int Nq)
{
 int i,j,j1,j2;

 for (j1=0;j1<Nq;j1++)
     for (j2=j1;j2<Nq;j2++)
         {
          j = Nq-j2+j1-1;

          for (i=0;i<NASubst;i++)
              CC->Y[i][j] -= CC->Y[i][j+1];
         }
#if defined(__IMITATE__) && !defined(WIN32)
 Check();
#endif
}

/* ВЫЧИСЛЕНИЕ КОЭФФИЦИЕНТОВ МЕТОДА ИНТЕГРИРОВАНИЯ */
/* СООТВЕТСТВУЮЩЕГО ПОРЯДКА */
/* COSET */
void CalcKoeffs(KineticContext * CC, int Nq)
{
 int i;

 CC->El[1]=1.0;
 if (CC->MaxDer>5) CC->MaxDer = 5;
 switch (Nq) {
   case 1: CC->El[0] = 1.0;
           break;
   case 2: CC->El[0] = 6.666666666667E-01;
           CC->El[2] = 3.333333333333E-01;
           break;
   case 3: CC->El[0] = 5.454545454546E-01;
           CC->El[2] = 1.0;
           CC->El[3] = 9.090909090909E-02;
           break;
   case 4: CC->El[0] = 0.48;
           CC->El[2] = 0.7;
           CC->El[3] = 0.2;
           CC->El[4] = 0.02;
           break;
   case 5: CC->El[0] = 4.379562043796E-01;
           CC->El[2] = 8.211678832117E-01;
           CC->El[3] = 3.102189781022E-01;
           CC->El[4] = 5.474452554745E-02;
           CC->El[5] = 3.649635036496E-03;
 }
 for (i=0;i<3;i++)
     CC->Tq[i] = Koeffs[Nq-1][i];
 CC->Tq[3]=0.5/(Nq+2.0);
}

/* ОСУЩЕСТВЛЯЕТ LU-РАЗЛОЖЕНИЕ МАТРИЦЫ CC->e-beta*h*A, */
/* ГДЕ A-ЯКОБИАН СИСТЕМЫ, CC->e-ЕДИНИЧНАЯ МАТРИЦА, */
/* h-ТЕКУЩИЙ ШАГ ИНТЕГРИРОВАНИЯ, beta-КОНСТАНТА */
/* MATDEC */
void GetLU(KineticContext * CC, int * iRow)
{
 int i,j,k;

 if (NASubst<=1)
    {
     CC->Df[0] = 1.0/CC->Df[0];
     iRow[0]  = 0;
    }
 else
   {
    for (i=0;i<NASubst;i++)
        iRow[i]=i;
    for (i=0;i<NASubst-1;i++)
      {
       double Big  = 0.0;
       int    iBig = 0;

       double Kf;
       int    ir;

       for (j=i;j<NASubst;j++)
           {
            /* По непонятным причинам компилятор Borland C++ 3.1
               отказался скомпилировать следующую строку :
            double size = fabs(Df[NASubst*iRow[j]+i]);
            */
            double * Deriv = CC->Df;
            double size    = fabs(Deriv[NASubst*iRow[j]+i]);

            if (size>Big)
               {
                Big  = size;
                iBig = j;
               }
           }
       if (iBig!=i)
          {
           int V = iRow[i];
           iRow[i] = iRow[iBig];
           iRow[iBig] = V;
          }
       ir = NASubst*iRow[i];

       Kf = 1.0/CC->Df[ir+i];

       CC->Df[ir+i] = Kf;
       for (j=i+1;j<NASubst;j++)
           {
            int jr = NASubst*iRow[j];

            double Fact = Kf*CC->Df[jr+i];

            CC->Df[jr+i] = Fact;
            for (k=i+1;k<NASubst;k++)
                CC->Df[jr+k] -= Fact*CC->Df[ir+k];
           }
      }
    i = iRow[NASubst-1];

    CC->Df[NASubst*(i+1)-1] = 1.0/CC->Df[NASubst*(i+1)-1];
   }
#if defined(__IMITATE__) && !defined(WIN32)
 Check();
#endif
}

/* РЕШЕНИЕ СИСТЕМЫ ЛИНЕЙНЫХ УРАВНЕНИЙ */
void SolveLU(KineticContext * CC, int * iRow, double * Y, double * X)
{
 int i,j,k;

 if (NASubst<=1) X[0] = Y[0]*CC->Df[0];
 else
   {
    X[0] = Y[iRow[0]];
    for (i=1;i<NASubst;i++)
        {
         int   ir = iRow[i];
         double V  = Y[ir];

         for (j=0;j<i;j++)
             V -= CC->Df[NASubst*ir+j]*X[j];
         X[i] = V;
        }

    X[NASubst-1]*=CC->Df[NASubst*iRow[NASubst-1]+NASubst-1];

    for (i=1,j=NASubst-2;i<NASubst;i++,j--)
        {
         int   jr = NASubst*iRow[j];
         double V = X[j];

         for (k=j+1;k<NASubst;k++)
             V -= CC->Df[jr+k]*X[k];
         X[j] = V*CC->Df[jr+j];
        }
   }
#if defined(__IMITATE__) && !defined(WIN32)
 Check();
#endif
}

#define SetNq(CC,NewNq)          \
{                                \
 CC->Nq = NewNq;                 \
 CC->l  = CC->Nq+1;              \
 CalcKoeffs(CC,CC->Nq);          \
                                 \
 CC->lmax   = CC->MaxDer+1;      \
 CC->oldl0  = CC->El[0];         \
 CC->edn    = CC->Tq[0]*CalcEps; \
 CC->edn   *= CC->edn;           \
 CC->e      = CC->Tq[1]*CalcEps; \
 CC->e     *= CC->e;             \
 CC->eup    = CC->Tq[2]*CalcEps; \
 CC->eup   *= CC->eup;           \
 CC->bnd    = CC->Tq[3]*CalcEps; \
 CC->bnd   *= CC->bnd;           \
}

int Stiff(KineticContext * CC)
{
 int   kFlag=0;
 double Told=CC->Time;

 double pr1,pr2,pr3;
 double d;
 int newq;

 int Break590;
 int Break670;
 int ExitLoop;

 double Savel0;

 int i,j;

 if (CC->jStart==0)
    {
     DifFun(CC,(double *) CC->Y,7,CC->FSave1);
     for (i=0;i<NASubst;i++)
         CC->Y[i][1]=CC->FSave1[i]*CC->H;

     SetNq(CC,1);

     CC->idoub  = CC->l+1;
     CC->rmax   = rmxini;
     CC->trend  = 1.0;
     CC->rc     = 0.0;
     CC->hold   = CC->H;
     CC->evalja = 1;
     CC->conver = 0;
    }
 else
   if (CC->jStart<0)
      if (CC->H!=CC->hold)
         {
          CC->rh = CC->H/CC->hold;
          CC->H  = CC->hold;
          Rescale(CC);
         }

 Break590 = 0;
 Break670 = 0;
 ExitLoop = 0;
 do /* 200 */
   {
    if (fabs(CC->rc-1.0)>rctest) CC->evalja = 1;
    CC->Time+=CC->H;
    Predictor(CC,CC->Nq);

    do
      {
       int iter = 0;

       double d1;

       for (i=0;i<NASubst;i++)
           CC->Err[i] = 0.0;
       DifFun(CC,(double *) CC->Y,7,CC->FSave2);

       if (CC->evalja)
          {
           double rf = -CC->El[0] * CC->H;

           Jacobian(CC);

           for (i=0;i<NASubst;i++)
             for (j=0;j<NASubst;j++)
                CC->Df[NASubst*i+j]*=rf;

           for (i=0;i<NASubst;i++)
               CC->Df[(NASubst+1)*i]+=1.0;

           CC->evalja = 0;
           CC->conver = 0;

           CC->rc = 1.0;
           GetLU(CC,CC->iRow);
          }
       ExitLoop = 0;
       Break590 = 0;
       do
         {/* 460 */
          for (i=0;i<NASubst;i++)
              CC->FSave2[i] = CC->FSave2[i]*CC->H - CC->Y[i][1] - CC->Err[i];
          SolveLU(CC,CC->iRow,CC->FSave2,CC->FSave1);

          d = 0.0;
          for (i=0;i<NASubst;i++)
              {
               double V = CC->FSave1[i]/CC->Ymax[i];

               CC->Err[i] += CC->FSave1[i];
               CC->FSave1[i] = CC->Y[i][0]+CC->El[0]*CC->Err[i];
               d += V*V;
              }
          if (iter!=0) CC->trend = amax2(0.9*CC->trend,d/d1);
          if ((d*amin2(1.0,2.0*CC->trend))<=CC->bnd) Break590 = 1;
          else
             {
              d1 = d;
              iter++;
              if (iter==maxite) ExitLoop = 1;
              else DifFun(CC,CC->FSave1,1,CC->FSave2);
             }
         }
       while ((!ExitLoop) && (!Break590));
       ExitLoop = Break590 || (!CC->conver);
       if (!ExitLoop) CC->evalja=1;
      }
    while (!ExitLoop);

    if (!Break590)
       {
        CC->Time = Told;
        CC->rmax = rmxfai;
        RePredictor(CC,CC->Nq);

        if (fabs(CC->H)<=(MinH*1.00001))
           {
            CC->hold = CC->H;
            CC->jStart = CC->Nq;

            return -2;
           }
        CC->rh = rhcorr;
        Rescale(CC);
       }
    else
      {
       /* 590 */
       d = 0.0;
       for (i=0;i<NASubst;i++)
           {
            double V = CC->Err[i]/CC->Ymax[i];
            d += V*V;
           }
       CC->conver = 1;

       if (d>CC->e)
          {
           kFlag--;
           CC->Time = Told;
           RePredictor(CC,CC->Nq);

           CC->rmax = rmxfai;
           if (fabs(CC->H)<=(MinH*1.00001))
              {
               CC->hold = CC->H;
               CC->jStart = CC->Nq;

               return -1;
              }
           if (kFlag<=-maxfai)
              { /* 650 */
               CC->rh = amax2(MinH/fabs(CC->H),rherr3);
               CC->H *= CC->rh;
               DifFun(CC,(double *) CC->Y,7,CC->FSave1);
               for (i=0;i<NASubst;i++)
                   CC->Y[i][1] = CC->H * CC->FSave1[i];

               CC->evalja = 1;
               CC->idoub  = idelay;
               if (CC->Nq!=1) SetNq(CC,1);
              }
           else
             {
              pr2 = 1.0/(pow(d/CC->e,0.5/CC->l)*bias2+dfltzr);
              if (CC->Nq!=1)
                 {
                  double sum = 0.0;
                  for (i=0;i<NASubst;i++)
                      {
                       double V = CC->Y[i][CC->l-1] / CC->Ymax[i];
                       sum += V*V;
                      }
                  pr1 = 1.0/(pow(sum/CC->edn,0.5/CC->Nq)*bias1+dfltzr);
                  if (pr1>pr2)
                     { /* 640 */
                      double Savel0 = CC->oldl0;

                      SetNq(CC,CC->Nq - 1);
                      CC->rh =  pr1;
                      CC->rc *= CC->El[0]/Savel0;
                      Rescale(CC);
                      continue; /* GOTO 200 */
                     }
                 }
              CC->rh = pr2;
              Rescale(CC);
             }
          }
       else Break670 = 1;
      }
   }
 while (!Break670);
 /* 670 */
 for (j=0;j<CC->l;j++)
     for (i=0;i<NASubst;i++)
       CC->Y[i][j] += CC->El[j]*CC->Err[i];

 if (CC->idoub!=1)
    {
     CC->idoub--;
     if ((CC->idoub<=1) && (CC->Nq!=CC->MaxDer))
        for (i=0;i<NASubst;i++)
            CC->Y[i][CC->lmax-1] = CC->Err[i];
     CC->hold = CC->H;
     CC->jStart = CC->Nq;

     return 0;
    }
 pr3 = dfltzr;
 if (CC->Nq!=CC->MaxDer)
    {
     double sum = 0.0;
     for (i=0;i<NASubst;i++)
         {
          double V = (CC->Err[i] - CC->Y[i][CC->lmax-1]) / CC->Ymax[i];
          sum += V*V;
         }
     pr3 = 1.0/(pow(sum/CC->eup,0.5/(CC->l+1))*bias3+dfltzr);
    }
 /* 720 */
 pr2 = 1.0/(pow(d/CC->e,0.5/CC->l)*bias2+dfltzr);
 pr1 = dfltzr;
 if (CC->Nq!=1)
    {
     double sum = 0.0;
     for (i=0;i<NASubst;i++)
         {
          double V = CC->Y[i][CC->l-1] / CC->Ymax[i];
          sum += V*V;
         }
     pr1 = 1.0/(pow(sum/CC->edn,0.5/CC->Nq)*bias1+dfltzr);
    }
 /* 740 */
 if ((pr3>pr1) && (pr3>pr2))
    {
     newq = CC->l;
     CC->rh = pr3;
     if (CC->rh<thrshl)
        {
         CC->idoub  = idelay;
         CC->hold   = CC->H;
         CC->jStart = CC->Nq;

         return 0;
        }
     for (i=0;i<NASubst;i++)
         CC->Y[i][newq] = CC->Err[i]*CC->El[CC->l-1]/CC->l;
    }
 else
   if (pr1>pr2)
      {
       newq = CC->Nq-1;
       CC->rh = pr1;
       if (CC->rh<thrshl)
          {
           CC->idoub  = idelay;
           CC->hold   = CC->H;
           CC->jStart = CC->Nq;

           return 0;
          }
      }
   else
     {
      CC->rh = pr2;
      if (CC->rh<thrshl) CC->idoub = idelay;
      else
        {
         Rescale(CC);
         CC->rmax = rmxnor;
        }
      CC->hold   = CC->H;
      CC->jStart = CC->Nq;

      return 0;
     }
 /* 780 */
 Savel0 = CC->oldl0;
 SetNq(CC,newq);
 CC->rc *= CC->El[0]/Savel0;
 /* 790 */
 Rescale(CC);
 CC->rmax = rmxnor;
 /* 840 */
 CC->hold = CC->H;
 CC->jStart = CC->Nq;

 return 0;
}

int AdamsOrder = 4;
int Fail = 0;

double AdamsTable[5][6] = {
   {1.0,  1.0,  0.0,  0.0,   0.0,  0.0  },
   {2.0,  1.0,  1.0,  0.0,   0.0,  0.0  },
   {12.0, 5.0,  8.0,  -1.0,  0.0,  0.0  },
   {24.0, 9.0,  19.0, -5.0,  1.0,  0.0  },
   {720.0,251.0,646.0,-264.0,106.0,-19.0}
};

void Adams_Rozhkov(KineticContext * CC)
{
 int    i,j,k;
 int    CurOrder;
 double Xp[MaxActSubst];

 double AdamsH = CC->H/AdamsOrder;

 for (i=0; i<NASubst; i++)
     Xp[i] = CC->Y[i][0];
 Fail = 0;
 for (CurOrder = 1; CurOrder<=AdamsOrder && !Fail && !CC->LossPrecision; CurOrder++)
   {
    double F1 = AdamsH/AdamsTable[CurOrder-1][0];

    for (i=0;i<NASubst;i++)
        {
         CC->FSave1[i] = 0;
         CC->FSave2[i] = 0;
        }
    /* Вычисляем значения правых частей системы ОДУ, учитывая,
       что для каждого j-вещества правая часть записывается
       как CC->FSave1[j]+Conc[j]*CC->FSave2[j] */
    for (i=0;i<NReact;i++)
        {
         int Left  = LR[i].NumLeft;
         int Right = LR[i].NumRight;

         long Lj = 0; /* Номера веществ в левой части */
         /* Здесь W, WT содержит значение скорости реакции.
            Делить на AScale */
         int      Flags[MaxActSubst];
         double * W  = CC->Ymax;
         double   WT = CC->C[i];

         for (j=0;j<Left;j++)
             {
              k    =  LR[i].Indexes[j];
              Lj   |= (1<<k);
              W[k] =  CC->C[i];
              if (LR[i].k[j]==0.0)
                 WT *= CC->Y[k][CurOrder-1];
              else
                 WT *= pow(CC->Y[k][CurOrder-1],LR[i].k[j]);
              Flags[k] = 0;
             }
         for (k=0;k<Left;k++)
             if (!Flags[LR[i].Indexes[k]])
                {
                 for (j=0;j<Left;j++)
                   if (LR[i].k[j]!=0.0)
                      W[LR[i].Indexes[k]] *= pow(CC->Y[LR[i].Indexes[j]][CurOrder-1],LR[i].k[j]-(j==k ? 1.0 : 0.0));
                   else if (j!=k)
                      W[LR[i].Indexes[k]] *= CC->Y[LR[i].Indexes[j]][CurOrder-1];
                 Flags[LR[i].Indexes[k]] = 1;
                }

         for (j=0;j<Right;j++)
             {
              k = LR[i].Indexes[Left+j];
              if (Lj & (1<<k))
                 CC->FSave2[k] += W[k];
              else
                 CC->FSave1[k] += WT;
             }
         for (j=0;j<Left;j++)
             {
              k = LR[i].Indexes[j];
              CC->FSave2[k] -= W[k];
             }
        }

    for (i=0; i<NASubst && !Fail; i++)
        {
         double HH = 1-F1*CC->FSave2[i]*AdamsTable[CurOrder-1][1];

         if (fabs(HH)>dfltzr)
            {
             CC->Y[i][CurOrder] = CC->FSave1[i]*AdamsTable[CurOrder-1][1];
             for (j=2; j<=CurOrder; j++)
                 CC->Y[i][CurOrder] += AdamsTable[CurOrder-1][j]*CC->Y[i][CurOrder-j];
             CC->Y[i][CurOrder] = (CC->Y[i][CurOrder-1]+F1*CC->Y[i][CurOrder])/HH;
             CC->Y[i][CurOrder-1] = CC->FSave1[i]+CC->Y[i][CurOrder]*CC->FSave2[i];
            }
         else
            {
             Fail = 1;
             CC->H /= 2.0;
            }
        }
    if ((CC->Time + CC->H) == CC->Time) CC->LossPrecision = 1;
   }

 if (!CC->LossPrecision)
    {
     if (!Fail)
        {
         for (i=0;i<NASubst;i++)
             CC->Y[i][0] = CC->Y[i][CurOrder-1];
         CC->Time += CC->H;
        }
     else
        for (i=0;i<NASubst;i++)
            CC->Y[i][0] = Xp[i];
     if (CC->_EndTime - CC->Time < CC->H)
        CC->H = CC->_EndTime - CC->Time;
     if ((CC->Time + CC->H) == CC->Time)
        CC->Time = CC->_EndTime;
    }
 CC->Iters += AdamsOrder-1;
}

double RosenTable[4][4] = {
   {-1.0,     0.0,       0.0,      0.0     },
   {1.0/8.0,  3.0/8.0,   1.0,      0.0     },
   {3.0/8.0,  19.0/24.0, -1.0/6.0, 0.0     },
   {13.0/6.0, 1.0/6.0,   -2.0,     2.0/3.0 }
};

void Rosenbrock(KineticContext * CC)
{
 int i,j,k;

 Fail = 0;

 Jacobian(CC);
 for (i=0;i<NASubst;i++)
   for (j=0;j<NASubst;j++)
      CC->Df[NASubst*i+j] *= -CC->H;
 for (i=0;i<NASubst;i++)
     CC->Df[(NASubst+1)*i] += 1.0;
 GetLU(CC,CC->iRow);

 for (k=1; k<=4; k++)
     {
      for (i=0;i<NASubst;i++)
          {
           CC->Y[i][k] = CC->Y[i][0];
           for (j=0; j<k-1; j++)
               CC->Y[i][k] += RosenTable[k-2][j] * CC->Y[i][j+1];
          }
          
      DifFun(CC,(double *) &CC->Y[0][k],7,CC->FSave2);
      for (i=0;i<NASubst;i++)
          CC->FSave2[i] *= CC->H;
      SolveLU(CC,CC->iRow,CC->FSave2,CC->FSave1);
      for (i=0;i<NASubst;i++)
          CC->Y[i][k] = CC->FSave1[i];
     }

 for (i=0;i<NASubst;i++)
     for (k=1; k<=4; k++)
         CC->Y[i][0] += RosenTable[3][k-1] * CC->Y[i][k];
 CC->Time += CC->H;
 if (CC->_EndTime - CC->Time < CC->H)
    CC->H = CC->_EndTime - CC->Time;
 if ((CC->Time + CC->H) == CC->Time)
    CC->Time = CC->_EndTime;
}

typedef void (* EvalFunction)(KineticContext *);

EvalFunction EvalF[2] = {Adams_Rozhkov,Rosenbrock};

#if defined(__UNIX__) || defined(__MVS__) || defined(__MAY_BE_MICROSOFTC__)
double pow10(int Arg)
{
 double result = 1;

 if (Arg>0)
   while((Arg--)>0) result*=10;
 else if (Arg<0)
   while((Arg++)<0) result/=10;

 return result;
}
#endif

void OneTaktKinetic(KineticContext * CC,int UseStiffGearMethod, int OtherMethod)
{
 int i,j;

 CC->Tem = R*CC->Tk;

 for (i=0;i<NASubst;i++)
     if (CC->Conc0[i]<0.0) CC->Conc0[i] = 0.0;

 CC->Iters = 0;
 if (CC->Tem>dfltzr)
    {
     /* Автоматический выбор масштабного множителя */
     double maxY = dfltzr;

     for (i=0;i<NASubst;i++)
         if (CC->Conc0[i]>maxY) maxY = CC->Conc0[i];

     if (maxY>dfltzr)
        {
         CC->LgAScale = (int) -log10(maxY);
         CC->AScale   = pow10(CC->LgAScale);

         /* Вычисление коэффициентов скорости */
         for (i=0;i<NReact;i++)
             {
              double PEi = (_A[i]-(LR[i].NumLeft-1)*CC->LgAScale)/M_LOG10E;
              double Degree = PEi-EA[i]/CC->Tem;

              if (Degree>-80.0) CC->C[i] = exp(Degree)*pow(CC->Tk,(float)Tn[i]);
              else CC->C[i] = 0.0;
             }

         if (UseStiffGearMethod)
            {
             double EpsCom = 1.0E-12/CalcEps;
             /* STIFFC */
             for (i=0;i<NASubst;i++)
                 {
                  CC->Ymax[i]=0.0;
                  CC->Y[i][0]=CC->Conc0[i]*CC->AScale;
                  for (j=1;j<6;j++)
                      CC->Y[i][j]=0.0;
                 }

             CC->Time = InitTime;

             if (CC->_EndTime<InitH) CC->H=CC->_EndTime;
             else CC->H = InitH;

             CC->LossPrecision = (CC->Time+CC->H)==CC->Time;

             CC->MaxDer = 5;

             CC->jStart = 0;

             for (CC->Iters = 0;CC->Time<=CC->_EndTime && CC->Iters<MaxIterations && !CC->LossPrecision;CC->Iters++)
              {
               double ys = 0.0;

               for (i=0;i<NASubst;i++)
                   ys += fabs(CC->Y[i][0]);
               ys *= EpsCom;
               for (i=0;i<NASubst;i++)
                   CC->Ymax[i] = amax3(CC->Y[i][0],ys,dfltzr);

               if (Stiff(CC))
                  {
                   CC->jStart=-1;
                   CC->H/=10.0;
                  }
               if ((CC->Time+CC->H)==CC->Time) CC->LossPrecision = 1;
              }
             if (CC->Iters==MaxIterations) CC->LossPrecision = 1;

             if (CC->LossPrecision)
                {
                 CC->KinErrorInfo.LossH    = CC->Iters!=MaxIterations;
                 CC->KinErrorInfo.LastH    = CC->H;
                 CC->KinErrorInfo.ReachTau = CC->Time;
                 CC->KinErrorInfo.Tmp      = (float) (CC->Tk-ZeroK);
                }
             else
                for (i=0;i<NASubst;i++)
                    {
                     double sum = 0.0;
                     for (j=0;j<=CC->jStart;j++)
                         sum += CC->Y[i][j]*pow((CC->_EndTime-CC->Time)/CC->H,j);
                     CC->Conc1[i] = (float)(sum/CC->AScale);
                    }
            }
         else
            {
             double Divider = 1.0;
             double Divisor = 2.0;
             double _EndTime = CC->_EndTime;
             double _H = amax2((CC->_EndTime-InitTime)/Divider,InitH);

             CC->Iters = 0;

             for (i=0;i<NASubst;i++)
                 CC->Y[i][0] = CC->Conc0[i]*CC->AScale;

             CC->Time = InitTime;
             if (_H>InitH)
                {
                 float StFull[MaxActSubst];
                 int   EndFlag = 0;

                 for (i=0;i<NASubst;i++)
                     CC->Y[i][6] = CC->Y[i][0];

                 do
                   {
                    CC->_EndTime = _H;
                    CC->H = _H;
                    EvalF[OtherMethod](CC);
                    CC->Iters++;
                    for (i=0;i<NASubst;i++)
                        {
                         StFull[i] = (float) CC->Y[i][0];
                         CC->Y[i][0] = CC->Y[i][6];
                        }
                    if (Fail)
                       {
                        _H /= Divisor;
                        Divider *= Divisor;
                       }
                   }
                 while (Fail);

                 while (_H>InitH && Divisor*Divider<MaxIterations && !EndFlag)
                   {
                    float StHalf[MaxActSubst];
                    double Ymax = 0.0;

                    for (i=0;i<NASubst;i++)
                        CC->Y[i][0] = CC->Y[i][6];
                    CC->Time = InitTime;
                    CC->_EndTime = _H/2.0;
                    CC->H = CC->_EndTime;
                    EvalF[OtherMethod](CC);
                    CC->Iters++;
                    for (i=0;i<NASubst;i++)
                        StHalf[i] = (float) CC->Y[i][0];
                    CC->Time = InitTime;
                    CC->H = CC->_EndTime;
                    EvalF[OtherMethod](CC);
                    CC->Iters++;
                    for (i=0;i<NASubst;i++)
                        if (fabs(CC->Y[i][0])>Ymax) Ymax = fabs(CC->Y[i][0]);
                    EndFlag = 1;
                    for (i=0;i<NASubst;i++)
                        {
                         if (fabs(CC->Y[i][0]-StFull[i])>CalcEps*Ymax) EndFlag = 0;

                         StFull[i] = StHalf[i];
                        }
                    if (!EndFlag)
                       {
                        _H /= Divisor;
                        Divider *= Divisor;
                       }
                   }
                 CC->Time = InitTime+_H;
                }

             CC->_EndTime = _EndTime;
             CC->H       = amax2(_H,InitH);
             CC->LossPrecision = 0;

             for (;CC->Time<CC->_EndTime && CC->Iters<MaxIterations && !CC->LossPrecision;CC->Iters++)
                 EvalF[OtherMethod](CC);

             if (CC->Iters==MaxIterations) CC->LossPrecision = 1;

             if (CC->LossPrecision)
                {
                 CC->KinErrorInfo.LossH    = CC->Iters!=MaxIterations;
                 CC->KinErrorInfo.LastH    = CC->H;
                 CC->KinErrorInfo.ReachTau = CC->Time;
                 CC->KinErrorInfo.Tmp      = (float)(CC->Tk-ZeroK);
                }
             else
                for (i=0;i<NASubst;i++)
                    CC->Conc1[i] = (float)(CC->Y[i][0]/CC->AScale);
            }
        }
     else
        for (i=0;i<NASubst;i++)
            CC->Conc1[i] = CC->Conc0[i];
    }
 else
   for (i=0;i<NASubst;i++)
       CC->Conc1[i] = CC->Conc0[i];
}

#endif
