// Все double заменены на float
// Все константы преобразованы либо в макросимволы, либо записаны с модификатором __constant
// Глобальные данные помещены в общую структуру со статическими массивами и теперь находятся в глобальной памяти.
// Контексты размещены в локальной памяти и модифицированы: все массивы стали статическими, в связи с чем снизилось
// количество возможных задействованных потоков из-за потенциальной нехватки локальной памяти. По той же причине
// максимальное количество веществ снижено до 10.
// В ПППВ передаются: ссылка на блок глобальной памяти, блок локальной памяти и некоторые настроечные данные.

#ifndef __KINETIC__
#define __KINETIC__

#include <stdlib.h>
#include <math.h>

#pragma plan common begin

/* Максимальное число веществ, участвующих в реакции */
#define MaxActSubst 15 /* 10 !!! for GPU */
/* Максимальное количество реакций */
#define MaxReact 17 /* 6 !!! for GPU */
/* Максимальная длина имени вещества */
#define SubstNameLength 21

#define MaxOrder 5

#define Adams_Rozhkov_Method 0
#define Rosenbrock_Method 1

#define OTHER_KIN_METHOD Adams_Rozhkov_Method
// #define OTHER_KIN_METHOD Adams_Rozhkov_Method

typedef struct {
  int    NumLeft;
  int    NumRight;
  int    Indexes[MaxActSubst];
  float  k[MaxActSubst];
} Reaction;

typedef struct {
 int      LossH;    /* Признак малости H */
 float    LastH;    /* Последнее H */
 float    ReachTau; /* Достигнутое значение TAU */
 float    Tmp;      /* Температура проведения реакций */
} ErrorKinetic;

#define LAST_DERIVS (2*MaxOrder)

typedef float Derivs[LAST_DERIVS + 1];

typedef char   SubstName[SubstNameLength];

#ifndef M_LOG10E
#define M_LOG10E 0.43429448f
#endif

/* Constants */
/* Константы STIFF */
#define dfltzr 1.0E-18f
#define maxite 3
#define maxfai 3
#define rmxini 1.0E+4f
#define rmxnor 10.0f
#define rmxfai 2.0f
#define idelay 10
#define rhcorr 0.25f
#define rherr3 0.1f
#define thrshl 0.0f
#define rctest 0.3f
#define bias1  1.3f
#define bias2  1.2f
#define bias3  1.4f

#define ZeroK 273.0f
#define _R_   1.987E-3f

__constant float Koeffs[5][3] =
  {
   {1.0f,2.0f,3.0f},
   {1.0f,4.5f,6.0f},
   {0.5f,7.3333333f,9.1666667f},
   {0.1666667f,10.4166667f,12.5f},
   {0.04166667f,13.7f,1.0f}
  };

#define KinetAlpha 0.0025f

typedef struct {
	int MaxIterations; /* Максимальное число итераций */

	float InitTime; /* Начальное время интегрирования */
	float CalcEps; /* Относительная точность интегрирования ST1/eps */
	float MinH;/* Минимальный шаг интегрирования ST1/hmin */
	/* End of Constants */

	float DefaultT;

	/* Shared variables */
	float InitH;  /* Начальный шаг интегрирования h */
	/* Имена веществ [KG->NSubst] */
	SubstName Names[MaxActSubst];
	/* Скорости витания веществ [KG->NSubst] */
	float     Uw[MaxActSubst];
	/* Карта участия веществ в реакциях */
	char      Activity[MaxActSubst];
	/* Индексы веществ, участвующих в химических реакциях */
	int       TranMap[MaxActSubst];
	/* Ei : Энергии активации KINET/ea [KG->NReact] */
	float     EA[MaxReact];
	/* Ai : Lg( предэкспоненты констант скорости ) [KG->NReact] */
	float     _A[MaxReact];
	/* Ni : Показатели степеней в температурных множителях KINET/tn [KG->NReact] */
	float     Tn[MaxReact];
	/* Реакции KINET/KG->LR [KG->NReact] */
	Reaction  LR[MaxReact];

	int NReact;  /* Число реакций KINET/m */
	int NSubst;  /* Число веществ KINET/n = ST1/n */
	int NASubst; /* Число веществ, участвующих в реакции */

	int AdamsOrder;
	int Fail;
	/* End of shared variables */
} KineticGlobal;

typedef struct {
  ErrorKinetic KinErrorInfo; /* Блок информации по ошибке */

  int LossPrecision;        /* Флаг потери точности при вычислениях */
  int Iters; /* Действительное число итераций */

  float  Tk;   /* Температура в Кельвинах KINET/tk */
  float  _EndTime; /* Конечное время интегрирования */

  /* Начальная концентрация веществ [KG->NASubst] */
  float Conc0[MaxActSubst];
  /* Конечная концентрация [KG->NASubst] */
  float Conc1[MaxActSubst];

  /* Рабочие переменные [KG->NASubst] */
  Derivs * Y;
  /* KINET */
  int    LgAScale; /* ml = Lg(aml) */
  float  AScale; /* Масштабный множитель aml */
  float  C[MaxReact];    /* Вычисляемые константы скоростей реакций c [KG->NReact] */
  float  Tem;    /* Температура в ккал/моль tem */
  /* ST1 */
  float  Time;
  float  H;
  int    MaxDer;
  int    jStart;
  /* ST2 */
  float  Ymax[MaxActSubst]; /* [KG->NASubst] */
  /* ST3 */
  float  Err[MaxActSubst]; /* [KG->NASubst] */
  /* ST4 */
  float  * Df; /* [KG->NASubst][KG->NASubst] */
  /* ST5 */
  float  * FSave1; /* [KG->NASubst] */
  float  * FSave2; /* [KG->NASubst] */
  /* STIFF */
  int    iRow[MaxActSubst]; /* [KG->NASubst] */

  float  El[6];
  float  Tq[4];
  float  edn;
  float  e;
  float  eup;
  float  bnd;

  int   Nq;
  int   l;
  int   lmax;
  int   idoub;
  float rmax;
  float rc;

  float rh;

  float trend;
  float oldl0;
  float hold;

  int   evalja;
  int   conver;
} KineticContext;

float amax2(float First, float Second)
{
 return First>Second ? First : Second;
}

float amax3(float First, float Second, float Third)
{
 return First>Second ? amax2(First,Third) : amax2(Second,Third);
}

float amin2(float First, float Second)
{
 return First<Second ? First : Second;
}

float amin3(float First, float Second, float Third)
{
 return First<Second ? amin2(First,Third) : amin2(Second,Third);
}

float fabsr(float a) {
	return fabs(a);
}

float _powr(float a, float b) {
	return pow(a,b) /* native_powr(a,b) */;
}

float expr(float a) {
	return exp(a) /* native_exp(a) */;
}

#pragma plan common end

void lmemset(void * mem, char what, unsigned int size) {
 memset(mem, what, size);
}

void rmemset(void * mem, char what, unsigned int size) {
 memset(mem, what, size);
}

void gmemset(void * mem, char what, unsigned int size) {
 memset(mem, what, size);
}

#pragma plan gpu begin

void lmemset(__local void * mem, char what, unsigned int size) {
	__local char * _mem = (__local char *) mem;
	int i;

	for (i = 0; i < size; i++)
		_mem[i] = what;
}

void rmemset(void * mem, char what, unsigned int size) {
	char * _mem = (char *) mem;
	int i;

	for (i = 0; i < size; i++)
		_mem[i] = what;
}

void gmemset(__global void * mem, char what, unsigned int size) {
	__global char * _mem = (__global char *) mem;
	int i;

	for (i = 0; i < size; i++)
		_mem[i] = what;
}

#pragma plan gpu end

void SetTk(KineticContext * C, float Tmp)
{
	C->Tk = (float)(ZeroK + Tmp);
}

/* Выводим из рассмотрения вещества, фактически не участвующие в реакциях */
void ReTranslate(KineticGlobal * KG)
{
 char Map[MaxActSubst];
 char Delta = 0;
 int  i,j;

 KG->NASubst = 0;
 for (i=0; i<KG->NSubst; i++)
     if (KG->Activity[i])
        {
         Map[i] = Delta;
         KG->TranMap[KG->NASubst++] = i;
        }
     else
        Delta++;
 for (i=0; i<KG->NReact; i++)
     for (j=0; j<KG->LR[i].NumLeft+KG->LR[i].NumRight; j++)
         KG->LR[i].Indexes[j] -= Map[KG->LR[i].Indexes[j]];
}

int Translate(KineticGlobal * KG, int NumEq, char * Equation)
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
    float  How  = 0.0;
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
        How = (float)atof(&Equation[BegCount]);
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

    while (Index<KG->NSubst && !Found)
       {
        if (strlen(KG->Names[Index])==CurLong)
           {
            int CurN = 0;
            int CurE = Count;

            while (CurN<CurLong && KG->Names[Index][CurN]==Equation[CurE])
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
       KG->LR[NumEq].k[All] = How;
       KG->LR[NumEq].Indexes[All++] = Index;
      }
    KG->Activity[Index] = 1;

    if (End<Long-1)
       if (Equation[End]=='-')
          if (Equation[++End]=='>' && !Left)
             Left = All;
          else return 0;
    Count = End+1;
   }
 KG->LR[NumEq].NumLeft  = Left;
 KG->LR[NumEq].NumRight = All-Left;

 return (Left);
}

#pragma plan common begin

/* ВЫЧИСЛЕНИЕ ЗНАЧЕНИЙ ЯКОБИАНА ИНТЕГРИРУЕМОЙ СИСТЕМЫ ОДУ ПО */
/* ЗНАЧЕНИЯМ ВРЕМЕНИ t И РЕШЕНИЙ y, СООТВЕТСТВУЮЩИХ ЭТОМУ МОМЕНТУ */
/* ВРЕМЕНИ. ЭЛЕМЕНТЫ ЯКОБИАНА ЗАСЫЛАЮТСЯ В МАССИВ Df[50,50) */
/* PEDERV */
#ifndef __PARALLEL__
$ void Jacobian(__global KineticGlobal * KG, __local KineticContext * C)
#else
void Jacobian(__global KineticGlobal * KG, __local KineticContext * C)
#endif
{
 int i,j,k;

 int PtrIJ = 0;

 for (i=0;i<MaxActSubst*MaxActSubst;i++)
     C->Df[i] = 0.0f;
/* ЦИКЛ ПО ЧИСЛУ УРАВНЕНИЙ */
 for (i=0;i<KG->NReact;i++)
     {
      int Left  = KG->LR[i].NumLeft;
      int Right = KG->LR[i].NumRight;

      for (j=0;j<Left;j++)
          {
           int Cur = KG->LR[i].Indexes[j];

           float  W  = C->C[i];

           if (KG->LR[i].k[j]!=0.0) {
              W*=KG->LR[i].k[j]*_powr(C->Y[Cur][0],KG->LR[i].k[j]-1.0f);
           }
           for (k=0;k<Left;k++)
               if (k!=j)
                  if (KG->LR[i].k[k]==0.0) {
                     W*=C->Y[KG->LR[i].Indexes[k]][0];
                  } else
                     W*=_powr(C->Y[KG->LR[i].Indexes[k]][0],KG->LR[i].k[k]);
           for (k=0;k<Left;k++)
               C->Df[KG->LR[i].Indexes[k]*MaxActSubst+Cur] -= W;
           for (k=0;k<Right;k++)
               C->Df[KG->LR[i].Indexes[Left+k]*MaxActSubst+Cur] += W;
          }
     }
}

/* ВЫЧИСЛЯЕТ ЗНАЧЕНИЯ ПРАВЫХ ЧАСТЕЙ СИСТЕМЫ ОДУ ПО ЗНАЧЕНИЯМ */
/* ВРЕМЕНИ t И РЕШЕНИЙ y0, СООТВЕТСТВУЮЩИХ ЭТОМУ МОМЕНТУ ВРЕМЕНИ, */
/* И ЗАСЫЛАЕТ ИХ В МАССИВ f */
/* DIFFUN */
#ifndef __PARALLEL__
$ void DifFun(__global KineticGlobal * KG, __local KineticContext * CC,float  * Y0,int Ystep,float  * F)
#else
void DifFun(__global KineticGlobal * KG, __local KineticContext * CC,float  * Y0,int Ystep,float  * F)
#endif
{
 int i,j;

 for (i=0;i<KG->NASubst;i++)
     {
      if (Y0[Ystep*i]<0.0) Y0[Ystep*i] = 0.0;
      F[i] = 0.0;
     }

 for (i=0;i<KG->NReact;i++)
     {
      int Left  = KG->LR[i].NumLeft;
      int Right = KG->LR[i].NumRight;

      float  W = CC->C[i];

      for (j=0;j<Left;j++)
          if (KG->LR[i].k[j]==0.0) {
             W*=Y0[Ystep*KG->LR[i].Indexes[j]];
          } else
             W*=_powr(Y0[Ystep*KG->LR[i].Indexes[j]],KG->LR[i].k[j]);
      /* Здесь W содержит значение скорости соответствующей реакции.
         Делить на AScale */
      for (j=0;j<Left;j++)
          F[KG->LR[i].Indexes[j]] -= W;
      for (j=0;j<Right;j++)
          F[KG->LR[i].Indexes[Left+j]] += W;
     }
}

/* ПЕРЕМАСШТАБИРОВАНИЕ РЕШЕНИЙ СИСТЕМЫ И ИХ ПРОИЗВОДНЫХ */
/* ПРИ ИЗМЕНЕНИИ ШАГА ИНТЕГРИРОВАНИЯ */
/* RESCAL */
#ifndef __PARALLEL__
$ void Rescale(__global KineticGlobal * KG, __local KineticContext * CC)
#else
void Rescale(__global KineticGlobal * KG, __local KineticContext * CC)
#endif
{
 int i,j;

 float  R1=1.0;

 CC->rh = amax2(CC->rh,KG->MinH/fabsr(CC->H));
 CC->rh = amin3(CC->rh,CC->_EndTime/fabsr(CC->H),CC->rmax);
 for (j=1;j<CC->l;j++)
     {
      R1 *= CC->rh;
      for (i=0;i<KG->NASubst;i++)
          CC->Y[i][j] *= R1;
     }
 CC->H *= CC->rh;
 CC->rc *= CC->rh;
 CC->idoub = CC->l+1;
}

/* ВЫЧИСЛЕНИЕ ЗНАЧЕНИЯ ПРЕДИКТОРА */
void Predictor(__global KineticGlobal * KG, __local KineticContext * CC, int Nq)
{
 int i,j,j1,j2;

 for (j1=0;j1<Nq;j1++)
     for (j2=j1;j2<Nq;j2++)
         {
          j = Nq-j2+j1-1;

          for (i=0;i<KG->NASubst;i++)
              CC->Y[i][j] += CC->Y[i][j+1];
         }
}

void RePredictor(__global KineticGlobal * KG, __local KineticContext * CC, int Nq)
{
 int i,j,j1,j2;

 for (j1=0;j1<Nq;j1++)
     for (j2=j1;j2<Nq;j2++)
         {
          j = Nq-j2+j1-1;

          for (i=0;i<KG->NASubst;i++)
              CC->Y[i][j] -= CC->Y[i][j+1];
         }
}

/* ВЫЧИСЛЕНИЕ КОЭФФИЦИЕНТОВ МЕТОДА ИНТЕГРИРОВАНИЯ */
/* СООТВЕТСТВУЮЩЕГО ПОРЯДКА */
/* COSET */
void CalcKoeffs(__global KineticGlobal * KG, __local KineticContext * CC, int Nq)
{
 int i;

 CC->El[1]=1.0;
 if (CC->MaxDer>5) CC->MaxDer = 5;
 switch (Nq) {
   case 1: CC->El[0] = 1.0f;
           break;
   case 2: CC->El[0] = 6.666666666667E-01f;
           CC->El[2] = 3.333333333333E-01f;
           break;
   case 3: CC->El[0] = 5.454545454546E-01f;
           CC->El[2] = 1.0f;
           CC->El[3] = 9.090909090909E-02f;
           break;
   case 4: CC->El[0] = 0.48f;
           CC->El[2] = 0.7f;
           CC->El[3] = 0.2f;
           CC->El[4] = 0.02f;
           break;
   case 5: CC->El[0] = 4.379562043796E-01f;
           CC->El[2] = 8.211678832117E-01f;
           CC->El[3] = 3.102189781022E-01f;
           CC->El[4] = 5.474452554745E-02f;
           CC->El[5] = 3.649635036496E-03f;
 }
 for (i=0;i<3;i++)
     CC->Tq[i] = Koeffs[Nq-1][i];
 CC->Tq[3]=0.5f/(Nq+2.0f);
}

/* ОСУЩЕСТВЛЯЕТ LU-РАЗЛОЖЕНИЕ МАТРИЦЫ CC->e-beta*h*A, */
/* ГДЕ A-ЯКОБИАН СИСТЕМЫ, CC->e-ЕДИНИЧНАЯ МАТРИЦА, */
/* h-ТЕКУЩИЙ ШАГ ИНТЕГРИРОВАНИЯ, beta-КОНСТАНТА */
/* MATDEC */
void GetLU(__global KineticGlobal * KG, __local KineticContext * CC, __local int * iRow)
{
 int i,j,k;

 if (KG->NASubst<=1)
    {
     CC->Df[0] = 1.0f/CC->Df[0];
     iRow[0]  = 0;
    }
 else
   {
    for (i=0;i<KG->NASubst;i++)
        iRow[i]=i;
    for (i=0;i<KG->NASubst-1;i++)
      {
       float  Big  = 0.0;
       int    iBig = 0;

       float  Kf;
       int    ir;

       for (j=i;j<KG->NASubst;j++)
           {
            /* По непонятным причинам компилятор Borland C++ 3.1
               отказался скомпилировать следующую строку :
            float  size = fabs(Df[KG->NASubst*iRow[j]+i]);
            */
            float  size    = fabsr(CC->Df[iRow[j]*MaxActSubst+i]);

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
       ir = iRow[i];

       Kf = 1.0f/CC->Df[ir*MaxActSubst+i];

       CC->Df[ir*MaxActSubst+i] = Kf;
       for (j=i+1;j<KG->NASubst;j++)
           {
            int jr = iRow[j];

            float  Fact = Kf*CC->Df[jr*MaxActSubst+i];

            CC->Df[jr*MaxActSubst+i] = Fact;
            for (k=i+1;k<KG->NASubst;k++)
                CC->Df[jr*MaxActSubst+k] -= Fact*CC->Df[ir*MaxActSubst+k];
           }
      }
    i = iRow[KG->NASubst-1];

    CC->Df[i*MaxActSubst+KG->NASubst-1] = 1.0f/CC->Df[i*MaxActSubst+KG->NASubst-1];
   }
}

/* РЕШЕНИЕ СИСТЕМЫ ЛИНЕЙНЫХ УРАВНЕНИЙ */
void SolveLU(__global KineticGlobal * KG, __local KineticContext * CC, __local int * iRow, float  * Y, float  * X)
{
 int i,j,k;

 if (KG->NASubst<=1) X[0] = Y[0]*CC->Df[0];
 else
   {
    X[0] = Y[iRow[0]];
    for (i=1;i<KG->NASubst;i++)
        {
         int   ir = iRow[i];
         float  V  = Y[ir];

         for (j=0;j<i;j++)
             V -= CC->Df[ir*MaxActSubst+j]*X[j];
         X[i] = V;
        }

    X[KG->NASubst-1]*=CC->Df[iRow[KG->NASubst-1]*MaxActSubst+KG->NASubst-1];

    for (i=1,j=KG->NASubst-2;i<KG->NASubst;i++,j--)
        {
         int   jr = iRow[j];
         float  V = X[j];

         for (k=j+1;k<KG->NASubst;k++)
             V -= CC->Df[jr*MaxActSubst+k]*X[k];
         X[j] = V*CC->Df[jr*MaxActSubst+j];
        }
   }
}

#define SetNq(KG,CC,NewNq)       \
{                                \
 CC->Nq = NewNq;                 \
 CC->l  = CC->Nq+1;              \
 CalcKoeffs(KG,CC,CC->Nq);       \
                                 \
 CC->lmax   = CC->MaxDer+1;      \
 CC->oldl0  = CC->El[0];         \
 CC->edn    = CC->Tq[0]*KG->CalcEps; \
 CC->edn   *= CC->edn;           \
 CC->e      = CC->Tq[1]*KG->CalcEps; \
 CC->e     *= CC->e;             \
 CC->eup    = CC->Tq[2]*KG->CalcEps; \
 CC->eup   *= CC->eup;           \
 CC->bnd    = CC->Tq[3]*KG->CalcEps; \
 CC->bnd   *= CC->bnd;           \
}

#ifndef __PARALLEL__
$ void Stiff(__global KineticGlobal * KG, __local KineticContext * CC, int * result)
#else
void Stiff(__global KineticGlobal * KG, __local KineticContext * CC, int * result)
#endif
{
 int   kFlag=0;
 float  Told=CC->Time;

 float  pr1,pr2,pr3;
 float  d;
 int newq;

 int Break590;
 int Break670;
 int ExitLoop;

 float  Savel0;

 int i,j;

 if (CC->jStart==0)
    {
     DifFun(KG,CC,(float  *) CC->Y,7,CC->FSave1);
     for (i=0;i<KG->NASubst;i++)
         CC->Y[i][1]=CC->FSave1[i]*CC->H;

     SetNq(KG,CC,1);

     CC->idoub  = CC->l+1;
     CC->rmax   = rmxini;
     CC->trend  = 1.0;
     CC->rc     = 0.0;
     CC->hold   = CC->H;
     CC->evalja = 1;
     CC->conver = 0;
    }
 else {
   if (CC->jStart<0) {
      if (CC->H!=CC->hold)
         {
          CC->rh = CC->H/CC->hold;
          CC->H  = CC->hold;
          Rescale(KG,CC);
         }
   }
 }

 Break590 = 0;
 Break670 = 0;
 ExitLoop = 0;
 do /* 200 */
   {
    if (fabs(CC->rc-1.0)>rctest) CC->evalja = 1;
    CC->Time+=CC->H;
    Predictor(KG,CC,CC->Nq);

    do
      {
       int iter = 0;

       float  d1;

       for (i=0;i<KG->NASubst;i++)
           CC->Err[i] = 0.0;
       DifFun(KG,CC,(float  *) CC->Y,7,CC->FSave2);

       if (CC->evalja)
          {
           float  rf = -CC->El[0] * CC->H;

           Jacobian(KG,CC);

           for (i=0;i<KG->NASubst;i++)
             for (j=0;j<KG->NASubst;j++)
                CC->Df[i*MaxActSubst+j]*=rf;

           for (i=0;i<KG->NASubst;i++)
               CC->Df[i*MaxActSubst+i]+=1.0;

           CC->evalja = 0;
           CC->conver = 0;

           CC->rc = 1.0;
           GetLU(KG,CC,CC->iRow);
          }
       ExitLoop = 0;
       Break590 = 0;
       do
         {/* 460 */
          for (i=0;i<KG->NASubst;i++)
              CC->FSave2[i] = CC->FSave2[i]*CC->H - CC->Y[i][1] - CC->Err[i];
          SolveLU(KG,CC,CC->iRow,CC->FSave2,CC->FSave1);

          d = 0.0;
          for (i=0;i<KG->NASubst;i++)
              {
               float  V = CC->FSave1[i]/CC->Ymax[i];

               CC->Err[i] += CC->FSave1[i];
               CC->FSave1[i] = CC->Y[i][0]+CC->El[0]*CC->Err[i];
               d += V*V;
              }
          if (iter!=0) {
             CC->trend = amax2(0.9f*CC->trend,d/d1);
          }
          if ((d*amin2(1.0f,2.0f*CC->trend))<=CC->bnd) Break590 = 1;
          else
             {
              d1 = d;
              iter++;
              if (iter==maxite) ExitLoop = 1;
              else DifFun(KG,CC,CC->FSave1,1,CC->FSave2);
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
        RePredictor(KG,CC,CC->Nq);

        if (fabs(CC->H)<=(KG->MinH*1.00001))
           {
            CC->hold = CC->H;
            CC->jStart = CC->Nq;

            *result = -2; return;
           }
        CC->rh = rhcorr;
        Rescale(KG,CC);
       }
    else
      {
       /* 590 */
       d = 0.0;
       for (i=0;i<KG->NASubst;i++)
           {
            float  V = CC->Err[i]/CC->Ymax[i];
            d += V*V;
           }
       CC->conver = 1;

       if (d>CC->e)
          {
           kFlag--;
           CC->Time = Told;
           RePredictor(KG,CC,CC->Nq);

           CC->rmax = rmxfai;
           if (fabs(CC->H)<=(KG->MinH*1.00001))
              {
               CC->hold = CC->H;
               CC->jStart = CC->Nq;

               *result = -1; return;
              }
           if (kFlag<=-maxfai)
              { /* 650 */
               CC->rh = amax2(KG->MinH/fabsr(CC->H),rherr3);
               CC->H *= CC->rh;
               DifFun(KG,CC,(float  *) CC->Y,7,CC->FSave1);
               for (i=0;i<KG->NASubst;i++)
                   CC->Y[i][1] = CC->H * CC->FSave1[i];

               CC->evalja = 1;
               CC->idoub  = idelay;
               if (CC->Nq!=1) SetNq(KG,CC,1);
              }
           else
             {
              pr2 = 1.0f/(_powr(d/CC->e,0.5f/CC->l)*bias2+dfltzr);
              if (CC->Nq!=1)
                 {
                  float  sum = 0.0;
                  for (i=0;i<KG->NASubst;i++)
                      {
                       float  V = CC->Y[i][CC->l-1] / CC->Ymax[i];
                       sum += V*V;
                      }
                  pr1 = 1.0f/(_powr(sum/CC->edn,0.5f/CC->Nq)*bias1+dfltzr);
                  if (pr1>pr2)
                     { /* 640 */
                      float  Savel0 = CC->oldl0;

                      SetNq(KG,CC,CC->Nq - 1);
                      CC->rh =  pr1;
                      CC->rc *= CC->El[0]/Savel0;
                      Rescale(KG,CC);
                      continue; /* GOTO 200 */
                     }
                 }
              CC->rh = pr2;
              Rescale(KG,CC);
             }
          }
       else Break670 = 1;
      }
   }
 while (!Break670);
 /* 670 */
 for (j=0;j<CC->l;j++)
     for (i=0;i<KG->NASubst;i++)
       CC->Y[i][j] += CC->El[j]*CC->Err[i];

 if (CC->idoub!=1)
    {
     CC->idoub--;
     if ((CC->idoub<=1) && (CC->Nq!=CC->MaxDer)) {
        for (i=0;i<KG->NASubst;i++)
            CC->Y[i][CC->lmax-1] = CC->Err[i];
     }
     CC->hold = CC->H;
     CC->jStart = CC->Nq;

     *result = 0; return;
    }
 pr3 = dfltzr;
 if (CC->Nq!=CC->MaxDer)
    {
     float  sum = 0.0;
     for (i=0;i<KG->NASubst;i++)
         {
          float  V = (CC->Err[i] - CC->Y[i][CC->lmax-1]) / CC->Ymax[i];
          sum += V*V;
         }
     pr3 = 1.0f/(_powr(sum/CC->eup,0.5f/(CC->l+1.0f))*bias3+dfltzr);
    }
 /* 720 */
 pr2 = 1.0f/(_powr(d/CC->e,0.5f/CC->l)*bias2+dfltzr);
 pr1 = dfltzr;
 if (CC->Nq!=1)
    {
     float  sum = 0.0;
     for (i=0;i<KG->NASubst;i++)
         {
          float  V = CC->Y[i][CC->l-1] / CC->Ymax[i];
          sum += V*V;
         }
     pr1 = 1.0f/(_powr(sum/CC->edn,0.5f/CC->Nq)*bias1+dfltzr);
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

         *result = 0; return;
        }
     for (i=0;i<KG->NASubst;i++)
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

           *result = 0; return;
          }
      }
   else
     {
      CC->rh = pr2;
      if (CC->rh<thrshl) CC->idoub = idelay;
      else
        {
         Rescale(KG,CC);
         CC->rmax = rmxnor;
        }
      CC->hold   = CC->H;
      CC->jStart = CC->Nq;

      *result = 0; return;
     }
 /* 780 */
 Savel0 = CC->oldl0;
 SetNq(KG,CC,newq);
 CC->rc *= CC->El[0]/Savel0;
 /* 790 */
 Rescale(KG,CC);
 CC->rmax = rmxnor;
 /* 840 */
 CC->hold = CC->H;
 CC->jStart = CC->Nq;

 *result = 0;
}

__constant float  AdamsTable[MaxOrder][6] = {
   {1.0,  1.0,  0.0,  0.0,   0.0,  0.0  },
   {2.0,  1.0,  1.0,  0.0,   0.0,  0.0  },
   {12.0, 5.0,  8.0,  -1.0,  0.0,  0.0  },
   {24.0, 9.0,  19.0, -5.0,  1.0,  0.0  },
   {720.0,251.0,646.0,-264.0,106.0,-19.0}
};

__constant float AdamsTableB[MaxOrder][6] = {
   {1.0,  1.0,   0.0,    0.0,   0.0,    0.0  },
   {2.0,  3.0,   -1.0,   0.0,   0.0,    0.0  },
   {12.0, 23.0,  -16.0,  5.0,   0.0,    0.0  },
   {24.0, 55.0,  -59.0,  37.0,  -9.0,   0.0  },
   {720.0,1901.0,-2774.0,2616.0,-1274.0,251.0}
};

#ifndef __PARALLEL__
$ void Adams_Rozhkov(__global KineticGlobal * KG, __local KineticContext * CC)
#else
void Adams_Rozhkov(__global KineticGlobal * KG, __local KineticContext * CC)
#endif
{
 int    i,j,k;
 int    CurOrder;
 float  Xp[MaxActSubst];

 int    NP = 2;
 int    NQ = 2;
 for (i = 1; i < KG->AdamsOrder; i++) {
     NP = i*NQ + (i+2)*NQ;
     NQ *= 2;
 }

 float  AdamsH = CC->H/NP;
 float  HH = 0.0;

 int    np;

 for (i=0; i<KG->NASubst; i++)
     Xp[i] = CC->Y[i][0];
 KG->Fail = 0;

 NP = 0;

 for (CurOrder = 1; CurOrder<=KG->AdamsOrder && !KG->Fail && !CC->LossPrecision; CurOrder++)
 {
  for (np = 0; np < NP+2 && !KG->Fail && !CC->LossPrecision; np++)
  {
   int BashfortMode;

   for (BashfortMode = 1; BashfortMode >= 0 && !KG->Fail && !CC->LossPrecision; BashfortMode--)
   {
    float  F1 = AdamsH/(!BashfortMode ? AdamsTable[CurOrder-1][0] : AdamsTableB[CurOrder-1][0]);

    for (i=0;i<KG->NASubst;i++)
        {
         CC->FSave1[i] = 0;
         CC->FSave2[i] = 0;
        }
    /* Вычисляем значения правых частей системы ОДУ, учитывая,
       что для каждого j-вещества правая часть записывается
       как CC->FSave1[j]+Conc[j]*CC->FSave2[j] */
    for (i=0;i<KG->NReact;i++)
        {
         int Left  = KG->LR[i].NumLeft;
         int Right = KG->LR[i].NumRight;

         long Lj = 0; /* Номера веществ в левой части */
         /* Здесь W, WT содержит значение скорости реакции.
            Делить на AScale */
         int      Flags[MaxActSubst];
         __local float  * W  = CC->Ymax;
         float    WT = CC->C[i];

         for (j=0;j<Left;j++)
             {
              k    =  KG->LR[i].Indexes[j];
              Lj   |= (1<<k);
              W[k] =  CC->C[i];
              if (KG->LR[i].k[j]==0.0) {
                 WT *= CC->Y[k][np+CurOrder-BashfortMode];
              } else
                 WT *= _powr(CC->Y[k][np+CurOrder-BashfortMode],KG->LR[i].k[j]);
              Flags[k] = 0;
             }
         for (k=0;k<Left;k++)
             if (!Flags[KG->LR[i].Indexes[k]])
                {
                 for (j=0;j<Left;j++)
                   if (KG->LR[i].k[j]!=0.0) {
                      W[KG->LR[i].Indexes[k]] *= _powr(CC->Y[KG->LR[i].Indexes[j]][np+CurOrder-BashfortMode],KG->LR[i].k[j]-(j==k ? 1.0f : 0.0f));
                   } else if (j!=k)
                      W[KG->LR[i].Indexes[k]] *= CC->Y[KG->LR[i].Indexes[j]][np+CurOrder-BashfortMode];
                 Flags[KG->LR[i].Indexes[k]] = 1;
                }

         for (j=0;j<Right;j++)
             {
              k = KG->LR[i].Indexes[Left+j];
              if (Lj & (1<<k)) {
                 CC->FSave2[k] += W[k];
              } else
                 CC->FSave1[k] += WT;
             }
         for (j=0;j<Left;j++)
             {
              k = KG->LR[i].Indexes[j];
              CC->FSave2[k] -= W[k];
             }
        }

    for (i=0; i<KG->NASubst && !KG->Fail; i++)
        if (!BashfortMode) {
           float  HH = 1-F1*CC->FSave2[i]*AdamsTable[CurOrder-1][1];

           if (fabs(HH)>dfltzr)
              {
               float cnew = CC->FSave1[i]*AdamsTable[CurOrder-1][1];
               for (j=2; j<=CurOrder; j++)
                   cnew += AdamsTable[CurOrder-1][j]*CC->Y[i][np+CurOrder-j+1];
               CC->Y[i][np+CurOrder] = (CC->Y[i][np+CurOrder-1]+F1*cnew)/HH;
              }
           else
              {
               KG->Fail = 1;
              }
        } else {
           float  HH = -F1*CC->FSave2[i]*AdamsTableB[CurOrder-1][1];

           CC->Y[i][np+CurOrder] = CC->FSave1[i]*AdamsTableB[CurOrder-1][1];
           for (j=2; j<=CurOrder; j++)
               CC->Y[i][np+CurOrder] += AdamsTableB[CurOrder-1][j]*CC->Y[i][np+CurOrder-j];
           CC->Y[i][np+CurOrder] = CC->Y[i][np+CurOrder-1]*(1-HH)+F1*CC->Y[i][np+CurOrder];
        }
    if ((CC->Time + CC->H) == CC->Time || CC->H <= dfltzr || KG->Fail) CC->LossPrecision = 1;
   }
   HH += AdamsH;
  }
  if (!KG->Fail && !CC->LossPrecision)
     {
      for (i=0;i<KG->NASubst;i++) {
          for (j=1;j<=CurOrder;j++)
              CC->Y[i][j] = CC->Y[i][j*2];
      }
     }
  NP++;
  AdamsH *= 2;
 }

 if (!CC->LossPrecision)
    {
     if (!KG->Fail)
        {
         for (i=0;i<KG->NASubst;i++)
             CC->Y[i][0] = CC->Y[i][CurOrder-1];
         CC->Time += CC->H;
        }
     else
        for (i=0;i<KG->NASubst;i++)
            CC->Y[i][0] = Xp[i];
    }
 else {
    for (i=0;i<KG->NASubst;i++)
        CC->Y[i][0] = Xp[i];
 }
 CC->Iters += KG->AdamsOrder-1;
}

__constant float RosenTable[4][4] = {
   {-1.0f,      0.0f,        0.0f,       0.0f      },
   {1.0f/8.0f,  3.0f/8.0f,   1.0f,       0.0f      },
   {3.0f/8.0f,  19.0f/24.0f, -1.0f/6.0f, 0.0f      },
   {13.0f/6.0f, 1.0f/6.0f,   -2.0f,      2.0f/3.0f }
};

#ifndef __PARALLEL__
$ void Rosenbrock(__global KineticGlobal * KG, __local KineticContext * CC)
#else
void Rosenbrock(__global KineticGlobal * KG, __local KineticContext * CC)
#endif
{
 int i,j,k;

 KG->Fail = 0;

 Jacobian(KG,CC);
 for (i=0;i<KG->NASubst;i++)
   for (j=0;j<KG->NASubst;j++)
      CC->Df[i*MaxActSubst+j] *= -CC->H;
 for (i=0;i<KG->NASubst;i++)
     CC->Df[i*MaxActSubst+i] += 1.0;
 GetLU(KG,CC,CC->iRow);

 for (k=1; k<=4; k++)
     {
      for (i=0;i<KG->NASubst;i++)
          {
           CC->Y[i][k] = CC->Y[i][0];
           for (j=0; j<k-1; j++)
               CC->Y[i][k] += RosenTable[k-2][j] * CC->Y[i][j+1];
          }
          
      DifFun(KG,CC,(float  *) &CC->Y[0][k],7,CC->FSave2);
      for (i=0;i<KG->NASubst;i++)
          CC->FSave2[i] *= CC->H;
      SolveLU(KG,CC,CC->iRow,CC->FSave2,CC->FSave1);
      for (i=0;i<KG->NASubst;i++)
          CC->Y[i][k] = CC->FSave1[i];
     }

 for (i=0;i<KG->NASubst;i++)
     for (k=1; k<=4; k++)
         CC->Y[i][0] += RosenTable[3][k-1] * CC->Y[i][k];
 CC->Time += CC->H;
}

#ifndef __PARALLEL__
$ void DispatchF(int OtherMethod, __global KineticGlobal * KG, __local KineticContext * CC) {
#else
void DispatchF(int OtherMethod, __global KineticGlobal * KG, __local KineticContext * CC) {
#endif
 switch (OtherMethod) {
   case 0: Adams_Rozhkov(KG,CC); break;
   case 1: Rosenbrock(KG,CC); break;
 }
}

float  pow10(int Arg)
{
 float  result = 1;

 if (Arg>0)
   while((Arg--)>0) result*=10;
 else if (Arg<0)
   while((Arg++)<0) result/=10;

 return result;
}

#ifndef __PARALLEL__
$ void OneTaktKinetic(__global KineticGlobal * KG, __local KineticContext * CC,int UseStiffGearMethod, int OtherMethod)
#else
void OneTaktKinetic(__global KineticGlobal * KG, __local KineticContext * CC,int UseStiffGearMethod, int OtherMethod)
#endif
{
 Derivs Y[MaxActSubst];
 float  FSave1[MaxActSubst]; /* [KG->NASubst] */
 float  FSave2[MaxActSubst]; /* [KG->NASubst] */
 float  Df[MaxActSubst*MaxActSubst];

 int i,j;

 CC->Y = Y;
 CC->FSave1 = FSave1;
 CC->FSave2 = FSave2;
 CC->Df = Df;

 CC->Tem = _R_*CC->Tk;

 for (i=0;i<KG->NASubst;i++)
     if (CC->Conc0[i]<0.0) CC->Conc0[i] = 0.0;

 CC->Iters = 0;
 if (CC->Tem>dfltzr)
    {
     /* Автоматический выбор масштабного множителя */
     float  maxY = dfltzr;

     for (i=0;i<KG->NASubst;i++)
         if (CC->Conc0[i]>maxY) maxY = CC->Conc0[i];

     if (maxY>dfltzr)
        {
         CC->LgAScale = (int) -log10(maxY);
         CC->AScale   = pow10(CC->LgAScale);

         /* Вычисление коэффициентов скорости */
         for (i=0;i<KG->NReact;i++)
             {
              float  PEi = (KG->_A[i]-(KG->LR[i].NumLeft-1)*CC->LgAScale)/M_LOG10E;
              float  Degree = PEi-KG->EA[i]/CC->Tem;

              if (Degree>-80.0f) CC->C[i] = expr(Degree)*_powr(CC->Tk,(float)KG->Tn[i]);
              else CC->C[i] = 0.0;
             }

         if (UseStiffGearMethod)
            {
             float  EpsCom = 1.0E-12f/KG->CalcEps;
             /* STIFFC */
             for (i=0;i<KG->NASubst;i++)
                 {
                  CC->Ymax[i]=0.0;
                  CC->Y[i][0]=CC->Conc0[i]*CC->AScale;
                  for (j=1;j<6;j++)
                      CC->Y[i][j]=0.0;
                 }

             CC->Time = KG->InitTime;

             if (CC->_EndTime<KG->InitH) CC->H=CC->_EndTime;
             else CC->H = KG->InitH;

             CC->LossPrecision = (CC->Time+CC->H)==CC->Time;

             CC->MaxDer = 5;

             CC->jStart = 0;

             for (CC->Iters = 0;CC->Time<=CC->_EndTime && CC->Iters<KG->MaxIterations && !CC->LossPrecision;CC->Iters++)
              {
               float  ys = 0.0;
               int    resStiff = 0;

               for (i=0;i<KG->NASubst;i++)
                   ys += fabsr(CC->Y[i][0]);
               ys *= EpsCom;
               for (i=0;i<KG->NASubst;i++)
                   CC->Ymax[i] = amax3(CC->Y[i][0],ys,dfltzr);

               Stiff(KG,CC,&resStiff);
               if (resStiff)
                  {
                   CC->jStart=-1;
                   CC->H/=10.0;
                  }
               if ((CC->Time+CC->H)==CC->Time) CC->LossPrecision = 1;
              }
             if (CC->Iters==KG->MaxIterations) CC->LossPrecision = 1;

             if (CC->LossPrecision)
                {
                 CC->KinErrorInfo.LossH    = CC->Iters!=KG->MaxIterations;
                 CC->KinErrorInfo.LastH    = CC->H;
                 CC->KinErrorInfo.ReachTau = CC->Time;
                 CC->KinErrorInfo.Tmp      = (float) (CC->Tk-ZeroK);
                }
             else
                for (i=0;i<KG->NASubst;i++)
                    {
                     float  sum = 0.0;
                     for (j=0;j<=CC->jStart;j++)
                         sum += CC->Y[i][j]*_powr((CC->_EndTime-CC->Time)/CC->H,j*1.0f);
                     CC->Conc1[i] = (float)(sum/CC->AScale);
                    }
            }
         else
            {
             float  _EndTime = CC->_EndTime;
             float  _H;

#define max_history 5
#define period_of_check 2

#define START_HISTORY(n) \
  n_of_history = n; \
  history[n]._start = CC->Time; \
  history[n].order = ORDER; \
  history[n].error = KG->CalcEps; \
  for (i=0;i<KG->NASubst;i++) \
      history[n].conc0[i] = CC->Y[i][0];

#define RESTORE_HISTORY(n) \
  n_of_history = n; \
  CC->Time = history[n]._start; \
  KG->AdamsOrder = history[n].order; \
  KG->CalcEps = history[n].error; \
  for (i=0;i<KG->NASubst;i++) \
      CC->Y[i][0] = history[n].conc0[i];

#define STORE_HISTORY(err) \
  n_of_history++; \
  if (n_of_history == max_history) { \
     int kr; \
     for (kr = 0; kr < n_of_history-1; kr++) \
         history[kr] = history[kr+1]; \
     n_of_history--; \
  } \
  history[n_of_history].error = err; \
  START_HISTORY(n_of_history);

#define CHECK_REORDER() \
  if (ORDER < MaxOrder && OtherMethod == Adams_Rozhkov_Method) { \
     ORDER++; \
     KG->AdamsOrder = ORDER; \
     CC->H *= 0.25; \
     if (CC->Time + CC->H == CC->Time || CC->H <= dfltzr) { \
        CC->LossPrecision = (CC->_EndTime - CC->Time) > 20*CC->H; \
        break; \
     } \
     RESTORE_HISTORY(0); \
     CC->LossPrecision = 0; \
     checker = 0; \
     continue; \
  }

             typedef struct {
                float _start;
                float conc0[MaxActSubst];
                float error;
                int    order;
             } history_item;

             CC->Iters = 0;

             for (i=0;i<KG->NASubst;i++)
                 CC->Y[i][0] = CC->Conc0[i]*CC->AScale;

             int ORDER = 1;
             switch (OtherMethod) {
               case Adams_Rozhkov_Method : ORDER = KG->AdamsOrder; break;
               case Rosenbrock_Method :    ORDER = 4; break;
             }

             _H = (CC->_EndTime-KG->InitTime)/(256/(1 << (ORDER-1)));

             CC->Time = KG->InitTime;

             CC->_EndTime = _EndTime;
             CC->H        = _H;
             CC->LossPrecision = 0;
             KG->Fail = 0;

             history_item history[max_history];
             int n_of_history;
             int checker = 0;

             START_HISTORY(0);

             while (CC->Time < CC->_EndTime && !KG->Fail && !CC->LossPrecision) {
                if (checker == 0) {
                   checker = period_of_check - 1;
                   // Checked Work. Loop with check until loss h. If loss -> back history!
                   float startings[MaxActSubst];
                   float error = 0.0;
                   for (i=0;i<KG->NASubst;i++) {
                       startings[i] = CC->Y[i][0];
                   }
                   CC->H *= 0.5;
                   DispatchF(OtherMethod, KG, CC);
                   DispatchF(OtherMethod, KG, CC);
                   CC->H *= 2;
                   if (CC->LossPrecision) {
                      CHECK_REORDER();
                      error = 8.0 * KG->CalcEps;
                   } else {
                      for (i=0;i<KG->NASubst;i++) {
                          float maxY = fabsr(CC->Y[i][0]) > fabsr(startings[i]) ? fabsr(CC->Y[i][0]) : fabsr(startings[i]);
                          error += fabsr(CC->Y[i][0] - startings[i]) / (1 << ORDER - 1) / maxY;
                      }
                      error /= KG->NASubst;
                   }
                   if (error >= KG->CalcEps) {
                      // BIG ERROR
                      int delta = n_of_history > MaxOrder ? MaxOrder : n_of_history; // n_of_history * 0.8;
                      // Back Extrapolation
                      int min_i = n_of_history == 0 ? 0 : delta;
                      float newh = CC->H * 0.75;

                      if (error >= 20*KG->CalcEps) {
                         CHECK_REORDER();
                      }

                      if (n_of_history) {
                         float min_err = error;

                         while (delta >= 1) {
                            float err = error * AdamsTableB[delta-1][0];
                            for (i=1; i<delta; i++)
                                err -= AdamsTableB[delta-1][delta-i+1]*history[n_of_history - i + 1].error;
                            err /= AdamsTableB[delta-1][1];
                            if (err > 0.0) {
                               float C0 = error / _powr(CC->H, ORDER);
                               newh = amin2(newh, _powr(err / C0, 1.0 / ORDER));
                            }
                            if (err < min_err) {
                               min_err = err;
                               min_i = delta;
                            }
                            delta--;
                         }
                      }
                      int r = n_of_history - min_i;
                      RESTORE_HISTORY(r);
                      CC->H = amax2(CC->H * 0.55, newh);
                      if (CC->Time + CC->H == CC->Time || CC->H <= dfltzr) {
                         CC->LossPrecision = (CC->_EndTime - CC->Time) > 20*CC->H;
                         break;
                      }
                      START_HISTORY(0);
                      CC->LossPrecision = 0;
                      checker = 0;
                      continue;
                   }
                   else if (2*error < KG->CalcEps) {
                      // BIG ERROR
                      int delta = n_of_history > MaxOrder ? MaxOrder : n_of_history; // n_of_history * 0.8;
                      // Back Extrapolation
                      float newh = 1.25 * CC->H;
                      if (n_of_history) {
                         while (delta >= 1) {
                            float err = error * AdamsTableB[delta-1][0];
                            for (i=1; i<delta; i++)
                                err -= AdamsTableB[delta-1][delta-i+1]*history[n_of_history - i + 1].error;
                            err /= AdamsTableB[delta-1][1];
                            if (err > 0.0) {
                               float C0 = error / _powr(CC->H, ORDER);
                               newh = amax2(newh, _powr(err / C0, 1.0 / ORDER));
                            }
                            delta--;
                         }
                         if (newh / CC->H > 10) newh = 10 * CC->H;
                      }
                      CC->H = newh;
                      if (CC->Time + CC->H == CC->Time || CC->H <= dfltzr) {
                         CC->LossPrecision = (CC->_EndTime - CC->Time) > 20*CC->H;
                         break;
                      }
                      START_HISTORY(0);
                      CC->LossPrecision = 0;
                      checker = 0;
                      continue;
                   }
                   else {
                      STORE_HISTORY(error);
                   }
                } else {
                   DispatchF(OtherMethod, KG, CC);
                   checker--;
                   if (CC->LossPrecision) {
                      int r = n_of_history ? n_of_history-1 : 0;
                      RESTORE_HISTORY(r);
                      CC->LossPrecision = 0;
                      checker = 0;
                      continue;
                   }
                }
                // NOT Fail, NOT LossPrecision
                if (CC->_EndTime - CC->Time < CC->H) {
                   CC->H = CC->_EndTime - CC->Time;

                   START_HISTORY(0);
                }
                if ((CC->Time + CC->H) == CC->Time || CC->H <= dfltzr)
                   CC->Time = CC->_EndTime;
             }

             if (CC->Iters==KG->MaxIterations) CC->LossPrecision = 1;

             if (CC->LossPrecision)
                {
                 CC->KinErrorInfo.LossH    = CC->Iters!=KG->MaxIterations;
                 CC->KinErrorInfo.LastH    = CC->H;
                 CC->KinErrorInfo.ReachTau = CC->Time;
                 CC->KinErrorInfo.Tmp      = (float)(CC->Tk-ZeroK);
                 for (i=0;i<KG->NASubst;i++)
                     CC->Conc1[i] = CC->Conc0[i];
                }
             else
                for (i=0;i<KG->NASubst;i++)
                    CC->Conc1[i] = (float)(CC->Y[i][0]/CC->AScale);
            }
        }
     else
        for (i=0;i<KG->NASubst;i++)
            CC->Conc1[i] = CC->Conc0[i];
    }
 else
   for (i=0;i<KG->NASubst;i++)
       CC->Conc1[i] = CC->Conc0[i];
}

#pragma plan common end

#endif
