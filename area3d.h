#ifndef __AREA3D__
#define __AREA3D__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef __USE_VCL__
#include <stdarg.h>
#include <Dialogs.hpp>
#endif

#define nDims 3

const char * Axis[3] = {"x","y","z"};

#define IBM_3D_SIGN 2    /* Сигнатура IBM PC */

#define MaxNX 512

typedef struct {
 unsigned char * Map;
 float * Vals;
 int     ValNum;
} MapItem;

int    NX;      /* Число узлов по оси X */
int    NY;      /* Число узлов по оси Y */
int    NZ;      /* Число узлов по оси Z */

double HXmin;/* Минимальный шаг по оси X  */
double HXmax;/* Максимальный шаг по оси X */
int    HXreg;
double HYmin;/* Минимальный шаг по оси Y  */
double HYmax;/* Максимальный шаг по оси Y */
int    HYreg;
double HZmin;/* Минимальный шаг по оси Z  */
double HZmax;/* Максимальный шаг по оси Z */
int    HZreg;

#undef _X
#undef _Y
#undef _Z

float * _X;
float * _Y;
float * _Z;

double ModelTime;
double TAU;  /* Шаг по времени */
double TAU0;
double FastTAU = 0.0; /* Шаг по времени для быстрой подсистемы */
int    FastTAUDivider = 1; /* Количество быстрых итераций в одной медленной */
int    FastMode;
int    SlowMode;
int    FastIteration;

#ifndef __USE_VCL__
_int takt = 0; /* Счетчик итераций */
#endif

int    nPoissons = 0;

int    nSMP = 1;

int    UseSpecialBounds = 1;

typedef struct {
 float * WX1;
 float * WY1;
 float * WZ1;
 float * WX2;
 float * WY2;
 float * WZ2;
} WKoeffs;

#define MaxMDLength 40

char   ModelDate[MaxMDLength] = "";
double ModelGMTHour;
double ModelAltitude;
double ModelLongitude;
double ModelAngle;

double SolarDeclination = 0.0; /* Склонение Солнца на дату моделирования */

double SolarX = 0.0; /* Компоненты вектора падения солнечных лучей */
double SolarY = 0.0;
double SolarZ = -1.0;

/* Граничные условия */
typedef struct {
   signed int X1;
   signed int X2;
} DescLineItem;

typedef struct {
   unsigned char  Num;
   unsigned char   Dir1;
   unsigned char   Dir2;
   DescLineItem * Items;
} DescLine;

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

/* We need this for some OpenMP compilers... */
int imax(int a, int b)
{
 return max(a,b);
}

#define MinIndexBound1 'A'
#define Zero           '0'
#define Block          '1'
#define Bound1         '1'
#define Bound2         '2'
#define Bound3         '3'
#define BoundClosed    '9'

#define AreaBlock      (Block-Zero)
#define AreaFree       ('8'-Zero)
/* Специальный тип граничного условия - необходимо предварительное   */
/* вычисление по схеме Головичева, а затем коррекция (прием значения */
/* из соседнего процессора)                                          */
#define ExchngBound 0x7F

#define msBoundary 0xFF000000
#define shBoundary 12

#define fld2RF 0x01

#define XP ((x+1)   % NX)
#define XM ((x+NXs) % NX)
#define YP ((y+1)   % NY)
#define YM ((y+NYs) % NY)
#define ZP ((z+1)   % NZ)
#define ZM ((z+NZs) % NZ)

#define ZYXP (x<NXs ? Ptr+1 : z*NY*NX+y*NX+XP)
#define ZYXM (x>0 ? Ptr-1 : z*NY*NX+y*NX+XM)
#define ZYPX (y<NYs ? Ptr+NX : z*NY*NX+YP*NX+x)
#define ZYMX (y>0 ? Ptr-NX : z*NY*NX+YM*NX+x)
#define ZPYX (z<NZs ? Ptr+NY*NX : ZP*NY*NX+y*NX+x)
#define ZMYX (z>0 ? Ptr-NY*NX : ZM*NY*NX+y*NX+x)

int BoardSize;
int ScKfSize;

int NXs;
int NXs1;
int NXs2;
int NYs;
int NYs1;
int NYs2;
int NZs;
int NZs1;
int NZs2;

typedef double (*_InitVal) (int i, int x, int y, int z, unsigned char Map);

#ifdef __USE_VCL__
int ErrorMessage(const char * fmt, ...)
{
 int     Result;
 char    Buf[200];
 va_list Args;

 va_start(Args,fmt);
 Result = vsprintf(Buf,fmt,Args);
 MessageDlg(Buf,mtError,TMsgDlgButtons()<<mbOK,0);
 va_end(Args);

 return Result;
}
#else
#define ErrorMessage printf
#endif

void AbortServer(int Code)
{
#ifdef __MPI__
 MPI_Abort(MPI_COMM_WORLD,Code);
#else
#ifdef __MVS__
 ErrorMessage("exit(%i)\n",Code);
#endif
 exit(Code);
#endif
}

void * SafeMalloc(size_t size)
{
 void * Result = malloc(size);

 if (size && !Result)
    {
     ErrorMessage("There is not enough memory for malloc\n");
     AbortServer(10);
    }

 return Result;
}

void * SafeRealloc(void * block, size_t size)
{
 void * Result = realloc(block, size);

 if (size && !Result)
    {
     ErrorMessage("There is not enough memory for realloc\n");
     AbortServer(10);
    }

 return Result;
}

#define __MAXPATH  81

void ChangeExt(char * DestName, const char * SrcName, const char * NewExt)
{
 char * BeginExt = strrchr(strcpy(DestName,SrcName),'.');

 if (BeginExt)
    strcpy(BeginExt,NewExt);
 else
    strcat(DestName,NewExt);
}

unsigned char * ReadMap(int NVals, int SubClass, const char * CfgName, const char * MapBase, float ** Vals, int * ValNum)
{
 unsigned char * Result = (unsigned char *) SafeMalloc(NX*NY*NZ*sizeof(char));

 char MapName[__MAXPATH];
 char ValName[__MAXPATH];

 unsigned char MaxAlfa = 0;
 int  z,y,x,Ptr;

 FILE * MapFile;

 ChangeExt(MapName,CfgName,MapBase);
 if (Vals)
    {
     strcat(MapName,"m");
     ChangeExt(ValName,CfgName,MapBase);
     if (!NVals) strcat(ValName,"v");
    }

 if (!(MapFile = fopen(MapName,"rt")))
    {
     ErrorMessage("File '%s' does not exist\n",MapName);
     AbortServer(1);
    }

 for (z=0; z<NZ; z++)
     for (y=0; y<NY; y++)
         {
          char Buf[MaxNX+1] = "";

          if (!fscanf(MapFile,"%s\n",Buf) || strlen(Buf)!=NX)
             {
              ErrorMessage("Error in '%s' line (z=%i):(y=%i)\n",MapName,z,y);
              AbortServer(1);
             }
          memmove(&Result[(NZs-z)*NY*NX+(NYs-y)*NX],Buf,NX*sizeof(char)); /* NYs-y,NZs-z !!!!! */
          for (x=0; x<NX; x++)
              if (Buf[x]>MaxAlfa)
                 MaxAlfa = Buf[x];
         }
 fclose(MapFile);

 if (!Vals)
    for (Ptr = 0; Ptr<NX*NY*NZ; Ptr++)
        Result[Ptr] -= Zero;
 
 if (MaxAlfa>=MinIndexBound1 && Vals)
    {
     unsigned char Count = (MaxAlfa-MinIndexBound1+1);
     int NN = NVals==0 ? 1 : NVals;
     int CurVal;

     (*Vals) = (float *) SafeMalloc(NN*Count*sizeof(float));

     if (ValNum)
        (*ValNum) = Count;

     for (CurVal = 0; CurVal<NN; CurVal++)
         {
          char   FName[127];
          FILE * ValFile;

          if (!NVals)
             if (SubClass<0)
                strcpy(FName,ValName);
             else
                sprintf(FName,"%s%i",ValName,SubClass);
          else
             sprintf(FName,"%s%i",ValName,CurVal);

          ValFile = fopen(FName,"rt");

          if (!ValFile)
             if (NVals)
                {
                 printf("Warning: file '%s' does not exist => Values = 0.0 by default\n",FName);
                 for (x=0; x<Count; x++)
                     (*Vals)[CurVal*Count+x] = 0.0;
                }
             else
             {
              ErrorMessage("File '%s' does not exist\n",FName);
              AbortServer(1);
             }
          else
            {
             for (x=0; x<Count; x++)
                 if (!fscanf(ValFile,"%f",&((*Vals)[CurVal*Count+x])))
                    {
                     ErrorMessage("Error in '%s' item %i\n",ValName,x);
                     AbortServer(1);
                    }
             fclose(ValFile);
            }
         }
    }
 else
   {
    if (ValNum)
       (*ValNum) = 0;
    if (Vals)
       (*Vals) = NULL;
   }

 return Result;
}

void InitBoard(float * H, unsigned char * HMap, float * HVals, float InitValue)
{
 int z,y,x;
 int Ptr = 0;

 for (z=0; z<NZ; z++)
     for (y=0; y<NY; y++)
         for (x=0; x<NX; Ptr++,x++)
             if (HMap[Ptr]>=MinIndexBound1)
                H[Ptr] = HVals[HMap[Ptr]-MinIndexBound1];
             else
                H[Ptr] = InitValue;
}

void InitBoardF(int SubClass, float * H, unsigned char * HMap, float * HVals, _InitVal InitValue)
{
 int z,y,x;
 int Ptr = 0;

 for (z=0; z<NZ; z++)
     for (y=0; y<NY; y++)
         for (x=0; x<NX; Ptr++,x++)
             if (HMap[Ptr]>=MinIndexBound1)
                H[Ptr] = HVals[HMap[Ptr]-MinIndexBound1];
             else
                H[Ptr] = (float) InitValue(SubClass,x,y,z,HMap[Ptr]);
}

void AddLine(DescLineItem * Buf, int * Count, int X1, int X2)
{
 Buf[*Count].X1 = X1;
 Buf[*Count].X2 = X2;
 (*Count)++;
}

#define ExchngBit 0x80
#define TopBit    0x40
#define BottomBit 0x20
#define LeftBit   0x10
#define RightBit  0x08
#define ForwBit   0x04
#define BackBit   0x02

enum {_back = 0, _forw, _right, _left, _bottom, _top};

DescLine * CreateDescriptor(unsigned char * Area, unsigned char * OrgnArea,
                                                  int NumZ, int NumY, int NumX,
                                                  int ZStep, int YStep, int XStep,
                                                  char Mask1, char Mask2,
                                                  unsigned char * Boundaries)
{
 DescLine * Result = (DescLine *) SafeMalloc(NumZ*NumY*sizeof(DescLine));

 int z,y;

 for (z=0; z<NumZ; z++)
     for (y=0; y<NumY; y++)
         {
          DescLineItem Buf[MaxNX];
          int          Count = 0;
          int          Base = z*ZStep+y*YStep;

          int          XL = 0;

          for ( ; XL<NumX && Area[Base+XL*XStep]==0 ; XL++);

          if (XL==NumX)
             {
              Buf[Count].X1 = -1;
              Buf[Count].X2 = -1;
              Count++;
             }
          else
            {
             int X0 = XL;

             do
               {
                for ( ; X0<NumX && Area[Base+X0*XStep]!=0 ; X0++);

                if (X0==NumX)
                   {
                    if (XL>0)
                       AddLine(Buf,&Count,0,XL-1);
                   }
                else
                   {
                    int X1 = X0;

                    for ( ; X1<NumX && Area[Base+X1*XStep]==0 ; X1++);

                    if (X1==NumX)
                       if (!XL)
                          AddLine(Buf,&Count,X0,NumX-1);
                       else
                          AddLine(Buf,&Count,X0,XL-1);
                    else
                       AddLine(Buf,&Count,X0,X1-1);

                    X0 = X1;
                   }
               }
             while (X0!=NumX);
            }

          Result[z*NumY+y].Num  = Count;
          Result[z*NumY+y].Dir1 = _back;
          if (Mask1)
             while (Mask1!=(1<<(Result[z*NumY+y].Dir1+1)))
                Result[z*NumY+y].Dir1++;
          Result[z*NumY+y].Dir2 = _back;
          if (Mask2)
             while (Mask2!=(1<<(Result[z*NumY+y].Dir2+1)))
                Result[z*NumY+y].Dir2++;
          if (Count)
             {
              Result[z*NumY+y].Items = (DescLineItem *) SafeMalloc(Count*sizeof(DescLineItem));
              memmove(Result[z*NumY+y].Items,Buf,Count*sizeof(DescLineItem));
              while (--Count >= 0)
                if (Buf[Count].X1>=0 && Buf[Count].X2>=0)
                   {
                    Area[Base+((Buf[Count].X1+NumX-1)%NumX)*XStep] |= Mask1;
                    if (Boundaries && OrgnArea && OrgnArea[Base+((Buf[Count].X1+NumX-1)%NumX)*XStep]!=AreaFree)
                       Boundaries[Base+Buf[Count].X1*XStep] = 1;
                    Area[Base+((Buf[Count].X2+1)%NumX)*XStep] |= Mask2;
                    if (Boundaries && OrgnArea && OrgnArea[Base+((Buf[Count].X2+1)%NumX)*XStep]!=AreaFree)
                       Boundaries[Base+Buf[Count].X2*XStep] = 1;
                   }
             }
          else
             Result[z*NumY+y].Items = NULL;
         }

 return Result;
}

#define Gk(k,l) (exp(-(1.0*k/(l-1)*30.0-15.0)*(1.0*k/(l-1)*30.0-15.0)/30.0))
#define HXfun(k,l) ((Gk(k,l)-G0)*GZ+Hmin)

int ArrangeGrid(unsigned char * OrgnArea,
                 float * HX, float Hmin, float Hmax,
                 int NumZ,   int NumY,   int NumX,
                 int ZStep,  int YStep,  int XStep,
                 DescLine * Desc,
                 double * _HX)
{
 int    Item;
 int    z,y,x;
 int    Result = 1;

 unsigned char * MacroMap = (unsigned char *) SafeMalloc(NumX*sizeof(char));
 DescLine * MacroDesc;

 for (x=0; x<NumX && _HX; x++)
     if (_HX[x]==0.0) _HX = NULL;
 for (x=0; x<NumX; x++)
     HX[x] = (float)(_HX ? _HX[x] : Hmin);
 if (!_HX)
    {
     /* Инициализируем макрокарту */
     memset(MacroMap,0,NumX);
     /* Создаем макрокарту характерных узлов сетки */
     for (z=0; z<NumZ; z++)
         for (y=0; y<NumY; y++)
             for (Item = 0; Item<Desc[z*NumY+y].Num; Item++)
                 {
                  int X1 = Desc[z*NumY+y].Items[Item].X1;
                  unsigned char B1 = OrgnArea[z*ZStep+y*YStep+((X1+NumX-1)%NumX)*XStep];
                  int X2 = Desc[z*NumY+y].Items[Item].X2;
                  unsigned char B2 = OrgnArea[z*ZStep+y*YStep+((X2+1)%NumX)*XStep];

                  if (X1>=0 && X2>=0)
                     {
                      MacroMap[(X1+NumX-1) % NumX] = B1==AreaFree ? 0 : 1;
                      MacroMap[(X2+1) % NumX]      = B2==AreaFree ? 0 : 1;
                     }
                 }
     MacroDesc = CreateDescriptor(MacroMap,NULL,1,1,NumX,0,0,1,0,0,NULL);
     for (Item = 0; Item<MacroDesc[0].Num; Item++)
         {
          int X1 = MacroDesc[0].Items[Item].X1;
          int X2 = MacroDesc[0].Items[Item].X2;

          if (X1>=0 && X2>=0)
             {
              int Long = (X1<=X2 ? X2-X1+1 : NumX-(X1-X2-1))+1;

              double G0   = Gk(0.0,(double) Long);
              double Gmax = Gk(1.0*(Long-1.0)/2.0,(double) Long);

              double GZ   = (double) (Hmax-Hmin)/(Gmax-G0);

              for (x=0; x<Long; x++)
                  HX[(X1+x+NumX-1) % NumX] = (float) HXfun((double) x,(double) Long);
             }
         }
     if (MacroDesc[0].Items)
        free(MacroDesc[0].Items);
     free(MacroDesc);
     free(MacroMap);
    }
 for (x=1; x<NumX && Result; x++)
     if (fabs(HX[x]-HX[x-1])>1E-7)
        Result = 0;
 return Result;
}

void CalculateDists(int IsCyclic, int X, int N, float * H, float * L)
{
 float Lp,Ln;
 int   i;

 for (i=0; i<N; i++)
     L[i] = 1E6;

 L[X] = 0.0;

 for (i=1, Lp = 0.0; i<N-(IsCyclic!=0 ? 0 : X); i++)
     {
      Ln = Lp+H[(X+i-1)%N];
      Lp = Ln;

      L[(X+i)%N] = Ln;
     }
 for (i=1+(IsCyclic!=0 ? 0 : N), Lp = 0.0; i<N+(IsCyclic!=0 ? 0 : X+1); i++)
     {
      Ln = Lp+H[(X+N-i)%N];
      Lp = Ln;

      if (Ln<L[(X+N-i)%N])
         L[(X+N-i)%N] = Ln;
     }
}

void CalculateLmin2(unsigned char  * Area,
                    float * HZ, float * HY, float * HX,
                    float * Lmin2)
{
 float * LZ = (float *) SafeMalloc(NZ*sizeof(float));
 float * LY = (float *) SafeMalloc(NY*sizeof(float));
 float * LX = (float *) SafeMalloc(NX*sizeof(float));

 float L;

 int Z,Y,X;
 int z,y,x;
 int PTR = 0;
 int Ptr;

 int IsCyclicX = 0;
 int IsCyclicY = 0;
 int IsCyclicZ = 0;

 Ptr = 0;
 for (z=0; z<NZ && !IsCyclicX; z++)
     for (y=0; y<NY && !IsCyclicX; y++)
         for (x=0, IsCyclicX=1; x<NX; x++, Ptr++)
             if (Area[Ptr])
                IsCyclicX = 0;

 for (z=0; z<NZ && !IsCyclicY; z++)
     for (x=0; x<NX && !IsCyclicY; x++)
         for (y=0, IsCyclicY=1; y<NY; y++)
             if (Area[z*NY*NX+y*NX+x])
                IsCyclicY = 0;

 for (y=0; y<NY && !IsCyclicZ; y++)
     for (x=0; x<NX && !IsCyclicZ; x++)
         for (z=0, IsCyclicZ=1; z<NZ; z++)
             if (Area[z*NY*NX+y*NX+x])
                IsCyclicZ = 0;

 for (Z=0; Z<NZ; Z++)
     {
#if !defined(__MVS__) && !defined(__UNIX__) && !defined(__USE_VCL__)
      printf("Creating Lmin Table : %i%%\r",(int) (((double) Z/NZs)*100));
#endif

      CalculateDists(IsCyclicZ,Z,NZ,HZ,LZ);
      #pragma omp for schedule(dynamic,imax(1,NY/10)) private(Y,X,z,y,x,Ptr,PTR,L)
      for (Y=0; Y<NY; Y++)
          {
           CalculateDists(IsCyclicY,Y,NY,HY,LY);
           for (X=0, Ptr=(Y+Z*NY)*NX; X<NX; Ptr++,X++)
               if (Area[Ptr]==AreaBlock)
                  Lmin2[Ptr] = 1.0;
               else
                  {
                   CalculateDists(IsCyclicX,X,NX,HX,LX);

                   Lmin2[Ptr] = 1E6;

                   PTR = 0;

                   for (z=0; z<NZ; z++)
                       for (y=0; y<NY; y++)
                           for (x=0; x<NX; PTR++, x++)
                               if (Area[PTR]==AreaBlock)
                                  {
                                   L = LZ[z]*LZ[z]+LY[y]*LY[y]+LX[x]*LX[x];
                                   if (L<Lmin2[Ptr])
                                      Lmin2[Ptr] = L;
                                  }
                  }
          }
     }
#ifndef __USE_VCL__
 printf("Creating Lmin Table : 100%%\n");
#endif

 free(LZ);
 free(LY);
 free(LX);
}

char _ResBase[__MAXPATH];
char DResBase[__MAXPATH];
char SavName[__MAXPATH];
char DscName[__MAXPATH];

void CreateNames(char * CfgName)
{
 ChangeExt(_ResBase,CfgName,".");
 ChangeExt(DResBase,CfgName,".d");
 ChangeExt(SavName,CfgName,".sav");
 ChangeExt(DscName,CfgName,".dsc");
}

void GetSolarAngle(double Hours)
{
 double Pi = 3.1415926535897932;
 double Phi = ModelAltitude/180.0*Pi;
 double Lambda = ModelLongitude/180.0*Pi;
 double Alpha = ModelAngle/180.0*Pi;
 double DELTA = (12.0-ModelGMTHour-Hours)*Pi/12.0-Lambda;
 double sinDC = sin(SolarDeclination);
 double cosDC = cos(SolarDeclination);
 double sinPhi = sin(Phi);
 double cosPhi = cos(Phi);
 double cosAlpha = cos(Alpha);
 double sinAlpha = sin(Alpha);
 double sinH = sinDC*sinPhi+cosDC*cosPhi*cos(DELTA);
 double nx = (sinDC-sinH*sinPhi)/cosPhi;
 double ny = sin(DELTA)*cosDC;
 SolarX = -(nx*cosAlpha-ny*sinAlpha);
 SolarY = nx*sinAlpha+ny*cosAlpha;
 SolarZ = -sinH;
}

void InitSolar(void)
{
 double Pi = 3.1415926535897932;
 char Days[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
 char Buf[MaxMDLength];
 char * Token = strtok(strcpy(Buf,ModelDate),"./\\");
 int Day = Token ? atoi(Token) : 1;
 int Month = (Token = strtok(NULL,"./\\")) ? atoi(Token) : 1;
 int Year = (Token = strtok(NULL,"./\\")) ? atoi(Token) : 1900;
 int PassedDays;
 int i;

 if (Day<1) Day = 1;
 if (Month<1) Month = 1;
 if (Year<1900) Year = 1900;

 if (Day<22 && Month<=3)
    if (Month==3)
       PassedDays = 365-(22-Day);
    else
      {
       PassedDays = 22+Days[Month]-Day;
       for (i = Month+1; i<3; i++)
           PassedDays += Days[i];
       PassedDays = 365-PassedDays;
      }
 else
    if (Month==3)
       PassedDays = Day-22;
    else
      {
       PassedDays = Days[3]-22+Day;
       for (i = 4; i<Month; i++)
           PassedDays += Days[i];
      }
 SolarDeclination = 23.5/180.0*Pi*sin(2*Pi*PassedDays/365.0);

 GetSolarAngle(0.0);
}

void InitParameters(void)
{
 int  NYX;

 BoardSize = NZ*NY*NX*sizeof(float);

 NYX       = (NY>NX) ? NY : NX;
 ScKfSize  = (( (NYX>NZ) ? NYX : NZ )+2)*sizeof(double);
#ifdef __OPENMP__
 ScKfSize *= nSMP;
#endif

 NXs   = NX-1;
 NXs1  = NX-2;
 NXs2  = NX-3;
 NYs   = NY-1;
 NYs1  = NY-2;
 NYs2  = NY-3;
 NZs   = NZ-1;
 NZs1  = NZ-2;
 NZs2  = NZ-3;

 InitSolar();
}

#define SpecVAL  1E25
#define FictZero 1E24

float * VAL(float ** Bounds, char Dir, float * H)
{
 static float VAL_ZERO = 0.0;

 if (*H==FictZero)
    return &VAL_ZERO;
 else if (fabs(*H)<SpecVAL)
    return H;
 else
    {
     int POW;
     float Mantissa = frexp(*H,&POW);
     POW -= 100;
     if (POW<0 || POW>9)
        return H;
     else
        {
         int LOW = Mantissa>0 ? ((int)(((double) (Mantissa*1000000.0))+0.5)) : ((int)(((double) (-Mantissa*1000000.0))+0.5))-500000;
         return &(Bounds[Dir][POW*1000000+LOW-1]);
        }
    }
}

/* Рабочее поле турбулентной вязкости */
float *  Dt = NULL;
/* Поле квадратов кратчайших расстояний до твердой стенки */
float *  Lmin2 = NULL;
float *  DIV = NULL;
float *  DivSolar = NULL;

float * HX  = NULL;
float * HY  = NULL;
float * HZ  = NULL;
float * HXX = NULL;
float * HYY = NULL;
float * HZZ = NULL;
float * HXS = NULL;
float * HYS = NULL;
float * HZS = NULL;
float * HX2 = NULL;
float * HY2 = NULL;
float * HZ2 = NULL;
float * HXP = NULL;
float * HYP = NULL;
float * HZP = NULL;

typedef struct {
  float R1p, R2p, R3p, R1m, R2m, R3m, h;
} HKoeffs;

HKoeffs * HXR = NULL;
HKoeffs * HYR = NULL;
HKoeffs * HZR = NULL;

#ifdef __PARALLEL__
  float HZbf  = 0.0;
  float HZ2bf = 0.0;

#define HZ_ZM  (z==0 ? HZbf  : HZ[z-1])
#define HZ2_ZM (z==0 ? HZ2bf : HZ2[z-1])
#else
#define HZ_ZM  HZ[(z+NZs) % NZ]
#define HZ2_ZM HZ2[(z+NZs) % NZ]
#endif

#define MAP(Dir,M,Ptr) ( \
 (((Area[Ptr] & 0x7E) & ~(2<<Dir))==0) ? M[Ptr] : _MAP(Dir,Area,M,Ptr,NX,NY,NZ) \
)

#define IsBound ((unsigned char) Area[Ptr]>1)
#define IsExchng (Area[Ptr] & ExchngBit)
#define IsTop (Area[Ptr] & TopBit)
#define IsBottom (Area[Ptr] & BottomBit)
#define IsLeft (Area[Ptr] & LeftBit)
#define IsRight (Area[Ptr] & RightBit)
#define IsForw (Area[Ptr] & ForwBit)
#define IsBack (Area[Ptr] & BackBit)
#define IsClosed(F) (Maps[_Num##F].Map[Ptr]==BoundClosed)
#define IsCond(F,C) (Maps[_Num##F].Map[Ptr]==C)
#define IsDirCond(D,F,C) (MAP(D,Maps[_Num##F].Map,Ptr)==C)

unsigned char _MAP(char Dir, unsigned char * Area, unsigned char * Map, int Ptr, int NX, int NY, int NZ)
{
 int NXY = NX*NY;
 int NXYZ = NXY*NZ;
 int dx = 0;
 int dy = 0;
 int dz = 0;
 int Ptr1;
 
 if (Dir==_right || Dir==_left)
    {
     dy = IsForw ? +1 : IsBack ? -1 : 0;
     dz = IsTop ? +1 : IsBottom ? -1 : 0;
    }
 else  if (Dir==_top || Dir==_bottom)
    {
     dx = IsRight ? +1 : IsLeft ? -1 : 0;
     dy = IsForw ? +1 : IsBack ? -1 : 0;
    }
 else
    {
     dx = IsRight ? +1 : IsLeft ? -1 : 0;
     dz = IsTop ? +1 : IsBottom ? -1 : 0;
    }

 Ptr1 = Ptr+dz*NXY+dy*NX+dx;

 if (Ptr1>=0 && Ptr1<NXYZ && Map[Ptr1]!=Zero && Map[Ptr1]!=ExchngBound)
    return Map[Ptr1];
 else
    return Map[Ptr];
}

#define FZYXP(F) (*VAL(Bounds,_right,&F[ZYXP]))
#define FZYXM(F) (*VAL(Bounds,_left,&F[ZYXM]))
#define FZYPX(F) (*VAL(Bounds,_forw,&F[ZYPX]))
#define FZYMX(F) (*VAL(Bounds,_back,&F[ZYMX]))
#define FZPYX(F) (*VAL(Bounds,_top,&UH##F[y*NX+x]))
#define FZMYX(F) (*VAL(Bounds,_bottom,&DH##F[y*NX+x]))

#define _rg(F) (*VAL(Bounds,_right,&F))
#define _lf(F) (*VAL(Bounds,_left,&F))
#define _fw(F) (*VAL(Bounds,_forw,&F))
#define _bw(F) (*VAL(Bounds,_back,&F))
#define _tp(F) (*VAL(Bounds,_top,&F))
#define _bt(F) (*VAL(Bounds,_bottom,&F))

#ifdef __USE_VCL__
#undef FZPYX
#define FZPYX(F) (*VAL(Bounds,_top,&F[ZPYX]))
#undef FZMYX
#define FZMYX(F) (*VAL(Bounds,_bottom,&F[ZMYX]))
#define dFdx(F) (FZYXP(F)*HX2[XM] + F[Ptr]*(HX2[x]-HX2[XM]) - FZYXM(F)*HX2[x])/HXP[x]
#define dFdy(F) (FZYPX(F)*HY2[YM] + F[Ptr]*(HY2[y]-HY2[YM]) - FZYMX(F)*HY2[y])/HYP[y]
#define dFdz(F) (FZPYX(F)*HZ2[ZM] + F[Ptr]*(HZ2[z]-HZ2[ZM]) - FZMYX(F)*HZ2[z])/HZP[z]
#define d2Fdx2(F) 2.0*(FZYXP(F)*HX[XM] - F[Ptr]*HXS[x] + FZYXM(F)*HX[x])/HXP[x]
#define d2Fdy2(F) 2.0*(FZYPX(F)*HY[YM] - F[Ptr]*HYS[y] + FZYMX(F)*HY[y])/HYP[y]
#define d2Fdz2(F) 2.0*(FZPYX(F)*HZ[ZM] - F[Ptr]*HZS[z] + FZMYX(F)*HZ[z])/HZP[z]
#else
#define dFdx(F) ((FZYXP(F)*HX2[XM] + F[Ptr]*(HX2[x]-HX2[XM]) - FZYXM(F)*HX2[x])/HXP[x])
#define dFdy(F) ((FZYPX(F)*HY2[YM] + F[Ptr]*(HY2[y]-HY2[YM]) - FZYMX(F)*HY2[y])/HYP[y])
#define dFdz(F) ((FZPYX(F)*HZ2_ZM  + F[Ptr]*(HZ2[z]-HZ2_ZM) - FZMYX(F)*HZ2[z])/HZP[z])
#define _dFdx(d,F) ((*VAL(Bounds,d,&F[ZYXP])*HX2[XM] + *VAL(Bounds,d,&F[Ptr])*(HX2[x]-HX2[XM]) - *VAL(Bounds,d,&F[ZYXM])*HX2[x])/HXP[x])
#define _dFdy(d,F) ((*VAL(Bounds,d,&F[ZYPX])*HY2[YM] + *VAL(Bounds,d,&F[Ptr])*(HY2[y]-HY2[YM]) - *VAL(Bounds,d,&F[ZYMX])*HY2[y])/HYP[y])
#define _dFdz(d,F) ((*VAL(Bounds,d,&UH##F[y*NX+x])*HZ2_ZM + *VAL(Bounds,d,&F[Ptr])*(HZ2[z]-HZ2_ZM) - *VAL(Bounds,d,&DH##F[y*NX+x])*HZ2[z])/HZP[z])

#define d2Fdx2(F) 2.0*(FZYXP(F)*HX[XM] - F[Ptr]*HXS[x] + FZYXM(F)*HX[x])/HXP[x]
#define d2Fdy2(F) 2.0*(FZYPX(F)*HY[YM] - F[Ptr]*HYS[y] + FZYMX(F)*HY[y])/HYP[y]
#define d2Fdz2(F) 2.0*(FZPYX(F)*HZ_ZM  - F[Ptr]*HZS[z] + FZMYX(F)*HZ[z])/HZP[z]
#define d2NFdx2(F,N) 2.0*(FZYXP(F)*FZYXP(N)*HX[XM] - F[Ptr]*N[Ptr]*HXS[x] + FZYXM(F)*FZYXM(N)*HX[x])/HXP[x]
#define d2NFdy2(F,N) 2.0*(FZYPX(F)*FZYPX(N)*HY[YM] - F[Ptr]*N[Ptr]*HYS[y] + FZYMX(F)*FZYMX(N)*HY[y])/HYP[y]
#define d2NFdz2(F,N) 2.0*(FZPYX(F)*FZPYX(N)*HZ_ZM  - F[Ptr]*N[Ptr]*HZS[z] + FZMYX(F)*FZMYX(N)*HZ[z])/HZP[z]
#define dRdFdx2(F,R) ((FZYXP(F)*HX[XM] - F[Ptr]*HX[x])*(FZYXP(R)+R[Ptr]) - (F[Ptr]*HX[XM] - FZYXM(F)*HX[x])*(FZYXM(R)+R[Ptr]))/HXP[x]
#define dRdFdy2(F,R) ((FZYPX(F)*HY[YM] - F[Ptr]*HY[y])*(FZYPX(R)+R[Ptr]) - (F[Ptr]*HY[YM] - FZYMX(F)*HY[y])*(FZYMX(R)+R[Ptr]))/HYP[y]
#define dRdFdz2(F,R) ((FZPYX(F)*HZ_ZM - F[Ptr]*HZ[z])*(FZPYX(R)+R[Ptr]) - (F[Ptr]*HZ_ZM - FZMYX(F)*HZ[z])*(FZMYX(R)+R[Ptr]))/HZP[z]
#define _d2Fdx2(d,F) 2.0*(*VAL(Bounds,d,&F[ZYXP])*HX[XM] - *VAL(Bounds,d,&F[Ptr])*HXS[x] + *VAL(Bounds,d,&F[ZYXM])*HX[x])/HXP[x]
#define _d2Fdy2(d,F) 2.0*(*VAL(Bounds,d,&F[ZYPX])*HY[YM] - *VAL(Bounds,d,&F[Ptr])*HYS[y] + *VAL(Bounds,d,&F[ZYMX])*HY[y])/HYP[y]
#define _d2Fdz2(d,F) 2.0*(*VAL(Bounds,d,&F[ZPYX])*HZ_ZM  - *VAL(Bounds,d,&F[Ptr])*HZS[z] + *VAL(Bounds,d,&F[ZMYX])*HZ[z])/HZP[z]
#define dFdxl(F) (FZYXP(F)-_lf(F[Ptr]))/HX[x]
#define dFdxr(F) (_rg(F[Ptr])-FZYXM(F))/HX[XM]
#define dFdyb(F) (FZYPX(F)-_bw(F[Ptr]))/HY[y]
#define dFdyf(F) (_fw(F[Ptr])-FZYMX(F))/HY[YM]
#define dFdzb(F) (FZPYX(F)-_bt(F[Ptr]))/HZ[z]
#define dFdzt(F) (_tp(F[Ptr])-FZMYX(F))/HZ_ZM
#define dFdxn(F) (IsLeft ? dFdxl(F) : IsRight ? dFdxr(F) : 0.0)
#define dFdyn(F) (IsBack ? dFdyb(F) : IsForw ? dFdyf(F) : 0.0)
#define dFdzn(F) (IsBottom ? dFdzb(F) : IsTop ? dFdzt(F) : 0.0)
#endif

#define DeclareStandard \
          int z,y,x; \
          int Ptr;

#endif
