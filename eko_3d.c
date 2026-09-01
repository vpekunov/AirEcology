#pragma plan vectorized

/* При работе с транслятором MPI => Planning C в Windows рекомендуется: */
/* mpi2reent.exe <num_procs> -followdefines -D_WIN32 -D_MSC_VER eko_3d.c <outfile.cpp> */
/* Reenterable <outfile.cpp> <translated.cpp> */
/* Далее компилировать <translated.cpp> */
/* При работе с транслятором MPI => Planning C в UNIX/LINUX рекомендуется: */
/* mpi2reent.exe <num_procs> -followdefines -D__UNIX__ eko_3d.c <outfile.cpp> */
/* Reenterable <outfile.cpp> <translated.cpp> */
/* Далее компилировать <translated.cpp> */

/* Преобразование в C++ без MPI под Linux: Reenterable.exe -followdefines -nosourcelines -D__UNIX__ -D__OPENMP__ -D__MVS__ eko_3d.c reent1.cpp */
/* Преобразование в C++ без MPI под Windows: Reenterable.exe -followdefines -nosourcelines -D_MSC_VER -D_WIN32 -D__OPENMP__ -D__MVS__ eko_3d.c reent1.cpp */

/* (C) 1997-2015 V.V.Pekunov */
/* Программа моделирования (3D модель) с турбулентностью */
/* Последний вариант -- с дополнительной балансировкой по К-Я */
/* F(T) = kT*T */
/* Учитывается распространение пыли. */
/* Экстраполяция значений на стыках. Подстройка коэффициентов регрессии. */
/* Специфические моменты вынесены в отдельный файл specific.h */
/* Поддерживается подключение дополнительных скалярных переменных с решателями */
/* -- Июль 2004 -- */
/* Теперь математическая модель проектируется визуально. */
/* specific.h генерируется автоматически специальной системой */
/* -- Июль 2005 -- */
/* Изменена базовая модель. Переход к системе с векторным потенциалом. */
/* Теперь архимедова сила F(T) задается производными, исчезла необходимость */
/*   задавать опорную температуру */
/* Развит механизм дополнительных решателей. */
/* Теперь можно задавать ГРАНИЧНЫЕ условия функцией */
/* Добавлен алгебраический решатель (прямое вычисление) */
/* Добавлен решатель для уравнений Пуассона (верхняя релаксация с */
/*   четно-нечетным (шахматным) порядком обхода, с чебышевским */
/*   ускорением. Оптимальное Theta ищется сканированием (оптимизация) */
/* Внесена большая неявность -- функция K разбита на две части по Головичеву */
/*   K = K'+H*S */
/* Теперь можно свободно вводить источниковые члены в уравнения для веществ */
/* ВАЖНОЕ ИЗМЕНЕНИЕ: теперь K, S везде задаются с естественным знаком ПЛЮС */
/* ВАЖНОЕ ИЗМЕНЕНИЕ: */
/* Введено разделение правых частей на быструю и медленную подсистемы */
/* Медленная подсистема работает при полном интегрировании всех уравнений */
/*   с главным шагом TAU -- с предвычислением, предсказанием, обменами. */
/* Быстрая подсистема интегрируется "внутри" одной итерации медленной */
/*   подсистемы, с шагом FastTAU. При этом работа идет локально, */
/*   нет обменов данными, нет предвычислений по Головичеву, отключается */
/*   предсказание. Выключаются все решатели, КРОМЕ Рожкова и ДИНАМИЧЕСКОГО */
/*   который переходит в режим вычислений по схеме Рожкова. */
/* Добавлен сброс в ноль значений переменных, по модулю меньших 1E-20 */
/* Теперь можно задавать НАЧАЛЬНЫЕ значения функцией */
/* Добавлена возможность указать границы значений переменных: +/-/любая. */
/* Добавлен расчет угла падения солнечных лучей */
/* Добавлен решатель для уравнения распространения прямых солнечных лучей */
/* Добавлено излучение: прямое и диффузное солнечное (видимое и инфракрасное), */
/*  также тепловое инфракрасное. */
/* -- Февраль 2006 -- */
/* Программа адаптирована для поддержки Router МВС-100 */
/* -- Апрель 2007 -- */
/* Введена частичная прямая поддержка блоков с общей памятью (гетерогенных систем) */
/* -- Июль 2007 -- */
/* Введена поддержка OpenMP */
/* -- Август 2007 -- */
/* Теперь у веществ в реакциях могут быть отрицательные и дробные коэффициенты */
/* -- Сентябрь 2007 -- */
/* Усовершенствована обработка граничных условий 2-го рода на углах. Теперь */
/*  в угол записывается среднее значение по соседним узлам по нормалям. */
/* Реализована более корректная схема для противоточных производных. */
/* Введена динамическая настройка шага по времени (временная сетка множителей) */
/* -- Ноябрь 2007 -- */
/* Изменена базовая модель. Вернулся к системе с давлением (уравнение Пуассона) */
/* Введена возможность отключения счета по базовой модели в заданное время */
/* -- Сентябрь 2015 */
/* Изменение параметров при запуске. Подключение графического интерфейса ввода параметров params.exe */
/* -- Октябрь 2015 */
/* Новый контроль точности */
/* -- Ноябрь 2015 - Май 2016 */
/* Теперь сложные граничные узлы с несколькими гранями сохраняются дополнительно в файл *.bnd */

/* For Microsoft Visual Studio */
#define _CRT_SECURE_NO_WARNINGS

#include "memoization.h"

/* Интерфейс распараллеливания */
// #define __MPI__          /* Устанавливается при работе с MPI            */
// #define __OPENMP__       /* Устанавливается при работе с OpenMP   */
/*#define __ROUTER__     */  /* Устанавливается при работе с Router+        */
/*#define __ROUTER_100__ */  /* Устанавливается при работе с Router MVS-100 */

/* Платформа */
#define __MVS__        /* Устанавливается при работе с МВС-1000         */
/* #define __UNIX__ */      /* Устанавливается при работе с UNIX/LINUX/PARIX */
/* По умолчанию -- DOS/WIN32 платформа */

/*#define __DEBUG__    */ /* Разрешить отладочный вывод в файл */

#ifdef __OPENMP__
int UseOpenMP = 1;
#else
int UseOpenMP = 0;
#endif

#define StopFileName "stop"

int OMP_KineticChunk = 50;

#ifdef __MVS__
#define MVSConfigFile "eko_3d.mvs"
#define MVSOutputBase "eko_3d."
#endif

#if defined(__MPI__) || defined(__ROUTER__) || defined(__ROUTER_100__)
#define __PARALLEL__
#endif

#if defined(__WIN32__)
#define __MAY_BE_BORLANDC__
#ifndef WIN32
#define WIN32
#endif
#elif defined(_WIN32)
#define __MAY_BE_MICROSOFTC__
#ifndef WIN32
#define WIN32
#endif
#endif

#if defined(__LARGE__) || defined(WIN32) || defined(__MVS__) || defined(__UNIX__)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if !defined(__MPI__) && !defined(__ROUTER__)
#if !defined(__UNIX__) && !defined(__MVS__)
#include <sys\timeb.h>
#else
#include <time.h>
#endif
#endif

#ifndef __ROUTER_100__

#if (defined(__MVS__) || defined(__UNIX__)) && !defined(WIN32)
#include <sys/shm.h>
#else
#include <winsock.h>
#endif

#endif

#if !defined(__MAY_BE_MICROSOFTC__) && !defined(__MVS__) && !defined(__UNIX__)
#include <conio.h>
#endif

#if defined(__MPI__)

#ifdef __OPENMP__

#include <omp.h>
#include <mpi.h>

#else

#include <mpi.h>

#endif

#elif defined(__ROUTER__)

#include "routelib.h"

#elif defined(__ROUTER_100__)

#include "Router\routelib.h"

#elif defined(__OPENMP__)

#include <omp.h>

#endif

#if defined(WIN32) || defined(__MVS__) || defined(__UNIX__)
typedef unsigned int  _int;
const char * _intScanf = "%u";
#else
typedef unsigned long _int;
const char * _intScanf = "%lu";
#endif

const char * _doubleScanf = "%lf";

#include "CFGFILE.H"
#include "area3d.h"
#include "KINETIC.H"

KineticContext * Contexts;
KineticGlobal KGlobal;
int KinetChunkSize = 120*32; // 3840

int NSubst;
int NASubst;
int NReact;

int SetSteadyKinetics = 0; /* !!! */
int InitOpenMPIters = 3; /* !!! */

#ifndef __PARALLEL__
TraceTypeHost * Traces = NULL;
int NKinets = 0;

double KinetTime = 0.0;
#endif

#include "specific.h"

#ifdef __OPENMP__
  typedef double TIME_STRUCT;
#define FTIME(Arg) *Arg=omp_get_wtime()
#define DIFFTIME(BeginTime,EndTime) (EndTime-BeginTime)
#elif defined(__MPI__)
  typedef double TIME_STRUCT;
#define FTIME(Arg) *Arg=MPI_Wtime()
#define DIFFTIME(BeginTime,EndTime) (EndTime-BeginTime)
#elif defined(__ROUTER__)
  typedef int TIME_STRUCT;
#define FTIME(Arg) *Arg=_cputime()
#define DIFFTIME(BeginTime,EndTime) ((double)(EndTime-BeginTime)/1000.0)
#elif defined(__UNIX__) || defined(__MVS__)
  typedef time_t TIME_STRUCT;
#define FTIME time
#define DIFFTIME(BeginTime,EndTime) (difftime(EndTime,BeginTime))
#elif defined(__MAY_BE_MICROSOFTC__)
  typedef struct _timeb TIME_STRUCT;
#define FTIME _ftime
#define DIFFTIME(BeginTime,EndTime) ((double) (EndTime.time-BeginTime.time) + \
                                       (double) (EndTime.millitm-BeginTime.millitm)/1000)
#else
  typedef struct timeb TIME_STRUCT;
#define FTIME ftime
#define DIFFTIME(BeginTime,EndTime) ((double) (EndTime.time-BeginTime.time) + \
                                       (double) (EndTime.millitm-BeginTime.millitm)/1000)
#endif

#ifdef __MVS__
 FILE * StdInput = NULL;
 FILE * StdOutput = NULL;
#else
 FILE * StdInput;
 FILE * StdOutput;
#endif

#define MIN_EPS 1E-8

double Alpha; /* Параметр для регуляризации по Булееву */

int    UseGear = 0;
int    MinGap  = 1;

int    UseEnhancedOpenMP = 1;

double EndTime = 0.0;

const double k1 = 50.0;
const double k2 = 5.0;
const double k4 = 0.001;

int    CheckNew = 1; /* 0 = IMPLICIT ONLY. 1 = IMPLICIT/SEMI-IMPLICIT */

int    CheckTau;
double MaxTau;
int    DecreaseOnly;
int    Period;
double Eps;
double EpsP;

#define MaxTimeN 50

int    TimeN = 0;
double TimeIters[MaxTimeN] = {0.0};
double TimeTAUk[MaxTimeN] = {0.0};
int    TimeStopBase[MaxTimeN] = {1};
int    Division[256] = {0};

int    CalcBase = 1;

double _HX[MaxNX] = {0.0};
double _HY[MaxNX] = {0.0};
double _HZ[MaxNX] = {0.0};

#define Lmin2FExt ".lm2"

double ZeroC[MaxActSubst];
int    MapSubsts[MaxActSubst];
int    NLightSubst = 0;

int    Saves[NumSaves];
int    SaveDn[MaxActSubst];
FILE * SaveFiles[NumSaves];

long SpecialBoundsNum = 0;

float ** OldRefAddrs = NULL;
float *  OldRefVals  = NULL;
long     OldRefNum   = 0;

char ReverseDirs[NumEqs+MaxActSubst] = {0};

typedef unsigned int word;

int UseSHMEM = 0;
int * ProcOnNode = NULL;
int SaveAllHistory = 1;

#ifdef __PARALLEL__
FILE * DebugFile = NULL;

int  nProcs   = 0;
int  MyProcID = 0;

int  TotalProcs = 0;
int  KinetProcs = 0;

int  CalculateAll = 0;

int  MasterFlag = 0;
int  FirstFlag  = 0;
int  LastFlag   = 0;

int  UpExchange   = 0;
int  SHMEMUp      = 0;
int  DownExchange = 0;
int  SHMEMDown    = 0;

#ifdef __MPI__
MPI_Comm  SlavesComm;
MPI_Group WorldGroup;
MPI_Group SlavesGroup;
#endif

int SlaveID = 0;

typedef unsigned char byte;
#define SlaveToProcID(i) ((i)+1)
#define ProcIDToSlave(i) ((i)-1)

#ifdef __DEBUG__

void PrepareDebug()
{
 if (DebugFile)
    fflush(DebugFile);
 else
    {
     char FName[30];

     sprintf(FName,"Log%i",MyProcID);
     DebugFile = fopen(FName,"wt");
    }
}

#define DebugPrintf PrepareDebug(); fprintf
#define DEBUG_FILE DebugFile

#else

#define DebugPrintf if (0) fprintf
#define DEBUG_FILE  stdout

#endif

int * FromGrid   = NULL;
int * LengthGrid = NULL;

unsigned int BuffNodes = 0;
float    *   HBuff;
float    *   eHBuff;

int RecvBoundSize;

#ifdef __ROUTER__
typedef struct {
 int Handle;
 int BufferNum;
} SendDsc;

typedef struct {
 int   Size;
 int * Buffer;
 int   Free;
} BufListItem;

int     * NumSends = NULL; /* TotalProcs */
int     * MaxSends = NULL; /* TotalProcs */
SendDsc * HndSends = NULL; /* TotalProcs*DefMaxSends */

int           MaxListItems = 0;
int           NumListItems = 0;
BufListItem * BufList      = NULL;

#define DefMaxSends 16

#define DeltaListSize 50

void * CreateBuffer(void * Buf, long Length, int * BufItem)
{
 int BestSize = 0;
 int Count;

 div_t LD = div(Length,sizeof(int));
 int   _Length = LD.quot + (LD.rem!=0);

 for (Count = 0; Count<NumListItems; Count++)
     if (BufList[Count].Size>=_Length && BufList[Count].Free)
        if (!BestSize || BestSize>BufList[Count].Size)
           {
            *BufItem = Count;
            BestSize = BufList[Count].Size;
           }
 if (BestSize)
    BufList[*BufItem].Free = 0;
 else
    {
     for (Count = 0; Count<NumListItems; Count++)
         if (BufList[Count].Size>BestSize && BufList[Count].Free)
            {
             *BufItem = Count;
             BestSize = BufList[Count].Size;
            }
     if (BestSize)
        BufList[*BufItem].Buffer = (int *) SafeRealloc(BufList[*BufItem].Buffer,_Length*sizeof(int));
     else
        {
         *BufItem = NumListItems;
         if (NumListItems==MaxListItems)
            {
             MaxListItems += DeltaListSize;
             BufList = (BufListItem *) SafeRealloc(BufList,MaxListItems*sizeof(BufListItem));
            }
         NumListItems++;
         BufList[*BufItem].Buffer = (int *) SafeMalloc(_Length*sizeof(int));
        }
     BufList[*BufItem].Free   = 0;
     BufList[*BufItem].Size   = _Length;
    }
 return memmove(BufList[*BufItem].Buffer,Buf,Length);
}

void FreeBuffer(int BufItem)
{
 BufList[BufItem].Free = 1;
}

void RouterSend(int Proc, void * Buf, long Length)
{
 int Base = Proc*DefMaxSends;

 if (NumSends[Proc] == MaxSends[Proc] ||
     (NumSends[Proc] && t_write(HndSends[Base].Handle)))
    {
     w_write(HndSends[Base].Handle);
     FreeBuffer(HndSends[Base].BufferNum);
     if (--NumSends[Proc])
        memmove(&HndSends[Base],&HndSends[Base+1],NumSends[Proc]*sizeof(SendDsc));
     RouterSend(Proc, Buf, Length);
    }
 else
   {
    int RetCode = r_write(Proc,
                          CreateBuffer(Buf,Length,&HndSends[Base+NumSends[Proc]].BufferNum),
                          Length);

    switch (RetCode) {
      case RUN_BUSY:
        FreeBuffer(HndSends[Base+NumSends[Proc]].BufferNum);
        MaxSends[Proc] = NumSends[Proc];
        DebugPrintf(DEBUG_FILE,"MaxSends[%i] = %i now\n",Proc,NumSends[Proc]);
        RouterSend(Proc, Buf, Length);
        break;
      case RUN_OK:
        HndSends[Base+NumSends[Proc]++].Handle = started_handle;
        break;
      default:
        DebugPrintf(DEBUG_FILE,"Error on r_write(%i,%x,%i) = %i\n",Proc,Buf,Length,RetCode);
#ifdef __MVS__
        fclose(StdOutput);
#endif
        AbortServer(-30);
	}
   }
}

void RouterIRecv(int Proc, void * Buf, long Length)
{
 int RetCode = r_read(Proc, Buf, Length);

 if (RetCode!=RUN_OK)
    {
     DebugPrintf(DEBUG_FILE,"Error on r_recv(%i,%x,%i) = %i\n",Proc,Buf,Length,RetCode);
#ifdef __MVS__
     fclose(StdOutput);
#endif
     AbortServer(-40);
    }
}
#endif

#ifdef __ROUTER_100__
void ROUTER_100_SEND(int ID, byte * Buf, int Length)
{
 if (Length%sizeof(int))
    {
     int    sLength = Length+sizeof(int)-Length%sizeof(int);
     byte * sBuf    = SafeMalloc(sLength);
     memmove(sBuf,Buf,Length);
     r_write(ID, sBuf, sLength);
     w_write(ID);
     free(sBuf);
    }
 else
    {
     r_write(ID, Buf, Length);
     w_write(ID);
    }
}

void ROUTER_100_RECV(int ID, byte * Buf, int Length)
{
 if (Length%sizeof(int))
    {
     int    rLength = Length+sizeof(int)-Length%sizeof(int);
     byte * rBuf    = SafeMalloc(rLength);
     r_read(ID, rBuf, rLength);
     w_read(ID);
     memmove(Buf,rBuf,Length);
     free(rBuf);
    }
 else
    {
     r_read(ID, Buf, Length);
     w_read(ID);
    }
}
#endif

#if !defined(__ROUTER_100__)

#if !((defined(__UNIX__) || !defined(__MVS__)) && !defined(WIN32))
#define MemMultiplier 100

int    MemCounter = 0;
HANDLE MemHandles[40];
HANDLE InHandles[2][40];
#endif

void * SHMEM_AT(int Handle, int Size)
{
#if (defined(__UNIX__) || defined(__MVS__)) && !defined(WIN32)
 return (void *) shmat(Handle,0,0);
#else
 int ProcID    = Handle/MemMultiplier;
 int MemNumber = Handle%MemMultiplier;

 if (ProcID==MyProcID)
    return (void *) MapViewOfFile(MemHandles[MemNumber],FILE_MAP_ALL_ACCESS,0,0,Size);
 else
    {
     int N = ProcID>MyProcID;
     char Buf[20];

     sprintf(Buf,"Ecology%i",Handle);
     if (!(InHandles[N][MemNumber] = OpenFileMapping(FILE_MAP_ALL_ACCESS,FALSE,Buf)))
        {
         fprintf(StdOutput,"Can't open Shared Memory Region. Error = %i\n",GetLastError());
#ifdef __MVS__
         fclose(StdOutput);
#endif
         AbortServer(-400);
        }

     return (void *) MapViewOfFile(InHandles[N][MemNumber],FILE_MAP_ALL_ACCESS,0,0,Size);
    }
#endif
}

void SHMEM_DT(void * Buf, int Handle)
{
#if (defined(__UNIX__) || defined(__MVS__)) && !defined(WIN32)
 shmdt(Buf);
#else
 int ProcID    = Handle/MemMultiplier;
 int MemNumber = Handle%MemMultiplier;

 UnmapViewOfFile(Buf);
 if (ProcID!=MyProcID)
    CloseHandle(InHandles[ProcID>MyProcID][MemNumber]);
#endif
}

int SHMEM_GET(int Size)
{
#if (defined(__UNIX__) || defined (__MVS__)) && !defined(WIN32)
 int result = shmget(IPC_PRIVATE,Size,IPC_CREAT | 0666);

 if (result < 0)
    printf("Shared memory allocation error: %i\n", errno);

 return result;
#else
 char Buf[20];
 int Result = MyProcID*MemMultiplier+MemCounter;

 sprintf(Buf,"Ecology%i",Result);
 if (!(MemHandles[MemCounter++] = CreateFileMapping(INVALID_HANDLE_VALUE, NULL,
       PAGE_READWRITE, 0, Size, Buf)))
    {
     fprintf(StdOutput,"Can't create Shared Memory Region. Error = %i\n",GetLastError());
#ifdef __MVS__
     fclose(StdOutput);
#endif
     AbortServer(-400);
    }

 return Result;
#endif
}

void SHMEM_RM(int Handle)
{
#if (defined(__UNIX__) || defined (__MVS__)) && !defined(WIN32)
 struct shmid_ds buf;

 shmctl(Handle,IPC_RMID,&buf);
#else
 CloseHandle((HANDLE)(Handle%MemMultiplier));
#endif
}

#endif

void SendMaster(byte * Buf, long Length)
{
#if defined(__ROUTER_100__)
 ROUTER_100_SEND(0,Buf,Length);
#elif defined(__ROUTER__)
 RouterSend(0, Buf, Length);
#else
 MPI_Send(Buf, Length, MPI_BYTE, 0, MyProcID, MPI_COMM_WORLD);
#endif
}

void SendMasterXXXTag(byte * Buf, long Length)
{
#ifdef __MPI__
 MPI_Send(Buf, Length, MPI_BYTE, 0, 123456, MPI_COMM_WORLD);
#else
 SendMaster(Buf,Length);
#endif
}

void RecvMaster(byte * Buf, long Length)
{
#if defined(__ROUTER_100__)
 ROUTER_100_RECV(0,Buf,Length);
#elif defined (__ROUTER__)
 RouterIRecv(0, Buf, Length);
 w_read(started_handle);
#else
 MPI_Status status;

 MPI_Recv(Buf, Length, MPI_BYTE, 0, MyProcID, MPI_COMM_WORLD, &status);
#endif
}

#if defined(__MPI__) || defined(__ROUTER__) || defined(__ROUTER_100__)

#ifdef __MPI__
typedef MPI_Request REQUEST;
#else
typedef int REQUEST;
#endif

void IRecvMaster(byte * Buf, long Length, REQUEST * Request)
{
#ifdef __MPI__
 MPI_Irecv(Buf, Length, MPI_BYTE, 0, MyProcID, MPI_COMM_WORLD, Request);
#elif defined(__ROUTER_100__)
 r_read(0, Buf, Length);
 *Request = 0;
#else
 RouterIRecv(0,Buf,Length);
 *Request = started_handle;
#endif
}

void IRecvNeighbour(int slave, byte * Buf, long Length, REQUEST * Request)
{
#ifdef __MPI__
 MPI_Irecv(Buf, Length, MPI_BYTE, slave, 1235, SlavesComm, Request);
#elif defined(__ROUTER_100__)
 r_read(SlaveToProcID(slave), Buf, Length);
 *Request = SlaveToProcID(slave);
#else
 RouterIRecv(SlaveToProcID(slave),Buf,Length);
 *Request = started_handle;
#endif
}

void Wait(REQUEST * Request)
{
#ifdef __MPI__
 MPI_Status Status;

 MPI_Wait(Request,&Status);
#else
 w_read(*Request);
#endif
}
#endif

void SendSlave(int slave, byte * Buf, long Length)
{
#if defined(__ROUTER_100__)
 ROUTER_100_SEND(SlaveToProcID(slave),Buf,Length);
#elif defined(__ROUTER__)
 RouterSend(SlaveToProcID(slave), Buf, Length);
#else
 MPI_Send(Buf, Length, MPI_BYTE, SlaveToProcID(slave), SlaveToProcID(slave), MPI_COMM_WORLD);
#endif
}

void RecvSlave(int slave, byte * Buf, long Length)
{
#if defined(__ROUTER_100__)
 ROUTER_100_RECV(SlaveToProcID(slave),Buf,Length);
#elif defined(__ROUTER__)
 RouterIRecv(SlaveToProcID(slave), Buf, Length);
 w_read(started_handle);
#else
 MPI_Status status;

 MPI_Recv(Buf, Length, MPI_BYTE, SlaveToProcID(slave), SlaveToProcID(slave), MPI_COMM_WORLD, &status);
#endif
}

void RecvSlaveXXXTag(int slave, byte * Buf, long Length)
{
#ifdef __MPI__
 MPI_Status status;

 MPI_Recv(Buf, Length, MPI_BYTE, SlaveToProcID(slave), 123456, MPI_COMM_WORLD, &status);
#else
 RecvSlave(slave,Buf,Length);
#endif
}

#endif

/* Коэффициенты,определяющие вхождение противоточных производных в решение */
WKoeffs WXYZ[NumPhases];

#ifdef __PARALLEL__
void CalculateWXY(float ** Bounds, unsigned char * Area, WKoeffs * W, float * Ux, float * Uy, int Delta)
{
#else
void CalculateWXY(float ** Bounds, unsigned char * Area, WKoeffs * W, float * Ux, float * Uy)
{
 int Delta = 0;
#endif
 DeclareStandard
 int zy;

 #pragma omp parallel if(UseOpenMP)
 #pragma omp for schedule(dynamic,imax(1,NY*NZ/imax(20,4*nSMP))) private(zy,z,y,x,Ptr)
 for (zy=Delta*NY; zy<(NZ-Delta)*NY; zy++)
    {
     z = zy/NY;
     y = zy%NY;
     for (x=0, Ptr=NX*zy; x<NX; Ptr++,x++)
      if (!Area[Ptr] || IsExchng)
         {
          float Uxs  = Ux[Ptr];
          float Uxm  = FZYXP(Ux)*HXR[x].R1m+Uxs*HXR[x].R2m+FZYXM(Ux)*HXR[x].R3m;
          float Uxp  = FZYXP(Ux)*HXR[x].R1p+Uxs*HXR[x].R2p+FZYXM(Ux)*HXR[x].R3p;
          float Uys  = Uy[Ptr];
          float Uym  = FZYPX(Uy)*HYR[y].R1m+Uys*HYR[y].R2m+FZYMX(Uy)*HYR[y].R3m;
          float Uyp  = FZYPX(Uy)*HYR[y].R1p+Uys*HYR[y].R2p+FZYMX(Uy)*HYR[y].R3p;

          W->WX1[Ptr]  = 0.5*(Uxs+Uxp)/HXR[x].h;
          W->WX2[Ptr]  = 0.5*(Uxs+Uxm)/HXR[x].h;
          W->WY1[Ptr]  = 0.5*(Uys+Uyp)/HYR[y].h;
          W->WY2[Ptr]  = 0.5*(Uys+Uym)/HYR[y].h;
         }
    }
}

#ifdef __PARALLEL__
void CalculateWZ(float ** Bounds, unsigned char * Area, WKoeffs * W, float * Uz, double Uw, int Delta, int UseOpenMP)
{
#else
void CalculateWZ(float ** Bounds, unsigned char * Area, WKoeffs * W, float * Uz, double Uw, int UseOpenMP)
{
 int Delta = 0;
#endif
 DeclareStandard
 int zy;

 #pragma omp parallel if(UseOpenMP)
 #pragma omp for schedule(dynamic,imax(1,NY*NZ/imax(20,4*nSMP))) private(zy,z,y,x,Ptr)
 for (zy=Delta*NY; zy<(NZ-Delta)*NY; zy++)
    {
     z = zy/NY;
     y = zy%NY;
     for (x=0, Ptr=NX*zy; x<NX; Ptr++,x++)
      if (!Area[Ptr] || IsExchng)
         {
          float Uzs  = Uz[Ptr]-Uw;
          float Uzm  = (*VAL(Bounds,_top,&Uz[ZPYX])-Uw)*HZR[z].R1m+Uzs*HZR[z].R2m+(*VAL(Bounds,_bottom,&Uz[ZMYX])-Uw)*HZR[z].R3m;
          float Uzp  = (*VAL(Bounds,_top,&Uz[ZPYX])-Uw)*HZR[z].R1p+Uzs*HZR[z].R2p+(*VAL(Bounds,_bottom,&Uz[ZMYX])-Uw)*HZR[z].R3p;

          W->WZ1[Ptr]  = 0.5*(Uzs+Uzp)/HZR[z].h;
          W->WZ2[Ptr]  = 0.5*(Uzs+Uzm)/HZR[z].h;
         }
    }
}

#ifdef __PARALLEL__
 #define OffsBuf(S) (((NZ+2)*_Num##S+1)*NX*NY)
 #define DeclareUpDown(S) \
           float * DH##S = z==0   ? &V[OffsBuf(S)-NY*NX] : &S[(z-1)*NY*NX]; \
           float * UH##S = z==NZs ? &V[OffsBuf(S)+NZ*NY*NX] : &S[(z+1)*NY*NX];
#else
 #define DeclareUpDown(S) \
           float * DH##S = z==0   ? &S[NZs*NY*NX] : &S[(z-1)*NY*NX]; \
           float * UH##S = z==NZs ? &S[0] : &S[(z+1)*NY*NX];
#endif

void CalculateDIV(float ** Bounds, unsigned char * Area, float * Ux, float * Uy, float * Uz)
{
#ifdef __PARALLEL__
 float * UX = &Ux[NY*NX];
 float * UY = &Uy[NY*NX];
 float * UZ = &Uz[NY*NX];
 int     Delta = 0;
#else
 float * UX = Ux;
 float * UY = Uy;
 float * UZ = Uz;
 int     Delta = 1;
#endif
 DeclareStandard
 int zy;

 memset(DIV,0,BoardSize);
 #pragma omp parallel if(UseOpenMP)
 #pragma omp for schedule(dynamic,imax(1,NY*(NZ-2*Delta)/imax(20,4*nSMP))) private(zy,z,y,x,Ptr)
 for (zy=Delta*NY; zy<(NZ-Delta)*NY; zy++)
    {
     z = zy/NY;
     y = zy%NY;
     {
      float * DHUX = &Ux[(z-Delta)*NY*NX];
      float * UHUX = &Ux[(z-Delta+2)*NY*NX];
      float * DHUY = &Uy[(z-Delta)*NY*NX];
      float * UHUY = &Uy[(z-Delta+2)*NY*NX];
      float * DHUZ = &Uz[(z-Delta)*NY*NX];
      float * UHUZ = &Uz[(z-Delta+2)*NY*NX];

      for (x=0, Ptr=NX*zy; x<NX; Ptr++,x++)
          if (!Area[Ptr] || IsExchng)
             DIV[Ptr] = dFdx(UX)+dFdy(UY)+dFdz(UZ);
     }
    }
}

#undef _S
#undef _L
#undef _M
#undef _G

void OneLineGo (unsigned char * Area,
                unsigned char * Boundaries,
                float   ** Bounds,
                float    * H,
                float    * K,
                float    * _S,
                float    * WX1,
                float    * WX2,
                HKoeffs  * HXR,
                float    * HX,
                float    * HXX,
                float    * HXS,
                float    * H1,
                double   * _L,
                double   * _M,
                double   * _G,
                float      R,
                float    * Rt,
                float      Kp,
                int        NumZ,
                int        NumY,
                int        NumX,
                int        ZStep,
                int        YStep,
                int        XStep,
                DescLine * Desc,
                char       Flags,
                char       Regular,
#ifdef __PARALLEL__
                unsigned char * Map,
                float           HXbf,
                int        UseOpenMP)
#else
                unsigned char * Map,
                int        UseOpenMP)
#endif
{
 memmove(H1,H,BoardSize);
 #pragma omp parallel if(UseOpenMP)
 {
#ifdef __OPENMP__
  int      NArray   = ScKfSize/nSMP/sizeof(double);
  int      ThreadID = omp_get_thread_num();
  double * L = &_L[ThreadID*NArray];
  double * G = &_G[ThreadID*NArray];
  double * M = &_M[ThreadID*NArray];
#else
  double * L = _L;
  double * G = _G;
  double * M = _M;
#endif
  int i,x,y,z,zy;
  int Item;
  char BND = Flags>>shBoundary;

  #pragma omp for schedule(dynamic,imax(1,NumY*NumZ/imax(20,4*nSMP))) private(zy,z,y,x,i,Item)
  for (zy=0;zy<NumY*NumZ;zy++)
      {
       z = zy/NumY;
       y = zy%NumY;
       for (Item = 0; Item<Desc[zy].Num; Item++)
           {
            int X1 = Desc[zy].Items[Item].X1;
            int X2 = Desc[zy].Items[Item].X2;
            int D1 = Desc[zy].Dir1;
            int D2 = Desc[zy].Dir2;
            double S,B,Z1,Z2,Z3;
            double aR1;
            double Rp,Rm,Rz;

            if (X1>=0 && X2>=0)
               {
                int Long = X1<=X2 ? X2-X1+1 : NumX-(X1-X2-1);

                int PrevX = (X1+NumX-1) % NumX;
                int PostX = (X2+1)      % NumX;

                int IndexGRN1 = z*ZStep + y*YStep + PrevX*XStep;
                int IndexGRN2 = z*ZStep + y*YStep + PostX*XStep;

                unsigned char TypeGRN1 = MAP(D1,Map,IndexGRN1);
                unsigned char TypeGRN2 = MAP(D2,Map,IndexGRN2);

                if (TypeGRN1==Bound2)
                   {
                    L[1] = 1;
                    M[1] = 0;
                   }
                else if (TypeGRN1==Bound3)
                   {
                    L[1] = 0.5;
                    M[1] = 0;
                   }
                else if (TypeGRN1>=MinIndexBound1 || TypeGRN1==BoundClosed)
                   {
                    L[1] = 0;
                    M[1] = *VAL(Bounds,D1,&H[IndexGRN1]);
                   }

                for (i=1; i<=Long ; i++)
                    {
                     x    = (X1+i-1) % NumX;

                     if (Boundaries[_ZYX] && Boundaries[_ZYX]<=BND)
                        {
                         L[i+1] = 0.0;
                         M[i+1] = K[_ZYX];
                        }
                     else
                        {
                         aR1  = TAU/HXR[x].h/HXR[x].h;

                         if (Flags & fld2RF)
                            {
                             Rm = R+Kp**VAL(Bounds,D1,&Rt[_ZYXM]);
                             Rp = R+Kp**VAL(Bounds,D2,&Rt[_ZYXP]);
                             Rz = 2.0*(R+Kp*Rt[_ZYX]);
                            }
                         else
                           {
                            Rp   = R+0.5*Kp*(Rt[_ZYX]+(*VAL(Bounds,D2,&Rt[_ZYXP])*HXR[x].R1p+Rt[_ZYX]*HXR[x].R2p+*VAL(Bounds,D1,&Rt[_ZYXM])*HXR[x].R3p));
                            Rm   = R+0.5*Kp*(Rt[_ZYX]+(*VAL(Bounds,D2,&Rt[_ZYXP])*HXR[x].R1m+Rt[_ZYX]*HXR[x].R2m+*VAL(Bounds,D1,&Rt[_ZYXM])*HXR[x].R3m));
                            Rz   = Rm+Rp;
                           }
                         Z1   = aR1*(Rm*HXR[x].R3m+Rp*HXR[x].R3p);
                         Z2   = aR1*(Rm*HXR[x].R1m+Rp*HXR[x].R1p);
                         Z3   = 1+aR1*(Rz-Rm*HXR[x].R2m-Rp*HXR[x].R2p)-(_S==NULL ? 0.0 : TAU/3.0*_S[_ZYX]);
                         if (WX1[_ZYX]>0.0)
                            Z3 += TAU*WX1[_ZYX];
                         else
                            {
                             Z1 -= TAU*WX1[_ZYX]*HXR[x].R3p;
                             Z2 -= TAU*WX1[_ZYX]*HXR[x].R1p;
                             Z3 += TAU*WX1[_ZYX]*HXR[x].R2p;
                            }
                         if (WX2[_ZYX]<=0.0)
                            Z3 -= TAU*WX2[_ZYX];
                         else
                            {
                             Z1 += TAU*WX2[_ZYX]*HXR[x].R3m;
                             Z2 += TAU*WX2[_ZYX]*HXR[x].R1m;
                             Z3 -= TAU*WX2[_ZYX]*HXR[x].R2m;
                            }
                         B    = H[_ZYX]+TAU*(K==NULL ? 0.0 : K[_ZYX])/3;
                         S    = Z3-Z1*L[i];
                         L[i+1] = Z2/S;
                         M[i+1] = (B+M[i]*Z1)/S;
                        }
                    }

                if (TypeGRN2==Bound2)
                    *VAL(Bounds,D2,&H1[IndexGRN2]) = M[Long+1]/(1-L[Long+1]);
                else if (TypeGRN2==Bound3)
                    *VAL(Bounds,D2,&H1[IndexGRN2]) = (M[Long+1]*(2-L[Long])-M[Long])/(1+L[Long+1]*(L[Long]-2));/*M[Long+1]/(2-L[Long+1]);*/
                else if (TypeGRN2>=MinIndexBound1 || TypeGRN2==BoundClosed)
                    *VAL(Bounds,D2,&H1[IndexGRN2]) = *VAL(Bounds,D2,&H[IndexGRN2]);

                for (i=Long; i>=0; i--)
                    {
                     x = (X1+NumX+i-1) % NumX;

                     if (i==Long)
                        H1[_ZYX] = *VAL(Bounds,D2,&H1[IndexGRN2])*L[i+1]+M[i+1];
                     else if (i==0)
                        *VAL(Bounds,D1,&H1[IndexGRN1]) = H1[_ZYXP]*L[i+1]+M[i+1];
                     else
                        H1[_ZYX] = H1[_ZYXP]*L[i+1]+M[i+1];
                    }
               }
            else
               {
                int IndexLast = z*ZStep+y*YStep+(NumX-1)*XStep;

                L[0] =  0;
                G[0] = -1;
                M[0] =  0;

                for (x=0;x<NumX;x++)
                    if (Boundaries[_ZYX] && Boundaries[_ZYX]<=BND)
                       {
                        L[i+1] = 0.0;
                        G[i+1] = 0.0;
                        M[i+1] = K[_ZYX];
                        Z1 = 0.0;
                        Z2 = 0.0;
                        Z3 = 1.0;
                        B = M[i+1];
                       }
                    else
                       {
                        aR1  = TAU/HXR[x].h/HXR[x].h;

                        if (Flags & fld2RF)
                           {
                            Rm = R+Kp**VAL(Bounds,D1,&Rt[_ZYXM]);
                            Rp = R+Kp**VAL(Bounds,D2,&Rt[_ZYXP]);
                            Rz = 2.0*(R+Kp*Rt[_ZYX]);
                           }
                        else
                          {
                           Rp   = R+0.5*Kp*(Rt[_ZYX]+(*VAL(Bounds,D2,&Rt[_ZYXP])*HXR[x].R1p+Rt[_ZYX]*HXR[x].R2p+*VAL(Bounds,D1,&Rt[_ZYXM])*HXR[x].R3p));
                           Rm   = R+0.5*Kp*(Rt[_ZYX]+(*VAL(Bounds,D2,&Rt[_ZYXP])*HXR[x].R1m+Rt[_ZYX]*HXR[x].R2m+*VAL(Bounds,D1,&Rt[_ZYXM])*HXR[x].R3m));
                           Rz   = Rm+Rp;
                          }
                        Z1 = aR1*(Rm*HXR[x].R3m+Rp*HXR[x].R3p);
                        Z2 = aR1*(Rm*HXR[x].R1m+Rp*HXR[x].R1p);
                        Z3 = 1+aR1*(Rz-Rm*HXR[x].R2m-Rp*HXR[x].R2p)-(_S==NULL ? 0.0 : TAU/3.0*_S[_ZYX]);
                        if (WX1[_ZYX]>0.0)
                           Z3 += TAU*WX1[_ZYX];
                        else
                           {
                            Z1 -= TAU*WX1[_ZYX]*HXR[x].R3p;
                            Z2 -= TAU*WX1[_ZYX]*HXR[x].R1p;
                            Z3 += TAU*WX1[_ZYX]*HXR[x].R2p;
                           }
                        if (WX2[_ZYX]<=0.0)
                           Z3 -= TAU*WX2[_ZYX];
                        else
                           {
                            Z1 += TAU*WX2[_ZYX]*HXR[x].R3m;
                            Z2 += TAU*WX2[_ZYX]*HXR[x].R1m;
                            Z3 -= TAU*WX2[_ZYX]*HXR[x].R2m;
                           }
                        B    = H[_ZYX]+TAU*(K==NULL ? 0.0 : K[_ZYX])/3;
                        if (x<NumX-1)
                           {
                            S    = Z3+Z1*L[x];
                            L[x+1] = -Z2/S;
                            G[x+1] = (Z1*G[x])/S;
                            M[x+1] = (B+M[x]*Z1)/S;
                           }
                       }

                Z2 = -Z2;

                for (x=0;x<NumX-1;x++)
                    {
                     Z3 = Z3-Z2*G[x+1];
                     B  = B-Z2*M[x+1];
                     Z2 = -Z2*L[x+1];
                    }

                H1[IndexLast] = (B+Z1*M[NumX-1])/(Z3+Z2+Z1*(G[NumX-1]+L[NumX-1]));

                for (x=NumX-2;x>=0;x--)
                    H1[_ZYX] = M[x+1]-H1[_ZYXP]*L[x+1]-G[x+1]*H1[IndexLast];
               }
           }
      }
 }
 memmove(H,H1,BoardSize);
}

void EulerCalculate(
                    float ** Bounds,
                    float   * H,
                    float   * K,
                    float   * S,
                    WKoeffs * W,
                    HKoeffs * HXR, HKoeffs * HYR, HKoeffs * HZR,
                    float   * HX,  float  * HY,  float * HZ,
                    float   * HXX, float  * HYY, float * HZZ,
                    float   * HXS, float  * HYS, float * HZS,
                    float   * H1,
                    float     R,
                    float   * Rt,
                    float     Kp,
                    char      Flags,
                    float *   ALPHA,
                    unsigned char * Map,
                    int UseOpenMP
                   ) {
 int x, zy;
 int Ptr;

 #pragma omp parallel if(UseOpenMP)
 #pragma omp for schedule(guided) private(zy,x,Ptr)
 for (zy = 0; zy < NY*NZ; zy++)
  {
   int z = zy/NY;
   int y = zy%NY;
   for (x=0, Ptr=zy*NX; x<NX; x++, Ptr++)
       if (Map[Ptr]==Zero && fabs(ALPHA[Ptr]-1.0) > 1E-5)
          {
           double aR1z  = TAU/HZR[z].h/HZR[z].h;
           double Z1z,Z2z,Z3z,Rpz,Rmz,Rzz;
           double aR1y  = TAU/HYR[y].h/HYR[y].h;
           double Z1y,Z2y,Z3y,Rpy,Rmy,Rzy;
           double aR1x  = TAU/HXR[x].h/HXR[x].h;
           double Z1x,Z2x,Z3x,Rpx,Rmx,Rzx;

           double RtZPYX = *VAL(Bounds,_top,&Rt[ZPYX]);
           double RtZMYX = *VAL(Bounds,_bottom,&Rt[ZMYX]);
           double RtZYPX = *VAL(Bounds,_forw,&Rt[ZYPX]);
           double RtZYMX = *VAL(Bounds,_back,&Rt[ZYMX]);
           double RtZYXP = *VAL(Bounds,_right,&Rt[ZYXP]);
           double RtZYXM = *VAL(Bounds,_left,&Rt[ZYXM]);

           if (Flags & fld2RF)
              {
               Rmz = R+Kp*RtZMYX;
               Rpz = R+Kp*RtZPYX;
               Rzz = 2.0*(R+Kp*Rt[Ptr]);
              }
           else
             {
              Rpz   = R+0.5*Kp*(Rt[Ptr]+(RtZPYX*HZR[z].R1p+Rt[Ptr]*HZR[z].R2p+RtZMYX*HZR[z].R3p));
              Rmz   = R+0.5*Kp*(Rt[Ptr]+(RtZPYX*HZR[z].R1m+Rt[Ptr]*HZR[z].R2m+RtZMYX*HZR[z].R3m));
              Rzz   = Rmz+Rpz;
             }
           Z1z = aR1z*(Rmz*HZR[z].R3m+Rpz*HZR[z].R3p);
           Z2z = aR1z*(Rmz*HZR[z].R1m+Rpz*HZR[z].R1p);
           Z3z = 1+aR1z*(Rzz-Rmz*HZR[z].R2m-Rpz*HZR[z].R2p)-(S==NULL ? 0.0 : TAU/3.0*S[Ptr]);
           
           if (W->WZ1[Ptr]>0.0)
              Z3z += TAU*W->WZ1[Ptr];
           else
              {
               Z1z -= TAU*W->WZ1[Ptr]*HZR[z].R3p;
               Z2z -= TAU*W->WZ1[Ptr]*HZR[z].R1p;
               Z3z += TAU*W->WZ1[Ptr]*HZR[z].R2p;
              }
           if (W->WZ2[Ptr]<=0.0)
              Z3z -= TAU*W->WZ2[Ptr];
           else
              {
               Z1z += TAU*W->WZ2[Ptr]*HZR[z].R3m;
               Z2z += TAU*W->WZ2[Ptr]*HZR[z].R1m;
               Z3z -= TAU*W->WZ2[Ptr]*HZR[z].R2m;
              }
           Z3z = 1 - Z3z;

           if (Flags & fld2RF)
              {
               Rmy = R+Kp*RtZYMX;
               Rpy = R+Kp*RtZYPX;
               Rzy = 2.0*(R+Kp*Rt[Ptr]);
              }
           else
             {
              Rpy   = R+0.5*Kp*(Rt[Ptr]+(RtZYPX*HYR[y].R1p+Rt[Ptr]*HYR[y].R2p+RtZYMX*HYR[y].R3p));
              Rmy   = R+0.5*Kp*(Rt[Ptr]+(RtZYPX*HYR[y].R1m+Rt[Ptr]*HYR[y].R2m+RtZYMX*HYR[y].R3m));
              Rzy   = Rmy+Rpy;
             }
           Z1y = aR1y*(Rmy*HYR[y].R3m+Rpy*HYR[y].R3p);
           Z2y = aR1y*(Rmy*HYR[y].R1m+Rpy*HYR[y].R1p);
           Z3y = 1+aR1y*(Rzy-Rmy*HYR[y].R2m-Rpy*HYR[y].R2p)-(S==NULL ? 0.0 : TAU/3.0*S[Ptr]);
           
           if (W->WY1[Ptr]>0.0)
              Z3y += TAU*W->WY1[Ptr];
           else
              {
               Z1y -= TAU*W->WY1[Ptr]*HYR[y].R3p;
               Z2y -= TAU*W->WY1[Ptr]*HYR[y].R1p;
               Z3y += TAU*W->WY1[Ptr]*HYR[y].R2p;
              }
           if (W->WY2[Ptr]<=0.0)
              Z3y -= TAU*W->WY2[Ptr];
           else
              {
               Z1y += TAU*W->WY2[Ptr]*HYR[y].R3m;
               Z2y += TAU*W->WY2[Ptr]*HYR[y].R1m;
               Z3y -= TAU*W->WY2[Ptr]*HYR[y].R2m;
              }
           Z3y = 1 - Z3y;

           if (Flags & fld2RF)
              {
               Rmx = R+Kp*RtZYXM;
               Rpx = R+Kp*RtZYXP;
               Rzx = 2.0*(R+Kp*Rt[Ptr]);
              }
           else
             {
              Rpx   = R+0.5*Kp*(Rt[Ptr]+(RtZYXP*HXR[x].R1p+Rt[Ptr]*HXR[x].R2p+RtZYXM*HXR[x].R3p));
              Rmx   = R+0.5*Kp*(Rt[Ptr]+(RtZYXP*HXR[x].R1m+Rt[Ptr]*HXR[x].R2m+RtZYXM*HXR[x].R3m));
              Rzx   = Rmx+Rpx;
             }
           Z1x = aR1x*(Rmx*HXR[x].R3m+Rpx*HXR[x].R3p);
           Z2x = aR1x*(Rmx*HXR[x].R1m+Rpx*HXR[x].R1p);
           Z3x = 1+aR1x*(Rzx-Rmx*HXR[x].R2m-Rpx*HXR[x].R2p)-(S==NULL ? 0.0 : TAU/3.0*S[Ptr]);
           
           if (W->WX1[Ptr]>0.0)
              Z3x += TAU*W->WX1[Ptr];
           else
              {
               Z1x -= TAU*W->WX1[Ptr]*HXR[x].R3p;
               Z2x -= TAU*W->WX1[Ptr]*HXR[x].R1p;
               Z3x += TAU*W->WX1[Ptr]*HXR[x].R2p;
              }
           if (W->WX2[Ptr]<=0.0)
              Z3x -= TAU*W->WX2[Ptr];
           else
              {
               Z1x += TAU*W->WX2[Ptr]*HXR[x].R3m;
               Z2x += TAU*W->WX2[Ptr]*HXR[x].R1m;
               Z3x -= TAU*W->WX2[Ptr]*HXR[x].R2m;
              }
           Z3x = 1 - Z3x;

           H1[Ptr] = H[Ptr]+
                       (*VAL(Bounds,_top,&H[ZPYX]))*Z2z+H[Ptr]*Z3z+(*VAL(Bounds,_bottom,&H[ZMYX]))*Z1z+
                       (*VAL(Bounds,_forw,&H[ZYPX]))*Z2y+H[Ptr]*Z3y+(*VAL(Bounds,_back,&H[ZYMX]))*Z1y+
                       (*VAL(Bounds,_right,&H[ZYXP]))*Z2x+H[Ptr]*Z3x+(*VAL(Bounds,_left,&H[ZYXM]))*Z1x+
                     TAU*(K==NULL ? 0.0 : K[Ptr]);
          }
  }
}

void CorrectBounds (float ** Bounds,
                    float         * H,
                    int             NumZ,
                    int             NumY,
                    int             NumX,
                    int             ZStep,
                    int             YStep,
                    int             XStep,
                    DescLine      * Desc,
                    unsigned char * Map,
                    unsigned char * Area,
                    unsigned char   Projection,
                    unsigned char   Direction,
                    int UseOpenMP)
{
 int z,y,zy;
 int Item;

 #pragma omp parallel if(UseOpenMP)
 #pragma omp for schedule(dynamic,imax(1,NumY*NumZ/imax(20,4*nSMP))) private(zy,z,y,Item)
 for (zy=0;zy<NumY*NumZ;zy++)
    {
     z = zy/NumY;
     y = zy%NumY;
     for (Item = 0; Item<Desc[zy].Num; Item++)
         {
          int X1 = Desc[zy].Items[Item].X1;
          int X2 = Desc[zy].Items[Item].X2;
          int D1 = Desc[zy].Dir1;
          int D2 = Desc[zy].Dir2;

          if (X1>=0 && X2>=0)
             {
              int PrevX = (X1+NumX-1) % NumX;
              int PostX = (X2+1)      % NumX;

              int IndexGRN1 = z*ZStep + y*YStep + PrevX*XStep;
              int IndexGRN2 = z*ZStep + y*YStep + PostX*XStep;

              unsigned char TypeGRN1 = MAP(D1,Map,IndexGRN1);
              unsigned char TypeGRN2 = MAP(D2,Map,IndexGRN2);

              if (TypeGRN1==Bound2)
                 *VAL(Bounds,D1,&H[IndexGRN1]) = H[z*ZStep+y*YStep+X1*XStep];
              else if (TypeGRN1==Bound3)
                 *VAL(Bounds,D1,&H[IndexGRN1]) = H[z*ZStep+y*YStep+X1*XStep]/2;

              if (TypeGRN2==Bound2)
                 *VAL(Bounds,D2,&H[IndexGRN2]) = H[z*ZStep+y*YStep+X2*XStep];
              else if (TypeGRN2==Bound3)
                 *VAL(Bounds,D2,&H[IndexGRN2]) = H[z*ZStep+y*YStep+X2*XStep]/2;
             }
         }
    }
}

void AnalyzeBoard (float * H, unsigned char * Map, char Restrict, int UseOpenMP)
{
 DeclareStandard
 int zy;

 #pragma omp parallel if(UseOpenMP)
 #pragma omp for schedule(dynamic,imax(1,NY*NZ/imax(20,4*nSMP))) private(zy,z,y,x,Ptr)
 for (zy=0; zy<NZ*NY; zy++)
    {
     z = zy/NY;
     y = zy%NY;
     for (x=0, Ptr=NX*zy; x<NX; Ptr++,x++)
         if (Map[Ptr]==Zero)
            if (fabs(H[Ptr])<1E-20 ||
                Restrict==rsPositive && H[Ptr]<0.0 && H[Ptr]>-SpecVAL ||
                Restrict==rsNegative && H[Ptr]>0.0 && H[Ptr]<SpecVAL)
               H[Ptr] = 0.0;
    }
}

void RemoveFZero (float * H, unsigned char * Map)
{
 DeclareStandard
 int zy;

 #pragma omp parallel if(UseOpenMP)
 #pragma omp for schedule(dynamic,imax(1,NY*NZ/imax(20,4*nSMP))) private(zy,z,y,x,Ptr)
 for (zy=0; zy<NZ*NY; zy++)
    {
     z = zy/NY;
     y = zy%NY;
     for (x=0, Ptr=NX*zy; x<NX; Ptr++,x++)
         if (Map[Ptr]==Zero)
            if (H[Ptr]==FictZero)
               H[Ptr] = 0.0;
    }
}

/* Процедура выполнения одной итерации над функцией H с граничными       */
/* условиями TypeGRN, полем коэффициентов K, коэффициентом при операторе */
/* Лапласа R, скоростями с предудущей итерации Ux,Uy.                    */
#ifndef __PARALLEL__
void Calculate (float ** Bounds,
                _Solver S,
                float   * H,
                float   * K,
                float   * _S,
                WKoeffs * W,
                HKoeffs * HXR, HKoeffs * HYR,HKoeffs  * HZR,
                float   * HX,  float  * HY,  float  * HZ,
                float   * HXX, float  * HYY, float  * HZZ,
                float   * HXS, float  * HYS, float  * HZS,
                float   * H1,  double * L,   double * M,   double * G,
                float     R,
                float   * Rt,
                float     Kp,
                DescLine * DescH, DescLine * DescP, DescLine * DescV,
                unsigned char * Map,
                unsigned char * Area,
                unsigned char * Boundaries,
                char Restrict,
                char Flags,
                char RegularX, char RegularY, char RegularZ,
                char Projection,
                int  UseOpenMP,
                char Reverse)
{
 if (FastMode && !SlowMode)
    if (!S) S = RozhkovSolver;
    else if (S != RozhkovSolver) S = NULL;
 
 if (S)
    {
     int Iteration = 1;
     int GO;

     do {
       GO = !S(Bounds,Iteration++,H,K,_S,H1,Map,Area,Boundaries,DescH,DescP,DescV,HX,HY,HZ,HXX,HYY,HZZ,HXS,HYS,HZS,L,M,G,Flags,Projection);
       AnalyzeBoard(H,Map,rsAny, UseOpenMP);
       CorrectBounds(Bounds,H, NZ,NX,NY, NY*NX,1,NX, DescP, Map, Area, Projection, 1, UseOpenMP);
       CorrectBounds(Bounds,H, NY,NX,NZ, NX,1,NY*NX, DescV, Map, Area, Projection, 2, UseOpenMP);
       CorrectBounds(Bounds,H, NZ,NY,NX, NY*NX,NX,1, DescH, Map, Area, Projection, 0, UseOpenMP);
     } while (GO);
	 AnalyzeBoard(H, Map, Restrict, UseOpenMP);
	 RemoveFZero(H, Map);
    }
 else if (SlowMode)
    {
     if (Reverse)
        {
         /* Прогонка по Z  */
         OneLineGo(Area,Boundaries,Bounds,H,K,_S,W->WZ1,W->WZ2,HZR,HZ,HZZ,HZS,
                   H1,L,M,G,R, Rt,Kp, NY,NX,NZ, NX,1,NY*NX, DescV, Flags, RegularZ,
                   Map, UseOpenMP);
         AnalyzeBoard(H,Map,Restrict, UseOpenMP);
         CorrectBounds(Bounds,H, NZ,NY,NX, NY*NX,NX,1, DescH,Map,Area, Projection, 0, UseOpenMP);
         CorrectBounds(Bounds,H, NZ,NX,NY, NY*NX,1,NX, DescP,Map,Area, Projection, 1, UseOpenMP);
        }
     else
        {
         /* Прогонка по Y  */
         OneLineGo(Area,Boundaries,Bounds,H,K,_S,W->WY1,W->WY2,HYR,HY,HYY,HYS,
                   H1,L,M,G,R, Rt,Kp, NZ,NX,NY, NY*NX,1,NX, DescP, Flags, RegularY,
                   Map, UseOpenMP);
         AnalyzeBoard(H,Map,Restrict, UseOpenMP);
         CorrectBounds(Bounds,H, NY,NX,NZ, NX,1,NY*NX, DescV,Map,Area, Projection, 2, UseOpenMP);
         CorrectBounds(Bounds,H, NZ,NY,NX, NY*NX,NX,1, DescH,Map,Area, Projection, 0, UseOpenMP);
        }
     /* Прогонка по X */
     OneLineGo(Area,Boundaries,Bounds,H,K,_S,W->WX1,W->WX2,HXR,HX,HXX,HXS,
               H1,L,M,G,R, Rt,Kp, NZ,NY,NX, NY*NX,NX,1, DescH, Flags, RegularX,
               Map, UseOpenMP);
     AnalyzeBoard(H,Map,Restrict, UseOpenMP);
     CorrectBounds(Bounds,H, NY,NX,NZ, NX,1,NY*NX, DescV,Map,Area, Projection, 2, UseOpenMP);
     CorrectBounds(Bounds,H, NZ,NX,NY, NY*NX,1,NX, DescP,Map,Area, Projection, 1, UseOpenMP);
     if (!Reverse)
        {
         /* Прогонка по Z  */
         OneLineGo(Area,Boundaries,Bounds,H,K,_S,W->WZ1,W->WZ2,HZR,HZ,HZZ,HZS,
                   H1,L,M,G,R, Rt,Kp, NY,NX,NZ, NX,1,NY*NX, DescV, Flags, RegularZ,
                   Map, UseOpenMP);
         AnalyzeBoard(H,Map,Restrict, UseOpenMP);
         CorrectBounds(Bounds,H, NZ,NY,NX, NY*NX,NX,1, DescH,Map,Area, Projection, 0, UseOpenMP);
         CorrectBounds(Bounds,H, NZ,NX,NY, NY*NX,1,NX, DescP,Map,Area, Projection, 1, UseOpenMP);
        }
     else
        {
         /* Прогонка по Y  */
         OneLineGo(Area,Boundaries,Bounds,H,K,_S,W->WY1,W->WY2,HYR,HY,HYY,HYS,
                   H1,L,M,G,R, Rt,Kp, NZ,NX,NY, NY*NX,1,NX, DescP, Flags, RegularY,
                   Map, UseOpenMP);
         AnalyzeBoard(H,Map,Restrict, UseOpenMP);
         CorrectBounds(Bounds,H, NY,NX,NZ, NX,1,NY*NX, DescV,Map,Area, Projection, 2, UseOpenMP);
         CorrectBounds(Bounds,H, NZ,NY,NX, NY*NX,NX,1, DescH,Map,Area, Projection, 0, UseOpenMP);
        }
    }
}
#else
void CalculateInXY
               (float  ** Bounds,
                float   * H,
                float   * K,
                float   * _S,
                WKoeffs * W,
                HKoeffs * HXR, HKoeffs * HYR,
                float   * HX,  float  * HY,
                float   * HXX, float  * HYY,
                float   * HXS, float  * HYS,
                float   * H1,  double * L, double * M, double * G,
                float     R,
                float   * Rt,
                float     Kp,
                DescLine      * DescH, DescLine * DescP, DescLine * DescV,
                unsigned char * Map, unsigned char * Area, unsigned char * Boundaries,
                char Restrict,
                char Flags,
                char RegularX, char RegularY,
                char Projection,
                int  UseOpenMP,
                char Reverse)
{
 if (SlowMode)
    {
     int    Mode;

     for (Mode=0; Mode<2; Mode++)
       if (Mode==Reverse)
          {
           /* Прогонка по X */
           OneLineGo(Area,Boundaries,Bounds,H,K,_S,W->WX1,W->WX2,HXR,HX,HXX,HXS,
                     H1,L,M,G,R, &Rt[NY*NX],Kp, NZ,NY,NX, NY*NX,NX,1, DescH, Flags, RegularX,
                     Map, HX[NXs], UseOpenMP);
           AnalyzeBoard(H,Map,Restrict, UseOpenMP);
           CorrectBounds(Bounds,H, NY,NX,NZ, NX,1,NY*NX, DescV,Map,Area, Projection, 2, UseOpenMP);
           CorrectBounds(Bounds,H, NZ,NX,NY, NY*NX,1,NX, DescP,Map,Area, Projection, 1, UseOpenMP);
          }
       else
          {
           /* Прогонка по Y  */
           OneLineGo(Area,Boundaries,Bounds,H,K,_S,W->WY1,W->WY2,HYR,HY,HYY,HYS,
                     H1,L,M,G,R, &Rt[NY*NX],Kp, NZ,NX,NY, NY*NX,1,NX, DescP, Flags, RegularY,
                     Map, HY[NYs], UseOpenMP);
           AnalyzeBoard(H,Map,Restrict, UseOpenMP);
           CorrectBounds(Bounds,H, NY,NX,NZ, NX,1,NY*NX, DescV,Map,Area, Projection, 2, UseOpenMP);
           CorrectBounds(Bounds,H, NZ,NY,NX, NY*NX,NX,1, DescH,Map,Area, Projection, 0, UseOpenMP);
          }
    }
}

/* Флаг разрешения экстраполяции */
int AllowPrediction = 0;
/* Число коэффициентов интерполируещего полинома */
#define NWB    4
/* Число точек, по которым идет интерполяция МНК */
#define NExp   5
/* Количество экстраполируемых точек */
#define NPred  1
/* Количество пропускаемых в начале итераций */
#define ByPass 50
/* Стадия: <0 = Запрет. 0..NExp-1 = сбор информации. NExp..NExp+NPred+1 = экстраполяция */
int Stage = -1;
/* Флаг использования взвешенных коэффициентов для МНК */
int AllowQW = 1;
/* Период контроля сеансов предикции */
#define ReCalcNPs 10
/* Количество экспериментальных точек для анализа при пересчете коэффициентов */
#define ReCalcNPoints 30
/* Максимальное значение коэффициента */
#define MaxKQ 490.0
/* Минимальное значение коэффициента */
#define MinKQ 10.0
/* Начальное значение коэффициента */
#define InitKQ 250.0
/* Шаг изменения коэффициента */
#define StepKQ 20.0

typedef double WBVect[NWB];  /* Вектор значений */
typedef int    WBInds[NWB];  /* Вектор индексов строк для решения СЛАУ */
typedef float  WBVals[NExp+NPred]; /* Вектор значений для контроля коэффициентов */

typedef union {
 /* Коэффициенты полинома */
 WBVect X;
 /* Вектор значений интерполируемой ф-ции для МНК */
 float F[NExp];
} TrVect;

double TAUs[NExp][2*NWB-1];
double sTAUs[NPred][NWB];
double KQ[NExp];
WBVect Koefs[NWB];
WBVect Bs;
WBInds LUInds;
WBVect LUMatr[NWB];

void PreparePrediction(double * KG)
{
 int i,j,k;
 double Q = 0.0;

 for (i=0; i<NExp; i++)
     TAUs[i][0] = KG==NULL ? 1.0 : KG[i];
 for (i=0; i<NExp; i++, Q+=TAU)
   for (j=1; j<2*NWB-1; j++)
     TAUs[i][j] = Q*TAUs[i][j-1]; 
 for (i=0; i<NWB; i++)
   for (j=0; j<NWB; j++)
     for (k=0, Koefs[i][j]=0.0; k<NExp; k++)
       Koefs[i][j]+=TAUs[k][i+j];
 _GetLU(NWB, LUInds, (double *) Koefs, (double *) LUMatr);
}

double CalcPredictErr(WBVals * Trace, double Min, double * KQ, int Index, double Incr)
{
 int i,j,k;
 double Err = 0.0;

 KQ[Index] += Incr; 
 if (KQ[Index]>=MinKQ && KQ[Index]<=MaxKQ)
    {
     PreparePrediction(KQ);
     /* МНК и предикция во всех точках. Суммируем отклонения в Err. */
     for (i=0; i<ReCalcNPoints; i++)
       {
        WBVect X;

        for (j=0; j<NWB; j++)
          for (k=0, Bs[j]=0.0; k<NExp; k++)
            Bs[j] += TAUs[k][j]*Trace[i][k];
        _SolveLU(NWB, LUInds, (double *) LUMatr, (double *) Bs, (double *) X);
        for (k=0; k<NPred; k++)
            {
             double V = 0.0;

             for (j=0; j<NWB; j++)
                 V += X[j]*sTAUs[k][j];
             Err += fabs(1-V/Trace[i][NExp+k]);
            }
       }
    }
 else
    Err = Min;
 KQ[Index] -= Incr;

 return min(Err,Min);
}

void PredictBoard(unsigned char * Map, TrVect * Trace, float * H)
{
 int x, y, j, k;
 int Ptr;
 
 for (y=0, Ptr=0; y<NY; y++)
     for (x=0; x<NX; x++, Ptr++)
         if (Map[Ptr]==ExchngBound)
            {
             if (Stage==NExp)
               {
                for (j=0; j<NWB; j++)
                  for (k=0, Bs[j]=0.0; k<NExp; k++)
                    Bs[j] += TAUs[k][j]*Trace[Ptr].F[k];
                _SolveLU(NWB, LUInds, (double *) LUMatr, (double *) Bs, (double *) Trace[Ptr].X);
               }
             for (j=0, H[Ptr]=0.0; j<NWB; j++)
                 H[Ptr] += Trace[Ptr].X[j]*sTAUs[Stage-NExp][j];
            }
}

void CalculateInZ
               (float   ** Bounds,
                _Solver S,
                float    * DH, float * H, float * UH,
                TrVect   * DTrace, TrVect * UTrace,
                float    * K,
                float    * _S,
                WKoeffs  * W,
                HKoeffs * HXR,  HKoeffs * HYR,HKoeffs * HZR,
                float    * HX,  float  * HXX, float  * HXS,
                float    * HY,  float  * HYY, float  * HYS,
                float    * HZ,  float  * HZZ, float  * HZS,
                float    * H1,  double * L,   double * M,   double * G,
                float      R,
                float    * Rt, /* (NZ+2)*NY*NX !!! */
                float      Kp,
                DescLine * DescH, DescLine * DescP, DescLine * DescV,
                unsigned char * Map, unsigned char * Area, unsigned char * Boundaries,
                char Restrict,
                char Flags,
                char RegularZ,
                char Projection,
                int  UseOpenMP)
{
 char BND = Flags>>shBoundary;
 double aR1;
 int    x, y;
 int    Ptr;

 if (FastMode && !SlowMode)
    if (!S) S = RozhkovSolver;
    else if (S != RozhkovSolver) S = NULL;

 memmove(H1,H,BoardSize);
 if (SlowMode)
    if (Stage<NExp)
       {
        if (DownExchange)
           {
            if (Stage>=0)
               for (y=0, Ptr=0; y<NY; y++)
                   for (x=0; x<NX; x++, Ptr++)
                       if (Map[Ptr]==ExchngBound)
                          DTrace[Ptr].F[Stage] = H1[Ptr];
            if (S==NULL)
               {
                /* Предвычисление по схеме Головичева на нижней границе */
                aR1  = TAU/HZR[0].h/HZR[0].h;
                #pragma omp parallel if(UseOpenMP)
                #pragma omp for schedule(dynamic,imax(4,NY/8)) private(y,x,Ptr)
                for (y=0; y<NY; y++)
                    for (x=0, Ptr=y*NX; x<NX; x++, Ptr++)
                        if (Map[Ptr]==ExchngBound)
                           {
                            double Z1,Z2,Z3,Rp,Rm,Rz;

                            if (Boundaries[Ptr] && Boundaries[Ptr]<=BND)
                               H[Ptr] = K[Ptr];
                            else
                               {
                                if (Flags & fld2RF)
                                   {
                                    Rm = R+Kp*Rt[Ptr];
                                    Rp = R+Kp*Rt[Ptr+2*NY*NX];
                                    Rz = 2.0*(R+Kp*Rt[Ptr+NY*NX]);
                                   }
                                else
                                  {
                                   Rp   = R+0.5*Kp*(Rt[Ptr+NY*NX]+(Rt[Ptr+2*NY*NX]*HZR[0].R1p+Rt[Ptr+NY*NX]*HZR[0].R2p+Rt[Ptr]*HZR[0].R3p));
                                   Rm   = R+0.5*Kp*(Rt[Ptr+NY*NX]+(Rt[Ptr+2*NY*NX]*HZR[0].R1m+Rt[Ptr+NY*NX]*HZR[0].R2m+Rt[Ptr]*HZR[0].R3m));
                                   Rz   = Rm+Rp;
                                  }
                                Z1 = aR1*(Rm*HZR[0].R3m+Rp*HZR[0].R3p);
                                Z2 = aR1*(Rm*HZR[0].R1m+Rp*HZR[0].R1p);
                                Z3 = 1+aR1*(Rz-Rm*HZR[0].R2m-Rp*HZR[0].R2p)-(_S==NULL ? 0.0 : TAU/3.0*_S[Ptr]);
                                
                                if (W->WZ1[Ptr]>0.0)
                                   Z3 += TAU*W->WZ1[Ptr];
                                else
                                   {
                                    Z1 -= TAU*W->WZ1[Ptr]*HZR[0].R3p;
                                    Z2 -= TAU*W->WZ1[Ptr]*HZR[0].R1p;
                                    Z3 += TAU*W->WZ1[Ptr]*HZR[0].R2p;
                                   }
                                if (W->WZ2[Ptr]<=0.0)
                                   Z3 -= TAU*W->WZ2[Ptr];
                                else
                                   {
                                    Z1 += TAU*W->WZ2[Ptr]*HZR[0].R3m;
                                    Z2 += TAU*W->WZ2[Ptr]*HZR[0].R1m;
                                    Z3 -= TAU*W->WZ2[Ptr]*HZR[0].R2m;
                                   }
                                H[Ptr] = (H1[Ptr]+H1[Ptr+NY*NX]*Z2+DH[Ptr]*Z1+
                                          (TAU/3.0)*(K==NULL ? 0.0 : K[Ptr])
                                         )/Z3;
                               }
                           }
               }
           }
        if (UpExchange)
           {
            if (Stage>=0)
               for (y=0, Ptr=NZs*NY*NX; y<NY; y++)
                   for (x=0; x<NX; x++, Ptr++)
                       if (Map[Ptr]==ExchngBound)
                          UTrace[Ptr-NZs*NY*NX].F[Stage] = H1[Ptr];
            if (S==NULL)
               {
                /* Предвычисление по схеме Головичева на верхней границе */
                aR1  = TAU/HZR[NZs].h/HZR[NZs].h;
                #pragma omp parallel if(UseOpenMP)
                #pragma omp for schedule(dynamic,imax(4,NY/8)) private(y,x,Ptr)
                for (y=0; y<NY; y++)
                    for (x=0, Ptr=(NZs*NY+y)*NX; x<NX; x++, Ptr++)
                        if (Map[Ptr]==ExchngBound)
                           {
                            double Z1,Z2,Z3,Rp,Rm,Rz;

                            if (Boundaries[Ptr] && Boundaries[Ptr]<=BND)
                               H[Ptr] = K[Ptr];
                            else
                               {
                                if (Flags & fld2RF)
                                   {
                                    Rm = R+Kp*Rt[Ptr];
                                    Rp = R+Kp*Rt[Ptr+2*NY*NX];
                                    Rz = 2.0*(R+Kp*Rt[Ptr+NY*NX]);
                                   }
                                else
                                  {
                                   Rp   = R+0.5*Kp*(Rt[Ptr+NY*NX]+(Rt[Ptr+2*NY*NX]*HZR[NZs].R1p+Rt[Ptr+NY*NX]*HZR[NZs].R2p+Rt[Ptr]*HZR[NZs].R3p));
                                   Rm   = R+0.5*Kp*(Rt[Ptr+NY*NX]+(Rt[Ptr+2*NY*NX]*HZR[NZs].R1m+Rt[Ptr+NY*NX]*HZR[NZs].R2m+Rt[Ptr]*HZR[NZs].R3m));
                                   Rz   = Rm+Rp;
                                  }
                                Z1 = aR1*(Rm*HZR[NZs].R3m+Rp*HZR[NZs].R3p);
                                Z2 = aR1*(Rm*HZR[NZs].R1m+Rp*HZR[NZs].R1p);
                                Z3 = 1+aR1*(Rz-Rm*HZR[NZs].R2m-Rp*HZR[NZs].R2p)-(_S==NULL ? 0.0 : TAU/3.0*_S[Ptr]);
                                if (W->WZ1[Ptr]>0.0)
                                   Z3 += TAU*W->WZ1[Ptr];
                                else
                                   {
                                    Z1 -= TAU*W->WZ1[Ptr]*HZR[NZs].R3p;
                                    Z2 -= TAU*W->WZ1[Ptr]*HZR[NZs].R1p;
                                    Z3 += TAU*W->WZ1[Ptr]*HZR[NZs].R2p;
                                   }
                                if (W->WZ2[Ptr]<=0.0)
                                   Z3 -= TAU*W->WZ2[Ptr];
                                else
                                   {
                                    Z1 += TAU*W->WZ2[Ptr]*HZR[NZs].R3m;
                                    Z2 += TAU*W->WZ2[Ptr]*HZR[NZs].R1m;
                                    Z3 -= TAU*W->WZ2[Ptr]*HZR[NZs].R2m;
                                   }
                                H[Ptr] = (H1[Ptr]+UH[Ptr-NZs*NY*NX]*Z2+H1[Ptr-NY*NX]*Z1+
                                          (TAU/3.0)*(K==NULL ? 0.0 : K[Ptr])
                                         )/Z3;
                               }
                           }
               }
           }
       }
    else
       {
        if (DownExchange) PredictBoard(Map, DTrace, H);
        if (UpExchange) PredictBoard(&Map[NZs*NY*NX], UTrace, &H[NZs*NY*NX]);
       }
 if (S)
    {
     int Iteration = 1;
     int GO;

     do {
       GO = !S(Bounds,Iteration++,H,K,_S,H1,Map,Area,Boundaries,DescH,DescP,DescV,HX,HY,HZ,HXX,HYY,HZZ,HXS,HYS,HZS,L,M,G,Flags,Projection);
       AnalyzeBoard(H,Map,rsAny,UseOpenMP);
       CorrectBounds(Bounds,H, NZ,NX,NY, NY*NX,1,NX, DescP, Map,Area, Projection, 1, UseOpenMP);
       CorrectBounds(Bounds,H, NY,NX,NZ, NX,1,NY*NX, DescV, Map,Area, Projection, 2, UseOpenMP);
       CorrectBounds(Bounds,H, NZ,NY,NX, NY*NX,NX,1, DescH, Map,Area, Projection, 0, UseOpenMP);
     } while (GO);
	 AnalyzeBoard(H,Map,Restrict,UseOpenMP);
	 RemoveFZero(H,Map);
    }
 else if (SlowMode)
    {
     /* Прогонка по Z */
     OneLineGo(Area,Boundaries,Bounds,H,K,_S,W->WZ1,W->WZ2,HZR,HZ,HZZ,HZS,
               H1,L,M,G,R, &Rt[NY*NX],Kp, NY,NX,NZ, NX,1,NY*NX, DescV, Flags, RegularZ,
               Map, HZbf, UseOpenMP);
     AnalyzeBoard(H,Map,Restrict, UseOpenMP);
     CorrectBounds(Bounds,H, NZ,NY,NX, NY*NX,NX,1, DescH,Map,Area, Projection, 0, UseOpenMP);
     CorrectBounds(Bounds,H, NZ,NX,NY, NY*NX,1,NX, DescP,Map,Area, Projection, 1, UseOpenMP);
    }
}
#endif

/* Коэффициенты для Ux, Uy, Uz, P, Nu, (T и Dn) */
float * Kfs[NumEqs];
float * KDn[MaxActSubst];
float * Sfs[NumEqs];
float * SDn[MaxActSubst];

/* Вспомогательные данные */
float * UzSave  = NULL;

int BaseZ = 0; /* Номер узла -- начальной координаты Z участка процессора */

int TotalNParts = 0; /* Общее число частиц */
int NParts = 0; /* Локальное число частиц */

/* Вспомогательные данные для проверки точности */
StoreStruct Templ;
StoreStruct Old;
/* Вычисляемые данные (текущая итерация) */
StoreStruct Vars;

/* Поле критичности для проверки точности */
float * GRAV = NULL;
/* Поле неявности для разностной схемы (по результатам проверки точности) */
float * ALPHA = NULL;
/* Поле максимальных ошибок */
float * ERRS = NULL;

/* Промежуточные значения H */
float * H1 = NULL;

/* Коэффициенты прогонки */
double * L  = NULL;
double * M  = NULL;
double * G  = NULL;

/* Карты */
MapItem Maps[NumEqs];

unsigned char  * Area  = NULL;
unsigned char  * Boundaries = NULL;
unsigned char  * CMap  = NULL;
float * CVals = NULL;

int     SubstValNum;

DescLine * DescH = NULL;
DescLine * DescP = NULL;
DescLine * DescV = NULL;

#ifdef __PARALLEL__
void SendCommand(int Command)
{
 int i;

 for (i=1; i<TotalProcs; i++)
     SendSlave(i-1,(byte *) &Command,sizeof(Command));
}

void RecvCommand(int * Command)
{
 RecvMaster((byte *) Command, sizeof(int));
}

typedef struct {
  int    _RecvCommand;
  int    RecvData;
  int    SendData;
  int    CheckCommand;
  double _TAU;
  int    _NReact;
  int    _takt;
  int    _CalcBase;
} HandlePacket;

void SendPacket(int _RecvCommand, int RecvData, int SendData,
                int CheckCommand, double _TAU, int _NReact)
{
 HandlePacket Buf;
 int i;

 Buf._RecvCommand = _RecvCommand;
 Buf.RecvData     = RecvData;
 Buf.SendData     = SendData;
 Buf.CheckCommand = CheckCommand;
 Buf._TAU         = _TAU;
 Buf._NReact      = _NReact;
 Buf._takt        = takt;
 Buf._CalcBase    = CalcBase;
 for (i=1; i<TotalProcs; i++)
     SendSlave(i-1,(byte *) &Buf,sizeof(Buf));
}

void RecvPacket(HandlePacket * Packet)
{
 RecvMaster((byte *) Packet, sizeof(HandlePacket));
}

void MoveToBuf(int * Offset, int slave, float * Board)
{
 memmove(&HBuff[*Offset], &Board[((FromGrid[slave]+NZs) % NZ)*NY*NX], NX*NY*sizeof(float));
 (*Offset) += NY*NX;
 memmove(&HBuff[*Offset], &Board[FromGrid[slave]*NY*NX], LengthGrid[slave]*NX*NY*sizeof(float));
 (*Offset) += LengthGrid[slave]*NX*NY;
 memmove(&HBuff[*Offset], &Board[((FromGrid[slave]+LengthGrid[slave]) % NZ)*NY*NX], NX*NY*sizeof(float));
 (*Offset) += NY*NX;
}

void MoveFromBuf(int * Offset, int slave, float * Board)
{
 if (slave==0)
    {
     memmove(&Board[FromGrid[slave]*NY*NX], &HBuff[*Offset], (LengthGrid[slave]-1)*NX*NY*sizeof(float));
     (*Offset) += (LengthGrid[slave]-1)*NX*NY;
    }
 else if (slave==nProcs-2)
    {
     memmove(&Board[(FromGrid[slave]+1)*NY*NX], &HBuff[*Offset], (LengthGrid[slave]-1)*NX*NY*sizeof(float));
     (*Offset) += (LengthGrid[slave]-1)*NX*NY;
    }
 else
    {
     memmove(&Board[(FromGrid[slave]+1)*NY*NX], &HBuff[*Offset], (LengthGrid[slave]-2)*NX*NY*sizeof(float));
     (*Offset) += (LengthGrid[slave]-2)*NX*NY;
    }
}
#endif

#ifndef __PARALLEL__
void MoveBlocks(StoreStruct * Dest, StoreStruct * Src)
{
 int i;

 for (i=0; i<NumEqs; i++)
     memmove(Dest->Name[i],Src->Name[i],BoardSize);
 for (i=0;i<NSubst;i++)
     memmove(Dest->Dn[i],Src->Dn[i],BoardSize);
 if (CheckTau)
    for (i = 0; i < 3; i++) {
        memmove(Dest->Parts[i], Src->Parts[i], NParts*sizeof(float));
        memmove(Dest->PartsV[i], Src->PartsV[i], NParts*sizeof(float));
    }
 for (i=_back;i<=_top;i++)
     memmove(Dest->Bounds[i],Src->Bounds[i],Src->BoundSize);
}
#endif

#define CHK_NoCheck      0
#define CHK_Phase1       1
#define CHK_BeginPhase2  2
#define CHK_Phase2       3
#define CHK_EndPhase2    4
#define CHK_Recalculate  5

int TranslateBoundaries(float * Val, int Offset) {
 if (*Val != FictZero && fabs(*Val) >= SpecVAL)
    {
     int POW;
     float Mantissa = frexp(*Val, &POW);
     POW -= 100;
     if (POW >= 0 && POW <= 9) {
         int LOW = Mantissa>0 ? ((int)(((double) (Mantissa*1000000.0))+0.5)) : ((int)(((double) (-Mantissa*1000000.0))+0.5))-500000;
         int locPtr = POW*1000000+LOW-1;
         int glbPtr = Offset + locPtr;

         int CODE = glbPtr+1; // 128 Eqs * 65536 points; Min = 1; Max = 8388480+1
         int FIRST = CODE/1000000; // 0..9
         int REM = CODE%1000000; // Min = 0; Max = 999999
         // 500000..999999 => 0.500000..0.999999
         // 1..499999 => -[0.500001..0.999999]
         Mantissa = REM>=500000 ? 0.000001*REM : -0.000001*(REM+500000);
         POW = 100+FIRST; // Max = 109
         // ABS[max] = max(Mantissa)*2^max(POW) = 6.4903645828*1E32
         // ABS[min] = ABS[min(Mantissa)*2^min(POW)] = 6.3382530011*1E29
         // Critical Number = +-1E25, for example...
         *Val = ldexp(Mantissa,POW); // rounding must present!!!!!
         return 1;
     }
    }
 return 0;
}

#ifndef __PARALLEL__

int TraceCmp(const void * A0, const void * B0)
{
	TraceTypeHost * A = (TraceTypeHost *)A0;
	TraceTypeHost * B = (TraceTypeHost *)B0;
	int d = B->Iters - A->Iters;
	int i;

	if (abs(d) > 0) return d;

	for (i = 0; i < markZ99; i++) {
		int p = B->KinetTraceWOrder[i];
		int q = (B->Tr[p] - A->Tr[p]);
		if (q != 0) return q;
	}

	return 0;
}

int TraceWCmp(const void * A0, const void * B0)
{
	TraceWSorter * A = (TraceWSorter *)A0;
	TraceWSorter * B = (TraceWSorter *)B0;
	float d = B->W - A->W;

	return d < 0 ? -1 : (d > 0 ? +1 : 0);
}

float determ(float a1, float b1, float a2, float b2) {
	return a1*b2 - a2*b1;
}

void GetPredictor(int p, int Iters3, int Iters2, int Iters1, int Iters, float * W1, float * Q) {
	float a1, b1, a2, b2, s1, s2;
	if (p == 3) {
		a1 = Iters3*Iters3 + Iters2*Iters2 + Iters1*Iters1;
		b1 = Iters3 + Iters2 + Iters1;
		a2 = b1;
		b2 = 3;
		s1 = Iters2*Iters3 + Iters1*Iters2 + Iters*Iters1;
		s2 = Iters2 + Iters1 + Iters;
	} else if (p == 2) {
		a1 = Iters2*Iters2 + Iters1*Iters1;
		b1 = Iters2 + Iters1;
		a2 = b1;
		b2 = 2;
		s1 = Iters1*Iters2 + Iters*Iters1;
		s2 = Iters1 + Iters;
	} else if (p == 1) {
		*W1 = 1.0;
		*Q = 0.0;
		return;
	}

	float d = determ(a1, b1, a2, b2);
	if (fabsf(d) < 1E-5) {
		*W1 = 1.0;
		*Q = 0.0;
	} else {
		*W1 = determ(s1, b1, s2, b2) / d;
		*Q = determ(a1, s1, a2, s2) / d;
	}
}

reenterable processKinetics(int init /* !!! 0 -- work, 1 -- GPUinit, 2 -- OpenMPinit */, int k, int group_size, int NKinets, int KinetChunkSize, TraceTypeHost * Traces,
	KineticContext * Contexts, _global(1) KineticGlobal * KGlobal,
	float * T, void * _Dn, int UseGear, int Method, _local(1) KineticContext * C) {
	if (init) {
#pragma pack(push, 1)
		typedef union { __global void * p1; __global float ** p2; } converter;
		converter c;
#pragma pack(pop)
	        __global float ** Dn;
		int i, j, s, p;

		c.p1 = _Dn;
	        Dn = c.p2;

		if (group_size != 0x7FFFFFFF) {
			int offs = k - group_size;
			for (j = 0; j < group_size; j++, offs++) {
				if (!Contexts[j].LossPrecision)
					for (i = 0; i < KGlobal->NASubst; i++)
						Dn[KGlobal->TranMap[i]][Traces[offs].tag] = Contexts[j].Conc1[i];
				Traces[offs].Iters = Contexts[j].Iters;
				for (s = HistorySize-1; s > 0; s--)
					Traces[offs].Timings[s] = Traces[offs].Timings[s-1];
				Traces[offs].Timings[0] = (init == 1 ? last_execution_time(NULL) : 0.0f);
				memmove(Traces[offs].Tr, Contexts[j].Trace, sizeof(TraceType));
			}
		} else {
			gmemset(KGlobal->aMNK, 0, sizeof(KGlobal->aMNK));
			gmemset(KGlobal->bMNK, 0, sizeof(KGlobal->bMNK));
		}

		group_size = 0;

		for (; k < NKinets; k++) {
			Contexts[group_size].LossPrecision = 0;
			if (T)
				SetTk(&Contexts[group_size], T[Traces[k].tag]);
			else
				SetTk(&Contexts[group_size], KGlobal->DefaultT);
			for (i = 0; i < KGlobal->NASubst; i++)
				Contexts[group_size].Conc0[i] = Dn[KGlobal->TranMap[i]][Traces[k].tag];
			plan_last(0 /* !!! */, 0, 0, NKinets, KinetChunkSize, Traces, Contexts, KGlobal, T, _Dn, UseGear, Method, Contexts+group_size);
			
			group_size++;
			if (k == NKinets-1 || group_size == KinetChunkSize) {
/* !!! */
				if (k == NKinets-1) {
					int p;
					for (p = group_size; p < KinetChunkSize; p++) {
						Contexts[p].LossPrecision = 0;
						if (T)
							SetTk(&Contexts[p], T[Traces[k].tag]);
						else
							SetTk(&Contexts[p], KGlobal->DefaultT);
						for (i = 0; i < KGlobal->NASubst; i++)
							Contexts[p].Conc0[i] = Dn[KGlobal->TranMap[i]][Traces[k].tag];
						plan_last(0 /* !!! */, 0, 0, NKinets, KinetChunkSize, Traces, Contexts, KGlobal, T, _Dn, UseGear, Method, Contexts+p);
					}
				}
				plan_group_last;
				plan_last(init /* !!! */, k+1, group_size, NKinets, KinetChunkSize, Traces, Contexts, KGlobal, T, _Dn, UseGear, Method, Contexts);

				if (init == 1) {
					plan_group_vectorize(NULL);
				} else {
					plan_group_parallelize;
				}

				return;
			}
		}

		if (KGlobal->KinetTraceWNotFound) { 
			for (j = 0; j < NKinets; j++) {
				for (s = 0; s < markZ99; s++) {
					KGlobal->bMNK[s] += Traces[j].Iters*Traces[j].Tr[s];
					for (p = s; p < markZ99; p++)
						KGlobal->aMNK[s][p] += Traces[j].Tr[p]*Traces[j].Tr[s];
				}
			}
			for (s = 0; s < markZ99; s++) {
				for (p = s; p < markZ99; p++)
					KGlobal->aMNK[p][s] = KGlobal->aMNK[s][p];
			}

			for (s = 0; s < markZ99; s++) {
				float pp = fabsf(KGlobal->bMNK[s]);
				for (p = 0; p < markZ99; p++) {
					float p1 = fabsf(KGlobal->aMNK[s][p]);
					if (p1 > pp) pp = p1;
				}
				if (pp > 1.0f) {
					for (p = 0; p < markZ99; p++)
						KGlobal->aMNK[s][p] /= pp;
					KGlobal->bMNK[s] /= pp;
				}
			}

			for (k = 0; k < 5; k++) {
				float alpha = 0.02f;
				float f = 0.0f;
				for (s = 0; s < markZ99; s++) {
					float it = -KGlobal->bMNK[s];
					for (p = 0; p < markZ99; p++)
						it += KGlobal->aMNK[s][p]*KGlobal->KinetTraceW[p];
					f += it*it;
				}
				for (j = 0; j < markZ99; j++) {
					float delta = 0.001f;
					float f1 = 0.0f;

					KGlobal->KinetTraceW[j] += delta;
					for (s = 0; s < markZ99; s++) {
						float it = -KGlobal->bMNK[s];
						for (p = 0; p < markZ99; p++)
							it += KGlobal->aMNK[s][p]*KGlobal->KinetTraceW[p];
						f1 += it*it;
					}
					KGlobal->KinetTraceW[j] -= delta;
					KGlobal->gradMNK[j] = (f1 - f)/delta;
				}
				for (s = 0; s < markZ99; s++)
					KGlobal->KinetTraceW[s] -= alpha*KGlobal->gradMNK[s];
			}
		}

		#pragma omp parallel for schedule(guided) private(j,k,s)
		for (j = 0; j < NKinets; j++) {
			__global unsigned short int * buf = Traces[j].Tr3;
			float f, maxt = 0.000000001f;
			int errs[3];
			short int ks;
			signed char nump;

			for (k = 0; k < HistorySize; k++)
        	            if (Traces[j].Timings[k] > maxt)
				maxt = Traces[j].Timings[k];

			for (ks = 0; ks < markZ99; ks++) {
				k = KGlobal->KinetTraceWOrder[ks];

				errs[0] = abs(Traces[j].Prognosed0[k] - Traces[j].Tr[k]);
				errs[1] = abs(Traces[j].Prognosed1[k] - Traces[j].Tr[k]);
				errs[2] = abs(Traces[j].Prognosed2[k] - Traces[j].Tr[k]);

				if (errs[0] < errs[1])
					if (errs[0] < errs[2])
						nump = 0;
					else
						nump = 2;
				else
					if (errs[2] < errs[1])
						nump = 2;
					else
						nump = 1;

				for (s = HistorySize-1; s > 0; s--)
					Traces[j].Besters[s][k] = Traces[j].Besters[s-1][k];
				Traces[j].Besters[0][k] = nump;

				if (ks >= (markZ99 / 4) || fabsf(KGlobal->KinetTraceW[k]) <= 1E-6f) {
					for (s = HistorySize-1; s > 0; s--)
						Traces[j].Winners[s][k] = Traces[j].Winners[s-1][k];
					Traces[j].Winners[0][k] = 0;

					Traces[j].Prognosed0[k] = Traces[j].Prognosed1[k] = Traces[j].Prognosed2[k] = Traces[j].Tr[k];

					buf[k] = Traces[j].Tr[k];
				} else {
					float Quality[3] = { 0.0f, 0.0f, 0.0f };
					float W1[3];
					float Q[3];

					for (s = HistorySize-1; s > 0; s--) {
						char w = Traces[j].Winners[s][k];
						float t = Traces[j].Timings[s];
						float d = t > 1E-9f ? maxt/t : maxt/1E-9f;
						if (w == Traces[j].Besters[s][k])
							Quality[w] += d;
						else
							Quality[w] -= d;
					}

					if (Quality[0] > Quality[1])
						if (Quality[0] > Quality[2])
							nump = 0;
						else
							nump = 2;
					else
						if (Quality[2] > Quality[1])
							nump = 2;
						else
							nump = 1;

					GetPredictor(
						1, Traces[j].Tr3[k], Traces[j].Tr2[k], Traces[j].Tr1[k], Traces[j].Tr[k],
						&W1[0], &Q[0]
					);
					Traces[j].KinetTraceW1[k] = W1[0];
					Traces[j].KinetTraceQ1[k] = Q[0];
					f = Traces[j].Tr[k]*W1[0] + Q[0];
					if (f < 0.0f) Traces[j].Prognosed0[k] = 0;
					else Traces[j].Prognosed0[k] = (unsigned short int) f;

					GetPredictor(
						2, Traces[j].Tr3[k], Traces[j].Tr2[k], Traces[j].Tr1[k], Traces[j].Tr[k],
						&W1[1], &Q[1]
					);
					Traces[j].KinetTraceW2[k] = W1[1];
					Traces[j].KinetTraceQ2[k] = Q[1];
					f = Traces[j].Tr[k]*W1[1] + Q[1];
					if (f < 0.0f) Traces[j].Prognosed1[k] = 0;
					else Traces[j].Prognosed1[k] = (unsigned short int) f;

					GetPredictor(
						3, Traces[j].Tr3[k], Traces[j].Tr2[k], Traces[j].Tr1[k], Traces[j].Tr[k],
						&W1[2], &Q[2]
					);
					Traces[j].KinetTraceW3[k] = W1[2];
					Traces[j].KinetTraceQ3[k] = Q[2];
					f = Traces[j].Tr[k]*W1[2] + Q[2];
					if (f < 0.0f) Traces[j].Prognosed2[k] = 0;
					else Traces[j].Prognosed2[k] = (unsigned short int) f;

					for (s = HistorySize-1; s > 0; s--)
						Traces[j].Winners[s][k] = Traces[j].Winners[s-1][k];
					Traces[j].Winners[0][k] = nump;

					buf[k] = nump == 0 ? Traces[j].Prognosed0[k] : (nump == 1 ? Traces[j].Prognosed1[k] : Traces[j].Prognosed2[k]);
				}
			}
			Traces[j].Tr3 = Traces[j].Tr2;
			Traces[j].Tr2 = Traces[j].Tr1;
			Traces[j].Tr1 = Traces[j].Tr;
			Traces[j].Tr = buf;
		}
	}
	else { /* !!! */
		OneTaktKinetic(KGlobal, C, UseGear, Method);
		if (C->LossPrecision && UseGear) {
			OneTaktKinetic(KGlobal, C, 0, Method);
		}
	}
}

#endif

#ifdef __PARALLEL__
void OneTakt(int RecvCommand, int RecvData, int SendData, int CheckCommand, double * MaxEps)
#else
void OneTakt()
#endif
{
 int        i,j,k,x,y,z;
 int        Ptr;
#ifdef __PARALLEL__
 int        Offset;
 int        slave;
 int        PartL;
 double     LocalEps;
#ifdef __MPI__
 MPI_Status Status;
#endif
#endif

#ifdef __PARALLEL__
 SendPacket(RecvCommand, RecvData, SendData, CheckCommand, TAU, NReact);
 /* Посылаем ВСЕ частичные поля, необходимые для обработки: Ux,Uy,Uz,P,T,Nu,Dn */
 if (RecvData)
    {
     for (i=1; i<nProcs; i++)
         {
          Offset = 0;
          for (j=0; j<NumEqs; j++)
              MoveToBuf(&Offset,i-1,Vars.Name[j]);
#if defined(__MPI__) || defined(__ROUTER__) || defined(__ROUTER_100__)
          SendSlave(i-1,(byte *) HBuff, (long) NumEqs*(LengthGrid[i-1]+2)*NX*NY*sizeof(float));
         }
     for (i=1; i<nProcs; i++)
         {
          Offset = 0;
#endif
          for (j=0; j<NSubst; j++)
              MoveToBuf(&Offset,i-1,Vars.Dn[j]);
#if defined(__MPI__) || defined(__ROUTER__) || defined(__ROUTER_100__)
          if (NSubst)
             SendSlave(i-1,(byte *) HBuff, (long) NSubst*(LengthGrid[i-1]+2)*NX*NY*sizeof(float));
#else
          SendSlave(i-1,(byte *) HBuff, (long) (NSubst+NumEqs)*(LengthGrid[i-1]+2)*NX*NY*sizeof(float));
#endif
         }
    }
#else
 if (NSubst)
    memmove(UzSave,Vars.Name[PhaseVars[CarrierPhase]._Uz],BoardSize);

 for (i=0; i<NumPhases; i++)
     {
      if (!PhaseVars[i].IsLight)
         CalculateWXY((float **) Vars.Bounds,Area,&WXYZ[i],Vars.Name[PhaseVars[i]._Ux],Vars.Name[PhaseVars[i]._Uy]);
      CalculateWZ((float **) Vars.Bounds,Area,&WXYZ[i],Vars.Name[PhaseVars[i]._Uz],*PhaseVars[i]._Uw, UseOpenMP);
     }
 /* Скопировать турбулентную вязкость в рабочий массив */
 if (PhaseVars[CarrierPhase]._Nu>=0)
    memmove(Dt,Vars.Name[PhaseVars[CarrierPhase]._Nu],BoardSize);
 else
    memset(Dt,0,BoardSize);

 SlowMode = 1;
 FastMode = FastTAUDivider==1;
 for (FastIteration = 0; FastIteration<FastTAUDivider+(FastTAUDivider>1); FastIteration++)
     {
      if (FastIteration==1)
        {
         FastMode = 1;
         SlowMode = 0;
         TAU /= FastTAUDivider;
        }
      CalculateDIV(Vars.Bounds, Area, Vars.Name[PhaseVars[CarrierPhase]._Ux], Vars.Name[PhaseVars[CarrierPhase]._Uy], Vars.Name[PhaseVars[CarrierPhase]._Uz]);
      CalculateK(&WXYZ[CarrierPhase],&Vars,Kfs,KDn,Sfs,SDn,Area,Boundaries,Maps);

      if (CheckTau && CheckNew && SlowMode) {
         const double h  = 1.0;
         const int    MP = 5;

         memset(ALPHA, 0, BoardSize);
         for (i = 0; i < NParts; i++) {
             x = (int)(Vars.Parts[0][i] / h);
             y = (int)(Vars.Parts[1][i] / h);
             z = (int)(Vars.Parts[2][i] / h);
             ALPHA[z*NY*NX+y*NX+x] += 1.0;
         }
         for (z = 0, Ptr = 0; z < NZ; z++)
             for (y = 0; y < NY; y++)
                 for (x = 0; x < NX; x++, Ptr++) {
                     int q = (int)(0.5 + ALPHA[Ptr]);
                     ALPHA[Ptr] = 1.0 - 0.5*(q <= 1 ? 0.0 : (q < MP ? 1.0*q/MP : 1.0));
                 }

         if (CalcBase) {
            for (i=0; i<NumEqs; i++)
                if (VDefs[i].Solver == NULL)
                   EulerCalculate((float **) Vars.Bounds,
                                  Vars.Name[i], Kfs[i], Sfs[i], &WXYZ[VDefs[i].Phase],
                                  HXR, HYR, HZR,
                                  HX, HY, HZ,
                                  HXX, HYY, HZZ,
                                  HXS, HYS, HZS,
                                  Vars.eName[i],
                                  *VDefs[i]._NuMol,Dt,*VDefs[i]._Kappa,
                                  VDefs[i].Flags,
                                  ALPHA, Maps[i].Map, UseOpenMP
                   );
            if (NSubst > 0) {
               int EnhanceOpenMP = UseOpenMP && UseEnhancedOpenMP && NLightSubst>0 && (NLightSubst%nSMP || NLightSubst>3*nSMP);
               #pragma omp parallel if (EnhanceOpenMP)
               {
                #pragma omp for schedule(dynamic) private(i,j)
                for (j=0;j<NLightSubst;j++)
                    {
                     i = MapSubsts[j];
                     EulerCalculate((float **) Vars.Bounds,
                                    Vars.Dn[i], KDn[i], SDn[i], &WXYZ[CarrierPhase],
                                    HXR, HYR, HZR,
                                    HX, HY, HZ,
                                    HXX, HYY, HZZ,
                                    HXS, HYS, HZS,
                                    Vars.eDn[i],
                                    D, Dt, 0.0,
                                    1,
                                    ALPHA, CMap, UseOpenMP && (1-EnhanceOpenMP)
                     );
                    }
               }
               for (j=0;j<NSubst-NLightSubst;j++)
                   {
                    i = MapSubsts[NLightSubst+j];
                    CalculateWZ((float **) Vars.Bounds, Area, &WXYZ[CarrierPhase], UzSave, (double)KGlobal.Uw[i], UseOpenMP);
                    EulerCalculate((float **) Vars.Bounds,
                                   Vars.Dn[i], KDn[i], SDn[i], &WXYZ[CarrierPhase],
                                   HXR, HYR, HZR,
                                   HX, HY, HZ,
                                   HXX, HYY, HZZ,
                                   HXS, HYS, HZS,
                                   Vars.eDn[i],
                                   D, Dt, 0.0,
                                   1,
                                   ALPHA, CMap, UseOpenMP
                    );
                   }
               if (NSubst-NLightSubst > 0)
                  CalculateWZ((float **) Vars.Bounds, Area, &WXYZ[CarrierPhase], UzSave, 0.0, UseOpenMP);
            }
         }  
      }

      if (CalcBase)
         for (i=0; i<NumEqs; i++)
             {
              Calculate((float **) Vars.Bounds,
                    VDefs[i].Solver,Vars.Name[i],Kfs[i],Sfs[i],&WXYZ[VDefs[i].Phase],
                    HXR,HYR,HZR,
                    HX,HY,HZ,HXX,HYY,HZZ,HXS,HYS,HZS,H1,L,M,G,
                    *VDefs[i]._NuMol,Dt,*VDefs[i]._Kappa,
                    DescH,DescP,DescV,Maps[i].Map,Area,Boundaries,VDefs[i].Restrict,
                    VDefs[i].Flags,
                    HXreg,HYreg,HZreg,
                    VDefs[i].Projection,
                    UseOpenMP, ReverseDirs[i]);
              ReverseDirs[i] = 1-ReverseDirs[i];
             }

      if (NSubst>0)
         {
          int EnhanceOpenMP = UseOpenMP && UseEnhancedOpenMP && NLightSubst>0 && (NLightSubst%nSMP || NLightSubst>3*nSMP);
          #pragma omp parallel if (EnhanceOpenMP)
          {
           #ifdef __OPENMP__
             int      NArray   = ScKfSize/sizeof(double);
             int      ThreadID = omp_get_thread_num();
             double * _L  = &L[ThreadID*NArray];
             double * _G  = &G[ThreadID*NArray];
             double * _M  = &M[ThreadID*NArray];
             float  * _H1 = &H1[ThreadID*BoardSize/sizeof(float)];
           #else
             double * _L  = L;
             double * _G  = G;
             double * _M  = M;
             float  * _H1 = H1;
           #endif
           #pragma omp for schedule(dynamic) private(i,j)
           for (j=0;j<NLightSubst;j++)
               {
                i = MapSubsts[j];
                Calculate((float **) Vars.Bounds,
                          NULL,Vars.Dn[i],KDn[i],SDn[i],&WXYZ[CarrierPhase],
                          HXR,HYR,HZR,
                          HX,HY,HZ,HXX,HYY,HZZ,HXS,HYS,HZS,_H1,_L,_M,_G,D,Dt,0.0,DescH,DescP,DescV,CMap,Area,Boundaries,rsPositive,
                          1,HXreg,HYreg,HZreg,
                          VDefs[i].Projection,
                          UseOpenMP && (1-EnhanceOpenMP),
                          ReverseDirs[NumEqs+i]);
                ReverseDirs[NumEqs+i] = 1-ReverseDirs[NumEqs+i];
               }
          }
          for (j=0;j<NSubst-NLightSubst;j++)
              {
               i = MapSubsts[NLightSubst+j];
               if (SlowMode) CalculateWZ((float **) Vars.Bounds,Area,&WXYZ[CarrierPhase],UzSave,(double)KGlobal.Uw[i], UseOpenMP);
               Calculate((float **) Vars.Bounds,
                         NULL,Vars.Dn[i],KDn[i],SDn[i],&WXYZ[CarrierPhase],
                         HXR,HYR,HZR,
                         HX,HY,HZ,HXX,HYY,HZZ,HXS,HYS,HZS,H1,L,M,G,D,Dt,0.0,DescH,DescP,DescV,CMap,Area,Boundaries,rsPositive,
                         1,HXreg,HYreg,HZreg,
                         VDefs[i].Projection,
                         UseOpenMP, ReverseDirs[NumEqs+i]);
               ReverseDirs[NumEqs+i] = 1-ReverseDirs[NumEqs+i];
              }
         }
     }
 if (FastTAUDivider>1) TAU *= FastTAUDivider;

 if (CheckTau && CheckNew) {
    if (CalcBase) {
       for (i=0; i<NumEqs; i++)
           if (VDefs[i].Solver == NULL)
              for (z = 0, Ptr = 0; z < NZ; z++)
                  for (y = 0; y < NY; y++)
                      for (x = 0; x < NX; x++, Ptr++)
                          if (Maps[i].Map[Ptr] == Zero && fabs(ALPHA[Ptr] - 1.0) > 1E-5)
                             Vars.Name[i][Ptr] = ALPHA[Ptr]*Vars.Name[i][Ptr] + (1.0-ALPHA[Ptr])*Vars.eName[i][Ptr];
       if (NSubst > 0)
           for (i = 0; i < NSubst; i++)
               for (z = 0, Ptr = 0; z < NZ; z++)
                   for (y = 0; y < NY; y++)
                       for (x = 0; x < NX; x++, Ptr++)
                           if (CMap[Ptr] == Zero && fabs(ALPHA[Ptr] - 1.0) > 1E-5)
                              Vars.Dn[i][Ptr] = ALPHA[Ptr]*Vars.Dn[i][Ptr] + (1.0-ALPHA[Ptr])*Vars.eDn[i][Ptr];
    }  
 }
#endif

#ifdef __PARALLEL__
 if (CheckCommand==CHK_EndPhase2)
    (*MaxEps) = 0.0;
 /* Принимаем ВСЕ обработанные частичные поля: (Ux,Uy,Uz), P, T, Nu, Dn */
 if (SendData)
    for (i=1; i<nProcs; i++)
        {
#ifdef __MPI__
          MPI_Probe(MPI_ANY_SOURCE, 123456, MPI_COMM_WORLD, &Status);
          slave = ProcIDToSlave(Status.MPI_SOURCE);
/* Changed, because in PowerXPlorer version we have no different links between same processors
#elif defined(__ROUTER__) || defined(__ROUTER_100__)
          slave = ProcIDToSlave(i);
#else
          slave = SelectList(nProcs-1,Options);
#endif
*/
#else
          slave = ProcIDToSlave(i);
#endif
          if (slave==0 || slave==nProcs-2)
             PartL = (LengthGrid[slave]-1)*NY*NX;
          else
             PartL = (LengthGrid[slave]-2)*NY*NX;
          RecvSlaveXXXTag(slave,(byte *) HBuff, (long) (NSubst+NumEqs)*PartL*sizeof(float));

          Offset = 0;
          for (j=0 ; j<NumEqs; j++)
              MoveFromBuf(&Offset,slave,Vars.Name[j]);
          for (j=0; j<NSubst; j++)
              MoveFromBuf(&Offset,slave,Vars.Dn[j]);
        }
 RecvBoundSize = 0;
 if (SendData) {
    static int BoundsAllocated = 0;
    for (i=1; i<nProcs; i++)
      {
       char * Buf = (char *) HBuff;
       int ToRecv = 0;

       RecvSlave(i-1,(byte *) &ToRecv,sizeof(int));
       memmove(&Buf[RecvBoundSize],&ToRecv,sizeof(int));
       if (!BoundsAllocated)
          Vars.BoundSize += ToRecv;
       RecvBoundSize+=sizeof(int);
       if (ToRecv)
          for (j=_back; j<=_top; j++)
              {
               RecvSlave(i-1,(byte *) &Buf[RecvBoundSize],ToRecv);
               RecvBoundSize+=ToRecv;
              }
      }
    if (!BoundsAllocated) {
       for (j=_back; j<=_top; j++)
           Vars.Bounds[j] = (float *) SafeMalloc(Vars.BoundSize);
       SpecialBoundsNum = Vars.BoundSize/sizeof(float);
       OldRefAddrs = (float **) SafeMalloc(SpecialBoundsNum*sizeof(float *));
       OldRefVals  = (float *) SafeMalloc(SpecialBoundsNum*sizeof(float));
       BoundsAllocated = 1;
    }
    OldRefNum = 0;
    /* Move from HBuf to Bounds */
    Ptr = 0; /* Read pointer */
    x = 0;
    z = 0;
    for (i = 1; i < nProcs; i++)
      {
       char * Buf = (char *) HBuff;
       int n = 0;

       memmove(&n, &Buf[Ptr], sizeof(int));
       Ptr += sizeof(int);
       y = n/sizeof(float);
       if (n) {
          for (j = _back; j <= _top; j++)
              {
               memmove(&Vars.Bounds[j][x], &Buf[Ptr], n);
               Ptr += n;
              }
          /* Now substitute a refs to local boundaries onto the refs to global boundaries */
          slave = i - 1;
          if (slave == 0 || slave == nProcs-2)
             PartL = (LengthGrid[slave]-1)*NY*NX;
          else
             PartL = (LengthGrid[slave]-2)*NY*NX;

          for (j = 0 ; j < NumEqs; j++) {
              int _Ptr = z*NY*NX;

              float * Addr = &Vars.Name[j][_Ptr];
              for (k = 0; k < PartL; k++, Addr++) {
                  float Val  = *Addr;
                  if (TranslateBoundaries(Addr, x)) {
                     OldRefAddrs[OldRefNum] = Addr;
                     OldRefVals[OldRefNum++] = Val;
                  }
              }
          }

          for (j = 0 ; j < NSubst; j++) {
              int _Ptr = z*NY*NX;

              float * Addr = &Vars.Dn[j][_Ptr];
              for (k = 0; k < PartL; k++, Addr++) {
                  float Val  = *Addr;
                  if (TranslateBoundaries(Addr, x)) {
                     OldRefAddrs[OldRefNum] = Addr;
                     OldRefVals[OldRefNum++] = Val;
                  }
              }
          }

          z += PartL/(NY*NX);
       }
       x += y;
      }
 }
 if (CheckCommand==CHK_EndPhase2)
    for (i=1; i<nProcs; i++)
        {
         RecvSlave(i-1,(byte *) &LocalEps, sizeof(LocalEps));
         if (LocalEps>(*MaxEps))
            (*MaxEps) = LocalEps;
        }
#endif

#ifndef __PARALLEL__
 /* Пересчет частиц-индикаторов */
 if (CheckTau && CheckNew) {
    #pragma omp parallel if (UseOpenMP && UseEnhancedOpenMP)
    #pragma omp for schedule(guided) private(i,x,y,z,Ptr)
    for (i = 0; i < NParts; i++) {
        double h = 1.0;
        float X = Vars.Parts[0][i];
        float Y = Vars.Parts[1][i];
        float Z = Vars.Parts[2][i];
        x = (int) (X/h);
        y = (int) (Y/h);
        z = (int) (Z/h);
        Ptr = z*NY*NX+y*NX+x;
        if (x > 0 && x < NXs && y > 0 && y < NYs && z > 0 && z < NZs) {
           double VX = Vars.PartsV[0][i];
           double VY = Vars.PartsV[0][i];
           double VZ = Vars.PartsV[0][i];
           double V = sqrt(VX*VX+VY*VY+VZ*VZ);
           double PX = k2*(GRAV[ZYXP]-GRAV[ZYXM])/2.0/h - k4*VX*V;
           double PY = k2*(GRAV[ZYPX]-GRAV[ZYMX])/2.0/h - k4*VY*V;
           double PZ = k2*(GRAV[ZPYX]-GRAV[ZMYX])/2.0/h - k4*VZ*V;
           float NewX, NewY, NewZ;
           Vars.PartsV[0][i] += TAU*PX;
           Vars.PartsV[1][i] += TAU*PY;
           Vars.PartsV[2][i] += TAU*PZ;
           NewX = X + TAU*Vars.PartsV[0][i];
           NewY = Y + TAU*Vars.PartsV[1][i];
           NewZ = Z + TAU*Vars.PartsV[2][i];
           if (NewX >= 1.0 && NewX < NXs && NewY >= 1.0 && NewY < NYs &&
               NewZ >= 1.0 && NewZ < NZs
              ) {
              Vars.Parts[0][i] = NewX;
              Vars.Parts[1][i] = NewY;
              Vars.Parts[2][i] = NewZ;
           }
        }
    }
 }

 if (NReact)
    {
     int EnhanceOpenMP = UseOpenMP && UseEnhancedOpenMP && NASubst>0 && (NASubst%nSMP || NASubst>3*nSMP);
     static int _Iters = 0;
     int init = _Iters < InitOpenMPIters;
     TraceWSorter WS[markZ99];
     TIME_STRUCT start, end;

     KGlobal.InitH = min(KGlobal.InitH,TAU*0.001f);
     for (i=0; i<KinetChunkSize; i++)
         Contexts[i]._EndTime = (float)TAU;

     FTIME(&start);
     qsort(Traces, NKinets, sizeof(Traces[0]), TraceCmp);

     if (KGlobal.KinetTraceWNotFound) {
        for (i = 0; i < markZ99; i++) {
            WS[i].Order = i;
            WS[i].W = KGlobal.KinetTraceW[i];
        }
        qsort(WS, markZ99, sizeof(WS[0]), TraceWCmp);
        for (i = 0; i < markZ99; i++)
            KGlobal.KinetTraceWOrder[i] = WS[i].Order;
        KGlobal.KinetTraceWNotFound--;
     }

     if (PhaseVars[CarrierPhase]._T >= 0)
	 processKinetics(init ? 2 : 1 /* !!! */, 0, 0x7FFFFFFF, NKinets, KinetChunkSize, Traces, Contexts, &KGlobal,
		Vars.Name[PhaseVars[CarrierPhase]._T], Vars.Dn, UseGear, Adams_Rozhkov_Method, NULL);
     else
	 processKinetics(init ? 2 : 1 /* !!! */, 0, 0x7FFFFFFF, NKinets, KinetChunkSize, Traces, Contexts, &KGlobal,
		NULL, Vars.Dn, UseGear, Adams_Rozhkov_Method, NULL);
     _Iters++; /* !!! */
     FTIME(&end);
     KinetTime += DIFFTIME(start, end) - last_compile_time(NULL);

#ifdef __DEBUG__
     for (i = 0; i < markZ99; i++)
         printf("%f ", KGlobal.KinetTraceW[KGlobal.KinetTraceWOrder[i]]);
     printf("\n");
#endif

     #pragma omp parallel if(EnhanceOpenMP)
     #pragma omp for schedule(dynamic) private(i)
     for (i=0;i<NASubst;i++)
         {
          CorrectBounds(Vars.Bounds,Vars.Dn[KGlobal.TranMap[i]], NZ,NY,NX, NY*NX,NX,1, DescH,CMap,Area, 0, 0, UseOpenMP && (1-EnhanceOpenMP));
          CorrectBounds(Vars.Bounds,Vars.Dn[KGlobal.TranMap[i]], NZ,NX,NY, NY*NX,1,NX, DescP,CMap,Area, 0, 1, UseOpenMP && (1-EnhanceOpenMP));
          CorrectBounds(Vars.Bounds,Vars.Dn[KGlobal.TranMap[i]], NY,NX,NZ, NX,1,NY*NX, DescV,CMap,Area, 0, 2, UseOpenMP && (1-EnhanceOpenMP));
         }
    }
#endif
}

void CalculateEps(int NumZ, double * MaxEps, float * Templ, float * Cur, float * ERRS)
{
 int x,y,z;
 int Ptr = 0;

 for (z=0; z<NumZ; z++)
     for (y=0; y<NY; y++)
         for (x=0; x<NX; Ptr++,x++)
			 if (!IsBound)
             {
              double CurEps = fabs(Templ[Ptr] - Cur[Ptr])/(2.0*2.0-1.0);

			  if (CurEps > (*MaxEps))
                 (*MaxEps) = CurEps;
              if (CurEps > ERRS[Ptr])
                 ERRS[Ptr] = CurEps;
             }
}

#ifdef __PARALLEL__

#define BegOffs(i)   ((i)*(NZ+2)*NX*NY)
#define BegBoard(i)  (((i)*(NZ+2)+1)*NX*NY)
#define OffsBufDn(i) (((NZ+2)*(NumEqs+(i))+1)*NX*NY)

TrVect * UpTrace   = NULL;
TrVect * DownTrace = NULL;

typedef struct {
  unsigned       int Index;
  unsigned short int Iters;
} SortInfo;

typedef struct {
  unsigned int  Index;
  long          Load;
} _ProcessorLoad;

#define MinLong   50
#define DeltaLong 50

void Process_XY(float ** Bounds,
                int VarNum,
                float * K, float * _S, double R, float Kp,
                WKoeffs * W,
                unsigned char * HMap, unsigned char * HArea, unsigned char * Boundaries,
                char Restrict,
                char Flags,
                char RegularX, char RegularY,
                char Projection,
                int  UseOpenMP,
                char Reverse)
{
 #ifdef __OPENMP__
  int      NArray   = ScKfSize/sizeof(double);
  int      ThreadID = omp_get_thread_num();
  double * _L  = &L[ThreadID*NArray];
  double * _G  = &G[ThreadID*NArray];
  double * _M  = &M[ThreadID*NArray];
  float  * _H1 = &H1[ThreadID*BoardSize/sizeof(float)];
 #else
  double * _L  = L;
  double * _G  = G;
  double * _M  = M;
  float  * _H1 = H1;
 #endif
 CalculateInXY(Bounds,&HBuff[BegBoard(VarNum)],K,_S,W,
               HXR,HYR,
               HX,HY,HXX,HYY,HXS,HYS,_H1,_L,_M,_G,
               R,Dt,Kp,DescH,DescP,DescV,HMap,HArea,Boundaries,Restrict,
               Flags,RegularX,RegularY,Projection, UseOpenMP, Reverse);
}

void Process_Z(float ** Bounds,
               _Solver S, int VarNum,
               float * K, float * _S, double R, float Kp,
               WKoeffs * W,
               unsigned char * HMap, unsigned char * HArea, unsigned char * Boundaries,
               char Restrict,
               char Flags,
               char RegularZ,
               char Projection,
               int  UseOpenMP)
{
 #ifdef __OPENMP__
  int      NArray   = ScKfSize/sizeof(double);
  int      ThreadID = omp_get_thread_num();
  double * _L  = &L[ThreadID*NArray];
  double * _G  = &G[ThreadID*NArray];
  double * _M  = &M[ThreadID*NArray];
  float  * _H1 = &H1[ThreadID*BoardSize/sizeof(float)];
 #else
  double * _L  = L;
  double * _G  = G;
  double * _M  = M;
  float  * _H1 = H1;
 #endif
 CalculateInZ(Bounds,
              S,&HBuff[BegBoard(VarNum)-NX*NY],&HBuff[BegBoard(VarNum)],&HBuff[BegBoard(VarNum)+NZ*NX*NY],
              (DownTrace==NULL ? NULL : &DownTrace[VarNum*NX*NY]),(UpTrace==NULL ? NULL : &UpTrace[VarNum*NX*NY]),
              K,_S,W,
              HXR,HYR,HZR,
              HX,HXX,HXS,HY,HYY,HYS,HZ,HZZ,HZS,_H1,_L,_M,_G,R,Dt,Kp,DescH,DescP,DescV,HMap,HArea,Boundaries,Restrict,
              Flags,RegularZ,Projection, UseOpenMP);
}

#undef _A
#undef _B

int SortInfoCmp(const void * _A, const void * _B)
{
 return (signed int) ((SortInfo *) _B)->Iters - (signed int) ((SortInfo *) _A)->Iters;
}

int ProcessorLoadCmp(const void * _A, const void * _B)
{
 return ((_ProcessorLoad *) _B)->Load - ((_ProcessorLoad *) _A)->Load;
}

int TransactionCmp(const void * _A, const void * _B)
{
 return (signed int) ((_ProcessorLoad *) _A)->Index - (signed int) ((_ProcessorLoad *) _B)->Index;
}

void PackForUpExchange(float * Buf, int GH)
{
 int i;
 float * ActualBuf = Buf;

#if !defined(__ROUTER_100__)
 if (SHMEMUp) ActualBuf = (float *) SHMEM_AT(GH,(NumEqs+NSubst)*2*NX*NY*sizeof(float));
#endif
 for (i=0; i<NumEqs+NSubst; i++)
     memmove(&ActualBuf[i*2*NX*NY], &HBuff[BegBoard(i)+NZs2*NY*NX], 2*NX*NY*sizeof(float));
#if !defined(__ROUTER_100__)
 if (SHMEMUp)
    {
     SHMEM_DT(ActualBuf,GH);
     memmove(Buf,&GH,sizeof(GH));
    }
#endif
}

void UnPackAfterUpExchange(float * Buf)
{
 int i;
 float * ActualBuf = Buf;

#if !defined(__ROUTER_100__)
 int * GHbuf = (int *) Buf;
 if (SHMEMUp) ActualBuf = (float *) SHMEM_AT(*GHbuf,(NumEqs+NSubst)*2*NX*NY*sizeof(float));
#endif
 for (i=0; i<NumEqs+NSubst; i++)
     {
      memmove(&HBuff[BegBoard(i)+NZ*NY*NX], &ActualBuf[i*2*NX*NY], NX*NY*sizeof(float));
      memmove(&HBuff[BegBoard(i)+NZs*NY*NX], &ActualBuf[(i*2+1)*NX*NY], NX*NY*sizeof(float));
     }
#if !defined(__ROUTER_100__)
 if (SHMEMUp) SHMEM_DT(ActualBuf,*GHbuf);
#endif
}

void PackForDownExchange(float * Buf, int GH)
{
 int i;
 float * ActualBuf = Buf;

#if !defined(__ROUTER_100__)
 if (SHMEMDown) ActualBuf = (float *) SHMEM_AT(GH,(NumEqs+NSubst)*2*NX*NY*sizeof(float));
#endif
 for (i=0; i<NumEqs+NSubst; i++)
     {
      memmove(&ActualBuf[i*2*NX*NY], &HBuff[BegBoard(i)+2*NY*NX], NX*NY*sizeof(float));
      memmove(&ActualBuf[(i*2+1)*NX*NY], &HBuff[BegBoard(i)+1*NY*NX], NX*NY*sizeof(float));
     }
#if !defined(__ROUTER_100__)
 if (SHMEMDown)
    {
     SHMEM_DT(ActualBuf,GH);
     memmove(Buf,&GH,sizeof(GH));
    }
#endif
}

void UnPackAfterDownExchange(float * Buf)
{
 int i;
 float * ActualBuf = Buf;

#if !defined(__ROUTER_100__)
 int * GHbuf = (int *) Buf;
 if (SHMEMDown) ActualBuf = (float *) SHMEM_AT(*GHbuf,(NumEqs+NSubst)*2*NX*NY*sizeof(float));
#endif
 for (i=0; i<NumEqs+NSubst; i++)
     memmove(&HBuff[BegOffs(i)], &ActualBuf[i*2*NX*NY], 2*NX*NY*sizeof(float));
#if !defined(__ROUTER_100__)
 if (SHMEMDown) SHMEM_DT(ActualBuf,*GHbuf);
#endif
}

void PackGRAVForUpExchange(float * Buf, int GH)
{
 int i;
 float * ActualBuf = Buf;

#if !defined(__ROUTER_100__)
 if (SHMEMUp) ActualBuf = (float *) SHMEM_AT(GH, NX*NY*sizeof(float));
#endif
 memmove(ActualBuf, &GRAV[NZs1*NY*NX], NX*NY*sizeof(float));
#if !defined(__ROUTER_100__)
 if (SHMEMUp)
    {
     SHMEM_DT(ActualBuf,GH);
     memmove(Buf,&GH,sizeof(GH));
    }
#endif
}

void UnPackGRAVAfterUpExchange(float * Buf)
{
 int i;
 float * ActualBuf = Buf;

#if !defined(__ROUTER_100__)
 int * GHbuf = (int *) Buf;
 if (SHMEMUp) ActualBuf = (float *) SHMEM_AT(*GHbuf, NX*NY*sizeof(float));
#endif
 memmove(&GRAV[NZs*NY*NX], ActualBuf, NX*NY*sizeof(float));
#if !defined(__ROUTER_100__)
 if (SHMEMUp) SHMEM_DT(ActualBuf,*GHbuf);
#endif
}

void PackGRAVForDownExchange(float * Buf, int GH)
{
 int i;
 float * ActualBuf = Buf;

#if !defined(__ROUTER_100__)
 if (SHMEMDown) ActualBuf = (float *) SHMEM_AT(GH, NX*NY*sizeof(float));
#endif
 memmove(ActualBuf, &GRAV[NY*NX], NX*NY*sizeof(float));
#if !defined(__ROUTER_100__)
 if (SHMEMDown)
    {
     SHMEM_DT(ActualBuf,GH);
     memmove(Buf,&GH,sizeof(GH));
    }
#endif
}

void UnPackGRAVAfterDownExchange(float * Buf)
{
 int i;
 float * ActualBuf = Buf;

#if !defined(__ROUTER_100__)
 int * GHbuf = (int *) Buf;
 if (SHMEMDown) ActualBuf = (float *) SHMEM_AT(*GHbuf, NX*NY*sizeof(float));
#endif
 memmove(GRAV, ActualBuf, NX*NY*sizeof(float));
#if !defined(__ROUTER_100__)
 if (SHMEMDown) SHMEM_DT(ActualBuf,*GHbuf);
#endif
}

void PackParts(float * ActualBuf, int up) {
 int maxNP = (NX*NY-1)/6;
 int Sended = 0;
 int StorePtr = 0;
 int BufPtr = 1;
 int i, j;

 for (i = 0; i < NParts; i++) {
     int bypass = 0;
     if (up && Vars.Parts[1][i] < 0 || !up && Vars.Parts[0][i] < 0) {
        if (up)
           Vars.Parts[1][i] = -Vars.Parts[1][i];
        else
           Vars.Parts[0][i] = -Vars.Parts[0][i];
        if (Sended < maxNP) {
           for (j = 0; j < 3; j++)
               ActualBuf[BufPtr++] = Vars.Parts[j][i];
           for (j = 0; j < 3; j++)
               ActualBuf[BufPtr++] = Vars.PartsV[j][i];
           Sended++;
           bypass = 1;
        }
     }
     if (StorePtr != i && !bypass) {
        for (j = 0; j < 3; j++)
            Vars.Parts[j][StorePtr] = Vars.Parts[j][i];
        for (j = 0; j < 3; j++)
            Vars.PartsV[j][StorePtr] = Vars.PartsV[j][i];
        StorePtr++;
     }
 }
 ActualBuf[0] = Sended;
 NParts -= Sended;
}

void UnPackParts(float * ActualBuf) {
 int Received = ActualBuf[0];
 int BufPtr = 1;
 int i, j;
 for (i = 0; i < Received; i++) {
     for (j = 0; j < 3; j++)
         Vars.Parts[j][NParts] = ActualBuf[BufPtr++];
     for (j = 0; j < 3; j++)
         Vars.PartsV[j][NParts] = ActualBuf[BufPtr++];
     NParts++;
 }
}

void PackPARTSForUpExchange(float * Buf, int GH)
{
 float * ActualBuf = Buf;

#if !defined(__ROUTER_100__)
 if (SHMEMUp) ActualBuf = (float *) SHMEM_AT(GH, NX*NY*sizeof(float));
#endif
 // memmove(ActualBuf, &DATA[UP], NX*NY*sizeof(float));
 PackParts(ActualBuf, 1);
#if !defined(__ROUTER_100__)
 if (SHMEMUp)
    {
     SHMEM_DT(ActualBuf,GH);
     memmove(Buf,&GH,sizeof(GH));
    }
#endif
}

void UnPackPARTSAfterUpExchange(float * Buf)
{
 int i, j;
 float * ActualBuf = Buf;

#if !defined(__ROUTER_100__)
 int * GHbuf = (int *) Buf;
 if (SHMEMUp) ActualBuf = (float *) SHMEM_AT(*GHbuf, NX*NY*sizeof(float));
#endif
 // memmove(&DATA[UP], ActualBuf, NX*NY*sizeof(float));
 UnPackParts(ActualBuf);
#if !defined(__ROUTER_100__)
 if (SHMEMUp) SHMEM_DT(ActualBuf,*GHbuf);
#endif
}

void PackPARTSForDownExchange(float * Buf, int GH)
{
 float * ActualBuf = Buf;

#if !defined(__ROUTER_100__)
 if (SHMEMDown) ActualBuf = (float *) SHMEM_AT(GH, NX*NY*sizeof(float));
#endif
 // memmove(ActualBuf, &DATA[DOWN], NX*NY*sizeof(float));
 PackParts(ActualBuf, 0);
#if !defined(__ROUTER_100__)
 if (SHMEMDown)
    {
     SHMEM_DT(ActualBuf,GH);
     memmove(Buf,&GH,sizeof(GH));
    }
#endif
}

void UnPackPARTSAfterDownExchange(float * Buf)
{
 int i, j;
 float * ActualBuf = Buf;

#if !defined(__ROUTER_100__)
 int * GHbuf = (int *) Buf;
 if (SHMEMDown) ActualBuf = (float *) SHMEM_AT(*GHbuf, NX*NY*sizeof(float));
#endif
 // memmove(DATA[DOWN], ActualBuf, NX*NY*sizeof(float));
 UnPackParts(ActualBuf);
#if !defined(__ROUTER_100__)
 if (SHMEMDown) SHMEM_DT(ActualBuf,*GHbuf);
#endif
}

void SendNeighbour(int NeighbourID, byte * Buf, long Length)
{
#if defined(__ROUTER_100__)
 r_write(SlaveToProcID(NeighbourID), Buf, Length);
 w_write(SlaveToProcID(NeighbourID));
#elif defined(__ROUTER__)
 RouterSend(SlaveToProcID(NeighbourID), Buf, Length);
#else
 MPI_Send(Buf, Length, MPI_BYTE, NeighbourID, 1235, SlavesComm);
#endif
}

void RecvNeighbour(int NeighbourID, byte * Buf, long Length)
{
#if defined(__ROUTER_100__)
 r_read(SlaveToProcID(NeighbourID), Buf, Length);
 w_read(SlaveToProcID(NeighbourID));
#elif defined(__ROUTER__)
 RouterIRecv(SlaveToProcID(NeighbourID), Buf, Length);
 w_read(started_handle);
#else
 MPI_Status status;

 MPI_Recv(Buf, Length, MPI_BYTE, NeighbourID, 1235, SlavesComm, &status);
#endif
}

#ifdef __MPI__
MPI_Request ExchRequests[4];
MPI_Request ExchGRAVRequests[4];
#endif

typedef void (* _packer)(float * Buf, int GH);
typedef void (* _unpacker)(float * Buf);

/* Обмен данными по кольцу процессоров */
void MakeExchange(const char * MSG1, const char * MSG2,
#ifdef __MPI__
                  MPI_Request * REQUESTS,
#endif
                  float * HU, float * HD, float * HUr, float * HDr,
                  _packer UPACK, _unpacker UUNPACK,
                  _packer DPACK, _unpacker DUNPACK,
                  int Usize, int Dsize,
                  int GHU, int GHD)
{
#ifdef __MPI__
 MPI_Status Statuses[4];
#elif defined(__ROUTER__)
 int US,UR,DS,DR;
#endif

 DebugPrintf(DEBUG_FILE, MSG1);
#if defined(__MPI__) || defined(__ROUTER__)
 if (UpExchange)   UPACK(HU,GHU);
 if (DownExchange) DPACK(HD,GHD);

 if (UpExchange && DownExchange)
    {
#ifdef __MPI__
     MPI_Startall(4,&REQUESTS[0]);
     MPI_Waitall(4,&REQUESTS[0],Statuses);
#else
     r_write(MyProcID+1, HU, Usize); US = started_handle;
     r_read(MyProcID+1, HUr, Usize); UR = started_handle;
     r_write(MyProcID-1, HD, Dsize); DS = started_handle;
     r_read(MyProcID-1, HDr, Dsize); DR = started_handle;
     w_write(DS); w_read(DR);
     w_read(UR); w_write(US);
#endif
    }
 else
    {
     if (UpExchange)
        {
#ifdef __MPI__
         MPI_Startall(2,&REQUESTS[0]);
         MPI_Waitall(2,&REQUESTS[0],Statuses);
#else
         r_write(MyProcID+1, HU, Usize); US = started_handle;
         r_read(MyProcID+1, HUr, Usize); UR = started_handle;
         w_read(UR); w_write(US);
#endif
        }
     if (DownExchange)
        {
#ifdef __MPI__
         MPI_Startall(2,&REQUESTS[2]);
         MPI_Waitall(2,&REQUESTS[2],Statuses);
#else
         r_write(MyProcID-1, HD, Dsize); DS = started_handle;
         r_read(MyProcID-1, HDr, Dsize); DR = started_handle;
         w_write(DS); w_read(DR);
#endif
        }
    }

 if (UpExchange)   UUNPACK(HUr);
 if (DownExchange) DUNPACK(HDr);
#else
 if (MyProcID%2)
    {
     if (UpExchange)
        {
         UPACK(HU,GHU);
#ifdef __ROUTER_100__
         r_read(MyProcID+1, HUr, Usize);
         r_write(MyProcID+1, HU, Usize);
         w_read(MyProcID+1); w_write(MyProcID+1);
#else
         SendNextRing(NextLink,(byte *) HU, Usize);
         RecvNextRing(NextLink,(byte *) HUr, Usize);
#endif
         UUNPACK(HUr);
        }
     if (DownExchange)
        {
         DPACK(HD,GHD);
#ifdef __ROUTER_100__
         r_read(MyProcID-1, HDr, Dsize);
         r_write(MyProcID-1, HD, Dsize);
         w_write(MyProcID-1); w_read(MyProcID-1);
#else
         SendPrevRing(PrevLink,(byte *) HD, Dsize);
         RecvPrevRing(PrevLink,(byte *) HDr, Dsize);
#endif
         DUNPACK(HDr);
        }
    }
 else
    {
     if (DownExchange)
        {
         DPACK(HD,GHD);
#ifdef __ROUTER_100__
         r_read(MyProcID-1, HDr, Dsize);
         r_write(MyProcID-1, HD, Dsize);
         w_write(MyProcID-1); w_read(MyProcID-1);
#else
         RecvPrevRing(PrevLink,(byte *) HDr, Dsize);
         SendPrevRing(PrevLink,(byte *) HD, Dsize);
#endif
         DUNPACK(HDr);
        }
     if (UpExchange)
        {
         UPACK(HU,GHU);
#ifdef __ROUTER_100__
         r_read(MyProcID+1, HUr, Usize);
         r_write(MyProcID+1, HU, Usize);
         w_read(MyProcID+1); w_write(MyProcID+1);
#else
         RecvNextRing(NextLink,(byte *) HUr, Usize);
         SendNextRing(NextLink,(byte *) HU, Usize);
#endif
         UUNPACK(HUr);
        }
    }
#endif
 DebugPrintf(DEBUG_FILE, MSG2);
}

SortInfo *  SortedIndex  = NULL;
SortInfo *  _SortedIndex = NULL;
int         ZeroNodes;

#ifdef __MPI__
typedef enum {psHasNodes,psRequesting,psStopped} psModes;

int IndexCmp(const void * _A, const void * _B)
{
 return *((signed int *) _A) - *((signed int *) _B);
}
#endif

int InitSpecialBounds(unsigned char * Area, int V1, int V2, int InitPtr)
{
 int i,j,z,y,x,Ptr;

 if (UseSpecialBounds)
    for (i=V1, j=InitPtr; i<=V2; i++)
      for (z=3, Ptr=3*NY*NX; z<NZ-3; z++)
        for (y=0; y<NY; y++)
          for (x=0; x<NX; x++,Ptr++)
            if (Area[Ptr] && !IsExchng && ((IsBack!=0)+(IsForw!=0)+(IsRight!=0)+(IsLeft!=0)+(IsBottom!=0)+(IsTop!=0))>1)
               {
                int CODE = j+1; // 128 Eqs * 65536 points; Min = 1; Max = 8388480+1
                int FIRST = CODE/1000000; // 0..9
                int REM = CODE%1000000; // Min = 0; Max = 999999
                // 500000..999999 => 0.500000..0.999999
                // 1..499999 => -[0.500001..0.999999]
                float Mantissa = REM>=500000 ? 0.000001*REM : -0.000001*(REM+500000);
                int   POW = 100+FIRST; // Max = 109
                // ABS[max] = max(Mantissa)*2^max(POW) = 6.4903645828*1E32
                // ABS[min] = ABS[min(Mantissa)*2^min(POW)] = 6.3382530011*1E29
                // Critical Number = +-1E25, for example...
                float Number = ldexp(Mantissa,POW); // rounding must present!!!!!
                int k;

                for (k=_back; k<=_top; k++)
                    Vars.Bounds[k][j] = HBuff[BegBoard(i)+Ptr];
                HBuff[BegBoard(i)+Ptr] = Number;
                j++;
               }
 return j;
}

unsigned int savedIters;

#pragma memoization(t,g,o) lin_extrapolator(150, 2) controlled(50, 1)
void prognose_load(int counter, unsigned int Index, double Iters[1]) {
   *Iters = savedIters;
}

void Slave()
{
 HandlePacket Packet;

 float * HU  = (float *) SafeMalloc(2*(NumEqs+NSubst)*NX*NY*sizeof(float));
 float * HD  = (float *) SafeMalloc(2*(NumEqs+NSubst)*NX*NY*sizeof(float));
 float * HUr = (float *) SafeMalloc(2*(NumEqs+NSubst)*NX*NY*sizeof(float));
 float * HDr = (float *) SafeMalloc(2*(NumEqs+NSubst)*NX*NY*sizeof(float));

 float * HUG  = CheckTau ? (float *) SafeMalloc(NX*NY*sizeof(float)) : NULL;
 float * HDG  = CheckTau ? (float *) SafeMalloc(NX*NY*sizeof(float)) : NULL;
 float * HUrG = CheckTau ? (float *) SafeMalloc(NX*NY*sizeof(float)) : NULL;
 float * HDrG = CheckTau ? (float *) SafeMalloc(NX*NY*sizeof(float)) : NULL;

 int     GHU  = 0;
 int     GHD  = 0;
 
 int     GHUG = 0;
 int     GHDG = 0;
 
 int     Usize = 0;
 int     Dsize = 0;

 int     UsizeG = 0;
 int     DsizeG = 0;

 WBVals * CheckTrace = NULL;
 int    * Ptrs = NULL;

 float        ** Buffers     = NULL;
 unsigned int ** ArrangeIndex     = NULL;
 unsigned int  * MaxLongs         = NULL;
 unsigned int  * Longs            = NULL;
 long          * WorkLoad         = NULL;
 _ProcessorLoad * Senders          = NULL;
 int             NumSenders;
 int             SenderCount;
 _ProcessorLoad * Receivers        = NULL;
 int             NumReceivers;
 int             ReceiverCount;
 _ProcessorLoad * Transactions     = NULL;
 int             NumTransactions;

#if defined(__MPI__) || defined(__ROUTER__) || defined(__ROUTER_100__)
 REQUEST   Request;
 REQUEST * GearRequests;
#endif
 long         TotalLoad;
 long         AverageLoad;
 long         RealLoad;
 signed   int DeltaLoad;
 signed   int Disbalance;
 unsigned int CurLoad;
 unsigned int CurLong;
 unsigned int CurProc;

 int GearSenderFlag;
 int GearReceiverFlag;

 long   Ptr, WPtr;
 double MaxEps;
 int    x,y,z;
 int    PartL, FromOffs;
 int    WCount;
 int    Count;
 int    Command = 1;
 int    i,j,k;
 int    PredictCount = 0;
 double Q;
 int    BoundsInitialized = 0;
 int    BoundsPtr = 0;

 static double SaveModelTime = 0.0;
 static int    SaveNParts = 0;

 if (CalculateAll)
    {
     int ToRecv;

     RecvMaster((byte *) &ToRecv,sizeof(int));
     if (ToRecv && ToRecv!=Vars.BoundSize)
        {
         fprintf(StdOutput,"Wow!!! My number of boundary points is not equal to number has received from Master. STOP! \n");
#ifdef __MVS__
         fclose(StdOutput);
#endif
         AbortServer(-80);
        }
     if (ToRecv)
        {
         for (i=_back; i<=_top; i++)
             RecvMaster((byte *) Vars.Bounds[i],ToRecv);
         if (!UseSpecialBounds)
            for (i=0; i<NumEqs+NSubst; i++)
             {
              float * Board = &HBuff[BegBoard(i)];
              float ** Bounds = Vars.Bounds;

              for (Ptr=0; Ptr<BoardSize; Ptr++)
                  if (fabs(Board[Ptr])>=SpecVAL)
                     {
                      int   Count = (IsBack!=0)+(IsForw!=0)+(IsRight!=0)+(IsLeft!=0)+(IsBottom!=0)+(IsTop!=0);
                      float Val = 0.0;

                      if (IsForw)   Val+=_fw(Board[Ptr]);
                      if (IsBack)   Val+=_bw(Board[Ptr]);
                      if (IsBottom) Val+=_bt(Board[Ptr]);
                      if (IsTop)    Val+=_tp(Board[Ptr]);
                      if (IsLeft)   Val+=_lf(Board[Ptr]);
                      if (IsRight)  Val+=_rg(Board[Ptr]);

                      Board[Ptr] = Val/Count;
                     }
             }
        }
     BoundsInitialized = ToRecv!=0;
    }

 if (AllowPrediction && !CheckTau)
    {
     if (UpExchange) UpTrace = (TrVect *) SafeMalloc((NumEqs+NSubst)*NX*NY*sizeof(TrVect));
     if (DownExchange) DownTrace = (TrVect *) SafeMalloc((NumEqs+NSubst)*NX*NY*sizeof(TrVect));

     for (i=0; i<NPred; i++)
         sTAUs[i][0] = 1.0;
     for (i=0, Q = NExp*TAU; i<NPred; i++, Q+=TAU)
       for (j=1; j<NWB; j++)
         sTAUs[i][j] = Q*sTAUs[i][j-1];
     if (AllowQW)
        {
         CheckTrace = (WBVals *) SafeMalloc(ReCalcNPoints*sizeof(WBVals));
         Ptrs       = (int *) SafeMalloc(ReCalcNPoints*sizeof(int));

         for (i=0; i<ReCalcNPoints; i++)
             do {
               x = 1+((double) rand()/RAND_MAX)*NXs2;
               y = 1+((double) rand()/RAND_MAX)*NYs2;
               z = 1+((double) rand()/RAND_MAX)*NZs2;
               Ptrs[i] = z*NY*NX+y*NX+x;
             } while (Area[Ptrs[i]]!=0);

         for (i=0; i<NExp; i++)
             KQ[i] = InitKQ;
         PreparePrediction(KQ);
        }
     else
        PreparePrediction(NULL);
    }

 if (NReact && UseGear)
    {
     ZeroNodes = 0;

#if defined(__MPI__) || defined(__ROUTER__) || defined(__ROUTER_100__)
     GearRequests  = (REQUEST *)       SafeMalloc((TotalProcs-1)*sizeof(REQUEST));
#endif
     Transactions  = (_ProcessorLoad *) SafeMalloc((TotalProcs-1)*sizeof(_ProcessorLoad));
     WorkLoad      = (long *)          SafeMalloc((TotalProcs-1)*sizeof(long));
     Senders       = (_ProcessorLoad *) SafeMalloc((TotalProcs-1)*sizeof(_ProcessorLoad));
     Receivers     = (_ProcessorLoad *) SafeMalloc((TotalProcs-1)*sizeof(_ProcessorLoad));

     Buffers       = (float **)        SafeMalloc((TotalProcs-1)*sizeof(float *));

     MaxLongs      = (unsigned int *)  SafeMalloc((TotalProcs-1)*sizeof(unsigned int));
     Longs         = (unsigned int *)  SafeMalloc((TotalProcs-1)*sizeof(unsigned int));

     if (CalculateAll)
        {
         Ptr = 0;
         for (z=0; z<NZ; z++)
             for (y=0; y<NY; y++)
                 for (x=0; x<NX; Ptr++,x++)
                     if (CMap[Ptr]==Zero)
                        ZeroNodes++;

         ArrangeIndex  = (unsigned int **) SafeMalloc((TotalProcs-1)*sizeof(unsigned int *));
         SortedIndex   = (SortInfo     *)  SafeMalloc(ZeroNodes*sizeof(SortInfo));
         _SortedIndex  = (SortInfo     *)  SafeMalloc(ZeroNodes*sizeof(SortInfo));

         Ptr = 0;
         for (Count=0; Count<ZeroNodes; Count++)
             {
              while (CMap[Ptr]!=Zero) Ptr++;
              SortedIndex[Count].Index = Ptr++;
              SortedIndex[Count].Iters = 1;
             }
        }

     CurLong = max(MinLong,ZeroNodes/(TotalProcs-1));
     for (i=0; i<TotalProcs-1; i++)
         {
          MaxLongs[i] = CurLong;
          Longs[i]    = 0;

          if (ArrangeIndex)
             ArrangeIndex[i] = (unsigned int *) SafeMalloc(MaxLongs[i]*sizeof(unsigned int));
          Buffers[i] = (float *) SafeMalloc(MaxLongs[i]*(NASubst+1)*sizeof(float));
         }
    }

 if (UpExchange) {
     Usize = SHMEMUp ? sizeof(int) : (NumEqs+NSubst)*2*NX*NY*sizeof(float);
     UsizeG = SHMEMUp ? sizeof(int) : NX*NY*sizeof(float);
 }
 if (DownExchange) {
     Dsize = SHMEMDown ? sizeof(int) : (NumEqs+NSubst)*2*NX*NY*sizeof(float);
     DsizeG = SHMEMDown ? sizeof(int) : NX*NY*sizeof(float);
 }

#if !defined(__ROUTER_100__)
 if (SHMEMUp) GHU = SHMEM_GET((NumEqs+NSubst)*2*NX*NY*sizeof(float));
 if (SHMEMDown) GHD = SHMEM_GET((NumEqs+NSubst)*2*NX*NY*sizeof(float));
 if (SHMEMUp) GHUG = CheckTau ? SHMEM_GET(NX*NY*sizeof(float)) : 0;
 if (SHMEMDown) GHDG = CheckTau ? SHMEM_GET(NX*NY*sizeof(float)) : 0;
#endif

#ifdef __MPI__
 if (UpExchange)
    {
     MPI_Send_init(HU, Usize, MPI_BYTE, MyProcID+1, 1000, MPI_COMM_WORLD,
                   &ExchRequests[0]);
     MPI_Recv_init(HUr, Usize, MPI_BYTE,
                   MyProcID+1, 1000, MPI_COMM_WORLD, &ExchRequests[1]);
     if (CheckTau) {
        MPI_Send_init(HUG, UsizeG, MPI_BYTE, MyProcID+1, 1101, MPI_COMM_WORLD,
                      &ExchGRAVRequests[0]);
        MPI_Recv_init(HUrG, UsizeG, MPI_BYTE,
                      MyProcID+1, 1101, MPI_COMM_WORLD, &ExchGRAVRequests[1]);
     }
    }
 if (DownExchange)
    {
     MPI_Send_init(HD, Dsize, MPI_BYTE, MyProcID-1, 1000, MPI_COMM_WORLD,
                   &ExchRequests[2]);
     MPI_Recv_init(HDr, Dsize, MPI_BYTE,
                   MyProcID-1, 1000, MPI_COMM_WORLD, &ExchRequests[3]);
     if (CheckTau) {
        MPI_Send_init(HDG, DsizeG, MPI_BYTE, MyProcID-1, 1101, MPI_COMM_WORLD,
                      &ExchGRAVRequests[2]);
        MPI_Recv_init(HDrG, DsizeG, MPI_BYTE,
                      MyProcID-1, 1101, MPI_COMM_WORLD, &ExchGRAVRequests[3]);
     }
    }
#endif

 if (FirstFlag)     { PartL = NZ-1; FromOffs = 1; }
 else if (LastFlag) { PartL = NZ-1; FromOffs = 2; }
 else               { PartL = NZ-2; FromOffs = 2; }

 do
   {
    RecvPacket(&Packet);
    TAU      = Packet._TAU;
    NReact   = Packet._NReact;
	KGlobal.NReact = NReact;
    takt     = Packet._takt;
    CalcBase = Packet._CalcBase;

    if (CalculateAll)
       {
        DebugPrintf(DEBUG_FILE,"Wait for work data\n");
#if defined(__MPI__) || defined(__ROUTER__) || defined(__ROUTER_100__)
        /* Получаем (Ux,Uy,Uz), P, Nu, T */
        if (Packet.RecvData)
           {
            RecvMaster((byte *) HBuff, (long) NumEqs*(NZ+2)*NX*NY*sizeof(float));
            if (!BoundsInitialized) BoundsPtr = InitSpecialBounds(Area,0,NumEqs-1,BoundsPtr);
           }
        else
           if (Stage<NExp)
	      {
	       MakeExchange("Start exchange\n","Stop exchange\n",
#ifdef __MPI__
                            ExchRequests,
#endif
                            HU,HD,HUr,HDr,
                            PackForUpExchange, UnPackAfterUpExchange,
                            PackForDownExchange, UnPackAfterDownExchange,
                            Usize,Dsize,GHU,GHD);
               #pragma omp parallel if(UseOpenMP && UseEnhancedOpenMP)
               #pragma omp for schedule(dynamic) private(i)
               for (i=0; i<NumEqs+NSubst; i++)
                   if (i<NumEqs)
                      CorrectBounds(Vars.Bounds,&HBuff[BegBoard(i)], NY,NX,NZ, NX,1,NY*NX, DescV, Maps[i].Map,Area,VDefs[i].Projection, 2, UseOpenMP && (1-UseEnhancedOpenMP));
                   else
                      CorrectBounds(Vars.Bounds,&HBuff[BegBoard(i)], NY,NX,NZ, NX,1,NY*NX, DescV, CMap,Area,0, 2, UseOpenMP && (1-UseEnhancedOpenMP));
	      }
        /* Инициируем получение Dn */
        if (Packet.RecvData && NSubst)
           IRecvMaster((byte *) &HBuff[BegOffs(NumEqs)], (long) NSubst*(NZ+2)*NX*NY*sizeof(float), &Request);
#else
        if (Packet.RecvData)
           {
            RecvMaster((byte *) HBuff, (long) (NSubst+NumEqs)*(NZ+2)*NX*NY*sizeof(float));
            if (!BoundsInitialized) BoundsPtr = InitSpecialBounds(Area,0,NSubst+NumEqs-1,BoundsPtr);
            BoundsInitialized = 1;
           }
        else
           if (Stage<NExp)
	      {
	       MakeExchange("Start exchange\n","Stop exchange\n",
#ifdef __MPI__
                            ExchRequests,
#endif
                            HU,HD,HUr,HDr,Usize,Dsize,GHU,GHD);
               #pragma omp parallel if(UseOpenMP && UseEnhancedOpenMP)
               #pragma omp for schedule(dynamic) private(i)
               for (i=0; i<NumEqs+NSubst; i++)
                   if (i<NumEqs)
                      CorrectBounds(Vars.Bounds,&HBuff[BegBoard(i)], NY,NX,NZ, NX,1,NY*NX, DescV, Maps[i].Map,Area,VDefs[i].Projection, 2, UseOpenMP && (1-UseEnhancedOpenMP));
                   else
                      CorrectBounds(Vars.Bounds,&HBuff[BegBoard(i)], NY,NX,NZ, NX,1,NY*NX, DescV, CMap,Area,0, 2, UseOpenMP && (1-UseEnhancedOpenMP));
	      }
#endif
        DebugPrintf(DEBUG_FILE,"Work data received\n");
        
        switch (Packet.CheckCommand) {
          case CHK_Phase1 :
               for (i=0; i<NumEqs; i++)
                   memmove(Old.Name[i],&HBuff[BegOffs(i)],(NZ+2)*NX*NY*sizeof(float));
               for (i=0; i<3; i++) {
                   memmove(Old.Parts[i],Vars.Parts[i],TotalNParts*sizeof(float));
                   memmove(Old.PartsV[i],Vars.PartsV[i],TotalNParts*sizeof(float));
               }
               for (i=_back; i<=_top; i++)
                   memmove(Old.Bounds[i],Vars.Bounds[i],NumEqs*Vars.BoundSize/(NumEqs+NSubst));
               break;
          case CHK_BeginPhase2 :
               for (i=0; i<NumEqs; i++)
                   {
                    memmove(Templ.Name[i],&HBuff[BegOffs(i)],(NZ+2)*NX*NY*sizeof(float));
                    memmove(&HBuff[BegOffs(i)],Old.Name[i],(NZ+2)*NX*NY*sizeof(float));
                   }
               for (i=0; i<3; i++) {
                   memmove(Templ.Parts[i],Vars.Parts[i],TotalNParts*sizeof(float));
                   memmove(Templ.PartsV[i],Vars.PartsV[i],TotalNParts*sizeof(float));
                   memmove(Vars.Parts[i],Old.Parts[i],TotalNParts*sizeof(float));
                   memmove(Vars.PartsV[i],Old.PartsV[i],TotalNParts*sizeof(float));
               }
               for (i=_back; i<=_top; i++) {
                   memmove(Templ.Bounds[i],Vars.Bounds[i],NumEqs*Vars.BoundSize/(NumEqs+NSubst));
                   memmove(Vars.Bounds[i],Old.Bounds[i],NumEqs*Vars.BoundSize/(NumEqs+NSubst));
               }
               break;
          case CHK_Recalculate :
               for (i=0; i<NumEqs; i++)
                   memmove(&HBuff[BegOffs(i)],Old.Name[i],(NZ+2)*NX*NY*sizeof(float));
               for (i=0; i<3; i++) {
                   memmove(Vars.Parts[i],Old.Parts[i],TotalNParts*sizeof(float));
                   memmove(Vars.PartsV[i],Old.PartsV[i],TotalNParts*sizeof(float));
               }
               for (i=_back; i<=_top; i++)
                   memmove(Vars.Bounds[i],Old.Bounds[i],NumEqs*Vars.BoundSize/(NumEqs+NSubst));
        }

        if (NSubst)
           memmove(UzSave,&HBuff[BegBoard(PhaseVars[CarrierPhase]._Uz)],BoardSize);

        for (i=0; i<NumPhases; i++)
            {
             if (!PhaseVars[i].IsLight) CalculateWXY(Vars.Bounds,Area,&WXYZ[i],&HBuff[BegBoard(PhaseVars[i]._Ux)],&HBuff[BegBoard(PhaseVars[i]._Uy)],Stage>=NExp);
             CalculateWZ(Vars.Bounds,Area,&WXYZ[i],&HBuff[BegBoard(PhaseVars[i]._Uz)],*PhaseVars[i]._Uw,Stage>=NExp, UseOpenMP);
            }
        /* Скопировать турбулентную вязкость в рабочий массив */
        if (PhaseVars[CarrierPhase]._Nu>=0)
           memmove(Dt,&HBuff[BegOffs(PhaseVars[CarrierPhase]._Nu)],(NZ+2)*NY*NX*sizeof(float));
        else
           memset(Dt,0,(NZ+2)*NY*NX*sizeof(float));

        SlowMode = 1;
        FastMode = FastTAUDivider==1;
        for (FastIteration = 0; FastIteration<FastTAUDivider+(FastTAUDivider>1); FastIteration++)
            {
             int EnhanceOpenMP = UseOpenMP && UseEnhancedOpenMP && NLightSubst>0 && (NLightSubst%nSMP || NLightSubst>3*nSMP);
             if (FastIteration==1)
               {
                FastMode = 1;
                SlowMode = 0;
                TAU /= FastTAUDivider;
               }

             CalculateDIV(Vars.Bounds, Area,
                &HBuff[BegOffs(PhaseVars[CarrierPhase]._Ux)],
                &HBuff[BegOffs(PhaseVars[CarrierPhase]._Uy)],
                &HBuff[BegOffs(PhaseVars[CarrierPhase]._Uz)]);
             CalculateK(&WXYZ[VDefs[i].Phase],Vars.Bounds,HBuff,Kfs,KDn,Sfs,SDn,Stage>=NExp,Area,Boundaries,Maps);

             if (CheckTau && CheckNew && SlowMode) {
                const double h  = 1.0;
                const int    MP = 5;

                memset(ALPHA, 0, BoardSize);
                for (i = 0; i < NParts; i++) {
                    x = (int)(Vars.Parts[0][i] / h);
                    y = (int)(Vars.Parts[1][i] / h);
                    z = (int)(Vars.Parts[2][i] / h) - BaseZ;
                    ALPHA[z*NY*NX+y*NX+x] += 1.0;
                }
                for (z = 0, Ptr = 0; z < NZ; z++)
                    for (y = 0; y < NY; y++)
                        for (x = 0; x < NX; x++, Ptr++) {
                            int q = (int)(0.5 + ALPHA[Ptr]);
                            ALPHA[Ptr] = 1.0 - 0.5*(q <= 1 ? 0.0 : (q < MP ? 1.0*q/MP : 1.0));
                        }

                if (CalcBase) {
                   for (i=0; i<NumEqs; i++)
                       if (VDefs[i].Solver == NULL)
                          EulerCalculate((float **) Vars.Bounds,
                                         &HBuff[BegBoard(i)], Kfs[i], Sfs[i], &WXYZ[VDefs[i].Phase],
                                         HXR, HYR, HZR,
                                         HX, HY, HZ,
                                         HXX, HYY, HZZ,
                                         HXS, HYS, HZS,
                                         &eHBuff[BegBoard(i)],
                                         *VDefs[i]._NuMol,&Dt[NY*NX],*VDefs[i]._Kappa,
                                         VDefs[i].Flags,
                                         ALPHA, Maps[i].Map, UseOpenMP
                          );
                   if (NSubst > 0) {
                      #pragma omp parallel if (EnhanceOpenMP)
                      {
                       #pragma omp for schedule(dynamic) private(i,j)
                       for (j=0;j<NLightSubst;j++)
                           {
                            i = MapSubsts[j];
                            EulerCalculate((float **) Vars.Bounds,
                                           &HBuff[BegBoard(NumEqs+i)], KDn[i], SDn[i], &WXYZ[CarrierPhase],
                                           HXR, HYR, HZR,
                                           HX, HY, HZ,
                                           HXX, HYY, HZZ,
                                           HXS, HYS, HZS,
                                           &eHBuff[BegBoard(NumEqs+i)],
                                           D, &Dt[NY*NX], 0.0,
                                           1,
                                           ALPHA, CMap, UseOpenMP && (1-EnhanceOpenMP)
                            );
                           }
                      }
                      for (j=0;j<NSubst-NLightSubst;j++)
                          {
                           i = MapSubsts[NLightSubst+j];
                           CalculateWZ((float **) Vars.Bounds, Area, &WXYZ[CarrierPhase], UzSave, (double)KGlobal.Uw[i], Stage>=NExp, UseOpenMP);
                           EulerCalculate((float **) Vars.Bounds,
                                          &HBuff[BegBoard(NumEqs+i)], KDn[i], SDn[i], &WXYZ[CarrierPhase],
                                          HXR, HYR, HZR,
                                          HX, HY, HZ,
                                          HXX, HYY, HZZ,
                                          HXS, HYS, HZS,
                                          &eHBuff[BegBoard(NumEqs+i)],
                                          D, &Dt[NY*NX], 0.0,
                                          1,
                                          ALPHA, CMap, UseOpenMP
                           );
                          }
                      if (NSubst-NLightSubst > 0)
                         CalculateWZ((float **) Vars.Bounds, Area, &WXYZ[CarrierPhase], UzSave, 0.0, Stage>=NExp, UseOpenMP);
                   }
                }  
             }

             if (CalcBase)
              {
               for (i=0; i<NumEqs; i++)
                   Process_Z(Vars.Bounds,VDefs[i].Solver,i,Kfs[i],Sfs[i],*VDefs[i]._NuMol,*VDefs[i]._Kappa,&WXYZ[VDefs[i].Phase],Maps[i].Map,Area,Boundaries,VDefs[i].Restrict,
                      VDefs[i].Flags,HZreg,VDefs[i].Projection, UseOpenMP);
               for (i=0; i<NumEqs; i++)
                   {
                    if (!VDefs[i].Solver)
                       Process_XY(Vars.Bounds,i,Kfs[i],Sfs[i],*VDefs[i]._NuMol,*VDefs[i]._Kappa,&WXYZ[VDefs[i].Phase],Maps[i].Map,Area,Boundaries,VDefs[i].Restrict,
                                  VDefs[i].Flags,HXreg,HYreg,VDefs[i].Projection, UseOpenMP, ReverseDirs[i]);
                    ReverseDirs[i] = 1-ReverseDirs[i];
                   }
              }
             if (SlowMode)
                {
#if defined(__MPI__) || defined(__ROUTER__) || defined(__ROUTER_100__)
                  /* Ожидаем Dn[i] */
                 DebugPrintf(DEBUG_FILE,"Wait for substances\n");
                 if (Packet.RecvData && NSubst)
                    {
                     Wait(&Request);
                     if (!BoundsInitialized)
                        BoundsPtr = InitSpecialBounds(Area,NumEqs,NSubst+NumEqs-1,BoundsPtr);
                    }
                 BoundsInitialized = 1;
                 DebugPrintf(DEBUG_FILE,"Substances received\n");
#endif
                 switch (Packet.CheckCommand) {
                   case CHK_Phase1 :
                        for (i=0; i<NSubst; i++)
                            memmove(Old.Dn[i],&HBuff[BegOffs(NumEqs+i)],(NZ+2)*NX*NY*sizeof(float));
                        j = NumEqs*Vars.BoundSize/((NumEqs+NSubst)*sizeof(float));
                        for (i=_back; i<=_top; i++)
                            memmove(&Old.Bounds[i][j],&Vars.Bounds[i][j],NSubst*Vars.BoundSize/(NumEqs+NSubst));
                        break;
                   case CHK_BeginPhase2 :
                        for (i=0; i<NSubst; i++)
                            {
                             memmove(Templ.Dn[i],&HBuff[BegOffs(NumEqs+i)],(NZ+2)*NX*NY*sizeof(float));
                             memmove(&HBuff[BegOffs(NumEqs+i)],Old.Dn[i],(NZ+2)*NX*NY*sizeof(float));
                            }
                        j = NumEqs*Vars.BoundSize/((NumEqs+NSubst)*sizeof(float));
                        for (i=_back; i<=_top; i++) {
                            memmove(&Templ.Bounds[i][j],&Vars.Bounds[i][j],NSubst*Vars.BoundSize/(NumEqs+NSubst));
                            memmove(&Vars.Bounds[i][j],&Old.Bounds[i][j],NSubst*Vars.BoundSize/(NumEqs+NSubst));
                        }
                        break;
                   case CHK_Recalculate :
                        for (i=0; i<NSubst; i++)
                            memmove(&HBuff[BegOffs(NumEqs+i)],Old.Dn[i],(NZ+2)*NX*NY*sizeof(float));
                        j = NumEqs*Vars.BoundSize/((NumEqs+NSubst)*sizeof(float));
                        for (i=_back; i<=_top; i++)
                            memmove(&Vars.Bounds[i][j],&Old.Bounds[i][j],NSubst*Vars.BoundSize/(NumEqs+NSubst));
                 }
                }

             #pragma omp parallel if (EnhanceOpenMP)
             #pragma omp for schedule(dynamic) private(i,j)
             for (j=0;j<NLightSubst;j++)
                 {
                  i = MapSubsts[j];
                  Process_Z(Vars.Bounds,NULL,NumEqs+i,KDn[i],SDn[i],D,0.0,&WXYZ[CarrierPhase],CMap,Area,Boundaries,rsPositive,1,HZreg,0, UseOpenMP && (1-EnhanceOpenMP));
                 }

             for (j=0;j<NSubst-NLightSubst;j++)
                 {
                  i = MapSubsts[NLightSubst+j];
                  if (SlowMode) CalculateWZ(Vars.Bounds,Area,&WXYZ[CarrierPhase],UzSave,(double)KGlobal.Uw[i],Stage>=NExp, UseOpenMP);
                  Process_Z(Vars.Bounds,NULL,NumEqs+i,KDn[i],SDn[i],D,0.0,&WXYZ[CarrierPhase],CMap,Area,Boundaries,rsPositive,1,HZreg,0, UseOpenMP);
                 }

             if (NSubst && SlowMode && NSubst!=NLightSubst) CalculateWZ(Vars.Bounds,Area,&WXYZ[CarrierPhase],UzSave,0.0,Stage>=NExp, UseOpenMP);
             #pragma omp parallel if (EnhanceOpenMP)
             #pragma omp for schedule(dynamic) private(i)
             for (i=0;i<NSubst;i++)
                 {
                  Process_XY(Vars.Bounds,NumEqs+i,KDn[i],SDn[i],D,0.0,&WXYZ[CarrierPhase],CMap,Area,Boundaries,rsPositive,1,HXreg,HYreg,0, UseOpenMP && (1-EnhanceOpenMP), ReverseDirs[NumEqs+i]);
                  ReverseDirs[NumEqs+i] = 1-ReverseDirs[NumEqs+i];
                 }
            }

        if (FastTAUDivider>1) TAU *= FastTAUDivider;

        if (CheckTau && CheckNew) {
           if (CalcBase) {
              for (i=0; i<NumEqs; i++)
                  if (VDefs[i].Solver == NULL) {
                     int _Ptr = BegBoard(i);
                     for (z = 0, Ptr = 0; z < NZ; z++)
                         for (y = 0; y < NY; y++)
                             for (x = 0; x < NX; x++, Ptr++, _Ptr++)
                                 if (Maps[i].Map[Ptr] == Zero && fabs(ALPHA[Ptr] - 1.0f) > 1E-5f)
                                    HBuff[_Ptr] = ALPHA[Ptr]*HBuff[_Ptr] + (1.0f-ALPHA[Ptr])*eHBuff[_Ptr];
                  }
              if (NSubst > 0)
                  for (i = 0; i < NSubst; i++) {
                      int _Ptr = BegBoard(NumEqs+i);
                      for (z = 0, Ptr = 0; z < NZ; z++)
                          for (y = 0; y < NY; y++)
                              for (x = 0; x < NX; x++, Ptr++, _Ptr++)
                                  if (CMap[Ptr] == Zero && fabs(ALPHA[Ptr] - 1.0f) > 1E-5f)
                                     HBuff[_Ptr] = ALPHA[Ptr]*HBuff[_Ptr] + (1.0f-ALPHA[Ptr])*eHBuff[_Ptr];
                  }
           }  
        }

        if (Packet.CheckCommand==CHK_Phase1) {
           SaveModelTime = ModelTime;
           SaveNParts = NParts;
        } else if (Packet.CheckCommand==CHK_BeginPhase2 || Packet.CheckCommand==CHK_Recalculate) {
           ModelTime = SaveModelTime;
           NParts = SaveNParts;
        }

        ModelTime += TAU;
        if (fmod(ModelTime,60.0)<2.0*TAU)
           GetSolarAngle(ModelTime/3600.0);

        if (CheckTau && CheckNew) {
           /* Пересчет частиц-индикаторов */
           #pragma omp parallel if (UseOpenMP && UseEnhancedOpenMP)
           #pragma omp for schedule(guided) private(i,x,y,z,Ptr)
           for (i = 0; i < NParts; i++) {
               double h = 1.0;
               float X = Vars.Parts[0][i];
               float Y = Vars.Parts[1][i];
               float Z = Vars.Parts[2][i] - BaseZ;
               x = (int) (X/h);
               y = (int) (Y/h);
               z = (int) (Z/h);
               Ptr = z*NY*NX+y*NX+x;
               if (x > 0 && x < NXs && y > 0 && y < NYs && z > 0 && z < NZs) {
                  double VX = Vars.PartsV[0][i];
                  double VY = Vars.PartsV[0][i];
                  double VZ = Vars.PartsV[0][i];
                  double V = sqrt(VX*VX+VY*VY+VZ*VZ);
                  double PX = k2*(GRAV[ZYXP]-GRAV[ZYXM])/2.0/h - k4*VX*V;
                  double PY = k2*(GRAV[ZYPX]-GRAV[ZYMX])/2.0/h - k4*VY*V;
                  double PZ = k2*(GRAV[ZPYX]-GRAV[ZMYX])/2.0/h - k4*VZ*V;
                  float NewX, NewY, NewZ;
                  Vars.PartsV[0][i] += (float)(TAU*PX);
                  Vars.PartsV[1][i] += (float)(TAU*PY);
                  Vars.PartsV[2][i] += (float)(TAU*PZ);
                  NewX = (float)(X + TAU*Vars.PartsV[0][i]);
                  NewY = (float)(Y + TAU*Vars.PartsV[1][i]);
                  NewZ = (float)(Z + TAU*Vars.PartsV[2][i]);
                  if (NewX >= 1.0 && NewX < NXs && NewY >= 1.0 && NewY < NYs &&
                      NewZ >= 0.0 && NewZ < NZ &&
                      !((FirstFlag || !DownExchange) && NewZ < 1.0) &&
                      !((LastFlag  || !UpExchange) && NewZ >= NZs)
                     ) {
                     Vars.Parts[0][i] = NewX;
                     Vars.Parts[1][i] = NewY;
                     Vars.Parts[2][i] = BaseZ + NewZ;
                     if (NewZ < 1.0)
                        Vars.Parts[0][i] = -Vars.Parts[0][i]; /* Это признак того, что частица улетит вниз из области этого процессора */
                     else if (NewZ >= NZs)
                        Vars.Parts[1][i] = -Vars.Parts[1][i]; /* Это признак того, что частица улетит вверх из области этого процессора */
                  }
               }
           }
           MakeExchange("Start PARTS exchange\n","Stop PARTS exchange\n",
   #ifdef __MPI__
                        ExchGRAVRequests,
   #endif
                        HUG,HDG,HUrG,HDrG,
                        PackPARTSForUpExchange, UnPackPARTSAfterUpExchange,
                        PackPARTSForDownExchange, UnPackPARTSAfterDownExchange,
                        UsizeG,DsizeG,GHUG,GHDG);
           DebugPrintf(DEBUG_FILE, "%i particles-indicators detected\n", NParts);
        }
       }

    if (NReact)
       {
        KGlobal.InitH   = min(KGlobal.InitH,TAU*0.001f);
        for (i=0; i<nSMP; i++)
            Contexts[i]._EndTime = (float)TAU;

        if (UseGear)
           {
            WorkLoad[SlaveID] = 0;
            /* Подсчитываем "свою" суммарную загрузку */
            if (CalculateAll)
               for (Count=0; Count<ZeroNodes; Count++)
                   WorkLoad[SlaveID] += SortedIndex[Count].Iters;
            /* Передаем значение "своей" суммарной загрузки всем остальным */
#ifdef __MPI__
            MPI_Allgather(MPI_IN_PLACE, sizeof(long), MPI_BYTE,
                          WorkLoad, sizeof(long), MPI_BYTE,
                          SlavesComm);
#else
            for (i=0; i<TotalProcs-1; i++)
                if (i==SlaveID)
                   {
                    for (j=0; j<TotalProcs-1; j++)
                        if (j!=SlaveID)
                           SendNeighbour(j, (byte *) &WorkLoad[i], sizeof(long));
                   }
                else
                   RecvNeighbour(i, (byte *) &WorkLoad[i], sizeof(long));
#endif
            /* Подсчитываем тотальную загрузку */
            TotalLoad = 0;
            for (i=0; i<TotalProcs-1; i++)
                TotalLoad += WorkLoad[i];
            AverageLoad = 1 + TotalLoad/(TotalProcs-1);
            DeltaLoad   = 1 + TotalLoad/1000;

            /* Определяем принимающие и передающие процессоры */
            NumSenders   = 0;
            NumReceivers = 0;
            GearSenderFlag   = 0;
            GearReceiverFlag = 0;
            for (i=0; i<TotalProcs-1; i++)
                {
                 Disbalance = WorkLoad[i] - AverageLoad;

                 if (Disbalance>=DeltaLoad)
                    {
                     Senders[NumSenders].Index  = i;
                     Senders[NumSenders++].Load = Disbalance;
                     if (i==SlaveID)
                        GearSenderFlag = 1;
                    }
                 else
                    if (Disbalance<=(-DeltaLoad))
                       {
                        Receivers[NumReceivers].Index  = i;
                        Receivers[NumReceivers++].Load = (-Disbalance);
                        if (i==SlaveID)
                           GearReceiverFlag = 1;
                       }
                }
            /* Сортируем передающие и принимающие процессоры по убыванию дисбаланса */
            if (NumSenders>1)
               qsort(Senders,NumSenders,sizeof(_ProcessorLoad),ProcessorLoadCmp);
            if (NumReceivers>1)
               qsort(Receivers,NumReceivers,sizeof(_ProcessorLoad),ProcessorLoadCmp);

            /* Создаем список транзакций */
            NumTransactions = 0;
            SenderCount     = 0;
            ReceiverCount   = 0;
            while (SenderCount<NumSenders && ReceiverCount<NumReceivers)
              if (Senders[SenderCount].Load==0)
                 SenderCount++;
              else
                if (Receivers[ReceiverCount].Load==0)
                   ReceiverCount++;
                else
                  {
                   if (Receivers[ReceiverCount].Load>Senders[SenderCount].Load)
                      CurLoad = Senders[SenderCount].Load;
                   else
                      CurLoad = Receivers[ReceiverCount].Load;
                   Senders[SenderCount].Load     -= CurLoad;
                   Receivers[ReceiverCount].Load -= CurLoad;
                   if (SlaveID==Senders[SenderCount].Index)
                      {
                       Transactions[NumTransactions].Index  = Receivers[ReceiverCount].Index;
                       Transactions[NumTransactions++].Load = CurLoad;
                      }
                   else
                     if (SlaveID==Receivers[ReceiverCount].Index)
                        {
                         Transactions[NumTransactions].Index  = Senders[SenderCount].Index;
                         Transactions[NumTransactions++].Load = CurLoad;
                        }
                  }
            /* Сортируем транзакции по возрастанию номеров процессоров */
            if (NumTransactions>1)
               qsort(Transactions,NumTransactions,sizeof(_ProcessorLoad),TransactionCmp);

            for (i=0; i<NumTransactions; i++)
                if (GearSenderFlag)
                   {
                    CurProc = Transactions[i].Index;
                    CurLoad = Transactions[i].Load;

                    Longs[CurProc] = 0;
                    Ptr            = 0;

                    /* Составляем списки узлов, которые будут отправлены на обработку */
                    for (Count = 0; CurLoad>0 && Count<ZeroNodes; Count++)
                        if (SortedIndex[Count].Iters)
                           if (SortedIndex[Count].Iters<=CurLoad)
                              {
                               CurLoad -= SortedIndex[Count].Iters;
                               if ((++Longs[CurProc])>MaxLongs[CurProc])
                                  {
                                   MaxLongs[CurProc] += DeltaLong;
                                   ArrangeIndex[CurProc] =
                                       (unsigned int *) SafeRealloc(ArrangeIndex[CurProc],
                                                                    MaxLongs[CurProc]*sizeof(unsigned int));
                                   Buffers[CurProc] =
                                       (float *) SafeRealloc(Buffers[CurProc],
                                                             MaxLongs[CurProc]*(NASubst+1)*sizeof(float));
                                  }
                               ArrangeIndex[CurProc][Longs[CurProc]-1] = SortedIndex[Count].Index;
                               for (j=0; j<NASubst; j++)
                                   Buffers[CurProc][Ptr++] = HBuff[OffsBufDn(KGlobal.TranMap[j])+SortedIndex[Count].Index];
                               if (PhaseVars[CarrierPhase]._T>=0)
                                  Buffers[CurProc][Ptr++] = HBuff[BegBoard(PhaseVars[CarrierPhase]._T)+SortedIndex[Count].Index];
                               else
                                  Buffers[CurProc][Ptr++] = (float)KGlobal.DefaultT;
                               SortedIndex[Count].Iters = 0;
                              }
                    /* Отправляем узлы списком */
                    SendNeighbour(CurProc,(byte *) &Longs[CurProc], sizeof(unsigned int));
                    if (Longs[CurProc])
                       SendNeighbour(CurProc,(byte *) Buffers[CurProc],Longs[CurProc]*(NASubst+1)*sizeof(float));
#if defined(__MPI__) || defined(__ROUTER__) || defined(__ROUTER_100__)
                    /* Сразу инициируем прием обработанного списка узлов */
                    if (Longs[CurProc])
                       IRecvNeighbour(CurProc,(byte *) Buffers[CurProc],
                           Longs[CurProc]*(NASubst+1)*sizeof(float),
                           &GearRequests[i]);
#endif
                   }
                else
                   if (GearReceiverFlag)
                      {
                       CurProc = Transactions[i].Index;

                       /* Принимаем узлы списком */
                       RecvNeighbour(CurProc,(byte *) &Longs[CurProc], sizeof(unsigned int));
                       if (Longs[CurProc]>MaxLongs[CurProc])
                          {
                           MaxLongs[CurProc] = Longs[CurProc]+DeltaLong;
                           ArrangeIndex[CurProc] =
                               (unsigned int *) SafeRealloc(ArrangeIndex[CurProc],
                                                            MaxLongs[CurProc]*sizeof(unsigned int));
                           Buffers[CurProc] =
                               (float *) SafeRealloc(Buffers[CurProc],
                                                     MaxLongs[CurProc]*(NASubst+1)*sizeof(float));
                          }
                       if (Longs[CurProc])
#if defined(__MPI__) || defined(__ROUTER__) || defined(__ROUTER_100__)
                          IRecvNeighbour(CurProc,(byte *) Buffers[CurProc],
                              Longs[CurProc]*(NASubst+1)*sizeof(float),
                              &GearRequests[i]);
#else
                          RecvNeighbour(CurProc,(byte *) Buffers[CurProc],Longs[CurProc]*(NASubst+1)*sizeof(float));
#endif
                      }

            RealLoad = 0;
            WCount = 0;
            if (CalculateAll)
              {
               /* Check defaults with order of handled nodes... */
               #pragma omp parallel if(UseOpenMP) private(Count,i,Ptr) shared(WCount)
               #pragma omp for schedule(dynamic,OMP_KineticChunk) reduction(+:RealLoad)
               for (Count = 0; Count<ZeroNodes; Count++)
                   if (SortedIndex[Count].Iters)
                      {
#ifdef __OPENMP__
                       KineticContext * CC = &Contexts[omp_get_thread_num()];
#else
                       KineticContext * CC = Contexts;
#endif
                       int _WCount;

                       Ptr = SortedIndex[Count].Index;

                       CC->LossPrecision = 0;
                       if (PhaseVars[CarrierPhase]._T>=0)
                          SetTk(CC,HBuff[BegBoard(PhaseVars[CarrierPhase]._T)+Ptr]);
                       else
                          SetTk(CC,(float)KGlobal.DefaultT);
                       for (i=0; i<NASubst; i++)
                           CC->Conc0[i] = HBuff[OffsBufDn(KGlobal.TranMap[i])+Ptr];

                       OneTaktKinetic(&KGlobal,CC,UseGear,Adams_Rozhkov_Method);
                       #pragma omp critical
                         _WCount = WCount++;
                       _SortedIndex[_WCount].Index = Ptr;
                       _SortedIndex[_WCount].Iters = CC->Iters+1;
                       RealLoad += CC->Iters+1;

                       if (!CC->LossPrecision)
                          for (i=0; i<NASubst; i++)
                              HBuff[OffsBufDn(KGlobal.TranMap[i])+Ptr] = CC->Conc1[i];
                       else {
                          fprintf(StdOutput,"Slave(%i) : Loss precision in node (%i) : %s\n",SlaveID,Ptr,
                                 CC->KinErrorInfo.LossH==0 ? "Too many iterations" : "Local Tau tends to zero");
                          fflush(StdOutput);
                       }
                      }
               memmove(SortedIndex,_SortedIndex,WCount*sizeof(SortInfo));
              }
            for (i=0; i<NumTransactions; i++)
                {
                 CurProc = Transactions[i].Index;

#if defined(__MPI__) || defined(__ROUTER__) || defined(__ROUTER_100__)
                 /* Ожидаем прихода сообщения: */
                 /*     GearSender   - обработанного списка узлов */
                 /*     GearReceiver - списка узлов для обработки */
                 if (Longs[CurProc])
                    Wait(&GearRequests[i]);
#endif
                 if (GearReceiverFlag)
                    {
                     double addRealLoad = 0.0;
		     int    MaxCount    = Longs[CurProc];

                     #pragma omp parallel if(UseOpenMP)
                     #pragma omp for schedule(dynamic,OMP_KineticChunk) reduction(+:addRealLoad) private(Count,j,Ptr,WPtr)
                     for (Count = 0; Count<MaxCount; Count++)
                         {
#ifdef __OPENMP__
                          KineticContext * CC = &Contexts[omp_get_thread_num()];
#else
                          KineticContext * CC = Contexts;
#endif
                          Ptr  = Count*(NASubst+1);
                          WPtr = Ptr;

                          CC->LossPrecision = 0;
                          for (j=0; j<NASubst; j++)
                              CC->Conc0[j] = Buffers[CurProc][Ptr++];
                          SetTk(CC,Buffers[CurProc][Ptr++]);

                          OneTaktKinetic(&KGlobal,CC,UseGear,Adams_Rozhkov_Method);
                          Buffers[CurProc][WPtr+NASubst] = CC->Iters+1.0f;
                          addRealLoad += CC->Iters+1;

                          if (!CC->LossPrecision)
                             for (j=0; j<NASubst; j++)
                                 Buffers[CurProc][WPtr++] = CC->Conc1[j];
                          else {
                             fprintf(StdOutput,"Slave(%i) : Loss precision in received node : %s\n",SlaveID,
                                    CC->KinErrorInfo.LossH==0 ? "Too many iterations" : "Local Tau tends to zero");
                             fflush(StdOutput);
                          }
                         }
                     RealLoad += (long) addRealLoad;
                     if (Longs[CurProc])
                        SendNeighbour(CurProc,(byte *) Buffers[CurProc],Longs[CurProc]*(NASubst+1)*sizeof(float));
                    }
                 else
                    if (GearSenderFlag)
                       {
#if !defined(__MPI__) && !defined(__ROUTER__)  && !defined(__ROUTER_100__)
                        /* При работе с MPI/Router обработанный список узлов уже принят */
                        if (Longs[CurProc])
                           RecvNeighbour(CurProc,(byte *) Buffers[CurProc],Longs[CurProc]*(NASubst+1)*sizeof(float));
#endif
                        for (Count=0, Ptr=0; (unsigned int)Count<Longs[CurProc]; Count++)
                            {
                             for (j=0; j<NASubst; j++)
                                 HBuff[OffsBufDn(KGlobal.TranMap[j])+ArrangeIndex[CurProc][Count]] = Buffers[CurProc][Ptr++];
                             SortedIndex[WCount].Index   = ArrangeIndex[CurProc][Count];
                             SortedIndex[WCount++].Iters = (unsigned short) Buffers[CurProc][Ptr++];
                            }
                       }
                }
            /* !!!!! */
            DebugPrintf(DEBUG_FILE,"Slave%i : PredictedLoad(%li) & RealLoad(%li)\n",SlaveID,AverageLoad,RealLoad);
            /* !!!!! */

            if (CalculateAll) {
               static int counter = 0;
               for (i = 0; i < ZeroNodes; i++) {
                   double p = SortedIndex[i].Iters;
                   savedIters = p;
                   prognose_load(counter, SortedIndex[i].Index, &p);
               }
               counter++;
               if (counter > 50)
                  for (i = 0; i < ZeroNodes; i++) {
                      double p = 100;
                      savedIters = SortedIndex[i].Iters;
                      predict_prognose_load(counter, SortedIndex[i].Index, &p);
                      if (p == p)
                         SortedIndex[i].Iters = p > 0 ? (int) p : 4;
                  }
               qsort(SortedIndex,ZeroNodes,sizeof(SortInfo),SortInfoCmp);
            }
           }
        else
          if (CalculateAll) {
             #pragma omp parallel if(UseOpenMP)
             #pragma omp for schedule(dynamic,OMP_KineticChunk) private(i,z,y,x,Ptr)
               for (Ptr=0; Ptr<NX*NY*NZ; Ptr++)
                   if (CMap[Ptr]==Zero)
                      {
#ifdef __OPENMP__
                       KineticContext * CC = &Contexts[omp_get_thread_num()];
#else
                       KineticContext * CC = Contexts;
#endif
                       CC->LossPrecision = 0;
                       if (PhaseVars[CarrierPhase]._T>=0)
                          SetTk(CC,HBuff[BegBoard(PhaseVars[CarrierPhase]._T)+Ptr]);
                       else
                          SetTk(CC,(float)KGlobal.DefaultT);
                       for (i=0;i<NASubst;i++)
                           CC->Conc0[i] = HBuff[OffsBufDn(KGlobal.TranMap[i])+Ptr];

                       OneTaktKinetic(&KGlobal,CC,UseGear,Adams_Rozhkov_Method);

                       if (!CC->LossPrecision)
                          for (i=0;i<NASubst;i++)
                              HBuff[OffsBufDn(KGlobal.TranMap[i])+Ptr] = CC->Conc1[i];
                       else
                          {
                           z = Ptr/(NX*NY);
                           y = (Ptr%(NX*NY))/NX;
                           x = Ptr%NX;
                           fprintf(StdOutput,"Slave(%i) : Loss precision in node (%i,%i,%i) : %s\n",SlaveID,x,y,z,
                                  CC->KinErrorInfo.LossH==0 ? "Too many iterations" : "Local Tau tends to zero");
                           fflush(StdOutput);
                          }
                      }
          }

        if (CalculateAll)
           {
            int EnhanceOpenMP = UseOpenMP && UseEnhancedOpenMP && NASubst>0 && (NASubst%nSMP || NASubst>3*nSMP);
            #pragma omp parallel if(EnhanceOpenMP)
            #pragma omp for schedule(dynamic) private(i)
            for (i=0;i<NASubst;i++)
                {
                 CorrectBounds(Vars.Bounds,&HBuff[OffsBufDn(KGlobal.TranMap[i])], NZ,NY,NX, NY*NX,NX,1, DescH,CMap,Area, 0, 0, UseOpenMP && (1-EnhanceOpenMP));
                 CorrectBounds(Vars.Bounds,&HBuff[OffsBufDn(KGlobal.TranMap[i])], NZ,NX,NY, NY*NX,1,NX, DescP,CMap,Area, 0, 1, UseOpenMP && (1-EnhanceOpenMP));
                 CorrectBounds(Vars.Bounds,&HBuff[OffsBufDn(KGlobal.TranMap[i])], NY,NX,NZ, NX,1,NY*NX, DescV,CMap,Area, 0, 2, UseOpenMP && (1-EnhanceOpenMP));
                }
           }
       }

    if (CalculateAll)
       {
        if (Packet.SendData)
           {
            for (i=0; i<NumEqs+NSubst; i++)
                {
                 memmove(&HD[(i*2)*NX*NY],&HBuff[BegOffs(i)],NX*NY*sizeof(float));
                 memmove(&HD[(i*2+1)*NX*NY],&HBuff[BegOffs(i+1)-NY*NX],NX*NY*sizeof(float));
                 memmove(&HBuff[i*PartL*NX*NY],&HBuff[(i*(NZ+2)+FromOffs)*NY*NX],PartL*NX*NY*sizeof(float));
                }
            SendMasterXXXTag((byte *) HBuff, (long) (NSubst+NumEqs)*PartL*NX*NY*sizeof(float));
            DebugPrintf(DEBUG_FILE,"Result sended to Master\n");
            for (i=NumEqs+NSubst-1; i>=0; i--)
                {
                 memmove(&HBuff[(i*(NZ+2)+FromOffs)*NY*NX],&HBuff[i*PartL*NX*NY],PartL*NX*NY*sizeof(float));
                 memmove(&HBuff[BegOffs(i)],&HD[(i*2)*NX*NY],NX*NY*sizeof(float));
                 memmove(&HBuff[BegOffs(i+1)-NY*NX],&HD[(i*2+1)*NX*NY],NX*NY*sizeof(float));
                }
            SendMaster((byte *) &Vars.BoundSize,sizeof(int));
            if (Vars.BoundSize)
               for (i=_back; i<=_top; i++)
                   SendMaster((byte *) Vars.Bounds[i],Vars.BoundSize);
           }
        if (Packet.CheckCommand==CHK_EndPhase2)
           {
            MaxEps = 0.0;
            memset(ERRS, 0, BoardSize);
            for (i=0; i<NumEqs; i++)
                if (VDefs[i].Solver == NULL)
                   CalculateEps(PartL,&MaxEps,&Templ.Name[i][FromOffs*NX*NY],&HBuff[(i*(NZ+2)+FromOffs)*NY*NX], ERRS);
            for (i=0;i<NSubst;i++)
                CalculateEps(PartL,&MaxEps,&Templ.Dn[i][FromOffs*NX*NY],&HBuff[((NumEqs+i)*(NZ+2)+FromOffs)*NY*NX], ERRS);

            if (CheckNew) {
               for (Ptr = 0; Ptr < NX*NY*NZ; Ptr++)
                   ERRS[Ptr] = (float)(ERRS[Ptr]<Eps ? +k1*ERRS[Ptr]/Eps : +k1);

               SolveCheckPoisson(GRAV, ERRS);

               MakeExchange("Start GRAV exchange\n","Stop GRAV exchange\n",
#ifdef __MPI__
                            ExchGRAVRequests,
#endif
                            HUG,HDG,HUrG,HDrG,
                            PackGRAVForUpExchange, UnPackGRAVAfterUpExchange,
                            PackGRAVForDownExchange, UnPackGRAVAfterDownExchange,
                            UsizeG,DsizeG,GHUG,GHDG);
            }
            SendMaster((byte *) &MaxEps, sizeof(MaxEps));
           }
       }

    if (CalculateAll && AllowPrediction && !CheckTau)
       {
        if (AllowQW && !PredictCount && Stage>=0)
           {
            for (i=0; i<ReCalcNPoints; i++)
                CheckTrace[i][Stage] = HBuff[BegBoard(PredictControlVar)+Ptrs[i]];
            if (Stage==NExp+NPred-1)
               {
                double AntiGrad[NExp];
                double Min;
                int    Flag = 1;

                for (i=0; i<NExp; i++)
                    KQ[i] = InitKQ;
                Min = CalcPredictErr(CheckTrace,1E+10,KQ,0,0.0);
                /* Вычисляем антиградиент */
                for (i=NExp-1; i>=0; i--)
                    {
                     double MinusV = CalcPredictErr(CheckTrace,Min,KQ,i,-StepKQ);
                     double PlusV  = CalcPredictErr(CheckTrace,Min,KQ,i,+StepKQ);

                     if (MinusV<PlusV)
                        {
                         AntiGrad[i] = -StepKQ;
                         Min = MinusV;
                        }
                     else if (PlusV<MinusV)
                        {
                         AntiGrad[i] = +StepKQ;
                         Min = PlusV;
                        }
                     else AntiGrad[i] = 0.0;
                     KQ[i] += AntiGrad[i];
                    }
                /* Спуск по антиградиенту */
                while (Flag)
                  for (i=NExp-1, Flag=0; i>=0; i--)
                      if (AntiGrad[i]!=0.0)
                         {
                          double NewMin = CalcPredictErr(CheckTrace,Min,KQ,i,AntiGrad[i]);
                          if (NewMin<Min)
                             {
                              Min = NewMin;
                              KQ[i] += AntiGrad[i];
                              Flag = 1;
                             }
                          else
                             AntiGrad[i] = 0.0;
                         }
                PreparePrediction(KQ);
                for (i=0; i<NExp; i++)
                    DebugPrintf(DEBUG_FILE,"%lf ",KQ[i]);
                DebugPrintf(DEBUG_FILE,"\n");
               }
           }
        if (AllowQW && Stage==NExp+NPred-1)
           if (++PredictCount==ReCalcNPs)
              PredictCount = 0;

        if (Stage<0)
           if (Stage == -ByPass) Stage = 0;
           else Stage--;
        else
           if (Stage == NExp+NPred-1) Stage = 0;
           else Stage++;
       }

    if (Packet._RecvCommand)
       RecvCommand(&Command);
    DebugPrintf(DEBUG_FILE,"End of current iteration\n");
   }
 while (Command);

#ifdef __MPI__
 if (UpExchange)
    {
     MPI_Request_free(&ExchRequests[0]);
     MPI_Request_free(&ExchRequests[1]);
     if (CheckTau) {
        MPI_Request_free(&ExchGRAVRequests[0]);
        MPI_Request_free(&ExchGRAVRequests[1]);
     }
    }
 if (DownExchange)
    {
     MPI_Request_free(&ExchRequests[2]);
     MPI_Request_free(&ExchRequests[3]);
     if (CheckTau) {
        MPI_Request_free(&ExchGRAVRequests[2]);
        MPI_Request_free(&ExchGRAVRequests[3]);
     }
    }
#endif

 free(HU);
 free(HD);
 free(HUr);
 free(HDr);

 if (CheckTau) {
    free(HUG);
    free(HDG);
    free(HUrG);
    free(HDrG);
 }

#if !defined(__ROUTER_100__)
 if (SHMEMUp) SHMEM_RM(GHU);
 if (SHMEMDown) SHMEM_RM(GHD);
 if (SHMEMUp && CheckTau) SHMEM_RM(GHUG);
 if (SHMEMDown && CheckTau) SHMEM_RM(GHDG);
#endif
 
 if (AllowPrediction && !CheckTau)
    {
     free(UpTrace); free(DownTrace);
     if (AllowQW)
        {
         free(CheckTrace); free(Ptrs);
        }
    }
 if (NReact && UseGear)
    {
#if defined(__MPI__) || defined(__ROUTER__) || defined(__ROUTER_100__)
     free(GearRequests);
#endif
     free(Transactions);
     free(WorkLoad);
     free(Senders);
     free(Receivers);
     for (i=0; i<TotalProcs-1; i++)
         free(Buffers[i]);
     free(Buffers);
     free(MaxLongs);
     free(Longs);
     if (CalculateAll)
        {
         for (i=0; i<TotalProcs-1; i++)
             free(ArrangeIndex[i]);
         free(ArrangeIndex);
         free(SortedIndex);
         free(_SortedIndex);
        }
    }
}
#endif

#define TauDiv 2

#ifdef __PARALLEL__
void OneCheckedTakt(int RecvData, int SendData)
#else
void OneCheckedTakt()
#endif
{
 int    SaveNR  = NReact;
 double SaveTau = TAU;

 double MaxEps = 0;
 double K;

 int Flag;
 int i;

// void OneTakt(int RecvCommand, int RecvData, int SendData, int CheckCommand, double * MaxEps)
// Old -- H[k]
// Templ -- H[k+1] : full TAU

#ifndef __PARALLEL__
 MoveBlocks(&Old,&Vars);
#endif
 if (NReact) {
	 NReact = 0;
	 KGlobal.NReact = 0;
 }
#ifdef __PARALLEL__
 OneTakt(0 /* No end syncrocommand */, RecvData, 0 /* No send data to master */, CHK_Phase1 /* Vars->Old */, NULL);
#else
 OneTakt();
 MoveBlocks(&Templ,&Vars);
 MoveBlocks(&Vars,&Old);
#endif

 TAU /= TauDiv;
 for (i=0;i<TauDiv;i++)
#ifdef __PARALLEL__
     OneTakt(i==(TauDiv-1) /* On the last stage will receive syncrocommand in the end */,
             0, /* No recv data from master */
             (i==(TauDiv-1) ? SendData : 0),
             (i==0 ? CHK_BeginPhase2 /* Vars->Templ, Old->Vars before calc */ :
                     (i==(TauDiv-1) ? CHK_EndPhase2 /* Master: collect EPS's. Slaves: calc EPS (Templ<->Vars), send to MASTER */:
                                      CHK_Phase2 /* Nothing special, calculate as in usual case */)),
             &MaxEps /* Master: union received (from slaves) MaxEps's */);
#else
     OneTakt();
 memset(ERRS, 0, BoardSize);
 for (i=0; i<NumEqs; i++)
     if (VDefs[i].Solver == NULL)
        CalculateEps(NZ,&MaxEps,Templ.Name[i],Vars.Name[i], ERRS);
 for (i=0;i<NSubst;i++)
     CalculateEps(NZ,&MaxEps,Templ.Dn[i],Vars.Dn[i], ERRS);
 if (CheckNew) {
    for (i = 0; i < NX*NY*NZ; i++)
        ERRS[i] = ERRS[i]<Eps ? +k1*ERRS[i]/Eps : +k1;
    SolveCheckPoisson(GRAV, ERRS);
 }
#endif

 fprintf(StdOutput,"Time = %lf : EPS = %lf\n", ModelTime, MaxEps);
 fflush(StdOutput);

 if (MaxEps<MIN_EPS && MaxEps<Eps)
    K = 2.0;
 else
    K = pow(Eps/max(1E-7, MaxEps), EpsP);

 if (K > 5.0) K = 5.0;
 if (K < 0.001) K = 0.001;

 NReact = SaveNR;
 KGlobal.NReact = NReact;
 if (Flag = (K<1.0 || (K>1.0 && !DecreaseOnly)))
    {
     TAU = SaveTau*K;
     if (TAU > MaxTau) TAU = MaxTau;
     fprintf(StdOutput,"ATTENTION : New TAU = %lf\n",TAU);
     fflush(StdOutput);
    }
 else
    TAU = SaveTau;

#ifdef __PARALLEL__
 SendCommand(1); // Send the syncrocommand
#endif
 if (Flag || NReact)
    {
#ifdef __PARALLEL__
     OneTakt(1, 0 /* No recv data from master */, SendData, CHK_Recalculate /* Old->Vars before calc */, NULL);
#else
     MoveBlocks(&Vars,&Old);
     OneTakt();
#endif
    }
}

#ifdef __PARALLEL__
void SendMap(int slave, int From, int Length, unsigned char * Map, int NumVals, float * Vals, int Power)
{
 SendSlave(slave,(byte *) &Map[From*NX*NY], Length*NX*NY*sizeof(char));
 SendSlave(slave,(byte *) &NumVals, sizeof(NumVals));
 if (NumVals)
    SendSlave(slave,(byte *) Vals, NumVals*Power*sizeof(float));
}

unsigned char * RecvMap(int FirstFlag, int LastFlag, unsigned char * AreaMap,
                        int * NumVals, float ** Vals, int Power)
{
 unsigned char * Result = (unsigned char *) SafeMalloc(NX*NY*NZ*sizeof(char));

 int Ptr;
 int x,y;

 RecvMaster((byte *) Result, NX*NY*NZ*sizeof(char));
 if (NumVals)
    {
     RecvMaster((byte *) NumVals, sizeof(int));
     if (*NumVals)
        {
         (*Vals) = (float *) SafeMalloc((*NumVals)*Power*sizeof(float));
         RecvMaster((byte *) (*Vals), (*NumVals)*Power*sizeof(float));
        }
     else
        (*Vals) = NULL;
    }
 if (AreaMap)
    {
     /* Корректируем карту, определяя, где необходим обмен */
     if (!FirstFlag)
        for (y=0, Ptr=0; y<NY; y++)
            for (x=0; x<NX; x++,Ptr++)
                if (AreaMap[Ptr]==0)
                   Result[Ptr] = ExchngBound;
     if (!LastFlag)
        for (y=0, Ptr=NZs*NY*NX; y<NY; y++)
            for (x=0; x<NX; x++,Ptr++)
                if (AreaMap[Ptr]==0)
                   Result[Ptr] = ExchngBound;
    }

 return Result;
}

#endif

FILE * OpenBoardFile(FILE * SavFile, int * Save, char * ResName, _int Cadres)
{
 FILE * Results = NULL;

 if (*Save)
    *Save =
       (
        SavFile && (Results = fopen(ResName,"r+b"))!=NULL
       ) ||
       (
        (Results = fopen(ResName,"w+b"))!=NULL
       );
 if (Results)
    {
     fseek(Results,0L,SEEK_END);
     if (ftell(Results)==0)
        {
         fprintf(Results,"%c%c%c%c",IBM_3D_SIGN,NX,NY,NZ);
         fwrite(&Cadres,sizeof(Cadres),1,Results);
        }
    }
 return Results;
}

void WriteBoardFile(float * Board, FILE * Results, _int Cadres)
{
 fwrite(Board,BoardSize,1,Results);
 fseek(Results,4L,SEEK_SET);
 fwrite(&Cadres,sizeof(Cadres),1,Results);
 fseek(Results,0L,SEEK_END);
 fflush(Results);
}

void CloseBoardFile(FILE * Results, _int Cadres)
{
 fseek(Results,4,SEEK_SET);
 fwrite(&Cadres,sizeof(Cadres),1,Results);
 fclose(Results);
}

enum EkoVars {
  varNX,varNY,varNZ,varHXmin,varHXmax,varHYmin,varHYmax,varHZmin,varHZmax,
        varTAU,varFastTAU,varTimeN,varTimeIters,varTimeTau,varTimeStopBase,
        varHX,varHY,varHZ,varSaveAllHistory,
  varNSubst,varNames,varUw,varZeroC,
  varNReact,varEA,varA,varTn,varReaction,varUseGear,
  varCheckTau,varMaxTau,varDecreaseOnly,varPeriod,varEps,varEpsP
};

void AfterNReact(void)
{
 KGlobal.NReact = NReact;
 KGlobal.NSubst = NSubst;
 if (NReact)
    if (Variables[varNSubst].VarUse && NSubst)
       {
		memset(KGlobal.Activity, 0, sizeof(KGlobal.Activity));;
        Variables[varEA].VarMap        = KGlobal.EA;
        Variables[varEA].VarMaxN       = NReact;
        Variables[varEA].VarMinN       = NReact;
        Variables[varA].VarMap         = KGlobal._A;
        Variables[varA].VarMaxN        = NReact;
        Variables[varA].VarMinN        = NReact;
        Variables[varTn].VarMap        = KGlobal.Tn;
        Variables[varTn].VarMaxN       = NReact;
        Variables[varTn].VarMinN       = NReact;
        Variables[varReaction].VarMaxN = 1;
        Variables[varReaction].VarMinN = 1;
       }
    else
      {
       fprintf(StdOutput,"NSubst>0 must be defined before NReact>0\n");
#ifdef __MVS__
       fclose(StdOutput);
#endif
       AbortServer(3);
      }
}

char * BeforeAnyRPrm(char * Msg)
{
 if (Variables[varNReact].VarUse)
    Msg[0] = 0;
 else
    strcpy(Msg,"NReact not defined");
 return Msg;
}

int HandleReaction(char * React)
{
 static int ReactCount = 0;

 return !NReact || Translate(&KGlobal,ReactCount++,React);
}

char * BeforeHX(char * Msg)
{
 if (Variables[varNX].VarUse)
    Msg[0] = 0;
 else
    strcpy(Msg,"NX not defined");
 return Msg;
}

char * BeforeHY(char * Msg)
{
 if (Variables[varNY].VarUse)
    Msg[0] = 0;
 else
    strcpy(Msg,"NX not defined");
 return Msg;
}

char * BeforeHZ(char * Msg)
{
 if (Variables[varNZ].VarUse)
    Msg[0] = 0;
 else
    strcpy(Msg,"NX not defined");
 return Msg;
}

unsigned char Eq2File[NumEqs];

int StandNSect = 0;

void AddVars()
{
 char VarName[20];
 int  Prev,i,j;

 AddSection("Area",1);
 AddSection("Substances",1);
 AddSection("Equations",1);
 AddSection("Accuracy",1);
 AddSection("Save",1);
 AddSection("Parameters",1);

 /* VarMap,VarName,VarType,VarMaxN,VarMinN,VarLong,VarBefore,VarHandle,VarAfter,VarUse,VarBcast,VarSect */
 AddVar(&NX,"NX",intT,1,1,1,NULL,NULL,NULL,0,0,0);
 AddVar(&NY,"NY",intT,1,1,1,NULL,NULL,NULL,0,0,0);
 AddVar(&NZ,"NZ",intT,1,1,1,NULL,NULL,NULL,0,0,0);
 AddVar(&HXmin,"HXmin",fltT,1,1,1,NULL,NULL,NULL,0,0,0);
 AddVar(&HXmax,"HXmax",fltT,1,1,1,NULL,NULL,NULL,0,0,0);
 AddVar(&HYmin,"HYmin",fltT,1,1,1,NULL,NULL,NULL,0,0,0);
 AddVar(&HYmax,"HYmax",fltT,1,1,1,NULL,NULL,NULL,0,0,0);
 AddVar(&HZmin,"HZmin",fltT,1,1,1,NULL,NULL,NULL,0,0,0);
 AddVar(&HZmax,"HZmax",fltT,1,1,1,NULL,NULL,NULL,0,0,0);
 AddVar(&TAU,"TAU",fltT,1,1,1,NULL,NULL,NULL,0,0,0);
 AddVar(&FastTAU,"FastTAU",fltT,1,1,1,NULL,NULL,NULL,1,0,0);
 AddVar(&TimeN,"TimeN",intT,1,1,1,NULL,NULL,NULL,1,1,0);
 AddVar(TimeIters,"TimeIters",fltT,MaxTimeN,0,1,NULL,NULL,NULL,1,1,0);
 AddVar(TimeTAUk,"TimeTAUk",fltT,MaxTimeN,0,1,NULL,NULL,NULL,1,1,0);
 AddVar(TimeStopBase,"TimeStopBase",intT,MaxTimeN,0,1,NULL,NULL,NULL,1,1,0);
 AddVar(_HX,"HX",fltT,MaxNX,0,1,BeforeHX,NULL,NULL,1,1,0);
 AddVar(_HY,"HY",fltT,MaxNX,0,1,BeforeHY,NULL,NULL,1,1,0);
 AddVar(_HZ,"HZ",fltT,MaxNX,0,1,BeforeHZ,NULL,NULL,1,1,0);
 AddVar(&SaveAllHistory,"SaveAllHistory",intT,1,1,1,NULL,NULL,NULL,1,1,0);
 AddVar(&NSubst,"NSubst",intT,1,1,1,NULL,NULL,NULL,0,0,1);
 AddVar(KGlobal.Names,"Names",chrT,MaxActSubst,0,SubstNameLength,NULL,NULL,NULL,0,0,1);
 AddVar(KGlobal.Uw,"Uw",ffltT,MaxActSubst,0,1,NULL,NULL,NULL,0,0,1);
 AddVar(ZeroC,"ZeroC",fltT,MaxActSubst,0,1,NULL,NULL,NULL,0,0,1);
 AddVar(&NReact,"NReact",intT,1,1,1,NULL,NULL,AfterNReact,0,0,2);
 AddVar(NULL,"EA",ffltT,0,0,1,BeforeAnyRPrm,NULL,NULL,0,0,2);
 AddVar(NULL,"A",ffltT,0,0,1,BeforeAnyRPrm,NULL,NULL,0,0,2);
 AddVar(NULL,"Tn",ffltT,0,0,1,BeforeAnyRPrm,NULL,NULL,0,0,2);
 AddVar(NULL,"Reaction",usrT,0,0,0,BeforeAnyRPrm,HandleReaction,NULL,0,0,2);
 AddVar(&UseGear,"UseGear",intT,1,1,1,NULL,NULL,NULL,0,0,2);
 AddVar(&SetSteadyKinetics,"SetSteady",intT,1,1,1,NULL,NULL,NULL,1,1,2); /* !!! */
 AddVar(&CheckTau,"CheckTau",intT,1,1,1,NULL,NULL,NULL,0,0,3);
 AddVar(&MaxTau,"MaxTau",fltT,1,1,1,NULL,NULL,NULL,0,0,3);
 AddVar(&DecreaseOnly,"DecreaseOnly",intT,1,1,1,NULL,NULL,NULL,0,0,3);
 AddVar(&Period,"Period",intT,1,1,1,NULL,NULL,NULL,0,0,3);
 AddVar(&Eps,"Eps",fltT,1,1,1,NULL,NULL,NULL,0,0,3);
 AddVar(&EpsP,"EpsP",fltT,1,1,1,NULL,NULL,NULL,0,0,3);

 for (i=0, j=-1, Prev=-1; i<NumEqs; i++)
     {
      if (VDefs[i].Projection<0 || Prev<0 || Prev>VDefs[i].Projection) j++;
      Eq2File[i] = j;
      if (VDefs[i].Projection<1)
         {
          if (VDefs[i].SubClass<0)
             sprintf(VarName,"Save%s",VDefs[i].FileBase);
          else
             sprintf(VarName,"Save%s%i",VDefs[i].FileBase,VDefs[i].SubClass);
          Saves[j] = 1;
          AddVar(&Saves[j],VarName,intT,1,1,1,NULL,NULL,NULL,1,0,4);
         }
      Prev = VDefs[i].Projection;
     }
 AddVar(SaveDn,"SaveDn",intT,MaxActSubst,0,1,NULL,NULL,NULL,0,0,4);

 AddVar(&Alpha,"Alpha",fltT,1,1,1,NULL,NULL,NULL,0,1,5);

 StandNSect = NSect;
 AddCustomVars();
 AddSection("Coordinates",1);
 AddVar(ModelDate,"Date",chrT,1,1,MaxMDLength,NULL,NULL,NULL,0,0,(char) (NSect-1));
 AddVar(&ModelGMTHour,"GMTHour",fltT,1,1,1,NULL,NULL,NULL,0,0,(char) (NSect-1));
 AddVar(&ModelAltitude,"Altitude",fltT,1,1,1,NULL,NULL,NULL,0,0,(char) (NSect-1));
 AddVar(&ModelLongitude,"Longitude",fltT,1,1,1,NULL,NULL,NULL,0,0,(char) (NSect-1));
 AddVar(&ModelAngle,"Angle",fltT,1,1,1,NULL,NULL,NULL,0,0,(char) (NSect-1));
 AddSection("Parallel",1);
 AddVar(&UseSHMEM,"UseSHMEM",intT,1,1,1,NULL,NULL,NULL,0,1,(char) (NSect-1));
 AddVar(&MinGap,"MinGap",intT,1,1,1,NULL,NULL,NULL,1,1,(char) (NSect-1));
 AddVar(Division,"Division",intT,sizeof(Division)/sizeof(Division[0]),1,1,NULL,NULL,NULL,1,1,(char) (NSect-1));
 AddSection("Default",1);
 AddVar(&KGlobal.DefaultT,"DefaultT",ffltT,1,1,1,NULL,NULL,NULL,1,1,(char) (NSect-1));
}

#ifdef __PARALLEL__
#define TempBoardSize ((NZ+2)*NX*NY*sizeof(float))
#else
#define TempBoardSize BoardSize
#endif

typedef char ProcessorName[128];

void RestoreLocalRefs() {
 int i;

 for (i = 0; i < OldRefNum; i++)
     *OldRefAddrs[i] = OldRefVals[i];
}

void CreateSavFile(FILE ** SavFile, _int maxtakt, int Quant, int takt, int SaveHisto)
{
 char Buf[256];
 int i, x;
#ifdef __PARALLEL__
 int BSize = RecvBoundSize;
#else
 int BSize = Vars.BoundSize;
#endif
 if (*SavFile)
    {
     fclose(*SavFile);
     *SavFile = NULL;
    }
 sprintf(Buf,"%s.%i",SavName,takt);
 if (!SaveHisto || rename(SavName,Buf))
    unlink(SavName);
 if (takt!=maxtakt && !(maxtakt == 0 && ModelTime >= EndTime))
    {
     *SavFile = fopen(SavName,"wb");
     fwrite(&maxtakt,sizeof(maxtakt),1,*SavFile);
     fwrite(&Quant,sizeof(Quant),1,*SavFile);
     fwrite(&takt,sizeof(takt),1,*SavFile);
#ifdef __PARALLEL__
     RestoreLocalRefs();
#endif
     for (i=0; i<NumEqs; i++)
         fwrite(Vars.Name[i],BoardSize,1,*SavFile);
     for (x=0;x<NSubst;x++)
         fwrite(Vars.Dn[x],BoardSize,1,*SavFile);
     if (BSize)
#ifdef __PARALLEL__
        fwrite(HBuff,BSize,1,*SavFile);
#else
        for (i=_back; i<=_top; i++)
            fwrite(Vars.Bounds[i],BSize,1,*SavFile);
#endif
     fwrite(&ModelTime,sizeof(ModelTime),1,*SavFile);
     fwrite(&EndTime,sizeof(EndTime),1,*SavFile);
     fprintf(StdOutput,"Experiment saved\n");
     fflush(StdOutput);
     fclose(*SavFile);
     *SavFile = NULL;
    }
}

/* Головная программа */
#if defined(__MPI__) || defined(__ROUTER_100__)
int main(int argc, char ** argv)
#else
int main(void)
#endif
{
 FILE * SavFile = NULL;
 FILE * DscFile = NULL;

 FILE * DResults[MaxActSubst];
 FILE * SaveBnds = NULL;

 char CfgFName[__MAXPATH];
 char Lmin2FName[__MAXPATH];
 ProcessorName * ProcessorNames = NULL;
 ProcessorName CurName = "";
 int ProcNameLength;

 FILE * Lmin2File;

 unsigned char * OrgnArea;
 int  i,j,x,y,z;
 int  Ptr;
 _int maxtakt; /* Число итераций */
 _int Cadres = 0;       /* Счетчик записываемых кадров */
 int  Quant;            /* Размер кадра в итерациях */
 int  QuantCount;
 int  TimeCounter = 0;
 int  Command;

 TIME_STRUCT _BeginTime;
 TIME_STRUCT _EndTime;

 int * Ranks;
 int Norm, NMod;
 int error;
 int START;
 int MaxDivision;
 int NN;

 unsigned int ChangeParams;

#ifdef __MVS__
 char Buf[256];

 sprintf(Buf,"%s%i",MVSOutputBase,rand());
 if (!(StdOutput = fopen(Buf,"w+t")))
    {
     printf("StdOutput '%s' creation error\n",Buf);
     AbortServer(100);
    }
#else
 StdInput  = stdin;
 StdOutput = stdout;
#endif

 KGlobal.MaxIterations = 1000;
 KGlobal.InitTime = 0.0f; /* Начальное время интегрирования */
 KGlobal.CalcEps = 1E-2f; /* Относительная точность интегрирования ST1/eps */
 KGlobal.MinH = 1.0E-18f;/* Минимальный шаг интегрирования ST1/hmin */

 KGlobal.DefaultT = 27.0f;

 KGlobal.InitH = 1E-7f;  /* Начальный шаг интегрирования h */

 KGlobal.AdamsOrder = 4;
 KGlobal.Fail = 0;

 KGlobal.KinetTraceWNotFound = 7;
 for (Ptr = 0; Ptr < markZ99; Ptr++) {
     KGlobal.KinetTraceW[Ptr] = 0.0f;
     KGlobal.KinetTraceWOrder[Ptr] = Ptr;
 }

#ifdef __PARALLEL__

#ifdef __MPI__
 if ((error = MPI_Init(&argc,&argv))!=MPI_SUCCESS)
    {
     fprintf(StdOutput,"MPI Initialization Error %i\n",error);
#ifdef __MVS__
     fclose(StdOutput);
#endif
     AbortServer(30);
    }
 if ((error = MPI_Barrier(MPI_COMM_WORLD))!=MPI_SUCCESS)
    {
     fprintf(StdOutput,"MPI Synchronization Error %i\n",error);
#ifdef __MVS__
     fclose(StdOutput);
#endif
     AbortServer(10);
    }
 MPI_Comm_size(MPI_COMM_WORLD,&nProcs);
 MPI_Comm_rank(MPI_COMM_WORLD,&MyProcID);
#elif defined(__ROUTER_100__)
 MyProcID = atoi(argv[1]);
 nProcs   = atoi(argv[2]);
#elif defined(__ROUTER__)
 if (!rf_create())
    {
     fprintf(StdOutput,"Router+ Initialization Error\n");
#ifdef __MVS__
     fclose(StdOutput);
#endif
     AbortServer(30);
    }
 nProcs   = sysProcTotalGet();
 MyProcID = sysProcNumGet();

 NumSends = (int *) SafeMalloc(nProcs*sizeof(int));
 MaxSends = (int *) SafeMalloc(nProcs*sizeof(int));
 HndSends = (SendDsc *) SafeMalloc(nProcs*DefMaxSends*sizeof(SendDsc));
 for (i=0; i<nProcs; i++)
     {
      NumSends[i] = 0;
      MaxSends[i] = DefMaxSends;
     }
#else
 nProcs   = GET_ROOT()->ProcRoot->nProcs;
 MyProcID = GET_ROOT()->ProcRoot->MyProcID;
#endif

 MasterFlag = MyProcID==0;

 if (nProcs<3)
    {
     if (MasterFlag)
        fprintf(StdOutput,"Must be more than 2 processors\n");
#ifdef __MVS__
     fclose(StdOutput);
#endif
     AbortServer(1);
    }
 DebugPrintf(DEBUG_FILE,"Initialised\n");
 /* Устанавливаем подтопологию "звезда" */
 #ifdef __MPI2REENT__
 sprintf(CurName,"MPI2REENT-%i",MyProcID);
 #elif !defined(__ROUTER_100__)
 gethostname(CurName,sizeof(CurName)-1);
 #else
 sprintf(CurName,"MVS-100-%i",MyProcID);
 #endif
 ProcNameLength = strlen(CurName);
 ProcOnNode = (int *) SafeMalloc(nProcs*sizeof(int));
 if (MasterFlag)
    {
     int CurProcID = 0;
     int ProcPerNode = 1;

     ProcessorNames = (ProcessorName *) SafeMalloc(nProcs*sizeof(ProcessorName));
     strcpy(ProcessorNames[0],CurName);
     fprintf(StdOutput,"Detected nodes:\n%s : ",ProcNameLength>0 ? CurName : "unknown node");
     fflush(StdOutput);
     ProcOnNode[0] = CurProcID;
     for (i=1; i<nProcs; i++)
         {
          RecvSlave(i-1,(byte *) &ProcNameLength, sizeof(ProcNameLength));
          RecvSlave(i-1,(byte *) (ProcessorNames+i), ProcNameLength+1);
          if (ProcNameLength==0 || strcmp(CurName,ProcessorNames[i])!=0)
             {
              strcpy(CurName,ProcessorNames[i]);
              fprintf(StdOutput,"%i\n%s : ",ProcPerNode,ProcNameLength>0 ? CurName : "unknown node");
              fflush(StdOutput);
              ProcPerNode = 1;
              CurProcID++;
             }
          else
             ProcPerNode++;
          ProcOnNode[i] = CurProcID;
         }
     fprintf(StdOutput,"%i\n",ProcPerNode);
     fflush(StdOutput);
     for (i=1; i<nProcs; i++)
         SendSlave(i-1,(byte *) ProcOnNode, nProcs*sizeof(int));
    }
 else
    {
     SendMaster((byte *) &ProcNameLength, sizeof(ProcNameLength));
     SendMaster((byte *) CurName, ProcNameLength+1);
     RecvMaster((byte *) ProcOnNode, nProcs*sizeof(int));
     if (MyProcID>1) SHMEMDown = ProcOnNode[MyProcID]==ProcOnNode[MyProcID-1];
     if (MyProcID<nProcs-1) SHMEMUp = ProcOnNode[MyProcID]==ProcOnNode[MyProcID+1];
    }

#endif

#ifdef __OPENMP__
#ifdef __OMP_MPI_H__
 nSMP = nSMPProcs;
#else
 nSMP = omp_get_num_procs();
#endif
#ifdef __PARALLEL__
 Contexts = (KineticContext *) malloc(nSMP*sizeof(KineticContext));
#else
 Contexts = (KineticContext *) malloc(imax(KinetChunkSize,nSMP)*sizeof(KineticContext));
#endif
#else
 KinetChunkSize = 1;
 Contexts = (KineticContext *) malloc(sizeof(KineticContext));
#endif

 InitVars();
 AddVars();

#ifdef __PARALLEL__
 TotalProcs = nProcs;
 if (MasterFlag)
    {
#endif
#ifndef __MVS__
     fprintf(StdOutput,"Input configuration file name = ");
#else
     fprintf(StdOutput,"Program started\n");
     fflush(StdOutput);
     if (!(StdInput = fopen(MVSConfigFile,"rt")))
        {
         fprintf(StdOutput,"StdInput '%s' not found\n",MVSConfigFile);
#ifdef __MVS__
         fclose(StdOutput);
#endif
         AbortServer(100);
        }
#endif
     if (StdInput==stdin)
        scanf("%s",CfgFName);
     else
        fscanf(StdInput,"%s",CfgFName);
     CreateNames(CfgFName);

#ifndef __MVS__
         fprintf(StdOutput,"Do you wish to change Params (1=Yes/0=No) : ");
#endif
         if (StdInput==stdin) {
            char _CRLF[51] = "";
	    i = scanf("%u",&ChangeParams);
            fgets(_CRLF, 50, stdin);
	 } else {
            char _CRLF[51] = "";
            i = fscanf(StdInput,"%u",&ChangeParams);
            fgets(_CRLF, 50, StdInput);
         }

         if (i < 1) {
            fprintf(StdOutput,"no input data\n");
#ifdef __MVS__
            fclose(StdOutput);
#endif
            AbortServer(-1300);
         }

     if (!ReadConfig(StdInput,StdOutput,ChangeParams,CfgFName,InitParameters))
        {
         fprintf(StdOutput,"%s\n",ErrorMsg);
#ifdef __MVS__
         fclose(StdOutput);
#endif
         AbortServer(-1);
        }
	 else {
		 KGlobal.NSubst = NSubst;
		 KGlobal.NReact = NReact;
		 KGlobal.NASubst = NASubst = 0;
		 if (NReact)
		 {
			 ReTranslate(&KGlobal);
			 NASubst = KGlobal.NASubst;
		 }
	 }
     if (FastTAU>0.0 && FastTAU<TAU) FastTAUDivider = (int) (TAU/FastTAU);
     TAU0 = TAU;
     for (i=0; i<NumPhases; i++)
         {
          PhaseLinks[i] = -1;
          if (PhaseVars[i].Source[0])
             {
              for (x=0; x<NSubst && PhaseLinks[i]<0; x++)
                  if (strcmp(KGlobal.Names[x],PhaseVars[i].Source)==0)
                     PhaseLinks[i] = x;
              if (PhaseLinks[i]<0)
                 {
                  fprintf(StdOutput,"Phase source '%s' is not found: Linking error\n",PhaseVars[i].Source);
#ifdef __MVS__
                  fclose(StdOutput);
#endif
                  AbortServer(-15);
                 }
             }
         }
     for (i=0; SubstRefs[i] && i<MaxActSubst; i++)
         {
          SLinks[i] = MaxActSubst;
          for (x=0; x<NSubst && SLinks[i]==MaxActSubst; x++)
              if (strcmp(KGlobal.Names[x],SubstRefs[i])==0)
                 SLinks[i] = x;
          if (SLinks[i]==MaxActSubst)
             {
              fprintf(StdOutput,"Substance '%s' is not found: Linking error\n",SubstRefs[i]);
#ifdef __MVS__
              fclose(StdOutput);
#endif
              AbortServer(-25);
             }
         }

     TotalNParts = NZ*NY*NX/4;
     NParts = TotalNParts;

#ifdef __PARALLEL__
 /* Определяем режимы работы процессоров */
 /* Должно быть не менее чем 2 ряда узлов на процессор */
     MaxDivision = NZ/(nProcs-1);
     if (Division[0])
        {
         for (i=1, KinetProcs = 0; i<TotalProcs && Division[i-1]; i++)
             {
              KinetProcs += Division[i-1];
              if (Division[i-1]>MaxDivision)
                 MaxDivision = Division[i-1];
             }
         if (i-1<2)
            {
             fprintf(StdOutput,"Config: 'Division' array: must be more than one non-zero contiguous blocks\n");
#ifdef __MVS__
             fclose(StdOutput);
#endif
             AbortServer(-900);
            }
         KinetProcs = TotalProcs-i;
        }
     else
        if (NZ<(MinGap+1)*(nProcs-1))
           KinetProcs = (nProcs-1)-(NZ/(MinGap+1));
        else
           KinetProcs = 0;
     if (KinetProcs) {
        fprintf(StdOutput,"Last (%i) processors will be used for Kinetic ONLY\n",KinetProcs);
        fflush(StdOutput);
     }
     nProcs -= KinetProcs;

     FromGrid   = (int *) SafeMalloc((TotalProcs-1)*sizeof(int));
     LengthGrid = (int *) SafeMalloc((TotalProcs-1)*sizeof(int));

     Norm = NZ / (nProcs-1);
     NMod = NZ % (nProcs-1);

     z = 0;
     for (i=1; i<TotalProcs; i++)
         {
          /* Отсылаем процессорам информацию о числе основных процессоров */
          SendSlave(i-1,(byte *) &nProcs,sizeof(nProcs));
          /* Отсылаем рабочим процессорам типы режимов работы */
          Command = (i<nProcs);
          SendSlave(i-1,(byte *) &Command,sizeof(Command));

          if (Command)
             {
              if (Division[0])
                 LengthGrid[i-1] = Division[i-1]+1;
              else
                 {
                  LengthGrid[i-1] = Norm+1;
                  if (i<=NMod)
                     LengthGrid[i-1]++;
                 }
              if (i>1 && i<nProcs-1)
                 LengthGrid[i-1]++;
              FromGrid[i-1] = z;
              z += LengthGrid[i-1]-2;
             }
          else
             {
              FromGrid[i-1]   = 0;
              LengthGrid[i-1] = 0;
             }
         }
     if (z+2!=NZ)
        {
         fprintf(StdOutput,"Config: 'Division' array: sum of block sizes is greater than area size\n");
         fflush(StdOutput);
        }
 /* Посылаем ВСЕМ рабочим процессорам прочитанные параметры эксперимента */
     for (i=1; i<TotalProcs; i++)
         {
          SendSlave(i-1,(byte *) &NX,sizeof(NX));
          SendSlave(i-1,(byte *) &NY,sizeof(NY));
          /* Вместо реального NZ отсылаем размер соответствующего участка */
          SendSlave(i-1,(byte *) &LengthGrid[i-1],sizeof(LengthGrid[i-1]));
          SendSlave(i-1,(byte *) &FromGrid[i-1],sizeof(int));
          SendSlave(i-1,(byte *) &TotalNParts,sizeof(TotalNParts));
          SendSlave(i-1,(byte *) &TAU,sizeof(TAU));
          SendSlave(i-1,(byte *) &FastTAU,sizeof(FastTAU));
          SendSlave(i-1,(byte *) &FastTAUDivider,sizeof(FastTAUDivider));
          SendSlave(i-1,(byte *) &CheckTau,sizeof(CheckTau));
          SendSlave(i-1,(byte *) &Eps,sizeof(Eps));
          SendSlave(i-1,(byte *) &NSubst,sizeof(NSubst));
          SendSlave(i-1,(byte *) &NReact,sizeof(NReact));
          SendSlave(i-1,(byte *) ModelDate,MaxMDLength*sizeof(char));
          SendSlave(i-1,(byte *) &ModelGMTHour,sizeof(ModelGMTHour));
          SendSlave(i-1,(byte *) &ModelAltitude,sizeof(ModelAltitude));
          SendSlave(i-1,(byte *) &ModelLongitude,sizeof(ModelLongitude));
          SendSlave(i-1,(byte *) &ModelAngle,sizeof(ModelAngle));
          if (NSubst)
             {
              SendSlave(i-1,(byte *) KGlobal.Uw,NSubst*sizeof(float));
              SendSlave(i-1,(byte *) PhaseLinks,NumPhases*sizeof(int));
              SendSlave(i-1,(byte *) SLinks,MaxActSubst*sizeof(char));
             }
          if (NReact)
             {
              SendSlave(i-1,(byte *) &NASubst,sizeof(NASubst));
              SendSlave(i-1,(byte *) KGlobal.TranMap,NASubst*sizeof(int));
              SendSlave(i-1,(byte *) &UseGear,sizeof(UseGear));
              SendSlave(i-1,(byte *) KGlobal.EA,NReact*sizeof(float));
              SendSlave(i-1,(byte *) KGlobal._A,NReact*sizeof(float));
              SendSlave(i-1,(byte *) KGlobal.Tn,NReact*sizeof(float));
              SendSlave(i-1,(byte *) KGlobal.LR,NReact*sizeof(Reaction));
             }
         }
     DebugPrintf(DEBUG_FILE,"Base information sended\n");
    }
 else
 /* Принимаем от мастера прочитанные параметры эксперимента */
    {
     RecvMaster((byte *) &nProcs,sizeof(nProcs));

     FirstFlag = MyProcID==1;
     LastFlag  = MyProcID==nProcs-1;

     RecvMaster((byte *) &CalculateAll,sizeof(CalculateAll));
     RecvMaster((byte *) &NX,sizeof(NX));
     RecvMaster((byte *) &NY,sizeof(NY));
     RecvMaster((byte *) &NZ,sizeof(NZ));
     RecvMaster((byte *) &BaseZ,sizeof(int));
     RecvMaster((byte *) &TotalNParts,sizeof(TotalNParts));
     RecvMaster((byte *) &TAU,sizeof(TAU));
     RecvMaster((byte *) &FastTAU,sizeof(FastTAU));
     RecvMaster((byte *) &FastTAUDivider,sizeof(FastTAUDivider));
     RecvMaster((byte *) &CheckTau,sizeof(CheckTau));
     RecvMaster((byte *) &Eps,sizeof(Eps));
     RecvMaster((byte *) &NSubst,sizeof(NSubst));
	 KGlobal.NSubst = NSubst;
     RecvMaster((byte *) &NReact,sizeof(NReact));
	 KGlobal.NReact = NReact;
     RecvMaster((byte *) ModelDate,MaxMDLength*sizeof(char));
     RecvMaster((byte *) &ModelGMTHour,sizeof(ModelGMTHour));
     RecvMaster((byte *) &ModelAltitude,sizeof(ModelAltitude));
     RecvMaster((byte *) &ModelLongitude,sizeof(ModelLongitude));
     RecvMaster((byte *) &ModelAngle,sizeof(ModelAngle));
     if (NSubst)
        {
         RecvMaster((byte *) KGlobal.Uw,NSubst*sizeof(float));
         RecvMaster((byte *) PhaseLinks,NumPhases*sizeof(int));
         RecvMaster((byte *) SLinks,MaxActSubst*sizeof(char));
        }
     if (NReact)
        {
		 memset(KGlobal.Activity,0,sizeof(KGlobal.Activity));;
         RecvMaster((byte *) &NASubst,sizeof(NASubst));
		 KGlobal.NASubst = NASubst;
         RecvMaster((byte *) KGlobal.TranMap,NASubst*sizeof(int));
         RecvMaster((byte *) &UseGear,sizeof(UseGear));
         RecvMaster((byte *) KGlobal.EA,NReact*sizeof(float));
         RecvMaster((byte *) KGlobal._A,NReact*sizeof(float));
         RecvMaster((byte *) KGlobal.Tn,NReact*sizeof(float));
         RecvMaster((byte *) KGlobal.LR,NReact*sizeof(Reaction));
        }
     DebugPrintf(DEBUG_FILE,"Base information received\n");
     if (CalculateAll)
        InitParameters();

     NParts = TotalNParts/(nProcs-1);
     i = TotalNParts%(nProcs-1);
     if (MyProcID <= i)
        NParts++;
    }
#endif

#ifdef __PARALLEL__

#ifdef __MPI__
 MPI_Comm_group(MPI_COMM_WORLD,&WorldGroup);

 Ranks = (int *) SafeMalloc((TotalProcs-1)*sizeof(int));
 for (i=0; i<TotalProcs-1; i++)
     Ranks[i] = i+1;
 MPI_Group_incl(WorldGroup,TotalProcs-1,Ranks,&SlavesGroup);
 free(Ranks);
 MPI_Comm_create(MPI_COMM_WORLD,SlavesGroup,&SlavesComm);
 if (!MasterFlag)
    MPI_Comm_rank(SlavesComm,&SlaveID);
#else
 SlaveID = MyProcID-1;
#endif

 if (MasterFlag)
    {
#endif
     /* Вычисляемые данные (текущая итерация) */
     for (i=0; i<NumEqs; i++) {
         Vars.Name[i]  = (float *) SafeMalloc(BoardSize);
         Vars.eName[i] = CheckTau == 0 ? NULL : (float *) SafeMalloc(BoardSize);
     }
     /* Вещества */
     Vars.Dn  = NSubst==0 ? NULL : (float **) SafeMalloc(NSubst*sizeof(float *));
     Vars.eDn = NSubst==0 ? NULL : (float **) SafeMalloc(NSubst*sizeof(float *));

     for (x=0;x<NSubst;x++) {
         Vars.Dn[x]  = (float *) SafeMalloc(BoardSize);
         Vars.eDn[x] = CheckTau == 0 ? NULL : (float *) SafeMalloc(BoardSize);
     }
     for (i = 0; i < 3; i++) {
         Vars.Parts[i]  = CheckTau == 0 ? NULL : (float *) SafeMalloc(TotalNParts*sizeof(float));
         Vars.PartsV[i] = CheckTau == 0 ? NULL : (float *) SafeMalloc(TotalNParts*sizeof(float));
     }
#ifdef __PARALLEL__
    }
 else
   if (CalculateAll)
      {
#endif
       Dt  = (float *) SafeMalloc(TempBoardSize); /* Турбулентная вязкость */

       /* Коэффициенты,определяющие вхождение противоточных производных в решение */
       for (i=0; i<NumPhases; i++)
           {
            WXYZ[i].WX1 = PhaseVars[i].IsLight ? WXYZ[CarrierPhase].WX1 : (float *) SafeMalloc(BoardSize);
            WXYZ[i].WY1 = PhaseVars[i].IsLight ? WXYZ[CarrierPhase].WY1 : (float *) SafeMalloc(BoardSize);
            WXYZ[i].WZ1 = (float *) SafeMalloc(BoardSize);
            WXYZ[i].WX2 = PhaseVars[i].IsLight ? WXYZ[CarrierPhase].WX2 : (float *) SafeMalloc(BoardSize);
            WXYZ[i].WY2 = PhaseVars[i].IsLight ? WXYZ[CarrierPhase].WY2 : (float *) SafeMalloc(BoardSize);
            WXYZ[i].WZ2 = (float *) SafeMalloc(BoardSize);
           }

       /* Коэффициенты */
       for (i=0; i<NumEqs; i++)
           {
            Kfs[i] = VDefs[i].NeedK!=0 ? (float *) SafeMalloc(BoardSize) : NULL;
            Sfs[i] = VDefs[i].NeedS!=0 ? (float *) SafeMalloc(BoardSize) : NULL;
           }
       for (i=0; i<NSubst; i++)
           {
            KDn[i] = NULL;
            SDn[i] = NULL;
           }
       for (i=0; i<NumPhases; i++)
           if (PhaseLinks[i]>=0)
              {
               KDn[PhaseLinks[i]] = (float *) SafeMalloc(BoardSize);
               SDn[PhaseLinks[i]] = (float *) SafeMalloc(BoardSize);
              }
       for (i=0; SubstRefs[i] && i<MaxActSubst; i++)
           {
            if (KDn[SLinks[i]]==NULL) KDn[SLinks[i]] = (float *) SafeMalloc(BoardSize);
            if (SDn[SLinks[i]]==NULL) SDn[SLinks[i]] = (float *) SafeMalloc(BoardSize);
           }
#ifdef __PARALLEL__
      }
#endif

#ifdef __PARALLEL__
 if (MasterFlag || CalculateAll)
    {
#endif
     UzSave = NSubst==0 ? NULL : (float *) SafeMalloc(BoardSize);

     HX  = (float *) SafeMalloc(NX*sizeof(float));
     HY  = (float *) SafeMalloc(NY*sizeof(float));
     HZ  = (float *) SafeMalloc(NZ*sizeof(float));
     _X  = (float *) SafeMalloc(NX*sizeof(float));
     _Y  = (float *) SafeMalloc(NY*sizeof(float));
     _Z  = (float *) SafeMalloc(NZ*sizeof(float));
     HXX = (float *) SafeMalloc(NX*sizeof(float));
     HYY = (float *) SafeMalloc(NY*sizeof(float));
     HZZ = (float *) SafeMalloc(NZ*sizeof(float));
     HXS = (float *) SafeMalloc(NX*sizeof(float));
     HYS = (float *) SafeMalloc(NY*sizeof(float));
     HZS = (float *) SafeMalloc(NZ*sizeof(float));
     HX2 = (float *) SafeMalloc(NX*sizeof(float));
     HY2 = (float *) SafeMalloc(NY*sizeof(float));
     HZ2 = (float *) SafeMalloc(NZ*sizeof(float));
     HXP = (float *) SafeMalloc(NX*sizeof(float));
     HYP = (float *) SafeMalloc(NY*sizeof(float));
     HZP = (float *) SafeMalloc(NZ*sizeof(float));
     HXR = (HKoeffs *) SafeMalloc(NX*sizeof(HKoeffs));
     HYR = (HKoeffs *) SafeMalloc(NY*sizeof(HKoeffs));
     HZR = (HKoeffs *) SafeMalloc(NZ*sizeof(HKoeffs));
#ifdef __PARALLEL__
    }
#endif

#ifdef __PARALLEL__
 if (!MasterFlag && CalculateAll)
    {
     H1 = (float *) SafeMalloc(nSMP*max(NZ,2*(NumEqs+NSubst))*NX*NY*sizeof(float));
#else
     H1 = (float *) SafeMalloc(nSMP*BoardSize);
#endif

     GRAV  = (float *) SafeMalloc(BoardSize);
     ALPHA = (float *) SafeMalloc(BoardSize);
     ERRS  = (float *) SafeMalloc(BoardSize);

     NN = BoardSize/sizeof(float);
     for (i = 0; i < NN; i++)
         ALPHA[i] = 0.5;
     memset(GRAV, 0, BoardSize);

     /* Коэффициенты прогонки */
     L  = (double *) SafeMalloc(nSMP*ScKfSize);
     M  = (double *) SafeMalloc(nSMP*ScKfSize);
     G  = (double *) SafeMalloc(nSMP*ScKfSize);
#ifdef __PARALLEL__
    }
#endif

 Lmin2    = (float *) SafeMalloc(BoardSize);
 DIV      = (float *) SafeMalloc(BoardSize);
 OrgnArea = (unsigned char *) SafeMalloc(NX*NY*NZ*sizeof(char));

 for (i=0; i<NSubst; i++)
     if (KGlobal.Uw[i]==0.0f)
        MapSubsts[NLightSubst++] = i;
 for (i=0, j=NLightSubst; i<NSubst; i++)
     if (KGlobal.Uw[i]!=0.0f)
        MapSubsts[j++] = i;

#ifdef __PARALLEL__
 DebugPrintf(DEBUG_FILE,"The main boards allocated\n");
 if (MasterFlag)
    {
#endif
     /* Читать информацию об области и граничных условиях */
     Area = ReadMap(0,-1,CfgFName,".map",NULL,NULL);

     memmove(OrgnArea,Area,NX*NY*NZ*sizeof(char));
     for (z=0, Ptr=0; z<NZ; z++)
       for (y=0; y<NY; y++)
         for (x=0; x<NX; x++,Ptr++)
             if (Area[Ptr])
                Area[Ptr] = AreaBlock;

     DscFile = fopen(DscName,"w+t");
     if (DscFile)
        fprintf(DscFile,"[Descriptor]\n");
     for (i=0; i<NumEqs; i++)
         {
          char MapExt[10] = ".";

          strcat(MapExt,VDefs[i].FileBase);
          if (VDefs[i].Projection>=0) strcat(MapExt,Axis[VDefs[i].Projection]);
          if (DscFile)
             if (VDefs[i].Projection<1)
                {
                 char NMBuf[90];
                 char ExtBase[30];

                 strcpy(ExtBase,VDefs[i].FileBase);
                 if (VDefs[i].Description[0])
                    strcpy(NMBuf,VDefs[i].Description);
                 else
                    strcpy(NMBuf,VDefs[i].FileBase);
                 if (VDefs[i].SubClass>=0)
                    sprintf(ExtBase,"%s%i",ExtBase,VDefs[i].SubClass);
                 for (x = i+1; x<NumEqs && VDefs[x].Projection>=1; x++);
                 if (Saves[Eq2File[i]])
                    fprintf(DscFile,"%s=%s %u\n",NMBuf,ExtBase,x-i);
                }
          Maps[i].Map = ReadMap(0,VDefs[i].SubClass,CfgFName,MapExt,
                                &Maps[i].Vals,&Maps[i].ValNum);
         }
     if (DscFile)
        {
         for (i=StandNSect; i<NSect; i++)
             fprintf(DscFile,"[%s]\n",Sections[i].Name);
         fclose(DscFile);
        }

     CMap = NSubst==0 ? NULL : ReadMap(NSubst,-1,CfgFName,".c",&CVals,&SubstValNum);

#ifndef __PARALLEL__
     NKinets = 0;
     if (CMap)
        for (Ptr = 0; Ptr < NX*NY*NZ; Ptr++)
       	    if (CMap[Ptr]==Zero)
               NKinets++;
     Traces = (TraceTypeHost *)malloc(NKinets*sizeof(TraceTypeHost));
     NKinets = 0;
     if (CMap)
        for (Ptr = 0; Ptr < NX*NY*NZ; Ptr++)
            if (CMap[Ptr]==Zero) {
                 int counter;

    		 Traces[NKinets].Tr = (unsigned short *) malloc(sizeof(TraceType));
    		 Traces[NKinets].Tr1 = (unsigned short *) malloc(sizeof(TraceType));
    		 Traces[NKinets].Tr2 = (unsigned short *) malloc(sizeof(TraceType));
    		 Traces[NKinets].Tr3 = (unsigned short *) malloc(sizeof(TraceType));
    		 memset(Traces[NKinets].Tr, 0, sizeof(TraceType));
    		 memset(Traces[NKinets].Tr1, 0, sizeof(TraceType));
    		 memset(Traces[NKinets].Tr2, 0, sizeof(TraceType));
    		 memset(Traces[NKinets].Tr3, 0, sizeof(TraceType));
    		 Traces[NKinets].KinetTraceQ1 = (float *) malloc(markZ99*sizeof(float));
    		 Traces[NKinets].KinetTraceW1 = (float *) malloc(markZ99*sizeof(float));
    		 Traces[NKinets].KinetTraceQ2 = (float *) malloc(markZ99*sizeof(float));
    		 Traces[NKinets].KinetTraceW2 = (float *) malloc(markZ99*sizeof(float));
    		 Traces[NKinets].KinetTraceQ3 = (float *) malloc(markZ99*sizeof(float));
    		 Traces[NKinets].KinetTraceW3 = (float *) malloc(markZ99*sizeof(float));
    		 Traces[NKinets].Prognosed0 = (unsigned short *) malloc(markZ99*sizeof(unsigned short));
    		 Traces[NKinets].Prognosed1 = (unsigned short *) malloc(markZ99*sizeof(unsigned short));
    		 Traces[NKinets].Prognosed2 = (unsigned short *) malloc(markZ99*sizeof(unsigned short));
                 for (counter = 0; counter < HistorySize; counter++) {
    		 	Traces[NKinets].Winners[counter] = (signed char *) malloc(markZ99*sizeof(signed char));
    		 	Traces[NKinets].Besters[counter] = (signed char *) malloc(markZ99*sizeof(signed char));
			Traces[NKinets].Timings[counter] = 1.0f;
			memset(Traces[NKinets].Winners[counter], 0, markZ99*sizeof(signed char));
			memset(Traces[NKinets].Besters[counter], 0, markZ99*sizeof(signed char));
		 }
                 for (counter = 0; counter < markZ99; counter++) {
                     Traces[NKinets].KinetTraceW1[counter] = 0.0;
                     Traces[NKinets].KinetTraceQ1[counter]  = 1.0;
                     Traces[NKinets].KinetTraceW2[counter] = 0.0;
                     Traces[NKinets].KinetTraceQ2[counter]  = 1.0;
                     Traces[NKinets].KinetTraceW3[counter] = 0.0;
                     Traces[NKinets].KinetTraceQ3[counter]  = 1.0;
                     Traces[NKinets].Prognosed0[counter] = 1;
                     Traces[NKinets].Prognosed1[counter] = 1;
                     Traces[NKinets].Prognosed2[counter] = 1;
                 }
                 Traces[NKinets].KinetTraceW = KGlobal.KinetTraceW;
                 Traces[NKinets].KinetTraceWOrder = KGlobal.KinetTraceWOrder;
                 Traces[NKinets].Iters = 1;
    		 Traces[NKinets++].tag = Ptr;
            }
#endif

#ifdef __PARALLEL__
 /* Работаем ТОЛЬКО в случае, если в верхней и нижней границах нет "щелей" */
     for (y=0, Ptr=0; y<NY; y++)
         for (x=0; x<NX; x++,Ptr++)
             if (Area[Ptr]==0 || Area[Ptr+NZs*NY*NX]==0)
                {
                 fprintf(StdOutput,"This program is not intended for processing such areas\n");
                 fprintf(StdOutput,"Area can not have a space at top/bottom\n");
#ifdef __MVS__
                 fclose(StdOutput);
#endif
                 AbortServer(4);
                }
 /* Посылаем КОЛЬЦУ процессоров информацию об области и граничных условиях */
     for (i=1; i<nProcs; i++)
         {
          SendSlave(i-1,(byte *) &Area[FromGrid[i-1]*NY*NX],LengthGrid[i-1]*NY*NX*sizeof(char));
          SendSlave(i-1,(byte *) &OrgnArea[FromGrid[i-1]*NY*NX],LengthGrid[i-1]*NY*NX*sizeof(char));
          DebugPrintf(DEBUG_FILE,"Area block sended to %i\n",i);
          for (x=0; x<NumEqs; x++)
              SendMap (i-1,FromGrid[i-1],LengthGrid[i-1],
                Maps[x].Map,Maps[x].ValNum,Maps[x].Vals,1);
          DebugPrintf(DEBUG_FILE,"Main maps sended to %i\n",i);
          if (NSubst)
             SendMap (i-1,FromGrid[i-1],LengthGrid[i-1],CMap,SubstValNum,CVals,NSubst);
         }
    }
 else
   if (CalculateAll)
      {
       /* Получаем информацию об области и граничных условиях и    */
       /* сразу корректируем карты, определяя, где необходим обмен */
       Area  = RecvMap(FirstFlag,LastFlag,NULL,NULL,NULL,1);
       RecvMaster((byte *) OrgnArea, NX*NY*NZ*sizeof(char));
       DebugPrintf(DEBUG_FILE,"Area block received\n");
       for (i=0; i<NumEqs; i++)
           Maps[i].Map = RecvMap(FirstFlag,LastFlag,Area,
                                 &Maps[i].ValNum,&Maps[i].Vals,1);
       DebugPrintf(DEBUG_FILE,"Main maps received\n");
       CMap  = NSubst==0 ? NULL : RecvMap(FirstFlag,LastFlag,Area,&SubstValNum,&CVals,NSubst);
       /* "Замыкаем" карту физической конфигурации сверху и снизу */
       if (!FirstFlag)
          for (y=0, Ptr=0; y<NY; y++)
              for (x=0; x<NX; x++,Ptr++)
                  if (Area[Ptr]==0)
                     {
                      DownExchange = 1;
                      Area[Ptr] = ExchngBit | 1;
                     }
       if (!LastFlag)
          for (y=0, Ptr=NZs*NY*NX; y<NY; y++)
              for (x=0; x<NX; x++,Ptr++)
                  if (Area[Ptr]==0)
                     {
                      UpExchange = 1;
                      Area[Ptr] = ExchngBit | 1;
                     }
      }
#endif

 Vars.BoundSize = 0;
 for (i=_back; i<=_top; i++)
     Vars.Bounds[i] = NULL;

#ifdef __PARALLEL__
 if (CalculateAll || MasterFlag)
    {
#endif
     Boundaries = (unsigned char *) SafeMalloc(BoardSize);
     memset(Boundaries,0,BoardSize);
     DescH = CreateDescriptor((unsigned char *) Area,(unsigned char *) OrgnArea,NZ,NY,NX,NY*NX,NX,1,LeftBit,RightBit,Boundaries);
     DescP = CreateDescriptor((unsigned char *) Area,(unsigned char *) OrgnArea,NZ,NX,NY,NY*NX,1,NX,BackBit,ForwBit,Boundaries);
     DescV = CreateDescriptor((unsigned char *) Area,(unsigned char *) OrgnArea,NY,NX,NZ,NX,1,NY*NX,BottomBit,TopBit,Boundaries);
#ifdef __PARALLEL__
    }
 if (!MasterFlag && CalculateAll)
    {
     for (z=3, Ptr=3*NY*NX; z<NZ-3; z++)
#else
     for (z=0, Ptr=0; z<NZ; z++)
#endif
       for (y=0; y<NY; y++)
         for (x=0; x<NX; x++,Ptr++)
           if (Area[Ptr] && !IsExchng && ((IsBack!=0)+(IsForw!=0)+(IsRight!=0)+(IsLeft!=0)+(IsBottom!=0)+(IsTop!=0))>1)
              Vars.BoundSize++;
     Vars.BoundSize *= (NumEqs+NSubst)*sizeof(float);
     for (i=_back; i<=_top; i++)
         {
          Vars.Bounds[i] = (float *) SafeMalloc(Vars.BoundSize);
          memset(Vars.Bounds[i],0,Vars.BoundSize);
         }

     for (i=0; i<NumEqs; i++)
         {
          Templ.Name[i] = CheckTau==0 ? NULL : (float *) SafeMalloc(TempBoardSize);
          Old.Name[i]   = CheckTau==0 ? NULL : (float *) SafeMalloc(TempBoardSize);
         }
     Templ.Dn = NSubst==0 ? NULL : (float **) SafeMalloc(NSubst*sizeof(float *));
     Templ.BoundSize = Vars.BoundSize;
     Old.Dn   = NSubst==0 ? NULL : (float **) SafeMalloc(NSubst*sizeof(float *));
     Old.BoundSize = Vars.BoundSize;
     for (i=_back; i<=_top; i++)
         {
          Templ.Bounds[i] = CheckTau==0 ? NULL : (float *) SafeMalloc(Templ.BoundSize);
          Old.Bounds[i] = CheckTau==0 ? NULL : (float *) SafeMalloc(Old.BoundSize);
         }

     for (x=0;x<NSubst;x++)
         {
          Templ.Dn[x] = CheckTau==0 ? NULL : (float *) SafeMalloc(TempBoardSize);
          Old.Dn[x]   = CheckTau==0 ? NULL : (float *) SafeMalloc(TempBoardSize);
         }
     for (i = 0; i < 3; i++) {
         Templ.Parts[i]  = CheckTau==0 ? NULL : (float *) SafeMalloc(TotalNParts*sizeof(float));
         Templ.PartsV[i] = CheckTau==0 ? NULL : (float *) SafeMalloc(TotalNParts*sizeof(float));
         Old.Parts[i]  = CheckTau==0 ? NULL : (float *) SafeMalloc(TotalNParts*sizeof(float));
         Old.PartsV[i] = CheckTau==0 ? NULL : (float *) SafeMalloc(TotalNParts*sizeof(float));
     }
#ifdef __PARALLEL__
    }
 DebugPrintf(DEBUG_FILE,"Prepare to the last stage\n");
 if (MasterFlag)
    {
#endif
     /* Инициализировать массивы данных о шагах сетки */
     HXreg = ArrangeGrid(OrgnArea,HX,(float)HXmin,(float)HXmax,NZ,NY,NX,NY*NX,NX,1,DescH,_HX[0]>0.0 ? _HX : NULL);
     HYreg = ArrangeGrid(OrgnArea,HY,(float)HYmin,(float)HYmax,NZ,NX,NY,NY*NX,1,NX,DescP,_HY[0]>0.0 ? _HY : NULL);
     HZreg = ArrangeGrid(OrgnArea,HZ,(float)HZmin,(float)HZmax,NY,NX,NZ,NX,1,NY*NX,DescV,_HZ[0]>0.0 ? _HZ : NULL);
     for (z=0,_Z[0]=0.0;z<NZs;z++)
         _Z[z+1] = _Z[z]+HZ[z];
     for (z=0;z<NZ;z++)
         {
          char S = HZ[z]<HZ[ZM];
          HZ2[z] = HZ[z]*HZ[z];
          HZZ[z] = HZ[z]*HZ[ZM];
          HZS[z] = HZ[z]+HZ[ZM];
          HZP[z] = HZZ[z]*HZS[z];
          HZR[z].h   = S ? HZ[z] : HZ[ZM];
          HZR[z].R1m = S ? (HZ[z]-HZ[ZM])/HZS[z] : 0.0f;
          HZR[z].R2m = S ? 2.0f*(1.0f-HZ[z]/HZ[ZM]) : 0.0f;
          HZR[z].R3m = S ? 2.0f*HZ2[z]/HZ[ZM]/HZS[z] : 1.0f;
          HZR[z].R1p = S ? 1.0f : 2.0f*HZ[ZM]*HZ[ZM]/HZ[z]/HZS[z];
          HZR[z].R2p = S ? 0.0f : 2.0f*(1.0f-HZ[ZM]/HZ[z]);
          HZR[z].R3p = S ? 0.0f : (HZ[ZM]-HZ[z])/HZS[z];
         }

     ChangeExt(Lmin2FName,CfgFName,Lmin2FExt);
     if ((Lmin2File = fopen(Lmin2FName,"rb")) && fseek(Lmin2File,0,SEEK_END)==0 && ftell(Lmin2File)==BoardSize)
        {
         fseek(Lmin2File,0,SEEK_SET);
         fread(Lmin2,1,BoardSize,Lmin2File);
         fprintf(StdOutput,"Lmin2 Table is loaded from cache\n");
         fflush(StdOutput);
        }
     else
        {
         Lmin2File = fopen(Lmin2FName,"wb");

         if (PhaseVars[CarrierPhase]._Nu>=0)
            CalculateLmin2(OrgnArea,HZ,HY,HX,Lmin2);
         else
            memset(Lmin2,0,BoardSize);
         fwrite(Lmin2,1,BoardSize,Lmin2File);
        }
     fclose(Lmin2File);
#ifdef __PARALLEL__
     for (i=1; i<nProcs; i++)
         {
          for (z=0; z<NVar; z++)
              if (Variables[z].VarBcast)
                 {
                  SendSlave(i-1,(byte *) &Variables[z].VarMaxN,sizeof(Variables[z].VarMaxN));
                  if (Variables[z].VarMaxN)
                     SendSlave(i-1,(byte *) Variables[z].VarMap,
                        Variables[z].VarLong*Variables[z].VarMaxN*cfgTypeSizes[Variables[z].VarType]);
                 }
          SendSlave(i-1,(byte *) HX,NX*sizeof(float));
          SendSlave(i-1,(byte *) HY,NY*sizeof(float));
          SendSlave(i-1,(byte *) &HZ[FromGrid[i-1]], LengthGrid[i-1]*sizeof(float));
          SendSlave(i-1,(byte *) &HZ[(FromGrid[i-1]+NZs) % NZ], sizeof(float));
          SendSlave(i-1,(byte *) &_Z[FromGrid[i-1]], LengthGrid[i-1]*sizeof(float));
          SendSlave(i-1,(byte *) &HZ2[FromGrid[i-1]],LengthGrid[i-1]*sizeof(float));
          SendSlave(i-1,(byte *) &HZ2[(FromGrid[i-1]+NZs) % NZ], sizeof(float));
          SendSlave(i-1,(byte *) &HZZ[FromGrid[i-1]],LengthGrid[i-1]*sizeof(float));
          SendSlave(i-1,(byte *) &HZS[FromGrid[i-1]],LengthGrid[i-1]*sizeof(float));
          SendSlave(i-1,(byte *) &HZP[FromGrid[i-1]],LengthGrid[i-1]*sizeof(float));
          SendSlave(i-1,(byte *) &HZR[FromGrid[i-1]],LengthGrid[i-1]*sizeof(HKoeffs));
          SendSlave(i-1,(byte *) &Lmin2[FromGrid[i-1]*NY*NX],LengthGrid[i-1]*NY*NX*sizeof(float));
          SendSlave(i-1,(byte *) &HXreg,sizeof(HXreg));
          SendSlave(i-1,(byte *) &HYreg,sizeof(HYreg));
          SendSlave(i-1,(byte *) &HZreg,sizeof(HZreg));
          DebugPrintf(DEBUG_FILE,"Main experiment data sended to %i\n",i);
         }
    }
 else
   {
    if (CalculateAll)
       {
        for (z=0; z<NVar; z++)
            if (Variables[z].VarBcast)
               {
                RecvMaster((byte *) &Variables[z].VarMaxN,sizeof(Variables[z].VarMaxN));
                if (Variables[z].VarMaxN)
                   RecvMaster((byte *) Variables[z].VarMap,
                      Variables[z].VarLong*Variables[z].VarMaxN*cfgTypeSizes[Variables[z].VarType]);
               }
        RecvMaster((byte *) HX,    NX*sizeof(float));
        RecvMaster((byte *) HY,    NY*sizeof(float));
        RecvMaster((byte *) HZ,    NZ*sizeof(float));
        RecvMaster((byte *) &HZbf, sizeof(float));
        RecvMaster((byte *) _Z,    NZ*sizeof(float));
        RecvMaster((byte *) HZ2,   NZ*sizeof(float));
        RecvMaster((byte *) &HZ2bf,sizeof(float));
        RecvMaster((byte *) HZZ,   NZ*sizeof(float));
        RecvMaster((byte *) HZS,   NZ*sizeof(float));
        RecvMaster((byte *) HZP,   NZ*sizeof(float));
        RecvMaster((byte *) HZR,   NZ*sizeof(HKoeffs));
        RecvMaster((byte *) Lmin2, BoardSize);
        RecvMaster((byte *) &HXreg, sizeof(HXreg));
        RecvMaster((byte *) &HYreg, sizeof(HYreg));
        RecvMaster((byte *) &HZreg, sizeof(HZreg));
        DebugPrintf(DEBUG_FILE,"Main experiment data received\n");

        SHMEMUp   = SHMEMUp && UseSHMEM;
        SHMEMDown = SHMEMDown && UseSHMEM;
       }
   }

 if (MasterFlag)
    BuffNodes = (NSubst+NumEqs)*(5+MaxDivision)*NY*NX;
 else {
    BuffNodes = (NSubst+NumEqs)*(2+NZ)*NY*NX;
    if (CalculateAll)
       for (i = 0; i < 3; i++) {
           Vars.Parts[i]  = CheckTau == 0 ? NULL : (float *) SafeMalloc(TotalNParts*sizeof(float));
           Vars.PartsV[i] = CheckTau == 0 ? NULL : (float *) SafeMalloc(TotalNParts*sizeof(float));
       }
 }
 BuffNodes = 1 + BuffNodes/(NSubst+1);
 HBuff = (float *) SafeMalloc((long) BuffNodes*(NSubst+1)*sizeof(float));
 eHBuff = CheckTau == 0 ? NULL : (float *) SafeMalloc((long) BuffNodes*(NSubst+1)*sizeof(float));
#endif

 free(OrgnArea);

#ifdef __PARALLEL__
 if (!MasterFlag && CalculateAll) {
#endif
    if (CheckTau) {
       for (i = 0; i < 3; i++)
           memset(Vars.PartsV[i], 0, TotalNParts*sizeof(float));
       for (i = 0; i < NParts; i++) {
           Vars.Parts[0][i] = 1.5f + 1.0f*(NX-3)*rand()/RAND_MAX;
           Vars.Parts[1][i] = 1.5f + 1.0f*(NY-3)*rand()/RAND_MAX;
           Vars.Parts[2][i] = 1.5f + 1.0f*(BaseZ + (NZ-3)*rand()/RAND_MAX);
       }
    }
#ifdef __PARALLEL__
 }

 if (MasterFlag || CalculateAll)
    {
#endif
     for (x=0,_X[0]=0.0;x<NXs;x++)
         _X[x+1] = _X[x]+HX[x];
     for (x=0;x<NX;x++)
         {
          char S = HX[x]<HX[XM];
          HX2[x] = HX[x]*HX[x];
          HXX[x] = HX[x]*HX[XM];
          HXS[x] = HX[x]+HX[XM];
          HXP[x] = HXX[x]*HXS[x];
          HXR[x].h   = S ? HX[x] : HX[XM];
          HXR[x].R1m = S ? (HX[x]-HX[XM])/HXS[x] : 0.0f;
          HXR[x].R2m = S ? 2.0f*(1.0f-HX[x]/HX[XM]) : 0.0f;
          HXR[x].R3m = S ? 2.0f*HX2[x]/HX[XM]/HXS[x] : 1.0f;
          HXR[x].R1p = S ? 1.0f : 2.0f*HX[XM]*HX[XM]/HX[x]/HXS[x];
          HXR[x].R2p = S ? 0.0f : 2.0f*(1.0f-HX[XM]/HX[x]);
          HXR[x].R3p = S ? 0.0f : (HX[XM]-HX[x])/HXS[x];
         }
     for (y=0,_Y[0]=0.0;y<NYs;y++)
         _Y[y+1] = _Y[y]+HY[y];
     for (y=0;y<NY;y++)
         {
          char S = HY[y]<HY[YM];
          HY2[y] = HY[y]*HY[y];
          HYY[y] = HY[y]*HY[YM];
          HYS[y] = HY[y]+HY[YM];
          HYP[y] = HYY[y]*HYS[y];
          HYR[y].h   = S ? HY[y] : HY[YM];
          HYR[y].R1m = S ? (HY[y]-HY[YM])/HYS[y] : 0.0f;
          HYR[y].R2m = S ? 2.0f*(1.0f-HY[y]/HY[YM]) : 0.0f;
          HYR[y].R3m = S ? 2.0f*HY2[y]/HY[YM]/HYS[y] : 1.0f;
          HYR[y].R1p = S ? 1.0f : 2.0f*HY[YM]*HY[YM]/HY[y]/HYS[y];
          HYR[y].R2p = S ? 0.0f : 2.0f*(1.0f-HY[YM]/HY[y]);
          HYR[y].R3p = S ? 0.0f : (HY[YM]-HY[y])/HYS[y];
         }
#ifdef __PARALLEL__
      }
#endif

#ifdef __PARALLEL__
 if (MasterFlag)
    {
#endif
     if (SavFile = fopen(SavName,"r+b"))
        {
         int RestFile;

         fread(&maxtakt,sizeof(maxtakt),1,SavFile);
         fread(&Quant,sizeof(Quant),1,SavFile);
         fread(&takt,sizeof(takt),1,SavFile);
         for (i=0; i<NumEqs; i++)
             fread(Vars.Name[i],1,BoardSize,SavFile);
         for (x=0;x<NSubst;x++)
             fread(Vars.Dn[x],1,BoardSize,SavFile);
         fprintf(StdOutput,"Experiment (%lu : %lu) continues\n",(long) maxtakt,(long) Quant);
         fflush(StdOutput);

#ifdef __PARALLEL__
         RestFile = ftell(SavFile);
         fseek(SavFile,0,SEEK_END);
         RestFile = ftell(SavFile)-RestFile;
         fseek(SavFile,-RestFile,SEEK_END);
         for (i=1; i<nProcs; i++)
             {
              if (RestFile)
                 {
                  int ToSend;
                  fread(&ToSend,sizeof(int),1,SavFile);
                  SendSlave(i-1,(byte *) &ToSend,sizeof(int));
                  if (ToSend)
                     for (j=_back; j<=_top; j++)
                         {
                          fread(HBuff,1,ToSend,SavFile);
                          SendSlave(i-1,(byte *) HBuff,ToSend);
                         }
                 }
              else
                 SendSlave(i-1,(byte *) &RestFile,sizeof(int));
             }
#else
         RestFile = ftell(SavFile);
         fseek(SavFile,0,SEEK_END);
         RestFile = ftell(SavFile)-sizeof(ModelTime)-sizeof(EndTime)-RestFile;
         if (RestFile!=(_top-_back+1)*Vars.BoundSize)
            {
             fprintf(StdOutput,"Can't read special boundary values from experiment file\n");
#ifdef __MVS__
             fclose(StdOutput);
#endif
             AbortServer(-30);
            }
         fseek(SavFile,-RestFile-sizeof(ModelTime)-sizeof(EndTime),SEEK_END);
         if (Vars.BoundSize)
            {
             for (i=_back; i<=_top; i++)
                 fread(Vars.Bounds[i],1,Vars.BoundSize,SavFile);
             SpecialBoundsNum = Vars.BoundSize/sizeof(float);
             if (!UseSpecialBounds)
                for (i=0; i<NumEqs+NSubst; i++)
                 {
                  float *  Board = i<NumEqs ? Vars.Name[i] : Vars.Dn[i-NumEqs];
                  float ** Bounds = Vars.Bounds;

                  for (Ptr=0; Ptr<BoardSize; Ptr++)
                      if (fabs(Board[Ptr])>=SpecVAL)
                         {
                          int   Count = (IsBack!=0)+(IsForw!=0)+(IsRight!=0)+(IsLeft!=0)+(IsBottom!=0)+(IsTop!=0);
                          float Val = 0.0;

                          if (IsForw)   Val+=_fw(Board[Ptr]);
                          if (IsBack)   Val+=_bw(Board[Ptr]);
                          if (IsBottom) Val+=_bt(Board[Ptr]);
                          if (IsTop)    Val+=_tp(Board[Ptr]);
                          if (IsLeft)   Val+=_lf(Board[Ptr]);
                          if (IsRight)  Val+=_rg(Board[Ptr]);

                          Board[Ptr] = Val/Count;
                         }
                 }
            }
#endif
         Cadres     = takt / Quant;
         QuantCount = takt % Quant;
         if (!fread(&ModelTime,sizeof(ModelTime),1,SavFile)) ModelTime = 0.0;
         if (!fread(&EndTime,sizeof(EndTime),1,SavFile)) EndTime = 0.0;
         fprintf(StdOutput,"Time = %lf\n",ModelTime);
         fflush(StdOutput);
        }
     else
        {
#ifndef __MVS__
         fprintf(StdOutput,"Input number of iterations : ");
#endif
         if (StdInput==stdin)
	    i = scanf(_intScanf,&maxtakt);
	 else
            i = fscanf(StdInput,_intScanf,&maxtakt);

         if (i < 1) {
            fprintf(StdOutput,"no input data\n");
#ifdef __MVS__
            fclose(StdOutput);
#endif
            AbortServer(-1100);
         }
         if (maxtakt == 0) {
#ifndef __MVS__
            fprintf(StdOutput,"Input end time : ");
#endif
            if (StdInput==stdin)
	       i = scanf(_doubleScanf,&EndTime);
	    else
               i = fscanf(StdInput,_doubleScanf,&EndTime);

            if (i < 1) {
               fprintf(StdOutput,"no input data\n");
#ifdef __MVS__
               fclose(StdOutput);
#endif
               AbortServer(-1100);
            }
         }
#ifndef __MVS__
         fprintf(StdOutput,"Input number of iterations in 1 cadre : ");
#endif
         if (StdInput==stdin)
	    i = scanf("%u",&Quant);
	 else
            i = fscanf(StdInput,"%u",&Quant);

         if (i < 1) {
            fprintf(StdOutput,"no input data\n");
#ifdef __MVS__
            fclose(StdOutput);
#endif
            AbortServer(-1200);
         }
         for (i=0; i<NumEqs; i++)
             InitBoardF(VDefs[i].SubClass,Vars.Name[i], Maps[i].Map, Maps[i].Vals, VDefs[i]._Zero);
         /* Установление кинетики, если нужно */ /* !!! */
         if (NReact && SetSteadyKinetics) {
            KineticContext CC = { 0 };
            float TimeQuant = TAU;
            float avr_dc = 0.0f;

            if (PhaseVars[CarrierPhase]._T>=0) {
               bool Flag = 0;
               for (z=0, Ptr=0; !Flag && z<NZ; z++)
                 for (y=0; !Flag && y<NY; y++)
                   for (x=0; !Flag && x<NX; x++,Ptr++)
                       if (!Area[Ptr]) {
                          SetTk(&CC,Vars.Name[PhaseVars[CarrierPhase]._T][Ptr]);
                          Flag = 1;
                       }
            } else
               SetTk(&CC,(float)KGlobal.DefaultT);
            for (i=0; i<NASubst; i++)
               CC.Conc0[i] = ZeroC[KGlobal.TranMap[i]];

            do {
                KGlobal.InitH = min(KGlobal.InitH,TimeQuant*0.001f);
                CC._EndTime = (float)TimeQuant;

                OneTaktKinetic(&KGlobal,&CC,UseGear,Adams_Rozhkov_Method);

                avr_dc = 1000.0f;

                if (CC.Iters >= KGlobal.MaxIterations)
                   TimeQuant /= 5.0;
                else if (CC.LossPrecision)
                   break;
                else {
                   avr_dc = 0.0f;
                   for (i=0; i<NASubst; i++) {
                       avr_dc += fabsf(CC.Conc1[i] - CC.Conc0[i]);
                       CC.Conc0[i] = CC.Conc1[i];
                   }
                   avr_dc /= NASubst;
                }
            } while (TimeQuant > 1E-10 && avr_dc > 1E-6*TimeQuant);
            if (TimeQuant <= 1E-10) {
               printf("Warning: Kinetics can't converge!\n");
               fprintf(StdOutput,"Warning: Kinetics can't converge!\n");
            }
            for (i=0; i<NASubst; i++)
               ZeroC[KGlobal.TranMap[i]] = CC.Conc0[i];
            fprintf(StdOutput,"Steady Concentrations = ");
            for (i = 0; i < NSubst; i++) {
                fprintf(StdOutput,"%le ", ZeroC[i]);
            }
            fprintf(StdOutput,"\n");
            fflush(StdOutput);
         }
         for (x=0;x<NSubst;x++)
             InitBoard(Vars.Dn[x],CMap,&(CVals[x*SubstValNum]),(float)ZeroC[x]);
#ifndef __PARALLEL__
         if (UseSpecialBounds) {
            for (i=0, j=0; i<NumEqs+NSubst; i++)
              for (z=0, Ptr=0; z<NZ; z++)
                for (y=0; y<NY; y++)
                  for (x=0; x<NX; x++,Ptr++)
                    if (Area[Ptr] && !IsExchng && ((IsBack!=0)+(IsForw!=0)+(IsRight!=0)+(IsLeft!=0)+(IsBottom!=0)+(IsTop!=0))>1)
                       {
                        int CODE = j+1; // 128 Eqs * 65536 points; Min = 1; Max = 8388480+1
                        int FIRST = CODE/1000000; // 0..9
                        int REM = CODE%1000000; // Min = 1; Max = 999999
                        // 500000..999999 => 0.500000..0.999999
                        // 1..499999 => -[0.500001..0.999999]
                        float Mantissa = REM>=500000 ? 0.000001*REM : -0.000001*(REM+500000);
                        int   POW = 100+FIRST; // Max = 109
                        // ABS[max] = max(Mantissa)*2^max(POW) = 6.4903645828*1E32
                        // ABS[min] = ABS[min(Mantissa)*2^min(POW)] = 6.3382530011*1E29
                        // Critical Number = +-1E25, for example...
                        float Number = ldexp(Mantissa,POW); // rounding must present!!!!!
                        int k;

                        for (k=_back; k<=_top; k++)
                            if (i<NumEqs)
                               Vars.Bounds[k][j] = Vars.Name[i][Ptr];
                            else
                               Vars.Bounds[k][j] = Vars.Dn[i-NumEqs][Ptr];
                        if (i<NumEqs)
                           Vars.Name[i][Ptr] = Number;
                        else
                           Vars.Dn[i-NumEqs][Ptr] = Number;
                        j++;
                       }
            SpecialBoundsNum = j;
         }
#else
         for (i=1; i<nProcs; i++)
             {
              int ToSend = 0;
              SendSlave(i-1,(byte *) &ToSend,sizeof(int));
             }
#endif
         QuantCount = 0;
         ModelTime  = 0.0;
        }

     if (UseSpecialBounds) {
        char Buf[__MAXPATH];

        sprintf(Buf,"%sbnd",_ResBase);

        if (SavFile)
           SaveBnds = fopen(Buf, "r+b");
        else
           SaveBnds = fopen(Buf, "w+b");
        if (SaveBnds) {
            fseek(SaveBnds,0L,SEEK_SET);
            fwrite(&SpecialBoundsNum,sizeof(SpecialBoundsNum),1,SaveBnds);
            fseek(SaveBnds,0L,SEEK_END);
        }
     }

     for (i=0;i<NumEqs;i++)
         if (VDefs[i].Projection<1 && Saves[Eq2File[i]])
            {
             char Buf[__MAXPATH];

             if (VDefs[i].SubClass<0)
                sprintf(Buf,"%s%s",_ResBase,VDefs[i].FileBase);
             else
                sprintf(Buf,"%s%s%i",_ResBase,VDefs[i].FileBase,VDefs[i].SubClass);
             SaveFiles[Eq2File[i]] = OpenBoardFile(SavFile,&Saves[Eq2File[i]],Buf,Cadres);
            }

     for (i=0;i<NSubst;i++)
         if (SaveDn[i])
            {
             char Buf[__MAXPATH];

             sprintf(Buf,"%s%i",DResBase,i);
             DResults[i] = OpenBoardFile(SavFile,&SaveDn[i],Buf,Cadres);
            }
#ifdef __MVS__
     if (StdInput) fclose(StdInput);
#endif

     while (TimeCounter<TimeN && TimeIters[TimeCounter] && takt>=TimeIters[TimeCounter])
       {
        if (TimeTAUk[TimeCounter]) TAU = TAU0*TimeTAUk[TimeCounter];
        TimeCounter++;
       }
     if (TimeCounter>0) CalcBase = TimeStopBase[TimeCounter-1];
     /* Цикл обработки */
     FTIME(&_BeginTime);
     START = 1;
     do
       {
        int StoppedByUser = 0;
        FILE * Stop = NULL;

        if (CheckTau && ((takt+1) % Period)==0)
#ifdef __PARALLEL__
           OneCheckedTakt(START, (QuantCount+1)==Quant);
        else
           OneTakt(1, START, (QuantCount+1)==Quant, CHK_NoCheck, NULL);
#else
           OneCheckedTakt();
        else
           OneTakt();
#endif

        ModelTime += TAU;
        if (TimeCounter<TimeN && TimeIters[TimeCounter])
           if ((takt+1)==TimeIters[TimeCounter])
              {
               CalcBase = TimeStopBase[TimeCounter];
               if (TimeTAUk[TimeCounter]) TAU = TAU0*TimeTAUk[TimeCounter];
               TimeCounter++;
              }
        if (fmod(ModelTime,60.0)<2.0*TAU)
           GetSolarAngle(ModelTime/3600.0);
        
        if ((++takt % 100) == 0)
#ifdef __PARALLEL__
           fprintf(StdOutput,"[%i]\n",takt);
#else
           fprintf(StdOutput,"Iteration number %i passed\n",takt);
#endif
        fflush(StdOutput);

        Stop = fopen(StopFileName,"rb");
        if (StoppedByUser = Stop!=NULL)
           {
            fclose(Stop);
            unlink(StopFileName);
            fprintf(StdOutput,"Stop - file detected : Experiment halted\n");
            fflush(StdOutput);
           }

        QuantCount++;
        if (QuantCount==Quant)
           {
            QuantCount = 0;
            Cadres++;
            /* Записать данные в файл результатов */
            for (i=0; i<NumEqs; i++)
                if (Saves[Eq2File[i]])
                   WriteBoardFile(Vars.Name[i],SaveFiles[Eq2File[i]],Cadres);
            for (i=0;i<NSubst;i++)
                if (SaveDn[i])
                   WriteBoardFile(Vars.Dn[i],DResults[i],Cadres);
            if (UseSpecialBounds && SaveBnds) {
               fseek(SaveBnds,0L,SEEK_SET);
               fwrite(&SpecialBoundsNum,sizeof(SpecialBoundsNum),1,SaveBnds);
               fseek(SaveBnds,0L,SEEK_END);
               for (i=_back; i<=_top; i++)
                   fwrite(Vars.Bounds[i], sizeof(float), SpecialBoundsNum, SaveBnds);
               fflush(SaveBnds);
            }
            CreateSavFile(&SavFile,maxtakt,Quant,takt,SaveAllHistory);
           }

        /* Закончить по достижении максимального числа итераций или если  */
        /* была нажата любая клавиша                                      */
#if !defined(__MAY_BE_MICROSOFTC__) && !defined(__MVS__) && !defined(__UNIX__)
        Command= (takt != maxtakt) && !(maxtakt == 0 && ModelTime >= EndTime) && (! kbhit()) && !StoppedByUser;
#else
        Command= (takt != maxtakt) && !(maxtakt == 0 && ModelTime >= EndTime) && !StoppedByUser;
#endif
#ifdef __PARALLEL__
        SendCommand(Command);
#endif
        START = 0;
       }
     while (Command);
     FTIME(&_EndTime);
     fprintf(StdOutput,"Number of iterations = %i.\n", takt);
     fprintf(StdOutput,"Model time = %lf sec.\n", ModelTime);
#ifndef __PARALLEL__
     fprintf(StdOutput,"Elapsed Kinetic time = %lf sec.\n", KinetTime);
#endif
     fprintf(StdOutput,"Elapsed time = %lf sec.\n", DIFFTIME(_BeginTime,_EndTime));
     CreateSavFile(&SavFile,maxtakt,Quant,takt,1);
#ifdef __PARALLEL__
    }
 else
    Slave();
#endif

#ifdef __PARALLEL__
 if (MasterFlag)
    {
#endif
     for (i=0; i<NumSaves; i++)
         if (Saves[i]) CloseBoardFile(SaveFiles[i], Cadres);

     for (i=0;i<NSubst;i++)
         if (SaveDn[i])
            CloseBoardFile(DResults[i], Cadres);

     if (UseSpecialBounds && SaveBnds)
        fclose(SaveBnds);

     for (i=0; i<NumEqs; i++) {
         free(Vars.Name[i]);
         if (CheckTau) free(Vars.eName[i]);
     }
     if (NSubst)
        {
         for (i=0; i<NSubst; i++) {
             free(Vars.Dn[i]);
             if (CheckTau) free(Vars.eDn[i]);
         }
         free(Vars.Dn);
         free(Vars.eDn);
        }

#ifdef __PARALLEL__
     free(FromGrid);
     free(LengthGrid);
     free(ProcessorNames);
    }
 free(ProcOnNode);

 free(HBuff);
 if (CheckTau) free(eHBuff);
#endif

#ifdef __PARALLEL__
 if (MasterFlag || CalculateAll)
    {
#endif
     if (CheckTau)
        for (i = 0; i < 3; i++) {
            free(Vars.Parts[i]);
            free(Vars.PartsV[i]);
        }
#ifdef __PARALLEL__
    }
#endif

 if (Vars.BoundSize) {
    for (i=_back; i<=_top; i++)
        free(Vars.Bounds[i]);
 }
 free(OldRefAddrs);
 free(OldRefVals);
 
#ifdef __PARALLEL__
 if (!MasterFlag && CalculateAll)
    {
#endif
     if (NSubst)
        {
         if (CheckTau) {
            for (i=0; i<NSubst; i++)
                {
                 free(Templ.Dn[i]);
                 free(Old.Dn[i]);
                }
            for (i = 0; i < 3; i++) {
                free(Templ.Parts[i]); free(Templ.PartsV[i]);
                free(Old.Parts[i]); free(Old.PartsV[i]);
            }
         }
         free(Templ.Dn);
         free(Old.Dn);
        }
     if (CheckTau)
        {
         for (i=0; i<NumEqs; i++)
             {
              free(Templ.Name[i]);
              free(Old.Name[i]);
             }
         for (i=_back; i<=_top; i++)
             {
              free(Templ.Bounds[i]);
              free(Old.Bounds[i]);
             }
        }

     free(H1);  free(L);   free(M);   free(G);

     free(ERRS);
     free(GRAV);
     free(ALPHA);

     for (i=0; i<NumPhases; i++)
         {
          if (!PhaseVars[i].IsLight)
             {
              free(WXYZ[i].WX1); free(WXYZ[i].WY1);
              free(WXYZ[i].WX2); free(WXYZ[i].WY2);
             }
          free(WXYZ[i].WZ1); free(WXYZ[i].WZ2);
         }
     for (i=0; i<NumEqs; i++)
         {
          free(Kfs[i]);
          free(Sfs[i]);
         }
     for (i=0; i<NSubst; i++)
         {
          free(KDn[i]);
          free(SDn[i]);
         }

     free(Dt);
#ifdef __PARALLEL__
    }
#endif

 free(Lmin2);
 free(DIV);

#ifndef __PARALLEL__
 for (Ptr = 0; Ptr < NKinets; Ptr++) {
	 int counter;

	 free(Traces[Ptr].Tr);
	 free(Traces[Ptr].Tr1);
	 free(Traces[Ptr].Tr2);
	 free(Traces[Ptr].Tr3);
	 free(Traces[Ptr].KinetTraceQ1);
	 free(Traces[Ptr].KinetTraceW1);
	 free(Traces[Ptr].KinetTraceQ2);
	 free(Traces[Ptr].KinetTraceW2);
	 free(Traces[Ptr].KinetTraceQ3);
	 free(Traces[Ptr].KinetTraceW3);
	 free(Traces[Ptr].Prognosed0);
	 free(Traces[Ptr].Prognosed1);
	 free(Traces[Ptr].Prognosed2);
         for (counter = 0; counter < HistorySize; counter++) {
		free(Traces[Ptr].Winners[counter]);
		free(Traces[Ptr].Besters[counter]);
	 }
 }
 free(Traces);
#endif

#ifdef __PARALLEL__
 if (MasterFlag || CalculateAll)
    {
#endif
     for (i=0; i<NZ*NY; i++)
         if (DescH[i].Items)
            free(DescH[i].Items);
     free(DescH);
     for (i=0; i<NZ*NX; i++)
         if (DescP[i].Items)
            free(DescP[i].Items);
     free(DescP);
     for (i=0; i<NY*NX; i++)
         if (DescV[i].Items)
            free(DescV[i].Items);
     free(DescV);
     free(Boundaries);

     if (NSubst)
        free(UzSave);
     free(Area);

     for (i=0; i<NumEqs; i++)
         {
          free(Maps[i].Map);
          if (Maps[i].Vals) free(Maps[i].Vals);
         }
     if (NSubst)
        {
         free(CMap);
         if (CVals) free(CVals);
        }

     free(HX);  free(HY);  free(HZ);
     free(_X);  free(_Y);  free(_Z);
     free(HXX); free(HYY); free(HZZ);
     free(HXS); free(HYS); free(HZS);
     free(HX2); free(HY2); free(HZ2);
     free(HXP); free(HYP); free(HZP);
     free(HXR); free(HYR); free(HZR);
#ifdef __PARALLEL__
    }
#endif

#ifdef __MPI__
 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Finalize();
#elif defined(__ROUTER__)
 for (i=0; i<TotalProcs; i++)
     for (x=0; x<NumSends[i]; x++)
         w_write(HndSends[i*DefMaxSends+x].Handle);
 free(NumSends); free(MaxSends); free(HndSends);
 for (i=0; i<NumListItems; i++)
     free(BufList[i].Buffer);
 free(BufList);
#endif

#ifdef __PARALLEL__
 DebugPrintf(DEBUG_FILE,"Finished\n");
 if (DebugFile) fclose(DebugFile);
#endif

#ifdef __MVS__
 fclose(StdOutput);
#endif

 free(Contexts);

 DoneVars();

 return 1;
}

#else

#error Must be a DOS/LARGE or WIN32 or MVS or UNIX model to compile this

#endif
