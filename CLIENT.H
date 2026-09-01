// client.h

/* (C) 1998-2000 В.В.Пекунов */
/* Модуль различных процедур, необходимых запускаемому приложению для
   взаимодействия с программой-имитатором */

/* Содержит функции AbortServer, printf, scanf и
   специальные функции :
     1) RecvParameters - прием блока параметров,
        передаваемого из функции StartThread.
     2) Check - обращение к Windows. Желательно периодически вызывать
        эту функцию из приложения.
   Прочие функции являются служебными и не должны вызываться
   непосредственно из приложения. */

/* Приложение может подключать этот модуль (#include) ТОЛЬКО ПОСЛЕ
   подключения всех стандартных библиотек BC++ */

/* Как имитатор, так и запускаемые с его помощью приложения должны
   компилироваться как Windows-приложения с моделью памяти Large */

#ifdef _Windows

#if defined(__LARGE__)

#ifndef __CLIENT__

#define __CLIENT__

#include <stdio.h>
#include <conio.h>
#include <stdarg.h>
#include <string.h>
#include <dos.h>
#include <dir.h>
#include <windows.h>
#include <ddeml.h>

#include "common.h"

void InitClient(void);

#include "root.h"

/* Выключаем выравнивание */
#pragma option -a-

/* Определения макросов ParallelList, ParallelRange, ParallelElse */
#define ParallelList if(PresentInList(
#define Of ,0xFF))

#define ParallelRange if(PresentInRanges(
#define To ,

#define ParallelElse else

/* "ParallelList" число ["," число] "Of"
       оператор
   [
   "ParallelElse"
       оператор
   ]

   "ParallelRange" ранг ["," ранг] "Of"
       оператор
   [
   "ParallelElse"
       оператор
   ]

   ранг = минимум "To" максимум
*/

/* Обработчик окна демонстрации связей */
static HWND LinksWindow = NULL;

/* Число виртуальных процессоров в системе */
static int NumProc = 1;
/* Идентификатор виртуального процессора, на котором запущен данный
   процесс (приложение) */
static int CurProc = 0;

/* Обработчик своего Windows-окна */
static HWND SelfHWND = 0;

/* Флаг отладочного режима (0/1) */
static unsigned char DebugMode = 0;

/* Идентификатор DDE-диалога с программой-сервером */
static HCONV ServerCNV = 0;

/* Признак того, что приложение активизировано как процесс, запущенный
   на имитаторе */
static unsigned char isInitialized = 0;

#pragma warn -par

/* Фиктивная функция обработки транзакций. */
HDDEDATA CALLBACK ClientCallback(UINT wType, UINT wFmt, HCONV hConv,
        HSZ hsz1, HSZ hsz2, HDDEDATA hData, DWORD dwData1, DWORD dwData2)
{
 return DDE_FNOTPROCESSED;
}

#pragma warn +par

/* Структура с информацией о соединении */
typedef struct LinkInfo_s {
  void              * LinkID;
  int                 ToID;
  struct LinkInfo_s * Next;
} LinkInfo_t;

/* Список соединений */
LinkInfo_t * LinkInfoList = NULL;

/* Добавление соединения в список */
LinkInfo_t * AddLink(void * LinkID, int ToID)
{
 LinkInfo_t * Result = (LinkInfo_t *) malloc(sizeof(LinkInfo_t));

 Result->Next   = LinkInfoList;
 LinkInfoList   = Result;
 Result->LinkID = LinkID;
 Result->ToID   = ToID;

 return Result;
}

/* Поиск соединения в списке */
LinkInfo_t * FindLink(void * LinkID)
{
 LinkInfo_t * Result = LinkInfoList;

 while (Result && Result->LinkID!=LinkID)
       Result = Result->Next;

 return Result;
}

/* Уничтожение списка соединений */
void DoneLinks(void)
{
 LinkInfo_t * Buf;

 while (LinkInfoList)
       {
        Buf = LinkInfoList->Next;
        free(LinkInfoList);
        LinkInfoList = Buf;
       }
}

/* Отослать информацию NewState об изменении состояния */
/* текущего процессора окну демонстрации связей */
void SendStateMessage(unsigned int NewState)
{
 if (LinksWindow)
    SendMessage(LinksWindow, WM_CHANGESTATE, CurProc, NewState);
}

/* Отослать информацию NewState об изменении состояния */
/* связи (CurProc,ToProc) окну демонстрации связей     */
void SendLinkMessage(int ToProc, unsigned int State)
{
 if (LinksWindow)
    SendMessage(LinksWindow, WM_CHANGELINK, (CurProc<<8) | ToProc, State);
}

/* Блок деинициализации приложения */
void NewExit(void)
{
 if (isInitialized)
    {
     SendStateMessage(PS_TERMINATED);
     DdeDisconnect(ServerHSZ);
     FreeHSZs();
     if (!DebugMode)
        DestroyWindow(SelfHWND);
     DdeUninitialize(Inst);
     DoneLinks();
    }
}

/* Полный аналог соответствующей функции для Parsytec.
   Прервать работу сервера (программы-имитатора) с кодом ошибки Code */
void AbortServer(int Code)
{
 char Out[80] = "AbortServer";
 /* Послать серверу требование прекратить работу */
 DdeClientTransaction(Out,strlen(Out)+1,
                      ServerCNV,NULL,CF_OWNERDISPLAY,
                      XTYP_EXECUTE,600000,
                      NULL);
 /* Вывод соответствующего сообщения */
 sprintf(Out,"AbortServer(%i) from ProcID%i",Code,CurProc);
 Error(Out);
}

/* Вывод сообщения об ошибке с указанием идентификатора */
/* текущего процессора. Работа приложения прекращается. */
void SlaveError(char * Msg)
{
 char Out[128];
 sprintf(Out,"ProcID%i : %s",CurProc,Msg);
 Error(Out);
}

/* Вывод форматированного сообщения об ошибке с  */
/* указанием идентификатора текущего процессора. */
/* Работа приложения прекращается. */
void SlaveErrorf(const char * format, ...)
{
 va_list args;
 char    Buffer[128];

 va_start(args,format);
 vsprintf(Buffer,format,args);
 va_end(args);
 SlaveError(Buffer);
}

/* Сообщения об ошибках DDE */
static char DdeErrors[DMLERR_LAST-DMLERR_FIRST+1][20] = {
        "AdvAckTimeOut",
        "Busy",
        "DataAckTimeOut",
        "DLL_Not_Initialized",
        "DLL_Usage",
        "ExecAckTimeOut",
        "InvalidParameter",
        "Low_Memory",
        "Memory_Error",
        "NotProcessed",
        "No_Conv_Established",
        "PokeAckTimeOut",
        "PostMsg_Failed",
        "ReEntrancy",
        "Server_Died",
        "Sys_Error",
        "UnAdvAckTimeOut",
        "UnFound_Queue_ID" };

/* Вызывается при ошибке DDE. Выдается сообщение Msg и сообщение о
   последней ошибке DDE. Работа приложения прекращается */
void DdeError(char * Msg)
{
 char Out[80];
 sprintf(Out,"%s (%s)",Msg,DdeErrors[DdeGetLastError(Inst)-DMLERR_FIRST]);
 SlaveError(Out);
}

/* Обращение к Windows (обработка сообщений) */
void Check(void)
{
 MSG Message;

 if (PeekMessage(&Message, NULL, 0, 0, PM_REMOVE))
    {
     TranslateMessage(&Message);
     DispatchMessage(&Message);
    }
}

/* Циклическое повторение POKE-транзакции (пока она не будет обработана
   сервером (имитатором)).
   Data - указатель на буфер с посылаемой информацией.
   Size - размер посылаемой информации в байтах
   Item - идентификатор "темы" транзакции
   Fmt  - Clipboard-формат передаваемых данных
   ErrMsg - Сообщение, которое будет выдано при ошибке */
void RepPokeTransaction(void * Data, long int Size,
                        HSZ Item, unsigned int Fmt,
                        char * ErrMsg)
{
 long Result;

 do
   {
    DdeClientTransaction(Data,Size,
                         ServerCNV,Item,Fmt,
                         XTYP_POKE,600000,
                         (DWORD *) &Result);
    Result &= 0xFFFF;

    if (Result!=DDE_FACK)
       if (Result==DDE_FBUSY)
          {
           /* Если сервер слишком занят, чтобы обработать транзакцию,
              то принудительно переключиться на следующее приложение */
            Yield();
            Check();
          }
       else
        DdeError(ErrMsg);
    }
 while (Result!=DDE_FACK);
}

/* Циклическое повторение REQUEST-транзакции (пока она не будет обработана
   сервером (имитатором)).
   Data - указатель на буфер, в который будет приниматься информация.
   Size - размер буфера в байтах
   Item - идентификатор "темы" транзакции */
void RepReqTransaction(void * Data, long int Size, HSZ Item)
{
 long     Result;
 HDDEDATA DHandle;

 do
   {
    DHandle = DdeClientTransaction(NULL,0,
                                   ServerCNV,Item,CF_OWNERDISPLAY,
                                   XTYP_REQUEST,600000,
                                   (DWORD *) &Result);
    if (!DHandle)
       {
        unsigned int Err = DdeGetLastError(Inst);

        if (Err==DMLERR_NOTPROCESSED || Err==DMLERR_BUSY)
           {
           /* Если сервер слишком занят или пока не может выдать требуемую
              информацию, то принудительно переключиться на следующее
              приложение */
             Yield();
             Check();
           }
        else
          {
           char Out[80];

           sprintf(Out,"Request ERROR (%s)",
                       DdeErrors[Err-DMLERR_FIRST]);
           SlaveError(Out);
          }
       }
   }
 while (!DHandle);

 DdeGetData(DHandle,Data,Size,0);
 DdeFreeDataHandle(DHandle);
}

/* Аналог обычной функции printf, осуществляющей вывод
   в окне данного приложения */
void AppPrintf(const char * format, ...)
{
 va_list args;

 va_start(args,format);
 vprintf(format,args);
 va_end(args);
}

/* Аналог обычной функции printf, осуществляющей вывод
   в окне данного приложения, если установлен отладочный режим */
void DebugMsg(const char * format, ...)
{
 va_list args;

 va_start(args,format);
 if (DebugMode)
    {
     char Buf[512];

     vsprintf(Buf,format,args);
     AppPrintf(Buf);
    }
 else
   Check();
 va_end(args);
}

/* Аналог обычной функции printf, осуществляющей вывод
   в окне программы-имитатора (Parsytec Server) */
int printf(const char * format, ...)
{
 va_list args;
 int     Result;
 char    Buffer[512];

 InitClient();

 va_start(args,format);
 /* Формирование выводимой строки */
 Result = vsprintf(Buffer,format,args);

 /* Передача сформированной строки программе-имитатору */
 RepPokeTransaction((void *) Buffer,strlen(Buffer)+1,
                    DataItemHSZ,CF_OWNERDISPLAY,
                    "Server (printf) ERROR");

 va_end(args);
 return Result;
}

/* Аналог обычной функции scanf, осуществляющей ввод с клавиатуры
   в окне программы-имитатора (Parsytec Server) */
int scanf(const char * format, ...)
{
 va_list args;
 int     Result;
 char    Buffer[512];

 InitClient();

 va_start(args,format);
 /* Запрос от программы-имитатора строки текста */
 RepReqTransaction((void *) Buffer, 512, DataItemHSZ);

 /* Извлечение данных из полученной строки */
 Result = vsscanf(Buffer,format,args);

 va_end(args);
 return Result;
}

/* Функция приема блока параметров в буфер Buf. Если точно известно, что
   приложение было запущено через функцию StartThread и
   ему при этом был передан блок параметров размером больше 0 байт, то
   приложение ОБЯЗАТЕЛЬНО должно вызвать эту функцию сразу же после
   запуска (до вызова любой иной функции, взаимодействующей с
   имитатором (printf,scanf,ConnectLink,SendLink,RecvLink,AbortServer,
   Select и т.п.).
   Функция возвращает размер принятого блока в байтах. */
int RecvParameters(void * Buf)
{
 HANDLE Memo;

 InitClient();

 DebugMsg("Receive parameters :");

 /* Сообщение имитатору о готовности принять блок параметров */
 RepPokeTransaction(&_hInstance,sizeof(unsigned int),
                    StartItemHSZ,CF_OWNERDISPLAY,
                    "RecvParameters (send ID) ERROR");

 /* Ожидание блока параметров */
 RepReqTransaction(&Memo,sizeof(Memo),ParamItemHSZ);

 if (Memo)
    {
     /* Копирование полученных данных в буфер */
     byte * MemoPtr = GlobalLock(Memo);
     int    ArgSize;

     memmove(&ArgSize,MemoPtr,sizeof(int));

     DebugMsg("Block [%i bytes] received\n",ArgSize);

     if (ArgSize)
        memmove(Buf,&MemoPtr[sizeof(int)],ArgSize);

     GlobalUnlock(Memo);
     GlobalFree(Memo);

     return ArgSize;
    }
 else
   {
    DebugMsg("No parameters\n");
    return 0;
   }
}

/* Вспомогательная функция реализации конструкции ParallelList.
   Неявно вызывается при использовании макроса ParallelList.
   Параметр : список целых чисел (по меньшей мере одно число),
     последним элементом должен быть код 0xFF.
   Возвращает 1, если в переданном списке содержится идентификатор
   виртуального процессора текущего процесса. Иначе возвращает 0. */
int PresentInList(int First, ...)
{
 va_list Arg;

 int Result;
 int CurArg;

 InitClient();

 Result = First==CurProc;

 va_start(Arg,First);
 while (!Result && (CurArg = va_arg(Arg,int)) != 0xFF)
    Result = CurArg==CurProc;
 va_end(ap);

 return Result;
}

/* Вспомогательная функция реализации конструкции ParallelRange.
   Неявно вызывается при использовании макроса ParallelRange.
   Параметр : список пар целых чисел (по меньшей мере два числа),
     последним элементом должен быть код 0xFF. Каждая пара чисел
     определяет соответственно минимальную и максимальную границу
     некоторого диапазона.
   Возвращает 1, если в переданном списке содержится пара, определяющая
   некоторый диапазон, в который попадает идентификатор виртуального
   процессора текущего процесса. Иначе возвращает 0. */
int PresentInRanges(int First, int Last, ...)
{
 va_list Arg;

 int Result;

 InitClient();

 Result = First<=CurProc && Last>=CurProc;

 va_start(Arg,Last);
 while (!Result && (First = va_arg(Arg,int))!=0xFF && (Last = va_arg(Arg,int))!=0xFF)
   Result = First<=CurProc && Last>=CurProc;
 va_end(ap);

 return Result;
}

/* Попытка инициализации приложения как процесса, запущенного под
   управлением имитатора */
void InitClient(void)
{
 char Debug;

 char FName[MAXPATH];
 char drive[MAXDRIVE];
 char dir[MAXDIR];
 char file[MAXFILE];
 char ext[MAXEXT];

 unsigned long LinkHWND;

 int DimX,DimY,DimZ;
 int MyX,MyY,MyZ;

 if (!isInitialized)
    if (GetWinFlags() & WF_PMODE)
       if (DdeInitialize(&Inst,ClientCallback,
                APPCLASS_STANDARD+APPCMD_CLIENTONLY,0) == DMLERR_NO_ERROR)
          /* Получение информации о конфигурации из командной строки
             и проверка ее на допустимость */
          if (_argc>=8                    &&
              sscanf(_argv[1],"%i",&MyX)  &&
              sscanf(_argv[2],"%i",&MyY)  &&
              sscanf(_argv[3],"%i",&MyZ)  &&
              sscanf(_argv[4],"%i",&DimX) &&
              sscanf(_argv[5],"%i",&DimY) &&
              sscanf(_argv[6],"%i",&DimZ) &&
              sscanf(_argv[7],"%c",&Debug))
             {
              CurProc = (MyZ*DimY+MyY)*DimX+MyX;
              NumProc = DimX*DimY*DimZ;
              if (CurProc>=0 && CurProc<NumProc &&
                  NumProc>=0 && NumProc<=MaxProc)
                 {
                  isInitialized = 1;
                  DebugMode = (Debug=='Y') || (Debug=='y');

                  if (_argc>=9 && sscanf(_argv[8],"%lu",&LinkHWND))
                     LinksWindow = (HWND) LinkHWND;
                  else
                     LinksWindow = (HWND) 0;

                  SendStateMessage(PS_ACTIVISED);
                  /* Подготовка данных, которые будет возвращать функция
                     GET_ROOT() */
                  RootProc.nProcs   = NumProc;
                  RootProc.MyProcID = CurProc;
                  RootProc.DimX     = DimX;
                  RootProc.DimY     = DimY;
                  RootProc.DimZ     = DimZ;
                  RootProc.MyX      = MyX;
                  RootProc.MyY      = MyY;
                  RootProc.MyZ      = MyZ;

                  RootData.ProcRoot = (RootProc_t *) &RootProc;
                  RootData.DummyRoot = NULL;

                  /* Инициализация прочих данных */
                  AppPrintf("Application = %s\n",_argv[0]);
                  AppPrintf("Processor ID = %i\n",CurProc);
                  AppPrintf("Processor Position = (%i,%i,%i)\n",MyX,MyY,MyZ);
                  AppPrintf("There are %i processors in system (%i,%i,%i)\n",
                            NumProc,DimX,DimY,DimZ);
                  SelfHWND = FindWindow(0,_argv[0]);

                  fnsplit(_argv[0],drive,dir,file,ext);

                  sprintf(FName,"%s%i",file,CurProc);
                  SetWindowText(SelfHWND,FName);

                  CreateHSZs();

                  /* Установление диалога с сервером (программой-имитатором) */
                  ServerCNV = DdeConnect(Inst,ServerHSZ,ServerIO_TopicHSZ,NULL);

                  atexit(NewExit);

                  if (!ServerCNV)
                     DdeError("Can not attach to Parsytec Server");
                 }
              else
                Error("Incorrect parameters");
             }
          else
            Error("Must be run under Parsytec Server");
       else
         Error("ERROR in DDE Initialisation");
    else
      Error("Must be a Windows Protected Mode");
 else
   Check();
}

#endif

#else

#error Must be a LARGE memory model

#endif

#else

#error Must be a Windows Application

#endif
