// thread.h

/* (C) 1998-2000 В.В.Пекунов */
/* Модуль содержит функцию StartThread, позволяющие запустить несколько
   процессов на одном виртуальном процессоре */

/* Приложение может подключать этот модуль (#include) ТОЛЬКО ПОСЛЕ
   подключения всех стандартных библиотек BC++ */

/* Как имитатор, так и запускаемые с его помощью приложения должны
   компилироваться как Windows-приложения с моделью памяти Large */

#ifdef _Windows

#if defined(__LARGE__)

#ifndef _THREAD_
#define _THREAD_

#include <dir.h>

#include "client.h"

/* Выключаем выравнивание */
#pragma option -a-

#define MIN_STACKSIZE 512

typedef void Thread_t;

/* В отличие от соответствующей функции для Parsytec, здесь запускается
   внешний процесс (скомпилированный EXE-файл).
   Function = имя файла запускаемого внешнего процесса (включая расширение)
   StackSize = введен для совместимости с Parsytec. Не используется.
   Error = указатель на переменную целого типа, которая при ошибке
            получит ненулевое значение.
   ArgSize = Размер передаваемых процессу параметров.
   ... = список передаваемых параметров. При ArgSize==0 игнорируется.
ВНИМАНИЕ : Если запускаемому процессу передаются параметры (ArgSize>0),
   то он ОБЯЗАТЕЛЬНО должен сразу же после запуска принять их с
   помощью функции RecvParameters.
   Если параметры не передаются, то функцию RecvParameters
   вызывать НЕЛЬЗЯ */
#pragma warn -par

Thread_t *StartThread (char * Function, int StackSize,
                       int* Error, int ArgSize, ...)
{
 ThreadData T;
 byte *     MemoPtr;

 va_list  args;
 long     result;
 char     ToRun[MAXPATH];

 InitClient();

 DebugMsg("Attempt to StartThread '%s'\n",Function);
 DebugMsg("Stack parameters will be ignored\n");

 /* Подготовка командной строки */
 sprintf(ToRun,"%s %i %i %i %i %i %i %c",
               Function,
               RootProc.MyX,RootProc.MyY,RootProc.MyZ,
               RootProc.DimX,RootProc.DimY,RootProc.DimZ,
               DebugMode ? 'Y' : 'N');

 va_start(args,ArgSize);
 if (ArgSize) /* Если передаются параметры */
    {
     /* Создать блок параметров */
     T.Data  = GlobalAlloc(GMEM_MOVEABLE | GMEM_SHARE,ArgSize+sizeof(int));
     MemoPtr = (byte *) GlobalLock(T.Data);
     memmove(MemoPtr,&ArgSize,sizeof(int));
     memmove(&MemoPtr[sizeof(int)],args,ArgSize);
     GlobalUnlock(T.Data);
     DebugMsg("Parameter Block [%i bytes] is sended\n",ArgSize);
    }
 else
   T.Data = 0;

 /* Запустить внешний процесс */
 T.hInst = WinExec(ToRun,SW_SHOWMINNOACTIVE);
 if (T.hInst<32)
    {
     DebugMsg("Can't execute Process : ERROR CODE %u\n",T.hInst);
     result = 0;
     *Error = T.hInst;
    }
 else
    {
     DebugMsg("Process is successfully initialized\n");
     if (T.Data) /* Если параметры передаются */
        {
         /* Послать блок параметров программе-имитатору */
         RepPokeTransaction(&T,sizeof(T),
                            ParamItemHSZ,CF_OWNERDISPLAY,
                            "Start data (send) ERROR");
         DebugMsg("Wait for RecvParameters confirmation = ");
         /* Ждать подтверждения (о получении блока порожденным
            процессом) от имитатора */
         RepPokeTransaction(&T.hInst,sizeof(T.hInst),
                            ConfirmPrmItemHSZ,CF_OWNERDISPLAY,
                            "Start data (confirm) ERROR");
         DebugMsg("Yes\n");
        }
     result = T.hInst;
     *Error = 0;
    }
 va_end(args);

 return (Thread_t *) result;
}

#pragma warn +par

#endif

#else

#error Must be a LARGE memory model

#endif

#else

#error Must be a Windows Application

#endif
