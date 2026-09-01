// common.h

/* (C) 1998-2000 В.В.Пекунов */
/* Имитатор многопроцессорной машины :
   Модуль общих (для имитатора и его приложений) определений и процедур */

/* Модуль является служебным. Описанные функции не должны вызываться
   непосредственно из приложения. */

/* Как имитатор, так и запускаемые с его помощью приложения должны
   компилироваться как Windows-приложения с моделью памяти Large */

#ifdef _Windows

#if defined(__LARGE__)

#ifndef __COMMON__

#define __COMMON__

#include <stdlib.h>
#include <stdarg.h>
#include <windows.h>
#include <ddeml.h>

/* Выключаем выравнивание */
#pragma option -a-

/* Ограничение на число виртуальных процессоров (выбрано произвольно) */
#define MaxProc 50

/* Команда инициализации окна демонстрации связей */
#define WM_INITIALIZE   WM_USER+1
/* Команда изменения состояния процессора */
#define WM_CHANGESTATE  WM_USER+2
/* Команда изменения состояния связи */
#define WM_CHANGELINK   WM_USER+3

/* Константы для отображения процессоров в окне демонстрации связей */
#define PS_NONACTIVISED 0
#define PS_ACTIVISED    1
#define PS_TERMINATED   2

/* Константы для отображения состояния связей в окне демонстрации связей */
#define LS_NOTHING      0x7F
#define LS_CONNECTED    0
#define LS_REQ_CONNECT  1
#define LS_RECEIVE      2
#define LS_RECEIVED     3
#define LS_SEND         4
#define LS_SENDED       5

typedef unsigned char byte;
typedef unsigned int  word;

/* Идентификационная запись запроса на соединение (ConnectLink) со стороны
   какого либо процесса */
typedef struct {
   byte SelfID; /* Идентификатор собственного виртуального процессора */
   byte ProcID; /* Идентификатор запрашиваемого виртуального процессора */
   long ReqID;  /* Идентификатор запроса. Процессы получат соединение,
                   если согласуются идентификаторы их виртуальных
                   процессоров и совпадут идентификаторы запроса. */
} LinkData;

/* Идентификационная запись процесса, запущенного через функцию
   StartThread или CreateThread */
typedef struct {
   unsigned int hInst; /* Windows-обработчик соответствующего приложения */
   HANDLE       Data;  /* Обработчик блока данных, содержащего информацию
                          о размере передаваемых процессу параметров и
                          собственно параметры */
} ThreadData;

/* Windows-обработчик процесса */
extern unsigned int _hInstance;

/* Строковые константы, используемые приложениями
   для идентификации транзакций */
static char * ServerName         = "SERVER";
static char * ServerIO_TopicName = "SERVER_IO_TOPIC";
static char * DataItemName       = "DATA_ITEM";
static char * ConnectItemName    = "CONNECT_ITEM";
static char * SendItemName       = "SEND_ITEM";
static char * RecvItemName       = "RECV_ITEM";
static char * ConfirmItemName    = "CONFIRM_ITEM";
static char * SelectItemName     = "SELECT_ITEM";
static char * StartItemName      = "START_ITEM";
static char * ParamItemName      = "PARAM_ITEM";
static char * ConfirmPrmItemName = "CONFIRM_PARAM_ITEM";

/* DDE-обработчик процесса */
static unsigned long Inst = 0;
/* Указатель на DDE-Callback функцию */
static FARPROC  DdeCallbackFunc = NULL;
/* Dde-обработчики перечисленных выше строковых констант */
static HSZ ServerHSZ = 0;
static HSZ ServerIO_TopicHSZ = 0;
static HSZ DataItemHSZ = 0;
static HSZ SendItemHSZ = 0;
static HSZ RecvItemHSZ = 0;
static HSZ ConfirmItemHSZ = 0;
static HSZ ConnectItemHSZ = 0;
static HSZ SelectItemHSZ = 0;
static HSZ StartItemHSZ = 0;
static HSZ ParamItemHSZ = 0;
static HSZ ConfirmPrmItemHSZ = 0;

/* Вывод сообщения об ошибке Msg и выход из программы */
void Error(char * Msg)
{
 MessageBox(0,Msg,NULL,MB_OK+MB_ICONSTOP);
 exit(1);
}

/* Создание DDE-обработчиков строковых констант */
void CreateHSZs(void)
{
 ServerHSZ = DdeCreateStringHandle(Inst,ServerName,0);
 ServerIO_TopicHSZ = DdeCreateStringHandle(Inst,ServerIO_TopicName,0);
 DataItemHSZ = DdeCreateStringHandle(Inst,DataItemName,0);
 SendItemHSZ = DdeCreateStringHandle(Inst,SendItemName,0);
 RecvItemHSZ = DdeCreateStringHandle(Inst,RecvItemName,0);
 ConfirmItemHSZ = DdeCreateStringHandle(Inst,ConfirmItemName,0);
 ConnectItemHSZ = DdeCreateStringHandle(Inst,ConnectItemName,0);
 SelectItemHSZ = DdeCreateStringHandle(Inst,SelectItemName,0);
 StartItemHSZ = DdeCreateStringHandle(Inst,StartItemName,0);
 ParamItemHSZ = DdeCreateStringHandle(Inst,ParamItemName,0);
 ConfirmPrmItemHSZ = DdeCreateStringHandle(Inst,ConfirmPrmItemName,0);
}

/* Высвобождение DDE-обработчиков строковых констант */
void FreeHSZs(void)
{
 DdeFreeStringHandle(Inst,ServerHSZ);
 DdeFreeStringHandle(Inst,ServerIO_TopicHSZ);
 DdeFreeStringHandle(Inst,DataItemHSZ);
 DdeFreeStringHandle(Inst,SendItemHSZ);
 DdeFreeStringHandle(Inst,RecvItemHSZ);
 DdeFreeStringHandle(Inst,ConfirmItemHSZ);
 DdeFreeStringHandle(Inst,ConnectItemHSZ);
 DdeFreeStringHandle(Inst,SelectItemHSZ);
 DdeFreeStringHandle(Inst,StartItemHSZ);
 DdeFreeStringHandle(Inst,ParamItemHSZ);
 DdeFreeStringHandle(Inst,ConfirmPrmItemHSZ);
}

#endif

#else

#error Must be a LARGE memory model

#endif

#else

#error Must be a Windows Application

#endif
