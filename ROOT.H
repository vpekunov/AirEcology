// root.h

/* (C) 1998-2000 В.В.Пекунов */
/* Модуль содержит функцию GET_ROOT запроса данных по конфигурации
   имитируемой системы, а также соответствующие типы данных */

/* Приложение может подключать этот модуль (#include) ТОЛЬКО ПОСЛЕ
   подключения всех стандартных библиотек BC++ */

/* Как имитатор, так и запускаемые с его помощью приложения должны
   компилироваться как Windows-приложения с моделью памяти Large */

#ifdef _Windows

#if defined(__LARGE__)

#ifndef __ROOT__

#define __ROOT__

#define MyProcID MyProcId

/* Выключаем выравнивание */
#pragma option -a-

typedef struct ProcInfo_s {
        int     nProcs;   /* Число виртуальных процессоров */
        int     MyProcId; /* Идентификатор собственного процессора */
        int     DimX;     /* Размер процессорной кладки */
        int     DimY;
        int     DimZ;
        int     MyX;      /* Позиция в процессорной кладке */
        int     MyY;
        int     MyZ;
} ProcInfo_t;

#define RootProc_t ProcInfo_t

typedef struct RootData_s {

    RootProc_t * ProcRoot;
    int        * DummyRoot;

} RootData_t;

static RootData_t RootData; /* Только для внутреннего пользования */
static RootProc_t RootProc; /* Только для внутреннего пользования */

#include "client.h"

/* Полный аналог соответствующей функции для Parsytec.
   Функция возвращает указатель на структуру, содержащую данные о
   конфигурации системы */
RootData_t * GET_ROOT(void)
{
 InitClient();
 return (RootData_t *) &RootData;
}

#endif

#else

#error Must be a LARGE memory model

#endif

#else

#error Must be a Windows Application

#endif
