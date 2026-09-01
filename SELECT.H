// select.h

/* (C) 1998-2000 В.В.Пекунов */
/* Модуль содержит вспомогательные функции, используемые при
   обмене данными через SendLink,RecvLink */

/* Содержит функции ReceiveOption,Select,SelectList */

/* Приложение может подключать этот модуль (#include) ТОЛЬКО ПОСЛЕ
   подключения всех стандартных библиотек BC++ */

/* Как имитатор, так и запускаемые с его помощью приложения должны
   компилироваться как Windows-приложения с моделью памяти Large */

#ifdef _Windows

#if defined(__LARGE__)

#ifndef _SELECT_

#define _SELECT_

#include <windows.h>
#include <ddeml.h>

#include "client.h"
#include "link.h"

/* Выключаем выравнивание */
#pragma option -a-

typedef LinkCB_t * Option_t;

/* Функция введена для совместимости с Parsytec.
   Возвращает идентификатор выбора соединения Link */
Option_t ReceiveOption (LinkCB_t *Link)
{
 InitClient();

 return Link;
}

/* Полный аналог соответствующей функции для Parsytec.
   На входе :
       nOptions = размер массива идентификаторов выборов соединения
       Options  = массив идентификаторов выбора соединения
   Функция ожидает, пока хотя бы на одном из соответствующих соединений
   не появится информация (передаваемая через SendLink) и вернет позицию
   в массиве Options соответствующего идентификатора выбора  */
int      SelectList (int nOptions, Option_t *Options)
{
 LinkInfo_t * LinkInfo;
 long         result;
 int          ret;
 int          i;

 InitClient();

 DebugMsg("Select from list of ProcID [");
 for (i=0; i<nOptions; i++)
     {
      LinkInfo = FindLink(Options[i]);
      if (LinkInfo)
         DebugMsg("%i%c",LinkInfo->ToID,(i==nOptions-1 ? ']' : ','));
      else
         SlaveErrorf("SelectList Parameter Error : Unknown Link ID = %lu in list",Options[i]);
     }
 DebugMsg(":%i = ",nOptions);
 do
   {
    ret = 0;
    /* Перебор массива Options */
    do
      {
       /* Запрос : присутствуют ли данные на соединении с идентификатором
          выбора Options[ret] ? */
       DdeClientTransaction(&Options[ret],sizeof(Option_t),
                            ServerCNV,SelectItemHSZ,CF_OWNERDISPLAY,
                            XTYP_POKE,600000,
                            (DWORD *) &result);
       result &= 0xFFFF;

       if (result!=DDE_FACK)
          if (result==DDE_FBUSY) ret++;
          else
            DdeError("Server ERROR : SelectList not processed");
      }
    while (ret<nOptions && result!=DDE_FACK);

    if (result!=DDE_FACK)
       {
        /* Ни на одном из соединений данных пока нет. Передаем управление
           следующему процессу */
        Yield();
        Check();
       }
   }
 while (result!=DDE_FACK);

 DebugMsg("Item%i\n",ret);

 return ret;
}

/* Полный аналог соответствующей функции для Parsytec.
   На входе :
nOptions = число(>0) переданных функции идентификаторов выборов соединения.
Option0  = первый из переданных идентификаторов.
   Функция ожидает, пока хотя бы на одном из соответствующих соединений
   не появится информация (передаваемая через SendLink) и вернет позицию
   соответствующего идентификатора выбора в списке параметров функции
   (идентификатор Option0 находится в позиции 0)  */
int      Select (int nOptions, Option_t Option0,...)
{
 /* Обработка осуществляется через обращение к SelectList */
 if (nOptions==1)
    return SelectList(1,(Option_t *) &Option0);
 else
   {
    va_list args;
    int     ret;

    va_start(args,nOptions);
    ret = SelectList(nOptions,(Option_t *) args);
    va_end(args);

    return ret;
   }
}

#endif

#else

#error Must be a LARGE memory model

#endif

#else

#error Must be a Windows Application

#endif
