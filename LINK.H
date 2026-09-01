// link.h

/* (C) 1998-2000 В.В.Пекунов */
/* Модуль процедур для установления соединений (Link) между процессами
   и передачи данных */

/* Содержит функции ConnectLink,SendLink,RecvLink */

/* Приложение может подключать этот модуль (#include) ТОЛЬКО ПОСЛЕ
   подключения всех стандартных библиотек BC++ */

/* Как имитатор, так и запускаемые с его помощью приложения должны
   компилироваться как Windows-приложения с моделью памяти Large */

#ifdef _Windows

#if defined(__LARGE__)

#ifndef _LINK_

#define _LINK_

#include "client.h"

/* Выключаем выравнивание */
#pragma option -a-

typedef void LinkCB_t;

/* Полный аналог соответствующей функции для Parsytec.
   Требование установления соединения текущего процесса с
   процессом, запущенным на виртуальном процессоре Processor.
   RequestId = код требования.
   Error = указатель на переменную целого типа, которая при ошибке
            получит ненулевое значение.
   Соединение будет установлено, если требуемый процесс тоже обратится
   к функции ConnectLink с соответствующим требованием, имеющим
   тот же код (RequestId).
   Функция возвращает идентификатор соединения, используемый в дальнейшем
   при работе с функциями SendLink, RecvLink . */
LinkCB_t * ConnectLink (int Processor, int RequestId, int *Error)
{
 LinkCB_t * Ret = NULL;
 char       Buf[80];
 LinkData   Data;

 InitClient();

 (*Error) = 0;
 DebugMsg("Attempt to ConnectLink with ProcID%i : TopicID = %i\n",
          Processor,RequestId);

 if (Processor>=0 && Processor<NumProc)
    {
     /* Подготовка идентификационных данных */
     Data.SelfID = CurProc;
     Data.ProcID = Processor;
     Data.ReqID  = RequestId;

     SendLinkMessage(Processor,LS_REQ_CONNECT);
     /* Пересылка имитатору идентификационных данных запрашиваемого соединения */
     RepPokeTransaction(&Data,sizeof(Data),
                        ConnectItemHSZ,CF_OWNERDISPLAY,
                        "Connect ERROR");

     DebugMsg("Wait for Link ID = ");
     /* Ожидание получения идентификатора соединения */
     RepReqTransaction((void *) &Ret,sizeof(Ret),ConnectItemHSZ);
     DebugMsg("%lu\n",Ret);
     /* Добавляем соединение в список */
     AddLink(Ret,Processor);
     SendLinkMessage(Processor,LS_CONNECTED);
    }
 else
    SlaveErrorf("ConnectLink Parameter Error : Processor ID%i does not exist",Processor);

 return Ret;
}

/* Полный аналог соответствующей функции для Parsytec.
   Синхронная передача данных по соединению Link из буфера Buf
   размером Size байт.
   Функция возвращает число переданных байт */
long SendLink (LinkCB_t * Link, byte * Buf, long Size)
{
 HANDLE       Memo;
 byte       * MemoPtr;
 LinkInfo_t * LinkInfo = FindLink(Link);

 InitClient();

 if (LinkInfo)
    {
     DebugMsg("Send to ProcID%i (LinkID = %lu) from adr [%p]:%li\n",LinkInfo->ToID,Link,Buf,Size);

     SendLinkMessage(LinkInfo->ToID,LS_SEND);
     /* Подготовка пакета данных */
     Memo = GlobalAlloc(GMEM_MOVEABLE | GMEM_SHARE,
                               sizeof(LinkCB_t *)+sizeof(long)+Size);
     MemoPtr = (byte *) GlobalLock(Memo);

     memmove(MemoPtr,&Link,sizeof(LinkCB_t *));
     memmove(&MemoPtr[sizeof(LinkCB_t *)],&Size,sizeof(long));
     memmove(&MemoPtr[sizeof(LinkCB_t *)+sizeof(long)],Buf,Size);
     GlobalUnlock(Memo);

     /* Пересылка пакета данных имитатору */
     RepPokeTransaction(&Memo,sizeof(Memo),
                        SendItemHSZ,CF_OWNERDISPLAY,
                        "Send ERROR");

     DebugMsg("Wait for Confirmation : ");
     /* Ожидание подтверждения о получении данных процессом с другой стороны
        соединения. */
     RepPokeTransaction(&Link,sizeof(Link),
                        ConfirmItemHSZ,CF_OWNERDISPLAY,
                        "Server ERROR : No Send confirmation");
     SendLinkMessage(LinkInfo->ToID,LS_SENDED);
     DebugMsg("Yes\n");
    }
 else
    SlaveErrorf("SendLink Parameter Error : Link ID = %lu does not exist",Link);

 return Size;
}

/* Полный аналог соответствующей функции для Parsytec.
   Синхронный прием данных по соединению Link в буфера Buf
   размером Size байт.
   Функция возвращает число фактически принятых байт */
long RecvLink (LinkCB_t * Link, byte * Buf, long Size)
{
 HANDLE       Memo;
 byte       * MemoPtr;
 LinkInfo_t * LinkInfo = FindLink(Link);

 LinkCB_t * RecvLink;
 long       RecvSize;

 InitClient();

 if (LinkInfo)
    {
     DebugMsg("Recv from ProcID%i (LinkID = %lu) to adr [%p]:%li\n",LinkInfo->ToID,Link,Buf,Size);

     SendLinkMessage(LinkInfo->ToID,LS_RECEIVE);
     /* Пересылка имитатору идентификатора соединения Link */
     RepPokeTransaction(&Link,sizeof(LinkCB_t *),
                        RecvItemHSZ,CF_OWNERDISPLAY,
                        "Receive ERROR");

     DebugMsg("Wait for Data = ");
     /* Ожидание получения пакета данных */
     RepReqTransaction(&Memo,sizeof(Memo),RecvItemHSZ);

     /* Анализ принятого пакета */
     MemoPtr = GlobalLock(Memo);
     memmove(&RecvLink,MemoPtr,sizeof(LinkCB_t *));
     if (Link!=RecvLink)
        SlaveErrorf("Unknown connection : Link ID = %lu",Link);

     memmove(&RecvSize,&MemoPtr[sizeof(LinkCB_t *)],sizeof(long));

     DebugMsg("%li bytes on Link\n",RecvSize);

     if (RecvSize>Size)
        RecvSize = Size; /* Можно принять не более Size байт */

     DebugMsg("Accepted %li bytes\n",RecvSize);

     SendLinkMessage(LinkInfo->ToID,LS_RECEIVED);
     /* Копирование данных в буфер */
     memmove(Buf,&MemoPtr[sizeof(LinkCB_t *)+sizeof(long)],RecvSize);

     GlobalUnlock(Memo);
     GlobalFree(Memo);
    }
 else
    SlaveErrorf("RecvLink Parameter Error : Link ID = %lu does not exist",Link);

 return RecvSize;
}

#endif

#else

#error Must be a LARGE memory model

#endif

#else

#error Must be a Windows Application

#endif
