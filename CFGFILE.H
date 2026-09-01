#ifndef __CFGFILE__
#define __CFGFILE__

#ifdef __UNIX__
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#else
#include <process.h>
#endif

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DESC_FILE          "eko_3d.params"
#define CHANGE_PARAMS_FILE "params.exe"
#define DATA_PARAMS_FILE   "params.out"

typedef enum {intT = 0,fltT,ffltT,chrT,usrT} cfgTypes;

int cfgTypeSizes[5] = {sizeof(int),sizeof(double),sizeof(float),sizeof(char),0};

typedef char * (*BeforeF)(char *);
typedef int    (*HandleF)(char *);
typedef void   (*AfterF) (void);

int  _NSect;
int  NSect;
int  _NVar;
int  NVar;

typedef struct {
 char * Name;
 char   Read;
} Section;

typedef struct {
 void    * VarMap;
 char    * VarName;
 char      VarType;
 int       VarMaxN;
 int       VarMinN;
 int       VarLong;
 BeforeF   VarBefore;
 HandleF   VarHandle;
 AfterF    VarAfter;
 char      VarUse;
 char      VarSect;
 char      VarBcast;
} Variable;

Section  * Sections;
Variable * Variables;

void InitVars()
{
 _NSect = NSect = 0;
 _NVar  = NVar  = 0;

 Sections = NULL;
 Variables = NULL;
}

#define REALLOC(Ptr,N,TypeItem) Ptr = (TypeItem *) realloc(Ptr,N*sizeof(TypeItem))

void AddSection(const char * Name, char Read)
{
 if (_NSect==NSect)
    {
     _NSect += 10;
     REALLOC(Sections,_NSect,Section);
    }
 Sections[NSect].Name = (char *) malloc((strlen(Name)+1)*sizeof(char));
 strcpy(Sections[NSect].Name, Name);
 Sections[NSect++].Read = Read;
}

void AddVar(void *  Ptr,    const char * Name,   char   Type,
            int     MaxN,   int          MinN,   int    Long,
            BeforeF Before, HandleF      Handle, AfterF After,
            char    Use,    char         Bcast,   char   Sect)
{
 if (_NVar==NVar)
    {
     _NVar += 20;
     REALLOC(Variables,_NVar,Variable);
    }
 Variables[NVar].VarMap = Ptr;
 Variables[NVar].VarName = (char *) malloc((strlen(Name)+1)*sizeof(char));
 strcpy(Variables[NVar].VarName, Name);
 Variables[NVar].VarType   = Type;   Variables[NVar].VarMaxN   = MaxN;
 Variables[NVar].VarMinN   = MinN;   Variables[NVar].VarLong   = Long;
 Variables[NVar].VarBefore = Before; Variables[NVar].VarHandle = Handle;
 Variables[NVar].VarAfter  = After;  Variables[NVar].VarUse    = Use;
 Variables[NVar].VarBcast  = Bcast;  Variables[NVar++].VarSect = Sect;
}

void DoneVars()
{
 int i;

 if (_NSect)
    {
     for (i=0; i<NSect; i++)
         free(Sections[i].Name);
     free(Sections);
    }
 if (_NVar)
    {
     for (i=0; i<NVar; i++)
         free(Variables[i].VarName);
     free(Variables);
    }
}

static char ErrorMsg[80];

char * fgetstr(char * s, int n, FILE * stream)
{
 char * Result = fgets(s,n,stream);

 if (Result)
    (*strchr(s,'\n')) = 0x0;

 return Result;
}

int compare_str(char * s1, char * s2)
{
 for (; *s1 && *s2 && tolower(*s1)==tolower(*s2); s1++, s2++);
 return (int) (tolower(*s1)-tolower(*s2));
}

char * DelSpaces(char * Buf)
{
 while ((*Buf) && isspace(*Buf)) Buf++;

 return Buf;
}

char * GetIdent(char * Buf, char * ID)
{
 Buf = DelSpaces(Buf);

 while ((*Buf) && ((*Buf)=='_' || isalnum(*Buf)))
       (*(ID++)) = (*(Buf++));
 (*ID) = 0x0;

 return Buf;
}

char * CheckSym(char * Buf, char Sym, int * Result)
{
 Buf = DelSpaces(Buf);
 if ((*Result) = (*Buf)==Sym)
    Buf++;
 return Buf;
}

#define ReadItems(Str,Templ,MaxN,Long,Type,Ptr) \
   for (Count = 0, Flag = 1; Count<MaxN && Flag; Count+=Flag) \
       { \
        WorkPtr = strtok(Count==0 ? Str : NULL," "); \
        Flag = WorkPtr && sscanf(WorkPtr,Templ,((Type) Ptr)+Long*Count); \
       }

int ReadSection(FILE * _IN, FILE * _OUT, FILE * _DSC, int ChangeParams, FILE * CfgFile, int * Line, char * Buf)
{
 char * Header = Buf;
 char * WorkPtr;
 char   ID[30];
 int    Section = 0;
 int    Flag;
 int    Count;

 if (!(*Buf)) /* Top of Configuration File */
    while (fgetstr(Buf,1000,CfgFile) && !(*DelSpaces(Buf))) (*Line)++;

 Header = CheckSym(DelSpaces(Header),'[',&Flag);
 if (!Flag) return !sprintf(ErrorMsg,"'[' expected in line %i\n",*Line);
 Header = GetIdent(Header,ID);
 while (Section<NSect && compare_str(ID,Sections[Section].Name))
    Section++;
 if (Section>=NSect) return !sprintf(ErrorMsg,"Unknown section identifier '%s' in line %i\n",ID,*Line);
 Header = CheckSym(DelSpaces(Header),']',&Flag);
 if (!Flag) return !sprintf(ErrorMsg,"']' expected in line %i\n",*Line);
 (*Line)++;

 while (fgetstr(Buf,1000,CfgFile))
   {
    if (*Buf)
       {
        char * Val = GetIdent(Buf,ID);

        if (*ID)
           {
            if (Sections[Section].Read)
               {
                Val = CheckSym(Val,'=',&Flag);

                if (Flag)
                   {
                    char BeforeMsg[1000];
                    int  VarIndex = 0;

                    while (VarIndex<NVar && compare_str(ID,Variables[VarIndex].VarName)) VarIndex++;

                    if (VarIndex<NVar && Variables[VarIndex].VarSect==Section)
                       if ( !(Variables[VarIndex].VarBefore) ||
                            _DSC || // if _DSC!=NULL then do not call Before
                            !(*(Variables[VarIndex].VarBefore(BeforeMsg))) )
                          {
                           char inpBuf[5001] = "";
                           Val = DelSpaces(Val);
                           if (*Val=='$' && *(Val+1)=='$') {
                              char * DescBeg = Val+2;
                              char * Default = NULL;
                              Val = DescBeg;
                              while (*Val && !(*Val=='$' && *(Val+1)=='$')) Val++;
                              if (*Val != '$')
                                 return !sprintf(ErrorMsg,"'%s': descriptor is incomplete. Terminal '$$' expected in line %i\n",ID,*Line);
                              *Val = 0x0;
                              if (_DSC)
                                 fprintf(_DSC, "%s\n%s\n%s\n", ID, DescBeg, Val+2);
                              if (ChangeParams) {
                                 Default = Val+2;
                                 Val = inpBuf;
                                 if (_OUT) fprintf(_OUT, "%s<%s> = ", ID, Default);
                                 if (fgetstr(Val, 5000, _IN)) {
                                    if (*Val) {
                                       if (_OUT) fprintf(_OUT, "'%s' is '%s' now\n", ID, Val);
                                    } else
                                       Val = Default;
                                 } else
                                    return !sprintf(ErrorMsg,"'%s': can't enter value in line %i\n",ID,*Line);
                              } else
                                 Val += 2;
                           }
                           if (!_DSC) { // if _DSC==NULL then read values, call Handle, After
                              Variables[VarIndex].VarUse = 1;
                              switch (Variables[VarIndex].VarType) {
                                case intT: ReadItems(Val,"%i",Variables[VarIndex].VarMaxN,Variables[VarIndex].VarLong,int *,Variables[VarIndex].VarMap)
                                           break;
                                case fltT: ReadItems(Val,"%lf",Variables[VarIndex].VarMaxN,Variables[VarIndex].VarLong,double *,Variables[VarIndex].VarMap)
                                           break;
								case ffltT:ReadItems(Val, "%f", Variables[VarIndex].VarMaxN, Variables[VarIndex].VarLong, float *, Variables[VarIndex].VarMap)
									break;
								case chrT: ReadItems(Val, "%s", Variables[VarIndex].VarMaxN, Variables[VarIndex].VarLong, char *, Variables[VarIndex].VarMap)
                                           break;
                                case usrT: if (Variables[VarIndex].VarHandle)
                                              Count = Variables[VarIndex].VarHandle(Val);
                                           else
                                              Count = 0;
                              }
                              if (Count<Variables[VarIndex].VarMinN)
                                 return !sprintf(ErrorMsg,"Value of '%s' is undefined in line %i\n",ID,*Line);
                              Variables[VarIndex].VarMaxN = Count;
                              if (Variables[VarIndex].VarAfter)
                                 Variables[VarIndex].VarAfter();
                           }
                          }
                       else
                         return !sprintf(ErrorMsg,"Initialization error of '%s' in line %i : %s\n",ID,*Line,BeforeMsg);
#ifndef __USE_VCL__
                    else
                      return !sprintf(ErrorMsg,"Unknown constant '%s' in line %i\n",ID,*Line);
#endif
                   }
                else
                   return !sprintf(ErrorMsg,"'=' expected in line %i\n",*Line);
               }
           }
        else
          if (*Val=='[')
             return 1;
          else
             return !sprintf(ErrorMsg,"'Identifier of constant' expected in line %i\n",*Line);
       }
    (*Line)++;
   }

 return 2;
}

typedef void (*InitFunction)();

int ReadConfig(FILE * _IN, FILE * _OUT, int ChangeParams, char * CfgName, InitFunction InitF)
{
 char Buf[1001] = "";
 int  RetCode;
 int  Count;
 int  Line = 1;

 FILE * Cfg = fopen(CfgName,"rt");
#if !defined(__MVS__)
 FILE * _DSC = ChangeParams ? fopen(DESC_FILE, "w+t") : NULL;
#else
 FILE * _DSC = NULL;
#endif

 if (!Cfg) return !sprintf(ErrorMsg,"Error opening configuration file\n");

 do
   if (!(RetCode = ReadSection(_IN, _OUT, _DSC, 0, Cfg, &Line, Buf))) return 0;
 while (RetCode!=2);

#if !defined(__MVS__)
 if (ChangeParams) {
    FILE * _DATA = NULL;
    int spawn_res = -1;

    long DSC_size = ftell(_DSC);
    fclose(_DSC);

    unlink(DATA_PARAMS_FILE);
    if (DSC_size) {
#ifdef __UNIX__
       pid_t child_pid;
       child_pid = fork();
       if (child_pid != 0) {
          int status;
		  wait(&status);
		  spawn_res = WEXITSTATUS(status) != 0;
	   } else {
          execl(CHANGE_PARAMS_FILE, CHANGE_PARAMS_FILE, DESC_FILE, DATA_PARAMS_FILE, NULL);
		  exit(-1);
       }
#else
       spawn_res = spawnl(P_WAIT, CHANGE_PARAMS_FILE, CHANGE_PARAMS_FILE, DESC_FILE, DATA_PARAMS_FILE, NULL);
#endif
    } else {
       _DATA = fopen(DATA_PARAMS_FILE,"w+t");
       fprintf(_DATA,"NO MODIFIABLE PARAMS DETECTED. THE MONSTER EKO_3D IS HAPPY\n");
       fclose(_DATA);
    }
    _DATA = fopen(DATA_PARAMS_FILE,"rt");

    if (_DATA) fseek(_DATA, 0L, SEEK_END);
    if (spawn_res != 0 || _DATA && ftell(_DATA) > 0) {
       if (_DATA) fseek(_DATA, 0L, SEEK_SET);
       fseek(Cfg, 0L, SEEK_SET);
       Line = 1;
       Buf[0] = 0x0;
       do
         if (!(RetCode = ReadSection(spawn_res == 0 ? _DATA : _IN, spawn_res == 0 ? NULL : _OUT, NULL, 1, Cfg, &Line, Buf))) return 0;
       while (RetCode!=2);
    } else {
       fprintf(_OUT, "Params entering is cancelled. Aborting the program\n");
       if (_DATA) fclose(_DATA);
       return 0;
    }
    if (_DATA) fclose(_DATA);
 }
#endif

 Count = 0;
 while (Count<NVar && Variables[Count].VarUse) Count++;

 if (Count<NVar)
    return !sprintf(ErrorMsg,"'%s' not defined in the section '%s'\n",
                    Variables[Count].VarName,Sections[Variables[Count].VarSect].Name);

 fclose(Cfg);

 InitF();

 return 1;
}

#endif
