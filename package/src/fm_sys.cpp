/* bls_sys.cpp  -	BLS System Object
 *
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */


#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fm_linux.h"
#include "fm_rlogger.h"
#include "fm_rsource.h"
#include "fm_err.h"
#include "fm_sys.h"

int CFmSys::g_nSysPid = 0;
int CFmSys::g_nTaskId = 1;

int R_sys_getpid()
{
    SEXP sys_pid, e1;
    int errorOccurred;
    PROTECT(e1 = lang1(install("Sys.getpid")  ) );
    PROTECT(sys_pid = R_tryEval(e1, R_GlobalEnv, &errorOccurred) );
    int* sys_pid_r = INTEGER(sys_pid);
    int ret =sys_pid_r[0];
    UNPROTECT(2);

    return(ret);
}

bool _change_fileext(char* filename, const char* newext, char* newfile)
{
    char* pDot = strrchr( filename, '.');
    if (pDot==NULL)
        sprintf(newfile, "%s%s",filename, newext);
    else
    {
        strcpy(newfile, filename);
        pDot = strrchr( newfile, '.');
        strcpy(pDot, newext);
    }

    return(true);
}

CFmSys::CFmSys()
{
}

CFmSys::~CFmSys()
{
    _log_debug(_HI_, "CFmSys is released successfully.");
}

int CFmSys::GetTempId(char* szGrpId)
{
    if (g_nSysPid<=0)
        g_nSysPid = R_sys_getpid();

    return(g_nSysPid*100+(g_nTaskId++));
}

int CFmSys::GetTempFile(char* szFile, const char*szExt, int nLen)
{
	SEXP sexpFile = run_script((char*)"tempfile()");
	const char* pszTmpfile = CHAR(STRING_ELT(sexpFile,0));

    if ( (unsigned int )nLen <= strlen(pszTmpfile) )
        return(-1);

	if(strlen(szExt)>0)
    	sprintf( szFile, "%s.%s",pszTmpfile, szExt);
    else
    	strcpy( szFile, pszTmpfile);

    return(0);
}

int CFmSys::GetSiblingFile(char* szFile, const char* szOrgFile, int nNum, bool bRBase)
{
    char szPath[MAX_PATH]={0};

    const char* pDot0 = strrchr( szOrgFile, '.');
    if (pDot0==NULL)
        sprintf(szPath, "%s-%d", szOrgFile, nNum );
    else
    {
        strcpy(szPath, szOrgFile);
        char* pDot = strrchr( szPath, '.');
        *pDot = '\0';
        sprintf(pDot, "-%d%s", nNum, pDot0);
    }

    strcpy(szFile, szPath);

    if (bRBase)
    {
        for(int i=0; i<strlen(szFile); i++)
            if (szFile[i]=='\\')
            {
                szFile[i]='/';
            }
    }

    return(0);
}


