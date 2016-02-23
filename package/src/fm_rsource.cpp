/* fm_rsource.cpp  -	Run R Script
 *
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */


#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>
#include <R_ext/Parse.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fm_rlogger.h"
#include "fm_err.h"
#include "fm_rsource.h"
#include "fm_matrix.h"
#include "fm_vector.h"
#include "fm_vector_str.h"


CFmRSource::CFmRSource()
{
}

CFmRSource::~CFmRSource()
{
    _log_debug(_HI_, "CFmSys is released successfully.");
}

int CFmRSource::Run2(char* srcFile, char* szCmd)
{
    //R Format:
    //sys.source(file, envir = baseenv(), chdir = FALSE,
    //           keep.source = getOption("keep.source.pkgs"))
    SEXP srcExp;
    PROTECT(srcExp = allocVector(STRSXP, 1));
    SET_STRING_ELT(srcExp, 0, mkChar(srcFile) );

    SEXP chdirExp;
    PROTECT(chdirExp = allocVector(LGLSXP, 1));
    LOGICAL(chdirExp)[0] = FALSE;

    SEXP OptionsExp;
    PROTECT(OptionsExp = allocVector(LGLSXP, 1));
    LOGICAL(OptionsExp)[0] = FALSE;

    SEXP e1;
    int errorOccurred;
    PROTECT(e1 = lang5(install("sys.source"), srcExp, R_GlobalEnv, chdirExp, OptionsExp ) );

    SEXP sys_ret;
    PROTECT(sys_ret = R_tryEval(e1, R_GlobalEnv, &errorOccurred) );
    //int* sys_ret_r = INTEGER(sys_ret);
    //int ret =sys_ret_r[0];
    if (errorOccurred!=0)
    {
        UNPROTECT(5);
        error("invalid call %s", szCmd);
        return(-1);
    }


    //R command
    SEXP cmdSexp, cmdExpr;
    ParseStatus status;

    PROTECT(cmdSexp = allocVector(STRSXP, 1));
    SET_STRING_ELT(cmdSexp, 0, mkChar(szCmd) );
    cmdExpr = PROTECT(R_ParseVector(cmdSexp, -1, &status, R_NilValue));
    if (status!=PARSE_OK)
    {
        UNPROTECT(7);
        error("invalid call %s", szCmd);
        return(-1);
    }

    SEXP ans = R_NilValue;
    for(int i=0; i<length(cmdExpr); i++)
        ans = eval(VECTOR_ELT(cmdExpr,i), R_GlobalEnv);

    UNPROTECT(7);

    return(0);
}


SEXP run_script(char* szSrc)
{
    SEXP cmdSexp, cmdExpr, ans = R_NilValue;
    int i;
    ParseStatus status;

    PROTECT(cmdSexp = allocVector(STRSXP, 1));
    SET_STRING_ELT(cmdSexp, 0, mkChar(szSrc) );
    cmdExpr = PROTECT(R_ParseVector(cmdSexp, -1, &status, R_NilValue));

    if (status!=PARSE_OK)
    {
        UNPROTECT(2);
        error("invalid call %s", szSrc);
    }

    for(i=0; i<length(cmdExpr); i++)
        ans = eval(VECTOR_ELT(cmdExpr,i), R_GlobalEnv);

    UNPROTECT(2);
    return(ans);
}


#define MAX_BUFFER_SIZE 65535*16
int CFmRSource::Run(char* srcFile)
{
    _log_debug(_HI_, "CFmRSource::Run:%s", srcFile);

    FILE* fp = fopen( srcFile, "rt");
    if (fp==NULL)
    {
        _log_error(_HI_, "The file can not be opened(%s)", srcFile);
        return(ERR_OPEN_FILE);
    }

    char szSrc[MAX_BUFFER_SIZE]={0};
    if( fgets( szSrc, MAX_BUFFER_SIZE, fp ) == NULL)
    {
        _log_error( _HI_, "Failed to read a line (%s)", srcFile);
        return(ERR_READ_FILE);
    }

    SEXP pRet = run_script(szSrc);
    return(0);
}

int CFmRSource::Run1(char* srcCmd)
{
    _log_debug(_HI_, "CFmRSource::Run:%s", srcCmd);

    SEXP pRet = run_script(srcCmd);
    return(0);
}

int CFmRSource::SetGlobalVariable(const char* szVar, CFmVector* pVct)
{
    SEXP expVct = GetSEXP(pVct);
    Rf_gsetVar( install(szVar), expVct, R_GlobalEnv );
    //UNPROTECT(1);

    return(0);
}

int CFmRSource::SetGlobalVariable(const char* szVar, CFmMatrix* pMat)
{
    SEXP expMat = GetSEXP( pMat );
    Rf_gsetVar( install(szVar), expMat, R_GlobalEnv );
    //UNPROTECT(1);

    return(0);
}

int CFmRSource::SetGlobalVariable(const char* szVar, CFmVectorStr* pVctStr)
{
    SEXP expVct = GetSEXP(pVctStr);
    Rf_gsetVar( install(szVar), expVct, R_GlobalEnv );
    //UNPROTECT(1);

    return(0);
}

int CFmRSource::GetGlobalVector(const char* szVar, CFmVector* pVct)
{
    SEXP expVct = Rf_findVar( install(szVar), R_GlobalEnv );
    if (expVct ==NULL)
    {
        _log_error( _HI_, "Can't find global variable(%s)", szVar);
        return(-1);
    }

    return( GetVector( expVct, pVct ) );
}

int CFmRSource::GetGlobalMatrix(const char* szVar, CFmMatrix* pMat )
{
    SEXP expMat = Rf_findVar( install(szVar), R_GlobalEnv );
    if (expMat == NULL)
    {
        _log_error( _HI_, "Can't find global variable(%s)", szVar);
        return(-1);
    }

    return( GetMatrix( expMat, pMat ) );
}
