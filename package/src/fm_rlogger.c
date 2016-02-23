/* logger.c  -	log functions
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include "fm_rlogger.h"

//typedef unsigned char BYTE;

#define FM_LOG_NONE		0
#define FM_LOG_WARN		1
#define FM_LOG_INFO		2
#define FM_LOG_DEBUG	3

static char szPid[256];
static char* pgm_name  = NULL;
static char* szLogName = NULL;
static int  nErrorCount= 0;
static FILE* pFLog     = NULL;
static int nDebug      = FM_LOG_NONE;

/****************
 * start the log file, if szLogFile is NULL, stdout or stderr is used to output log.
 * the Fd where logoutputs should go.
 */

int start_log(int nCmdDebug)
{
    nDebug = nCmdDebug;

#ifdef DEBUG
    if( nDebug )
       	pFLog = stdout;
	else
   		pFLog = stderr;
#else
	//pFLog = fopen("gwas.lasso.so.log", "w+");
	pFLog = NULL;
#endif

    return(0);
}

int stop_log()
{
    return(0);
}

void log_set_pid( int pid )
{
    if( pid )
        sprintf(szPid,"[%u]", (unsigned)pid );
    else
        *szPid = 0;
}

FILE* get_log_stream()
{
#ifdef DEBUG
    if( !pFLog )
        pFLog = stderr;
#endif
    return pFLog;
}

const char* get_log_name(void)
{
    return szLogName? szLogName : "";
}

int get_log_errorcount( int clear)
{
    int n = nErrorCount;
    if( clear )
        nErrorCount = 0;

    return n;
}

void log_print_prefix(const char *text)
{
	if(!pFLog) return;

	if( pgm_name )
		fprintf(pFLog, "%s%s: %s", pgm_name, szPid, text );
	else
		fprintf(pFLog, "?%s:", text );
}

static void log_print_prefix_f(const char *text, const char *fname)
{
	if(!pFLog) return;

	if( pgm_name )
		fprintf(pFLog, "%s%s:%s: %s", pgm_name, szPid, fname, text );
	else
		fprintf(pFLog, "%s:%s ", text, fname );
}

void _log_debug( const char* szSrc, int nSrcLine, const char* fmt, ... )
{
	if(!pFLog) return;

    va_list arg_ptr ;
    if ( nDebug >= FM_LOG_DEBUG )
    {
        log_print_prefix_f("DBG: ", szSrc);
        va_start( arg_ptr, fmt ) ;
        vfprintf(pFLog,fmt,arg_ptr) ;
        va_end(arg_ptr);
        fprintf(pFLog, "\n");
        fflush(pFLog);
    }
}

void _log_info( const char* szSrc, int nSrcLine, const char* fmt, ... )
{
	if(!pFLog) return;

    if ( nDebug >= FM_LOG_DEBUG )
    {
        va_list arg_ptr ;
        log_print_prefix_f("[INF]", szSrc);
        va_start( arg_ptr, fmt ) ;
        vfprintf(pFLog,fmt,arg_ptr) ;
        va_end(arg_ptr);

	    fprintf(pFLog, "\n");
	    fflush(pFLog);
    }
}


void _log_prompt( const char* szSrc, int nSrcLine, const char* fmt, ... )
{
	if(!pFLog) return;

    if ( nDebug >= FM_LOG_WARN )
    {
        va_list arg_ptr ;
        log_print_prefix_f("[*]", "");
        va_start( arg_ptr, fmt ) ;
        vfprintf(pFLog,fmt,arg_ptr) ;
        va_end(arg_ptr);

	    fprintf(pFLog, "\n");
	    fflush(pFLog);
    }
}

void _log_error( const char* szSrc, int nSrcLine, const char*  fmt, ... )
{
    va_list arg_ptr ;
    va_start( arg_ptr, fmt ) ;
	if(!pFLog)
	{
		static char szMsg[1024]={0};
		vsprintf( szMsg, fmt, arg_ptr);
		Rprintf(szMsg);
		Rprintf("\n");
	}
	else
	{
#ifdef DEBUG
	    log_print_prefix_f("[!]", "");
		vfprintf( stderr, fmt,arg_ptr) ;
		fprintf( stderr, "\n");
#endif
	}
	va_end(arg_ptr);
}

void _log_fatal( const char* szSrc, int nSrcLine, const char*  fmt, ... )
{
    va_list arg_ptr ;
    va_start( arg_ptr, fmt ) ;

	if(!pFLog)
	{
		static char szMsg[1024]={0};
		vsprintf( szMsg, fmt, arg_ptr);
		Rprintf("!!!!! ");
		Rprintf(szMsg);
		Rprintf("\n");
	}
	else
	{
#ifdef DEBUG
		log_print_prefix_f("[!]", "");
		vfprintf( stderr, fmt,arg_ptr) ;
		fprintf( stderr, "\n");
#endif
	}
    va_end(arg_ptr);

    stop_log();
    error("! An Exception was caugth, please see log file.");
}

void _log_hexdump( const char *text, const char *buf, size_t len )
{
 	if(!pFLog) return;

 	int i;

    log_print_prefix(text);
    for(i=0; i < len; i++ )
    {
        char c = buf[i];
        fprintf(pFLog, " %02X", (unsigned int)c );
    }

    fputc('\n', pFLog);
}
