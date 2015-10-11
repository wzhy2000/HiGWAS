/* fm_pcf.cpp  -	CFmPcf System Object
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

#ifdef WIN32
#include <share.h>
#endif

#include "fm_linux.h"
#include "fm_pcf.h"
#include "fm_rlogger.h"
#include "fm_err.h"

char CFmPcf::g_szPcfFile[];
PCFDAT CFmPcf::g_PCF;

double R_proc_time()
{
    SEXP proc_time, e1;
    int errorOccurred;
    PROTECT(e1 = lang1(install("proc.time")  ) );
    PROTECT(proc_time = R_tryEval(e1, R_GlobalEnv, &errorOccurred) );
    double* proc_time_r = REAL(proc_time);
    double ret = proc_time_r[2];
    UNPROTECT(2);

    return(ret);
}

CFmPcf::CFmPcf(char* szPcfFile)
{
    if (szPcfFile)
    {
        _log_debug( _HI_, "PCF File: %s", szPcfFile);
        strncpy( g_szPcfFile, szPcfFile, MAX_PATH );
        memset( &g_PCF, 0, sizeof(PCFDAT) );
    }
}

CFmPcf::~CFmPcf()
{
    _log_debug(_HI_, "CFmPcf is released successfully.");
}


PCFDAT* CFmPcf::GetPcfBlock()
{
    return &g_PCF;
}

void CFmPcf::SetPcfFile(char* szPcfFile)
{
    strncpy(g_szPcfFile, szPcfFile, MAX_PATH);
}

void CFmPcf::UpdatePcfFile( int nStatus, int nSnpSect,
        int nTotalSect,  int nMcmcRound, int nTotalRound, int exitcode)
{
    if (nStatus != LONGPROC_GOING)
        g_PCF.bGoingWell = nStatus;

    g_PCF.nError = exitcode;
    if (nSnpSect!=-1) g_PCF.nSnpSect = nSnpSect;
    if (nTotalSect!=-1) g_PCF.nTotalSnpSect = nTotalSect;
    if (nMcmcRound!=-1) g_PCF.nMcmcRound = nMcmcRound;
    if (nTotalRound!=-1) g_PCF.nTotalMcmcRound = nTotalRound;

    if (nStatus == LONGPROC_1 || nStatus == LONGPROC_2)
    {
        g_PCF.start_timer = R_proc_time();
        g_PCF.fELapsedSeconds = 0;
        g_PCF.fEstimatedSeconds=0;
        g_PCF.fProgress = 0;
    }

    if (nStatus == LONGPROC_END)
    {
        g_PCF.fELapsedSeconds = R_proc_time() - g_PCF.start_timer;
        g_PCF.fEstimatedSeconds=0;
        g_PCF.fProgress = 1;
    }

    if (nStatus == LONGPROC_GOING)
    {
        if (g_PCF.nTotalMcmcRound!=0 && g_PCF.nTotalSnpSect!=0)
        {
            g_PCF.fProgress = (1.0/g_PCF.nTotalSnpSect)*((double)g_PCF.nSnpSect + (double)(g_PCF.nMcmcRound)/(double)g_PCF.nTotalMcmcRound);

            double nt = R_proc_time() - g_PCF.start_timer;
            g_PCF.fELapsedSeconds = nt;
            if ( g_PCF.fProgress>0 )
                g_PCF.fEstimatedSeconds = nt*( 1 - g_PCF.fProgress)/g_PCF.fProgress;
        }
    }

    char sztemp[256] = {0};
    sprintf(sztemp, "Status:%d, Progress:%.3f, Total:%.2f, Elapsed:%.2f, Estimated:%.2f",
                (int)g_PCF.bGoingWell,
                g_PCF.fProgress*100,
                g_PCF.fELapsedSeconds+g_PCF.fEstimatedSeconds,
                g_PCF.fELapsedSeconds,
                g_PCF.fEstimatedSeconds);

    _log_info(_HI_, sztemp);

    if (strlen(g_szPcfFile)==0)
        return;

    FILE *stream;

    // Open output file for writing. Using _fsopen allows us to
    // ensure that no one else writes to the file while we are
    // writing to it.
#ifdef WIN32
    if( (stream = _fsopen( g_szPcfFile, "wb", _SH_DENYWR )) != NULL )
#else
    if( (stream = fopen( g_szPcfFile, "wb" )) != NULL )
#endif
    {
        fwrite( &g_PCF, sizeof(PCFDAT), 1, stream );
        fclose( stream );
    }
    else
    {
        _log_error(_HI_, "Failed to open PCF file(%s)", g_szPcfFile );
    }
}

