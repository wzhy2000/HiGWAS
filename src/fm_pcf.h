//fm_pcf.h: interface for the CFmPcf class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_PCF_H_)
#define _FM_PCF_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "fm_linux.h"

#define PCF_EXCEPT      -1
#define PCF_START       1
#define PCF_PAR_READ    2
#define PCF_DAT_LOAD    3
#define PCF_DAT_MERGE   4
#define PCF_SIMULATE    5
#define PCF_VARSEL      10
#define PCF_REFIT       11
#define PCF_QTLSCAN     10
#define PCF_PERMU       11
#define PCF_GOING       12
#define PCF_END         13

#define LONGPROC_1      10
#define LONGPROC_2      11
#define LONGPROC_GOING  12
#define LONGPROC_END    13


typedef struct pcf_dat{
    long bGoingWell;
    long nSnpSect;
    long nTotalSnpSect;
    long nMcmcRound;
    long nTotalMcmcRound;
    double fELapsedSeconds;
    double fEstimatedSeconds;
    double fProgress;
    double start_timer;
    int  nError;
    char szInfo[256];
} PCFDAT;

class CFmPcf
{
public:
    CFmPcf(char* szCfgFile=NULL);
    virtual ~CFmPcf();

    PCFDAT* GetPcfBlock();
    void SetPcfFile(char* szPcfFile);
    void UpdatePcfFile( int nStatus, int nSnpSect=0,
                int nTotalSect=0,  int nMcmcRound=0,
                int nTotalRound=0, int exitcode=0);

protected:
    static PCFDAT g_PCF;
    static char g_szPcfFile[MAX_PATH+1];
};

#endif
