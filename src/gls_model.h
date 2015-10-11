// gls_model.h: interface for the GLS class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GLS_H__NCLUDED_)
#define AFX_GLS_H__NCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "fm_matrix.h"
#include "fm_vector.h"
#include "fm_snpmat.h"
#include "gls_R.h"

class GLS_par;
class GLS_dat;
class GLS_res;
class GLS_cfg;
class CFmFileMatrix;

#define RUNMODE_UNK  0
#define RUNMODE_SIM  1
#define RUNMODE_PLINK  2
#define RUNMODE_SIMPLE 3
#define RUNMODE_MERGE  4

class GLS
{
public:
    GLS();
    virtual ~GLS();
public:
    int LoadSimulate( CMDOPTIONS *pCmd, GLS_par *pPar );
    int LoadPlink( CMDOPTIONS *pCmd );
    int LoadSimple( CMDOPTIONS *pCmd );
    int LoadSnpmat( CFmMatrix* pFmPhe, CFmMatrix* pFmSnp, CMDOPTIONS *pCmd );

    int Varsel(GLS_cfg* pCfg);
    int Refit(GLS_cfg* pCfg);

    SEXP GetRObj();

private:
    double func_invGau(double theta, double chi);
    void UpdatePcfFile( int nStatus, int nSnpSect, int nTotalSect,  int nMcmcRound, int nTotalRound);
    int proc_mcmc( CFmMatrix& Y, CFmMatrix& Z, CFmMatrix& Z0, CFmMatrix& X, CFmSnpMat& gen, CFmFileMatrix* pMatRet);
    int ExportFig(char* szFigFile, bool bRefit );

    int m_nDataType;
    bool m_bRefit;

    GLS_dat* m_pDat;
    GLS_res* m_pRes;
    GLS_cfg* m_pCfg;

    CMDOPTIONS *m_pCmd;
    char m_szTraceTag[128];
};

#endif

