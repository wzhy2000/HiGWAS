// bls_model.h: interface for the BLS class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BLS_H__NCLUDED_)
#define AFX_BLS_H__NCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "bls_R.h"
#include "bls_par.h"
#include "bls_dat.h"
#include "bls_res.h"
#include "bls_cfg.h"
#include "fm_filematrix.h"
#include "fm_matrix.h"
#include "fm_vector.h"
#include "fm_snpmat.h"


#define RUNMODE_UNK  0
#define RUNMODE_SIM  1
#define RUNMODE_PLINK  2
#define RUNMODE_SIMPLE 3
#define RUNMODE_MERGE  4

class BLS
{
public:
    BLS();
    virtual ~BLS();
public:
    int LoadSimulate( CMDOPTIONS *pCmd, BLS_par * pPar);
    int LoadPlink( CMDOPTIONS *pCmd );
    int LoadSimple( CMDOPTIONS *pCmd );
    int LoadSnpmat( CFmMatrix* pFmPhe, CFmMatrix* pFmSnp, CMDOPTIONS *pCmd );

    int Varsel(BLS_cfg* pCfg);
    int Refit(BLS_cfg* pCfg);

    SEXP GetRObj();

private:
    CFmMatrix* GetPrePcaAd();
    int GetH2(CFmMatrix* gen_pA, CFmVector &pA, CFmVector &pD, CFmVector& pH2, double var_y );
    int ExportFig(char* szFigFile, bool bRefit );
    double func_invGau(double theta, double chi);
    void UpdatePcfFile( int nStatus, int nSnpSect, int nTotalSect,  int nMcmcRound, int nTotalRound);
    int proc_prec( int nOrder);
    int proc_mcmc( CFmVector& Y0, CFmMatrix& Covs, CFmSnpMat& gen, CFmFileMatrix* r_pFileMat );

    int m_nDataType;
    bool m_bRefit;

    BLS_dat* m_pDat;
    BLS_res* m_pRes;
    BLS_cfg* m_pCfg;

    CMDOPTIONS *m_pCmd;
    char m_szFileMat[MAX_PATH];
    char m_szTraceTag[128];

};

void destroy( BLS* p );

#endif

