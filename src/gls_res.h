// GLS_res.h: interface for the result class
//
//////////////////////////////////////////////////////////////////////

#if !defined(GLS_RES_H__INCLUDED_)
#define GLS_RES_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdlib.h>
#include "fm_linux.h"
#include "gls_R.h"

#define TMP_FILE_RA     1
#define TMP_FILE_RD     2
#define TMP_FILE_RALPHA 3
#define TMP_FILE_RMU    4

#define LG	4

class CFmVectorStr;
class CFmMatrix;
class CFmVector;
class GLS_cfg;

class GLS_res
{
public:
    GLS_res(CMDOPTIONS *pCmd, GLS_cfg* pCfg);
    virtual ~GLS_res();

    int InitVarsel(int mRound, int nSnpP, int SubjN, int mcmcIter);
    int SetMcmcResults(bool bRefit, CFmVectorStr* pVctSnpName, CFmVector* pVctChr, CFmVector* pVctPos,CFmFileMatrix* pMatRet, int nCov);
    int InitRefit( int mcmcIter );

    char* GetVarselSnpName(int idx);
    int Summary(char* szOutFile = NULL);
    int GetRefitSnp(CFmVector** pVct);
    int GetSigSNP();

    int SaveRData( char* szRdataFile );
    int SaveCsv4Fig( char* szCsvFile, bool bRefit );
    int SaveSigFile( char* szCsvFile, bool bRefit );
    SEXP GetRObj();

private:
    int SortoutMcmc(CFmFileMatrix* pMatRet, CFmMatrix* pMu, CFmMatrix* pAlpha, CFmMatrix* pRa, CFmMatrix* pRd, int nSnp, int nCov);
    int GetMcmcInfo(CFmFileMatrix* pMatRet, int idx, CFmVector* pFmInfo, CFmVector* pFmModel, double fQval );
    int SortoutBestQ(CFmFileMatrix* pMatRet, CFmMatrix* pBestQ, int nSnp, int nCov);
    int SortoutPSRF(CFmFileMatrix* pMatRet, CFmMatrix* pMatPSRF, int nSnp, int nCov);
    int GetBestQInfo(CFmFileMatrix* pMatRet, int idx, CFmVector& fmQBest, CFmVector& fmQPosMean, CFmVector& fmQPosMin, CFmVector& fmQPosMax);
    int GetPSRFInfo(CFmFileMatrix* pMatRet, int idx, CFmMatrix& fmPSRF_R);

    int m_nSimuRound;
    int m_nSnpP;
    int m_nSubjN;
    int m_nMcmcIter;
    int m_nRefitSnp;

    GLS_cfg* m_pCfg;
    CMDOPTIONS *m_pCmd;

    CFmVectorStr* m_pVarsel_SnpName;
    CFmVector* m_pVarsel_SnpChr;
    CFmVector* m_pVarsel_SnpPos;
    CFmMatrix* m_pVarsel_Mu;	//[1,19], 1-4: Model; 5-9:L2/Median; 10-14: L2/Min; 15-19:L2/Max;
    CFmMatrix* m_pVarsel_Alpha; //[nCovar, 19], 1-4: Model; 5-9:L2/Median; 10-14: L2/Min; 15-19:L2/Max;
    CFmMatrix* m_pVarsel_Ra;	//[P, 19], 1-4: Model; 5-9:L2/Median; 10-14: L2/Min; 15-19:L2/Max;
    CFmMatrix* m_pVarsel_Rd;    //[P, 19], 1-4: Model; 5-9:L2/Median; 10-14: L2/Min; 15-19:L2/Max;
    CFmMatrix* m_pVarsel_BestQ; //[P, 19], 1-4: Q; 5-9:L2/Median; 10-14: L2/Min; 15-19:L2/Max;
    CFmMatrix* m_pVarsel_PSRF;

    CFmVectorStr* m_pRefit_SnpName;
    CFmVector* m_pRefit_SnpChr;
    CFmVector* m_pRefit_SnpPos;
    CFmMatrix* m_pRefit_Mu;	//[1,19], 1-4: Model; 5-9:L2/Median; 10-14: L2/Min; 15-19:L2/Max;
    CFmMatrix* m_pRefit_Alpha;  //[nCovar, 19], 1-4: Model; 5-9:L2/Median; 10-14: L2/Min; 15-19:L2/Max;
    CFmMatrix* m_pRefit_Ra;	//[P, 19], 1-4: Model; 5-9:L2/Median; 10-14: L2/Min; 15-19:L2/Max;
    CFmMatrix* m_pRefit_Rd;     //[P, 19], 1-4: Model; 5-9:L2/Median; 10-14: L2/Min; 15-19:L2/Max;
    CFmMatrix* m_pRefit_BestQ;  //[P, 19], 1-4: Q; 5-9:L2/Median; 10-14: L2/Min; 15-19:L2/Max;
    CFmMatrix* m_pRefit_PSRF;

    CFmVector* m_pRefitSnps;
    CFmVector* m_pSigSnps;
    CFmVector* m_pSigAddSnps;
    CFmVector* m_pSigDomSnps;
    int m_nTotalSnp;
};

void destroy( GLS_res* p);

#endif
