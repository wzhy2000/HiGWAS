// BLS_res.h: interface for the result class
//
//////////////////////////////////////////////////////////////////////

#if !defined(BLS_RES_H__INCLUDED_)
#define BLS_RES_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdlib.h>


#define TMP_FILE_COEFS  0
#define TMP_FILE_RA     1
#define TMP_FILE_RD     2
#define TMP_FILE_RH2    3

class BLS_cfg;
class CFmVectorStr;
class CFmMatrix;
class CFmVector;
class CFmFileMatrix;

class BLS_res
{
public:
    BLS_res(CMDOPTIONS *pCmd, BLS_cfg* pCfg);
    virtual ~BLS_res();

    int Summary(char* szOutFile = NULL);
    int SetResult(bool bRefit, CFmMatrix* gen_pA, double y_var, int nCov,
		CFmVectorStr* pVctSnp, CFmVector* pVctChr, CFmVector* pVctPos, CFmFileMatrix* pFileMat);
    int InitVarsel(int mRound, int nSnpP, int SubjN, int mcmcSect, int mcmcSnps, int mcmcIter, CFmVectorStr* pCovNames);
    int InitRefit( int mcmcIter, CFmVectorStr* pCovNames );

    int GetSigSNP();
    int GetRefitSnp(CFmVector* pVct);
    char* GetVarselSnpName(int idx);

    int SaveRData( char* szRdataFile );
    int SaveCsv4Fig( char* szRdataFile, bool bRefit );
    int SaveSigFile( char* szCsvFile, bool bRefit );
    SEXP GetRObj();

private:
    int GetH2(CFmMatrix* gen_pA, CFmVector& pA, CFmVector& pD, CFmVector* pH2, double y_var );
    int SortoutResult( CFmFileMatrix* pFileMat, CFmMatrix* gen_pA, double y_var, int nSnp, int nCov,
		 CFmMatrix* pmatCoef, CFmMatrix* pmatA, CFmMatrix* pmatD, CFmVector* pvctH2 );
    int GetMcmcResult( CFmFileMatrix* pFileMat, long int idx, double fQval, CFmVector& fmVctRet );
    int SortoutQBest( CFmFileMatrix* pFileMat, CFmMatrix* pBestQ, int nSnp, int nCov );
    int GetQBestInfo(CFmFileMatrix* pMatRet, int idx, CFmVector& fmQBest);

private:
    int m_nSimuRound;
    int m_nSnpP;
    int m_nSubjN;
    int m_nMcmcSect;
    int m_nMcmcSnp;
    int m_nMcmcIter;
    int m_nRefitSnp;

    BLS_cfg* m_pCfg;
    CMDOPTIONS *m_pCmd;

    CFmVectorStr* m_pVarsel_SnpName;
    CFmVector* m_pVarsel_SnpChr;
    CFmVector* m_pVarsel_SnpPos;
    CFmMatrix* m_pVarsel_Ra;
    CFmMatrix* m_pVarsel_Rd;
    CFmVector* m_pVarsel_Rh2;
    CFmMatrix* m_pVarsel_Coefs;
    CFmMatrix* m_pVarsel_QBest;

    CFmVectorStr* m_pRefit_SnpName;
    CFmVector* m_pRefit_SnpChr;
    CFmVector* m_pRefit_SnpPos;
    CFmMatrix* m_pRefit_Ra;
    CFmMatrix* m_pRefit_Rd;
    CFmVector* m_pRefit_Rh2;
    CFmMatrix* m_pRefit_Coefs;
    CFmMatrix* m_pRefit_QBest;

    CFmVector* m_pRefitSnps;
    CFmVector* m_pSigSnps;
    CFmVector* m_pSigAddSnps;
    CFmVector* m_pSigDomSnps;
    CFmVectorStr* m_pCovNames;
    int m_nTotalSnp;

};

void destroy( BLS_res* p );

#endif
