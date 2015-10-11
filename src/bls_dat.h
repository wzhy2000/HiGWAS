// BLS_dat.h: interface for the BLS_dat class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(BLS_DAT_H__INCLUDED_)
#define BLS_DAT_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "fm_dat.h"
#include "bls_R.h"

class BLS_par;
class CFmSnpMat;
class CFmSimulate;
class CFmMatrix;
class CFmVector;
class CFmVectorStr;
class CFmPackedSNP;

class BLS_dat
{
public:
    BLS_dat(CMDOPTIONS *pCmd);
    virtual ~BLS_dat();

    int LoadPlink( char* szFile_tped, char* szFile_tfam, char* szFile_pheno, bool bZnorm, char* szFile_presig);
    int LoadSimple( char* szFile_snp, char* szFile_pheno, bool bZnorm, char* szFile_presig);
    int AttachSnpmat( CFmMatrix* pFmPhe, CFmMatrix* pFmSnp, bool bZnorm);
    int LoadPhenoOnly(char* szFile_pheno, bool bZnorm );
    int Simulate( BLS_par* par, CMDOPTIONS* pCmd );

    char* GetSnpName(int idx);
    int GetPartialSNP( CFmSnpMat* pMat, int nSnpStart, int nSnpCnt );
    int AppendSnp( CFmSnpMat* pMat );
    int SelectRefitGenos(CFmVector* pSelSnp, CFmSnpMat* pMat, bool bGenNorm);
    int Summary(char* szOutFile=NULL);

public:
    CFmPackedSNP* m_pPackedSNP;
    CFmMatrix* m_pSimuSnps;   //Matrix[P,N] for simulation data

    CFmVector* m_pPhenoY;     //Matrix[N,1];-->CFmVector
    CFmMatrix* m_pCovars;     //Matrix[N,1];-->CFmVector
    CFmVectorStr* m_pSnpNames;

    bool m_bSimulate;
    int m_nRound;
    int m_nSubjN;
    int m_nSnpP;
    int m_nMesuTime;

private:
    void ResetAllObjects();

    CFmDat_Plink*  m_pPlink;
    CFmDat_Simple* m_pSimple;
    CFmSimulate*   m_pSimulate;

    CMDOPTIONS *m_pCmd;
    char* m_sPheno_file;
    char* m_sGeno_file;
};

void destroy( BLS_dat* p);

#endif


