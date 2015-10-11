// gls_dat.h: interface for the GLS_dat class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(GLS_DAT_H__INCLUDED_)
#define GLS_DAT_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "gls_R.h"
#include "fm_dat.h"

class GLS_par;
class CFmSnpMat;
class CFmSimulate;
class CFmMatrix;
class CFmVector;
class CFmPackedSNP;
class CFmVectorStr;

class GLS_dat
{
public:
    GLS_dat(CMDOPTIONS* pCmd);
    virtual ~GLS_dat();

    int AppendSnp( CFmSnpMat* pMat );
    int LoadPlink( char* szFile_tped, char* szFile_tfam, char* szFile_pheno, bool bZnorm );
    int LoadSimple( char* szFile_snp, char* szFile_pheno, bool bZnorm );
    int LoadPhenoOnly(char* szFile_pheno, bool bZnorm );

    int AttachSnpmat( CFmMatrix* pFmPhe, CFmMatrix* pFmSnp, bool bZnorm );

    int Simulate( GLS_par* par, CMDOPTIONS* pCmd );
    int Simulate_test( GLS_par* par, CMDOPTIONS* pCmd );

    int GetPartialSNP( CFmSnpMat* pMat, int nSnpStart, int nSnpCnt );
    int SelectRefitGenos(CFmVector* pSelSnp, CFmSnpMat* pMat, bool bGenNorm);
    int Summary(char* szOutFile=NULL);
    char* GetSnpName(int idx);

public:
    CFmPackedSNP* m_pPackedSNP;
    CFmMatrix* m_pSimuSnps;   //Matrix[P,N] for simulation data
    CFmMatrix* m_pPhenoY;     //Matrix[N,Q];-->CFmVector
    CFmMatrix* m_pPhenoZ;     //Matrix[N,Q];
    CFmMatrix* m_pPhenoZ0;     //Matrix[N,Q];
    CFmMatrix* m_pCovars;     //VECTOR[N];
    CFmVector* m_pZRange;     //VECTOR[Q];
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
    CMDOPTIONS*    m_pCmd;

    char* m_sPheno_file;
    char* m_sGeno_file;
};

void destroy( GLS_dat* p );

#endif
