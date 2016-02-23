// fm_dat.h: interface for the CFmDat class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_DAT_H_)
#define _FM_DAT_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "fm_snpmat.h"
#include "fm_matrix.h"
#include "fm_vector.h"
#include "fm_packedsnp.h"

#define MAX_TIME_INPHENO   128
#define MAX_SUBJ_INDATA    4096
#define MAX_SNPS_INTPED    16*256*256  //104,8576 SNPs*4096 Subjs/4 = 1, 073,741,824 (1GB memory)

class CFmSnpMat;

class CFmDat_Plink
{
public:
    CFmDat_Plink(char* szFile_tped, char* szFile_tfam );
    virtual ~CFmDat_Plink();

    int Load(char* szPresigFile=NULL);

    int m_nSnpP;
    CFmPackedSNP m_PackedSNP;
    CFmVectorStr* m_pSnpNames;
    CFmVectorStr* m_pSubIds;

private:
    bool Extract_snpinfo( char* szLine, CFmVectorStr* pSigSnp );
    int  LoadTped( char* szFile_tped, char* szPresig );
    int  LoadTfam( char* szFile_tfam );
    char* m_sTped_file;
    char* m_sTfam_file;
};

class CFmDat_Simple
{
public:
    CFmDat_Simple(char* szFile_snp );
    virtual ~CFmDat_Simple();

    int Load(char* szPresigFile=NULL);

    CFmPackedSNP m_PackedSNP;
    CFmVectorStr* m_pSnpNames;
    CFmVectorStr* m_pSubIds;

private:
    char* m_szFile_snp;
};

class CFmDat_Pheno
{
public:
    CFmDat_Pheno(char* szFile_pheno, bool bzNorm, char* szYname, char* szZname, char* szXname );
    CFmDat_Pheno(CFmMatrix* pFmPhe, bool bzNorm, char* szYname, char* szZname, char* szXname );
    virtual ~CFmDat_Pheno();

    void Init();
    int LoadNonlongdt(  CFmPackedSNP* pPackedSNP, CFmVectorStr* pFamSubjs  );
    int LoadLongdt(  CFmPackedSNP* pPackedSNP, CFmVectorStr* pFamSubjs );

    CFmVectorStr* m_pSubjs;   //Vector[ N ];
    CFmMatrix* m_pPhenoY;     //Matrix[ N, Q ];
    CFmMatrix* m_pPhenoZ;     //Matrix[ N, Q ];
    CFmMatrix* m_pPhenoZ0;    //Matrix[ N, Q ];
    CFmMatrix* m_pCovars;     //Matrix[ N, Model.par ];

    int m_nSubjN;
    int m_nMesuQ;
    int m_nCovCount;

    char* m_pszXname;
    char* m_pszYname;
    char* m_pszZname;

    CFmVectorStr* m_pXnames;

protected:
    int RemoveMissing( CFmPackedSNP* pPackedSNP );

private:
    char* m_szFile_pheno;
    bool m_bZnorm;
    CFmMatrix* m_pAttachedPhe;
};

void destroy( CFmDat_Plink* p );
void destroy( CFmDat_Simple* p );
void destroy( CFmDat_Pheno* p );

#endif
