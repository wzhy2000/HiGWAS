// fm_packedsnp.h: interface for the CFmPackedSNP class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_PACKEDSNP_H_)
#define _FM_PACKEDSNP_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "fm_snpmat.h"
#include "fm_matrix.h"
#include "fm_vector.h"

class CFmSnpMat;

typedef struct _snp_info
{
    int  nSnpIdx;
    char szChrId[32];
    char szSnpId[32];
    char szGenDist[32];
    char szBasePair[32];
} SNPINFO;

class CFmPackedSNP
{
public:
    CFmPackedSNP();
    virtual ~CFmPackedSNP();

    /* snp_type: 0, aa, 1:Aa 2:AA 3:NA*/
    inline void Set(long nOutSubj, long nOutSnp, unsigned short snp_type)
    {
        if(nOutSubj >= m_nSubjN || nOutSubj<0 || nOutSnp >= m_nSnpP || nOutSnp<0 )
            throw("invalid subject id or snp id.");

        int nSnp = (int)(m_pSnpIdx->Get(nOutSnp));
        int nSubj = (int)(m_pSubjIdx->Get(nOutSubj));

        //HERE INSIDE index must be provided
        unsigned char* pSnpStart = m_pData + (m_nBaseSubjN/4+1)*nSnp;
        pSnpStart += (nSubj/4);

        unsigned char cSnp4 = *pSnpStart;
        unsigned char newBit = (~( 0x3 << ((3 - nSubj%4)*2) ) & cSnp4) | (snp_type<<((3 - nSubj%4)*2));

        *pSnpStart = newBit;
        return;
    };

    inline unsigned short Get( long nOutSubj, long nOutSnp )
    {
        if(nOutSubj >= m_nSubjN || nOutSubj<0 || nOutSnp >= m_nSnpP || nOutSnp<0 )
            throw("invalid subject id or snp id.");

        int nSnp = (int)(m_pSnpIdx->Get(nOutSnp));
        int nSubj = (int)(m_pSubjIdx->Get(nOutSubj));

        //HERE INSIDE index must be provided
        unsigned char* pSnpStart = m_pData + (m_nBaseSubjN/4+1)*nSnp;
        pSnpStart += (nSubj/4);

        unsigned char cSnp4 = *pSnpStart;
        unsigned char snp = cSnp4 >> ((3 - nSubj%4)*2) & 0x03;

        if ((int)snp>3 )
            throw("invalid operator >>");

        return (unsigned short)snp;
    }

    inline int  GetNumSnps()
    {
        return m_nSnpP;
    };

    inline int  GetNumSubjs()
    {
        return m_nSubjN;
    };

    void InitData(long nSubj, long nSnp);
    int  AppendSnp();
    bool GetSnpRow(int nOutSnp, CFmVector* pVct);
    bool GetSnpRow(int nOutSnp, double* pValue);
    bool SetSnpRow(int nInSnp, double* pValue);

    bool SetSnpInfo( int nInSnp, char* szSnpId, char* szChrId, char* szGenDist, char* szBasePair, char sBigA );
    bool GetSnpInfo( int nOutSnp, SNPINFO** ppInfo );
    bool RemoveSNPs(CFmVector& nOutRemoves);
    bool RemoveSubjs(CFmVector& nOutRemoves);
    int  RemoveRareSNP( double fRareMaf );

    bool OrderByChrPos();
    int WriteAsCSVFile(const char* filename);

private:
    void quickSort(SNPINFO** items, int left, int right);

private:
    long m_nBaseSubjN;
    long m_nSubjN;
    long m_nSnpP;
    long m_nInsideSnpP;
    long m_nInsideSubjN;
    long m_nMaxSnps;
    CFmVector* m_pSnpIdx;
    CFmVector* m_pSubjIdx;
    SNPINFO** m_ppSnpInfo;
    unsigned char* m_pData;
};

void destroy(CFmPackedSNP* p);

#endif
