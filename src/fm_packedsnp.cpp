/* fm_packedsnp.cpp  -	Packed SNP Object
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <stdio.h>
#include <math.h>

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>
#include <Rmath.h>

#include "fm_packedsnp.h"
#include "fm_rlogger.h"
#include "fm_err.h"
#include "fm_new.h"

#define MAX_TIME_INPHENO   128
#define MAX_SUBJ_INDATA    4096
#define MAX_SNPS_INTPED    16*256*256  //104,8576 SNPs*4096 Subjs/4 = 1, 073,741,824 (1GB memory)

CFmPackedSNP::CFmPackedSNP()
{
    m_nSubjN = 0;
    m_nSnpP = 0;
	m_nBaseSubjN = 0;
    m_nInsideSnpP = 0;
    m_nInsideSubjN = 0;
    m_nMaxSnps = 0;
    m_pData = NULL;
    m_ppSnpInfo = NULL;
    m_pSnpIdx = NULL;
    m_pSubjIdx = NULL;
}

//HERE INSIDE index must be provided
CFmPackedSNP::~CFmPackedSNP()
{
    if(m_pData)
        Free(m_pData);

    if (m_ppSnpInfo)
    {
        for (int i=0;i<m_nMaxSnps;i++)
            if (m_ppSnpInfo[i]) delete (m_ppSnpInfo[i]);

        Free(m_ppSnpInfo);
    }

    if(m_pSnpIdx)
        destroy( m_pSnpIdx );

    if(m_pSubjIdx)
        destroy(  m_pSubjIdx );

    _log_debug(_HI_, "CFmPackedSNP is released successfully.");
}

void CFmPackedSNP::InitData(long nSubj, long nSnp)
{
    m_nSnpP = nSnp;
    m_nSubjN = nSubj;
    m_nBaseSubjN = nSubj;

    m_nInsideSubjN = nSubj;
    m_nInsideSnpP = nSnp;

    if (nSnp < 1000)
        m_nMaxSnps = 1000;
    else
        m_nMaxSnps = nSnp;

    //** initialize m_pData
    m_pData = (unsigned char*)Calloc( (m_nBaseSubjN/4+1)*m_nMaxSnps, char );
    if (!m_pData)
    {
        _log_fatal(_HI_, "Failed to alloc memory for the packed SNP structure:(%d bytes)",(m_nBaseSubjN/4+1)*m_nMaxSnps);
    }
    memset(m_pData, 0, (m_nBaseSubjN/4+1)*m_nMaxSnps );

    //** initialize SNP info pointer
    m_ppSnpInfo = (SNPINFO**)Calloc( m_nMaxSnps, SNPINFO* );
    if (!m_ppSnpInfo)
    {
        _log_fatal(_HI_, "Failed to alloc memory for the packed SNP info:(%d bytes)",sizeof(SNPINFO*)*m_nMaxSnps);
    }
    memset(m_ppSnpInfo, 0, sizeof(SNPINFO*)*m_nMaxSnps );

    //** initialize outside<->inside index table
    CFmNewTemp refNew;
    m_pSnpIdx = new (refNew) CFmVector(m_nMaxSnps, 0.0);
    for (int i=0; i<m_nMaxSnps; i++)
        m_pSnpIdx->Set(i, i);

    m_pSubjIdx  = new (refNew) CFmVector(m_nBaseSubjN, 0.0);
    for (int i=0; i<m_nBaseSubjN; i++)
        m_pSubjIdx->Set(i, i);
}

int CFmPackedSNP::AppendSnp()
{
    if ( m_nMaxSnps > m_nSnpP+1 )
    {
        m_nSnpP++;
        m_nInsideSnpP = m_nSnpP;
        return(m_nSnpP-1);
    }

    int newMax = m_nMaxSnps + 1000;

    // enlarge m_pData buffer
    unsigned char* pNewData = (unsigned char*)Calloc( (m_nBaseSubjN/4+1)*newMax, char  );
    if (!pNewData)
    {
        _log_fatal(_HI_, "Failed to alloc memory for the packed SNP structure:(%d bytes)", (m_nBaseSubjN/4+1)*newMax);
        return( ERR_MEM_ALLOC );
    }

    memset(pNewData, 0, (m_nBaseSubjN/4+1)*newMax );
    memcpy(pNewData, m_pData, (m_nBaseSubjN/4+1)*m_nMaxSnps );
    Free(m_pData);
    m_pData = pNewData;

    // enlarge m_ppSnpInfo buffer
    SNPINFO** ppNewSnpInfo = (SNPINFO**)Calloc( newMax, SNPINFO* );
    if (!ppNewSnpInfo)
    {
        _log_fatal(_HI_, "Failed to alloc memory for the packed SNP structure:(%d bytes)", sizeof(SNPINFO*)*newMax);
        return( ERR_MEM_ALLOC );
    }

    memset(ppNewSnpInfo, 0, sizeof(SNPINFO*)*newMax );
    memcpy(ppNewSnpInfo, m_ppSnpInfo, sizeof(SNPINFO*)*m_nMaxSnps );
    Free(m_ppSnpInfo);
    m_ppSnpInfo = ppNewSnpInfo;

	CFmNewTemp refNew;

    //** enlarge outside<->inside index table
    destroy( m_pSnpIdx );
    m_pSnpIdx = new (refNew) CFmVector(newMax, 0.0);
    for (int i=0; i<newMax; i++)
        m_pSnpIdx->Set(i, i);

    m_nMaxSnps = newMax;
    m_nSnpP++;
    m_nInsideSnpP = m_nSnpP;

    return(m_nSnpP-1);
}

bool CFmPackedSNP::GetSnpRow(int nOutSnpIdx, CFmVector* pVct)
{
    if (nOutSnpIdx<0 || nOutSnpIdx>=m_nSnpP )
    {
        _log_error(_HI_, "CFmPackedSNP::GetSnpRow: invalid row(%d).", nOutSnpIdx);
        return false;
    }

    pVct->Resize(0);
    for (int i=0;i<m_nSubjN;i++)
        pVct->Put( Get(i, nOutSnpIdx) );

    return (true);
}

bool CFmPackedSNP::GetSnpRow(int nOutSnpIdx, double* pValue)
{
    if (nOutSnpIdx<0 || nOutSnpIdx>=m_nSnpP )
    {
        _log_error(_HI_, "CFmPackedSNP::GetSnpRow: invalid row(%d).", nOutSnpIdx);
        return false;
    }

    for (int i=0;i<m_nSubjN;i++)
        pValue[i] = Get(i, nOutSnpIdx);

    return (true);
}

bool CFmPackedSNP::SetSnpRow(int nOutSnpIdx, double* pValue)
{
    if (nOutSnpIdx<0 || nOutSnpIdx>=m_nSnpP)
    {
        _log_error(_HI_, "CFmPackedSNP::SetSnpRow: invalid row(%d).", nOutSnpIdx);
        return false;
    }

    for (int i=0;i<m_nSubjN;i++)
         Set(i, nOutSnpIdx, (short unsigned int)pValue[i]);

    return (true);
}

bool CFmPackedSNP::SetSnpInfo( int nOutSnpIdx, char* szSnpId, char* szChrId, char* szbpPos, char* szBasePair, char sBigA )
{
    if (nOutSnpIdx<0 || nOutSnpIdx>=m_nSnpP)
    {
        _log_error(_HI_, "CFmPackedSNP::SetSnpInfo: invalid row(%d).", nOutSnpIdx);
        return false;
    }

    int nSnpIdx = (int)m_pSnpIdx->Get( nOutSnpIdx );

    //HERE INSIDE index must be provided
    if (m_ppSnpInfo[nSnpIdx]==NULL)
        m_ppSnpInfo[nSnpIdx] = new SNPINFO;

    m_ppSnpInfo[nSnpIdx]->nSnpIdx = nSnpIdx;

    strncpy( m_ppSnpInfo[nSnpIdx]->szChrId   , szChrId   ,31);
    strncpy( m_ppSnpInfo[nSnpIdx]->szSnpId   , szSnpId   ,31);
    strncpy( m_ppSnpInfo[nSnpIdx]->szGenDist , szbpPos ,31);
    strncpy( m_ppSnpInfo[nSnpIdx]->szBasePair, szBasePair,31);

    return(true);
}

bool CFmPackedSNP::GetSnpInfo( int nOutSnpIdx, SNPINFO** ppInfo )
{
    if (nOutSnpIdx<0 || nOutSnpIdx>=m_nSnpP)
    {
        _log_error(_HI_, "CFmPackedSNP::GetSnpInfo: invalid row(%d).", nOutSnpIdx);
        return false;
    }

    int nSnpIdx = (int)m_pSnpIdx->Get( nOutSnpIdx );

    //HERE INSIDE index must be provided
    *ppInfo = m_ppSnpInfo[nSnpIdx];

    return(true);
}

int CFmPackedSNP::RemoveRareSNP(double fRareMaf)
{
    _log_debug( _HI_, " RemoveRareSNP() Start, SnpP=%d, SubjN=%d", m_nSnpP, m_nSubjN);

    double* pSNPs = new double[m_nSubjN];
    CFmVector nRemoves(0, 0.0);

	R_set_seed( 100 );
	GetRNGstate();

    for (int i=0; i<m_nSnpP; i++)
    {
        GetSnpRow(i, pSNPs);

        //test code
        //for(int j=0;j<m_nSubjN; j++)
        //{
        //    if (pSNPs[j] != Get(j, i))
        //        Rprintf("%d", pSNPs[j]);
        //    Rprintf("\n");
        //}

        int n0=0, n1=0, n2=0, n3=0;
        for(int j=0; j<m_nSubjN; j++)
        {
            switch((int)(pSNPs[j]))
            {
                case 0: n0++;break;
                case 1: n1++;break;
                case 2: n2++;break;
                case 3: n3++;break;
            }
        }

        double p1=(2.0*n0+n1)/(2.0*n0+2.0*n1+2.0*n2);
        double p2=(2.0*n2+n1)/(2.0*n0+2.0*n1+2.0*n2);

        if ( fmin(p1,p2) < fRareMaf )
            nRemoves.Put(i);
		else
	    if (n3>0)
        {
            double prob[3]={0};
            prob[1] = 1.0*n0/(n0+n1+n2);
            prob[2] = 1.0*(n0+n1)/(n0+n1+n2);

            for(int j=0; j<m_nSubjN; j++)
            {
                if (pSNPs[j]==3)
                {
                    double missing_snp = runif(0,1);

                    if ( prob[2]<=missing_snp && missing_snp<1 )
                        pSNPs[j] = 2.0;
                    if ( prob[1]<=missing_snp && missing_snp<prob[2] )
                        pSNPs[j] = 1.0;
                    if ( 0<=missing_snp && missing_snp<prob[1] )
                        pSNPs[j] = 0.0;
                }
            }

            SetSnpRow(i, pSNPs);
        }
    }

    if (nRemoves.GetLength()>0)
        RemoveSNPs(nRemoves);

	PutRNGstate();

	_log_debug(_HI_, "RemoveRareSNP Total:%d, removed:%d", m_nSnpP, nRemoves.GetLength());

    return(0);
}

bool CFmPackedSNP::RemoveSNPs(CFmVector& nRemoves)
{
    for (int i=0;i<nRemoves.GetLength();i++)
        if ( nRemoves.Get(i)<0 || nRemoves.Get(i)>=m_nSnpP )
        {
            _log_error(_HI_, "CFmPackedSNP::RemoveSNPs: invalid snp(%d).", nRemoves.Get(i));
            return false;
        }

    m_pSnpIdx->RemoveElements( nRemoves );
    m_nSnpP -= nRemoves.GetLength();

    return(true);
}

bool CFmPackedSNP::RemoveSubjs(CFmVector& nRemoves)
{
    for (int i=0;i<nRemoves.GetLength();i++)
        if ( nRemoves.Get(i)<0 || nRemoves.Get(i)>=m_nSubjN )
        {
            _log_error(_HI_, "CFmPackedSNP::RemoveSubjs: invalid subjs(%d).", nRemoves.Get(i));
            return false;
        }

    m_pSubjIdx->RemoveElements( nRemoves );
    m_nSubjN -= nRemoves.GetLength();
    return(true);
}

bool IsLess(SNPINFO* x, SNPINFO* y)
{
    int nChrX = atoi(x->szChrId);
    if ( strcmp (x->szChrId, "M")==0)
        nChrX = 0;
    else if ( strcmp (x->szChrId, "X")==0)
        nChrX = 23;
    else if ( strcmp (x->szChrId, "Y")==0)
        nChrX = 24;
    else if ( strcmp (x->szChrId, "XY")==0)
        nChrX = 25;

    int nChrY = atoi(y->szChrId);
    if ( strcmp (y->szChrId, "M")==0)
        nChrY = 0;
    else if ( strcmp (y->szChrId, "X")==0)
        nChrY = 23;
    else if ( strcmp (y->szChrId, "Y")==0)
        nChrY = 24;
    else if ( strcmp (y->szChrId, "XY")==0)
        nChrY = 25;

    if (nChrX==nChrY)
    {
        int nPosX = atoi(x->szGenDist);
        int nPosY = atoi(y->szGenDist);

        return(nPosX<nPosY);
    }
    else
        return(nChrX<nChrY);
}

void CFmPackedSNP::quickSort(SNPINFO** items, int left, int right)
{
    int i = left;
    int j = right;
    SNPINFO* x = items[(left+right)/2];

    do
    {
        while(( IsLess(items[i], x) ) && (i < right))
            i++;
        while(( IsLess(x, items[j] ) ) && (j > left))
            j--;

        if(i <= j)
        {
            SNPINFO* y = items[i];
            items[i] = items[j];
            items[j] = y;

            double f = m_pSnpIdx->Get(i);
            m_pSnpIdx->Set( i, m_pSnpIdx->Get(j) );
            m_pSnpIdx->Set( j, f );

            i++;
            j--;
        }

    } while(i <= j);

    if(i < right)
        quickSort(items, i, right);

    if(left < j)
        quickSort(items, left, j);
}

bool CFmPackedSNP::OrderByChrPos()
{
    _log_info(_HI_, "CFmPackedSNP::OrderByChrPos: m_nSnpP=%d.", m_nSnpP);

    SNPINFO** ppSnpInfo = (SNPINFO**)Calloc( m_nSnpP, SNPINFO* );
    memcpy(ppSnpInfo, m_ppSnpInfo, sizeof(SNPINFO*)*m_nSnpP );

    quickSort(ppSnpInfo, 0, m_nSnpP-1);

    Free(ppSnpInfo);

    /*SNPINFO* pSnpInfo=NULL;
    for (int i=0; i<300; i++)
    {
        GetSnpInfo( i, &pSnpInfo );
        Rprintf("%s %s %s\n", pSnpInfo->szChrId, pSnpInfo->szGenDist, pSnpInfo->szChrId);
    }*/

    return(true);
}

int CFmPackedSNP::WriteAsCSVFile(const char* filename)
{
    int nSubCol = m_nSubjN;
    int nSnpRow = m_nSnpP;
	if (nSubCol>500) nSubCol=500;
	if (nSnpRow>500) nSnpRow=500;

    CFmMatrix mat(nSnpRow, nSubCol);
	double* pFval = new double[m_nSubjN];

	for( int i=0;i<	nSnpRow; i++)
	{
		GetSnpRow(i, pFval);
		SNPINFO* pInfo=NULL;
		GetSnpInfo( i, &pInfo );

		mat.SetRow(i, pFval, pInfo->szSnpId);
    }

    int ret = mat.WriteAsCSVFile(filename);

    return(ret);
}

void destroy(CFmPackedSNP* p)
{
	CFmNewTemp  fmRef;
	p->~CFmPackedSNP();
	operator delete(p, fmRef);
}
