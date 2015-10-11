/* matrix.cpp  -	CFmSnpMat Class
 *
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>

#include "fm_rlogger.h"
#include "fm_snpmat.h"
#include "fm_vector_str.h"
#include "fm_vector.h"
#include "fm_matrix.h"
#include "fm_err.h"
#include "fm_new.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
//#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction CFmSnpMat
//////////////////////////////////////////////////////////////////////

#ifdef _DEBUG
int CFmSnpMat::m_NextObjectNumber = 1 ;
#endif

//Normal:   [SnpP, SubjN]
//Transpose:[SnpN, SubjP]

CFmSnpMat::CFmSnpMat(int nSnps, int nSubjs, int nMaxSnps, int nMaxSubjs)
{
    m_nBaseSnps= nSnps;
    m_nNumSubjs = nSubjs;
    m_nNumSnps = nSnps;

    if (nMaxSnps>nSnps)
        m_nMaxSnps = nMaxSnps;
    else
        m_nMaxSnps = nSnps;

    if (nMaxSubjs>nSubjs)
        m_nMaxSubjs = nMaxSubjs;
    else
        m_nMaxSubjs = nSubjs;

    m_pData = AllocateChar(m_nMaxSnps, m_nMaxSubjs) ;

    m_pSnpInfos = NULL;
    m_pSubjInfos = NULL;

	m_bSD = FALSE;
	m_pSD = NULL;
}

CFmSnpMat::~CFmSnpMat()
{
    FreeMemory();
}

void CFmSnpMat::FreeMemory()
{
    Free( m_pData );
    m_pData = NULL ;

    if (m_pSnpInfos) destroy( m_pSnpInfos );
    if (m_pSubjInfos) destroy( m_pSubjInfos );
}

char* CFmSnpMat::AllocateChar( int nSnps, int nSubjs )
{
//**MEMTEST: Rprintf("ALLOC MEMORY: %d\n", (nSubjs * nSnps + 1)*sizeof(double));
    if (nSubjs * nSnps>1024*1024)
        _log_debug( _HI_ , "MEMORY: Try to allocate %.3fM bytes", (nSubjs * nSnps + 1)/1024.0/1024.0);

    char* pData = Calloc( nSubjs * nSnps + 1, char ) ;
	if (pData==NULL)
	{
        _log_fatal( _HI_ , "MEMORY: failed to allocate %d bytes to CFmSnpMat[%d,%d].", (nSubjs * nSnps + 1)*4, nSnps, nSubjs );
	}

    memset(pData, 0, sizeof(char) * (nSubjs * nSnps + 1)) ;
	return pData ;
}

//R: Me = Mo;
void CFmSnpMat::operator=( CFmSnpMat &other)
{
    FreeMemory();

    m_nNumSubjs = other.m_nNumSubjs;
    m_nNumSnps = other.m_nNumSnps;
    m_nBaseSnps = other.m_nBaseSnps;
    m_nMaxSnps = other.m_nMaxSnps;
    m_nMaxSubjs = other.m_nMaxSubjs;

    m_pData = AllocateChar(m_nMaxSnps, m_nMaxSubjs);
    for (int i = 0 ; i < m_nNumSnps ; i++)
        for (int j = 0 ; j <m_nNumSubjs ; j++)
            Set(i, j, other.Get_a(i, j)) ;

}

bool CFmSnpMat::SetSnpInfo(int nSnp, char* szSnpName, char* szChr, char* szPos)
{
    if (nSnp<0 || nSnp>= m_nNumSnps)
    {
        _log_error(_HI_, "MATRIX: failed to set Snpname for %d Snp", nSnp);
        return(false);
    }

	CFmNewTemp refNew;
    if (szSnpName)
    {
        if (m_pSnpInfos==NULL)
            m_pSnpInfos = new (refNew) CFmVectorStr( m_nNumSnps, m_nMaxSnps );

        if (m_nNumSnps>m_pSnpInfos->GetLength())
            m_pSnpInfos->Resize(m_nNumSnps);

        char szInfo[256];
        sprintf(szInfo, "%s/%s/%s", szSnpName, szChr, szPos);

        m_pSnpInfos->Set( nSnp, szInfo );
    }
    return(true);
}

bool CFmSnpMat::SetSubjName(int nSubj, char* szSubjName)
{
    if (nSubj<0 || nSubj>= m_nNumSubjs)
    {
        _log_error(_HI_, "MATRIX: failed to set Subjname for %d Subj", nSubj);
        return(false);
    }

	CFmNewTemp refNew;
    if (szSubjName)
    {
        if (m_pSubjInfos==NULL)
            m_pSubjInfos = new (refNew) CFmVectorStr( m_nNumSubjs, m_nMaxSubjs );

        if (m_nNumSubjs > m_pSubjInfos->GetLength())
            m_pSubjInfos->Resize(m_nNumSubjs);

        m_pSubjInfos->Set( nSubj, szSubjName);
    }

    return(true);
}

bool CFmSnpMat::GetSnpInfo(CFmVectorStr* pVctSnp, CFmVector* pVctChr, CFmVector* pVctPos )
{
    if (m_pSnpInfos==NULL)
        return false;

    for (int i=0; i<m_nNumSnps; i++)
    {
        char szChr[256]={""};
        char szPos[256]={""};
        char szSnpName[256]={""};

        char* szSnpInfo = m_pSnpInfos->Get( i );
        char* pos0 = strchr(szSnpInfo, '/');
        if ( pos0==NULL )
            return(false);
        if ( pos0 )
        {
			strncpy( szSnpName, szSnpInfo, (pos0-szSnpInfo) );
			szSnpName[pos0 - szSnpInfo] = '\0';
		}

        char* pos1 = strchr(pos0+1, '/');
        if ( pos1 )
        {
            strncpy( szChr, pos0+1, (pos1-pos0) );
			szChr[ pos1-pos0 ] = '\0';

            strcpy( szPos, pos1+1);
        }

        pVctSnp->Put( szSnpName );
        pVctChr->Put(atof(szChr));
        pVctPos->Put(atof(szPos));
    }

    return(true);
}

bool CFmSnpMat::GetSnpInfo(int nSnp, char* szSnpName, char* szChr, char* szPos)
{
    if (nSnp<0 || nSnp>= m_nNumSnps)
        return(false);
    else
    {
        if (m_pSnpInfos==NULL)
            return false;

        char* szSnpInfo = m_pSnpInfos->Get(nSnp);
        char* pos0 = strchr(szSnpInfo, '/');
        if (pos0==NULL)
            return(false);
        if (pos0 && szSnpName)
            strncpy( szSnpName, szSnpInfo, (pos0-szSnpInfo) );

        char* pos1 = strchr(pos0+1, '/');
        if ( pos1 )
        {
            if (szChr)
                strncpy( szChr, pos0+1, (pos1-pos0) );
            if (szPos)
                strcpy( szPos, pos1+1);
        };

        return(true);
    }
}

char* CFmSnpMat::GetSubjName(int nSubj)
{
    if (nSubj<0 || nSubj>= m_nNumSubjs)
        return(NULL);
    else
    {
        if (m_pSubjInfos==NULL)
            return NULL;
        return(m_pSubjInfos->Get(nSubj));
    }
}

CFmVector& CFmSnpMat::GetSnp(int nSnp)
{
    CFmVector* pTmp = CFmVector::FindReuseVector( m_nNumSubjs );

    for(int i=0; i<m_nNumSubjs; i++)
        pTmp->Set(i, Get_a(nSnp, i));

	return *pTmp;
}

CFmVector& CFmSnpMat::GetSubj(int nSubj)
{
    CFmVector* pTmp = CFmVector::FindReuseVector( m_nNumSnps );

    for(int i=0; i<m_nNumSnps; i++)
        pTmp->Set(i, Get_a(i, nSubj));

	return *pTmp;
}

void CFmSnpMat::SetSnp( int nSnp, CFmVector& vct, char* szSnpName, char* szChr, char* szPos )
{
    if ( nSnp >= m_nNumSnps)
        throw("over maximum Snps' size." );

    if (m_nNumSubjs != vct.GetLength() )
        throw( "Cannot concatenate matrices, not same size" );

    // now add the other matrix
    for (int j = 0 ; j <m_nNumSubjs; ++j)
        Set(nSnp, j, vct.Get(j)) ;

    SetSnpInfo(nSnp, szSnpName, szChr, szPos );

    if ( vct.IsReusable()) vct.Release();
}

void CFmSnpMat::SetSnp( int nSnp, double* fBuf , int nLen, char* szSnpName, char* szChr, char* szPos)
{
    if ( nSnp >= m_nNumSnps)
        throw("over maximum Snps' size." );

    // now add the other matrix
    for (int j = 0 ; j <m_nNumSubjs; ++j)
        Set(nSnp, j, fBuf[j] ) ;

    SetSnpInfo(nSnp, szSnpName, szChr, szPos );
}

void CFmSnpMat::SetSubj( int nSubj, CFmVector& vct , char* szSubjName)
{
    if ( nSubj >= m_nNumSubjs)
        throw("over maximum Subjumns' size." );

    if (m_nNumSnps != vct.GetLength() )
        throw( "Cannot concatenate matrices, not same size" );

    // now add the other matrix
    for (int i = 0 ; i <m_nNumSnps; i++)
        Set(i, nSubj, vct.Get(i)) ;

    SetSubjName(nSubj, szSubjName);

    if ( vct.IsReusable()) vct.Release();
}

void CFmSnpMat::SetSubj( int nSubj, double* fBuf, int nLen,char* szSubjName )
{
    if ( nSubj >= m_nNumSubjs)
        throw("over maximum Subjumns' size." );

    // now add the other matrix
    for (int i = 0 ; i <m_nNumSnps; i++)
        Set(i, nSubj, fBuf[i]) ;

    SetSubjName(nSubj, szSubjName);
}

int CFmSnpMat::WriteAsCSVFile(const char* filename)
{
    CFmMatrix* pMat = this->ToMatrix();
    int ret = pMat->WriteAsCSVFile(filename);
    destroy(  pMat );

    return(ret);
}

int CFmSnpMat::ReadFromCSVFile(const char* filename,  bool bSubjName, bool bSnpName)
{
	CFmNewTemp refNew;
    CFmMatrix* pMat = new (refNew) CFmMatrix(0,0);
    int ret = pMat->ReadFromCSVFile(filename, bSubjName, bSnpName);
    if (ret!=0)
        return(ret);

    this->FromMatrix(pMat);
    destroy(  pMat );

    return 0;
}

CFmMatrix* CFmSnpMat::ToMatrix()
{
	CFmNewTemp refNew;
    CFmMatrix* pMat = new (refNew) CFmMatrix(m_nNumSnps, m_nNumSubjs);
    for (int i=0;i<m_nNumSnps;i++)
        for (int j=0;j<m_nNumSubjs;j++)
            pMat->Set(i,j, Get_a(i,j)*1.0);

    if (m_pSnpInfos) pMat->SetRowNames(m_pSnpInfos);
    if (m_pSubjInfos) pMat->SetColNames(m_pSubjInfos);

    return(pMat);
}

CFmMatrix* CFmSnpMat::GetAdd()
{
    return(ToMatrix());
}

CFmMatrix* CFmSnpMat::GetDom()
{
	CFmNewTemp refNew;
    CFmMatrix* pMat = new (refNew) CFmMatrix(m_nNumSnps, m_nNumSubjs);
    for (int i=0;i<m_nNumSnps;i++)
        for (int j=0;j<m_nNumSubjs;j++)
            pMat->Set(i,j, 1-(int)abs(Get_a(i,j))*1.0);

    return(pMat);
}

void CFmSnpMat::FromMatrix(CFmMatrix* pMat)
{
    FreeMemory();

    m_nNumSubjs = pMat->GetNumCols();
    m_nMaxSubjs = pMat->GetNumCols();
    m_nNumSnps = pMat->GetNumRows();
    m_nBaseSnps = pMat->GetNumRows();
    m_nMaxSnps = pMat->GetNumRows();

    m_pData = AllocateChar(m_nMaxSnps, m_nMaxSubjs);
    for (int i = 0 ; i < m_nNumSnps ; i++)
        for (int j = 0 ; j <m_nNumSubjs ; j++)
            Set(i, j, pMat->Get(i, j)) ;

	CFmNewTemp refNew;
    m_pSnpInfos = new (refNew) CFmVectorStr(pMat->GetRowNames());
    m_pSubjInfos = new (refNew) CFmVectorStr(pMat->GetColNames());
}

void CFmSnpMat::sd()
{
	CFmNewTemp refNew;
	m_pSD = new (refNew) CFmVector( m_nNumSnps, 0 );
	CFmVector snpVct(m_nNumSnps);
	for (int i=0; i<m_nNumSnps; i++)
	{
		snpVct = GetSnp( i );

		float sd = sqrt( snpVct.GetVar() );
		m_pSD->Set(i, sd);
	}

	m_bSD = TRUE;
}

void destroy(CFmSnpMat* p)
{
	CFmNewTemp  fmRef;
	p->~CFmSnpMat();
	operator delete(p, fmRef);
}
