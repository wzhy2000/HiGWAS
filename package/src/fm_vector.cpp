/* fm_vector.cpp  -	CFmVector Class
 *
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>
#include <Rdefines.h>

#include "fm_matrix.h"
#include "fm_vector.h"
#include "fm_rlogger.h"
#include "fm_err.h"
#include "fm_new.h"

int CFmVector::_DV=0;
CFmVector** CFmVector::g_pReused =NULL;
int CFmVector::g_nObjCount=0;

CFmVector::CFmVector(bool bReused)
{
	if (bReused)
		m_nObjId = g_nObjCount++;
	else
		m_nObjId = -1;

    m_pNames = NULL;
    m_pszBuf = NULL;
	m_bReuse = 	bReused;
	m_nRefer = 0;
	m_nActLen= 1;
	m_nMaxLen= 100;
	m_pData = AllocateMemory( m_nMaxLen );
	AddRef();
}

CFmVector::CFmVector(CFmMatrix* pMat, int nRow_or_Cols, bool bRow, bool bReused)
{
	if (bReused)
		m_nObjId = g_nObjCount++;
	else
		m_nObjId = -1;

    m_pNames = NULL;
    m_pszBuf = NULL;
    m_bReuse = 	bReused;
	m_nRefer = 0;
	if (bRow)
		m_nActLen= pMat->GetNumCols();
	else
		m_nActLen= pMat->GetNumRows();

	m_nMaxLen= m_nActLen>100?m_nActLen:100;
	m_pData = AllocateMemory( m_nMaxLen );

	for(int i=0; i<m_nActLen; i++)
		if (bRow)
			m_pData[i] = pMat->Get( nRow_or_Cols, i);
		else
			m_pData[i] = pMat->Get( i, nRow_or_Cols);
	AddRef();
}

CFmVector::CFmVector(int nSize, double fInit, int nMaxSize, bool bReused)
{
	if (bReused)
		m_nObjId = g_nObjCount++;
	else
		m_nObjId = -1;

    m_pNames = NULL;
    m_pszBuf = NULL;
    m_bReuse = 	bReused;
	m_nRefer = 0;
	m_nActLen= nSize;
	if (nMaxSize<=0)
		m_nMaxLen= m_nActLen>100?m_nActLen:100;
	else
		m_nMaxLen= nMaxSize;

	m_pData = AllocateMemory( m_nMaxLen );

	for(int i=0; i<m_nActLen; i++)
		m_pData[i]=	fInit;

	AddRef();
}

CFmVector::CFmVector(CFmVector* pVct, bool bReused )
{
    if (bReused)
        m_nObjId = g_nObjCount++;
    else
        m_nObjId = -1;

    m_pNames = NULL;
    m_pszBuf = NULL;
    m_bReuse = 	bReused;
    m_nRefer = 0;
    m_nActLen= pVct->m_nActLen;
    m_nMaxLen= pVct->m_nMaxLen;
    m_pData = AllocateMemory( m_nMaxLen );

    for(int i=0; i<m_nActLen; i++)
        m_pData[i] = pVct->m_pData[i];

    AddRef();
}

CFmVector::CFmVector(const CFmVector& pOther)
{
    Rprintf("Copy & construct matrix.............\n");
    throw ("Copy constructor is not good for CFmMatrix");
}

CFmVector::~CFmVector()
{
    if (m_pszBuf) Free(m_pszBuf);
    m_pszBuf = NULL;

 //Rprintf("Destructor Vector(%X)\n", this);

	Release();
}

bool CFmVector::SetLength(int nLen)
{
	if (nLen <= m_nMaxLen && nLen>=0)
	{
		m_nActLen = nLen;
		return (true);
	}
	else
		return false;
}

bool CFmVector::operator==(CFmVector &other)
{
	if (_DV) Rprintf(" OP  b = Ve == Vo..\n");
	if (m_nActLen != other.m_nActLen)
        throw( "Vectors do not have common size in op(==)" );

	int ret = memcmp(m_pData, other.m_pData, sizeof(double) * m_nActLen ) ;

	if (IsReusable()) Release();
	if (other.IsReusable()) other.Release();

	if (ret== 0)
		return true ;
	else
		return false ;
}

CFmVector& CFmVector::operator+(double other)
{
	if (_DV) Rprintf(" OP  Ve + v..\n");

    CFmVector* pTmp = FindReuseVector( m_nActLen );

	int nLen = m_nActLen;
	if (m_nActLen==0)
		nLen = m_nMaxLen;

	for(int i=0; i<nLen; i++)
		pTmp->m_pData[i] = m_pData[i] + other;

	if (IsReusable()) Release();

	return(*pTmp);
}

CFmVector& CFmVector::operator-(double other)
{
	if (_DV) Rprintf(" OP  Ve - v..\n");

    CFmVector* pTmp = FindReuseVector( m_nActLen );

	int nLen = m_nActLen;
	if (m_nActLen==0)
		nLen = m_nMaxLen;

	for(int i=0; i<nLen; i++)
		pTmp->m_pData[i] = m_pData[i] - other;

	if (IsReusable()) Release();

	return(*pTmp);
}

CFmVector& CFmVector::operator*(double other)
{
	if (_DV) Rprintf(" OP  Ve * v..\n");

    CFmVector* pTmp = FindReuseVector( m_nActLen );

	int nLen = m_nActLen;
	if (m_nActLen==0)
		nLen = m_nMaxLen;

	for(int i=0; i<nLen; i++)
		pTmp->m_pData[i] = m_pData[i] * other;

	if (IsReusable()) Release();

	return(*pTmp);
}

CFmVector& CFmVector::operator/(double other)
{
	if (_DV) Rprintf(" OP  Ve / v..\n");

    CFmVector* pTmp = FindReuseVector( m_nActLen );

	int nLen = m_nActLen;
	if (m_nActLen==0)
		nLen = m_nMaxLen;

	for(int i=0; i<nLen; i++)
		pTmp->m_pData[i] = m_pData[i] / other;

	if (IsReusable()) Release();

	return(*pTmp);
}

CFmVector& CFmVector::operator^(double other)
{
	if (_DV) Rprintf(" OP  Ve ^ v..\n");

    CFmVector* pTmp = FindReuseVector( m_nActLen );

	int nLen = m_nActLen;
	if (m_nActLen==0)
		nLen = m_nMaxLen;

	for(int i=0; i<nLen; i++)
		pTmp->m_pData[i] = pow( m_pData[i] , other);

	if (IsReusable()) Release();

	return(*pTmp);
}

CFmVector& CFmVector::operator+(CFmVector& other)
{
	if (_DV) Rprintf(" OP  Ve + Vo..\n");

	int nLen=0;
	if (m_nActLen!=0)
	{
		if (m_nActLen != other.m_nActLen)
		{
            throw( "Vectors do not have common size in op(+)" );
		}

		nLen = m_nActLen;
	}
	else
		nLen = other.m_nActLen;

    CFmVector* pTmp = FindReuseVector( nLen );

	for(int i=0; i<nLen; i++)
		pTmp->m_pData[i] = m_pData[i] + other.m_pData[i];

	if (IsReusable()) Release();
	if (other.IsReusable()) other.Release();

	return(*pTmp);

}

CFmVector& CFmVector::operator-(CFmVector& other)
{
	if (_DV) Rprintf(" OP  Ve - Vo..\n");

	int nLen=0;
	if (m_nActLen!=0)
	{
		if (m_nActLen != other.m_nActLen)
		{
            throw( "Vectors do not have common size in op(-)" );
		}

		nLen = m_nActLen;
	}
	else
		nLen = other.m_nActLen;

    CFmVector* pTmp = FindReuseVector( nLen );

	for(int i=0; i<nLen; i++)
		pTmp->m_pData[i] = m_pData[i] - other.m_pData[i];

	if (IsReusable()) Release();
	if (other.IsReusable()) other.Release();

	return(*pTmp);
}

CFmVector& CFmVector::operator*(CFmVector& other)
{
	if (_DV) Rprintf(" OP  Ve * Vo..\n");


	int nLen=0;
	if (m_nActLen!=0)
	{
		if (m_nActLen != other.m_nActLen)
		{
            throw( "Vectors do not have common size in op(*)" );
		}

		nLen = m_nActLen;
	}
	else
		nLen = other.m_nActLen;

    CFmVector* pTmp = FindReuseVector( nLen );

	for(int i=0; i<nLen; i++)
		pTmp->m_pData[i] = m_pData[i] * other.m_pData[i];

	if (IsReusable()) Release();
	if (other.IsReusable()) other.Release();

	return(*pTmp);
}

CFmVector& CFmVector::operator/(CFmVector& other)
{
	if (_DV) Rprintf(" OP  Ve / Vo..\n");

	int nLen=0;
	if (m_nActLen!=0)
	{
		if (m_nActLen != other.m_nActLen)
		{
            throw( "Vectors do not have common size in op(/)" );
		}

		nLen = m_nActLen;
	}
	else
		nLen = other.m_nActLen;

    CFmVector* pTmp = FindReuseVector( nLen );

	for(int i=0; i<nLen; i++)
		pTmp->m_pData[i] = m_pData[i] / other.m_pData[i];

	if (IsReusable()) Release();
	if (other.IsReusable()) other.Release();

	return(*pTmp);
}

CFmMatrix& CFmVector::operator*(CFmMatrix& other)
{
	if (_DV) Rprintf(" OP  Ve * Mo..\n");

	// first check for a valid multiplication operation
	if (m_nActLen != other.GetNumRows() )
        throw( "Matrices do not have common size");

    CFmMatrix* pTmp = CFmMatrix::FindReuseMatrix( 1, other.GetNumCols() );

	double	 value ;
	for (int i = 0 ; i < other.GetNumCols() ; ++i)
	{
		value = 0.0 ;
		for (int k = 0 ; k < other.GetNumRows(); ++k)
			value += Get(k) * other.Get(k,i) ;

		pTmp->Set(0, i, value) ;
	}

	if ( other.IsReusable()) other.Release();
	if ( this->IsReusable()) this->Release();

	return *pTmp ;
}

void CFmVector::operator=(double other)
{
	if (_DV) Rprintf(" OP  Ve = v.\n");

    int nLen = m_nActLen;
    if (m_nActLen==0)
        nLen = m_nMaxLen;

    for(int i=0; i<nLen; i++)
		m_pData[i] = other;
}

void CFmVector::operator=(double* other)
{
	if (_DV) Rprintf(" OP  Ve = *pv.\n");

	int nLen = m_nActLen;
	if (m_nActLen==0)
		nLen = m_nMaxLen;

	for(int i=0; i<nLen; i++)
		m_pData[i] = other[i];
}

void CFmVector::operator=(CFmVector& other)
{
	if (_DV) Rprintf(" OP  Ve = Vo.\n");

	if (m_nMaxLen<other.GetLength())
    {
        double* pData = AllocateMemory(other.GetLength() + 100);
        FreeMemory(m_pData);
        m_pData = pData;
        m_nMaxLen = other.GetLength() + 100;
    }

	m_nActLen = other.GetLength();
	for(int i=0; i<m_nActLen; i++)
		m_pData[i] = other.m_pData[i];

	if (other.IsReusable()) other.Release();
}

void CFmVector::Set(CFmVector& other)
{
    if (_DV) Rprintf(" OP  Ve = Vo.\n");


    if (m_nMaxLen<other.GetLength())
    {
        Rprintf(" %x, OP  Ve = Vo%d, %d\n",this, m_nMaxLen, other.GetLength() );

        double* pData = AllocateMemory(other.GetLength() + 100);
        FreeMemory( m_pData );
        m_pData = pData;
        m_nMaxLen = other.GetLength() + 100;
    }

    m_nActLen = other.GetLength();
    for(int i=0; i<m_nActLen; i++)
        m_pData[i] = other.m_pData[i];

    if (other.IsReusable()) other.Release();
}

CFmVector& CFmVector::GetTransposed()
{
    if (_DV) Rprintf(" OP  t(Ve).\n");

    CFmVector* pTmp = FindReuseVector(m_nActLen );
    for (int i=0; i<m_nActLen; i++)
        pTmp->Set( i, m_pData[i]);

    if (IsReusable()) Release();

    return(*pTmp);
}

CFmVector& CFmVector::abs()
{
    if (_DV) Rprintf(" OP  abs(Ve).\n");

    CFmVector* pTmp = FindReuseVector(m_nActLen );
    for (int i=0; i<m_nActLen; i++)
        pTmp->Set( i, fabs(m_pData[i]) );

    if (IsReusable()) Release();

    return(*pTmp);
}

double CFmVector::Sum()
{
	double sum = 0;
	for(int i=0; i<m_nActLen; i++)
		sum += m_pData[i];

	if (IsReusable()) Release();

	return(sum);
}

double CFmVector::GetMax()
{
    double max = m_pData[0];
    for(int i=1; i<m_nActLen; i++)
        if ( m_pData[i] > max )
            max = m_pData[i];

    if (IsReusable()) Release();

    return(max);
}

double CFmVector::GetMin()
{
    double min = m_pData[0];
    for(int i=1; i<m_nActLen; i++)
        if ( m_pData[i] < min )
            min = m_pData[i];

    if (IsReusable()) Release();

    return(min);
}


double CFmVector::GetMean()
{
    double sum = 0;
    for(int i=0; i<m_nActLen; i++)
        sum += m_pData[i];

    if (IsReusable()) Release();

    return(sum/m_nActLen);
}

double CFmVector::GetVar()
{
    double mean = GetMean();
    double sum = 0;
    for(int i=0; i<m_nActLen; i++)
        sum += (m_pData[i]-mean)*(m_pData[i]-mean);

    if (IsReusable()) Release();

    if (m_nActLen>1)
        return(sum/(m_nActLen-1));
    else
        return(0);
}

double CFmVector::Prod(CFmVector& other)
{
    if (_DV) Rprintf(" OP  Ve %*% Vo.\n");

    if (m_nActLen<other.GetLength())
        throw( "unequal length \n");

    double sum=0;
    for(int i=0; i<m_nActLen; i++)
		sum += (m_pData[i]*other.m_pData[i]);

    if (other.IsReusable()) other.Release();
    if (this != &other)
        if (IsReusable()) Release();
    return(sum);
}

CFmVector& CFmVector::GetReciprocal()
{
    if (_DV) Rprintf(" OP  recip.\n");

    CFmVector* pTmp = FindReuseVector(m_nActLen );
    for (int i=0; i<m_nActLen; i++)
        pTmp->Set( i, 1/m_pData[i]);

    if (IsReusable()) Release();

    return(*pTmp);
}

void CFmVector::UniquePut(double fValue)
{
    for(int i=0; i<m_nActLen; i++)
    {
        if ( m_pData[i] == fValue )
            return;
    }

    m_nActLen += 1;
    m_pData[m_nActLen-1] =fValue;
}

void CFmVector::Put(double fValue)
{
    if (m_nActLen+ 1>m_nMaxLen)
    {
        m_nMaxLen += 100;
        double* pData = AllocateMemory( m_nMaxLen );
        memcpy( pData, m_pData, m_nActLen*sizeof(double));
        FreeMemory( m_pData );
        m_pData=pData;
    }

    m_nActLen += 1;
    m_pData[m_nActLen-1] =fValue;
}

void CFmVector::Append(CFmVector& other)
{
    if (m_nActLen+ other.GetLength() > m_nMaxLen)
    {
        m_nMaxLen += (other.GetLength()+100);
        double* pData = AllocateMemory( m_nMaxLen );
        memcpy( pData, m_pData, m_nActLen*sizeof(double));
        FreeMemory( m_pData );
        m_pData=pData;
    }

    for (int i=0; i<other.GetLength(); i++)
    {
        m_nActLen += 1;
        m_pData[m_nActLen-1] =other[i];
    }

    if (other.IsReusable())
        other.Release();
}

void CFmVector::Remove( int nPos )
{
	if( nPos<0 && nPos>=m_nActLen)
		return;

	m_nActLen--;
	if( nPos == m_nActLen)
		return;

	for(int i=nPos; i<m_nActLen; i++)
	{
		m_pData[i] = m_pData[i+1];
	}
}

void CFmVector::RemoveNan()
{
        int f_cnt=0;
        for(int i=0, nLast = 0; i<m_nActLen; i++)
        {
                if ( isnan( m_pData[i]) )
                {
                        //skip
                        f_cnt++;
                }
                else
                {
                        if (nLast <i)
                                m_pData[nLast] = m_pData[i];

                        nLast++;
                }
        }

        m_nActLen -= f_cnt;
}

int CFmVector::GetLengthNonNan()
{
    int cnt = 0;
    for(int i=0; i<m_nActLen; i++)
        if ( !isnan(m_pData[i]) )
            cnt++;

    if (IsReusable()) Release();

    return(cnt);
}

int CFmVector::GetLengthNonvalue(double fvalue)
{
    int cnt = 0;
    for(int i=0; i<m_nActLen; i++)
        if ( m_pData[i]!=fvalue ) cnt++;

    if (IsReusable()) Release();

    return(cnt);
}

int CFmVector::Find(double fVal)
{
    int nPos = -1;
    for(int i=0; i<m_nActLen; i++)
    {
        if (m_pData[i] == fVal )
        {
            nPos=i;
            break;
        }
    }

    if (IsReusable()) Release();

    return( nPos );
}


void CFmVector::AddRef()
{
	m_nRefer++;
}

void CFmVector::Release()
{
	m_nRefer--;
	if ( m_nRefer == 0)
	{
		if (!m_bReuse)
		{
//Rprintf("DOUBLEs Freed:%X!\n", m_pData);
			FreeMemory(m_pData);
			m_pData = NULL ;
		}
		else
		{
			if (_DV) Rprintf("***RECYCLE VECTOR(%d):\n", m_nObjId);
		}
	}
}

void CFmVector::FreeMemory( double* p )
{
	Free(p);
}

double* CFmVector::AllocateMemory( int nLen )
{
    double	*pData = Calloc( nLen + 1, double);
//Rprintf("%d DOUBLEs allocated(%X)!\n", nLen + 1, this);

//**MEMTEST:fprintf(stderr,"ALLOC MEMORY(CFmVector): %d[%d], %s\n", (nLen + 1)*sizeof(double),nLen, m_bReuse?"YES":"NO" );

	if (pData==NULL)
	{
        _log_fatal( _HI_ , "MEMORY: failed to allocate %d bytes to CFmVector.", (nLen+1)*4);
	}

	memset(pData, 0, sizeof(double) * (nLen + 1)) ;
	return pData ;
}

bool CFmVector::IsFitted(int nLen)
{
	return( m_nMaxLen >= nLen );
}

void CFmVector::Resize(int nLen, bool reset)
{
	if (m_nMaxLen < nLen )
    {
        FreeMemory( m_pData );
        m_nMaxLen = nLen+100;
        m_pData = AllocateMemory(m_nMaxLen);
    }

	m_nActLen = nLen;
    if(reset)
        memset( m_pData, 0, sizeof(double)*m_nMaxLen);
}


#define REUSE_VECTOR_COUNT	1000
CFmVector* CFmVector::FindReuseVector( int nLen)
{
	if (g_pReused==NULL)
	{
        g_pReused = (CFmVector**)Calloc( REUSE_VECTOR_COUNT,  CFmVector* );
        memset(g_pReused, 0, sizeof(CFmVector*)*REUSE_VECTOR_COUNT );
	}

	for(int i=0; i<REUSE_VECTOR_COUNT; i++)
	{
        CFmVector* pTmp = (CFmVector* )(g_pReused[i]);
		if ( pTmp!=NULL &&
			 pTmp->GetRefer()==0 &&
			 pTmp->IsFitted(nLen))
		{
			pTmp->Resize(nLen);
			pTmp->AddRef();

			if (_DV) Rprintf("***REUSE VECTOR(%d)\n", pTmp->m_nObjId);
			return(pTmp);
		}
	}

	CFmNewTemp fmRef;
    CFmVector* pVct = new (fmRef) CFmVector(nLen, 0.0, -1, true);
	if (_DV) Rprintf("***NEW VECTOR(%d): Len=%d\n", pVct->m_nObjId, nLen);
	if ( pVct->m_nObjId > 5000 )
		_DV=1;
	if ( pVct->m_nObjId > 5000*2 )
        throw( "Resue vectors have big problem\n...");

	for(int i=0; i<REUSE_VECTOR_COUNT; i++)
	{
        CFmVector* pTmp = (CFmVector* )(g_pReused[i]);
		if ( pTmp==NULL)
		{
			g_pReused[i] = pVct;
			break;
		}
	}

	return(pVct);
}

bool CFmVector::RemoveElements(CFmVector& nRows)
{
    for (int i=0; i<nRows.GetLength(); i++)
    {
        if ( nRows.Get(i)<0 || nRows.Get(i)>=m_nActLen )
            return(false);
    }

    int k=0;
    for (int i=0; i<m_nActLen; i++)
    {
        if (nRows.Find(i)>=0)
        {
        }
        else
        {
            m_pData[k] = m_pData[i];
            k++;
        }
    }

    for(;k<m_nActLen;k++)
        m_pData[k] = 0;

    m_nActLen -= nRows.GetLength();
    return(true);
}

int CFmVector::WriteAsCSVFile(const char* szCsvFile, bool bAppend, const char* szTag)
{
    FILE* fp = NULL;
    if (!bAppend)
        fp = fopen( szCsvFile, "wt");
    else
        fp = fopen( szCsvFile, "a+");

    if (fp==NULL)
    {
        _log_error( _HI_ , "The file can not be created(%s)", szCsvFile);
        return(ERR_CREATE_FILE);
    }

    if (szTag!=NULL)
    {
        fprintf(fp, "\n>>-----------------%s--------------------<<\n", szTag);
    }

    for (int i = 0; i < m_nActLen ; i++)
    {
        fprintf(fp, "%.8f", Get(i)  );
        fprintf(fp, "\n");
    }

    fclose(fp);
    return(0);
}

void CFmVector::Sort(bool bDecreasing)
{
    bool bMoved = true;
    while(bMoved)
    {
        bMoved = false;
        if (bDecreasing)
        {
            for(int i=0;i<m_nActLen-1;i++)
            {
                if ( m_pData[ i ] < m_pData[ i+1 ] )
                {
                    double d = m_pData[ i ];
                    m_pData[ i ] = m_pData[i+1];
                    m_pData[i+1] = d;
                    bMoved = true;
                }
            }
        }
        else
        {
            for(int i=0;i<m_nActLen-1;i++)
            {
                if ( m_pData[ i ] > m_pData[ i+1 ] )
                {
                    double d = m_pData[ i ];
                    m_pData[ i ] = m_pData[i+1];
                    m_pData[i+1] = d;
                    bMoved = true;
                }
            }
        }
    }
}

bool CFmVector::Order(CFmVector& order)
{
    CFmVector* pValue = FindReuseVector( m_nActLen );
    CFmVector* pIdx = FindReuseVector( m_nActLen );
    for(int i=0; i<m_nActLen; i++)
    {
        pValue->m_pData[i] = m_pData[i];
        pIdx->m_pData[i] = i;
    }

    bool bMoved = true;
    while(bMoved)
    {
        bMoved = false;
        for(int i=0;i<m_nActLen-1;i++)
        {
            if (pValue->Get(i) < pValue->Get(i+1) )
            {
                double d = pValue->Get(i);
                pValue->Set( i, pValue->Get(i+1) );
                pValue->Set( i+1, d );
                double n = pIdx->Get(i);
                pIdx->Set( i, pIdx->Get(i+1) );
                pIdx->Set( i+1, n );
                bMoved = true;
            }
        }
    }

    //pIdx->Show("pIdx");
    //pValue->Show("pValue");

    order.Resize(m_nActLen);
    for(int i=0; i<m_nActLen; i++)
        order[ (int)(pIdx->Get(i))] = (i+1);

    pValue->Release();
    pIdx->Release();

    return(true);
}

// order element>=1 not >=0!!!
bool CFmVector::Rearrange(CFmVector& order)
{
    CFmVector* pValue = FindReuseVector( m_nActLen );
    for(int i=0; i<m_nActLen; i++)
    {
        pValue->m_pData[i] = m_pData[i];
    }

    for(int i=0; i<m_nActLen; i++)
    {
        Set( (int)( order.Get(i) )-1 ,pValue->m_pData[i] );
    }

    pValue->Release();
    return(true);
}


double CFmVector::GetMedian()
{
    return GetMean();
}


bool CFmVector::IsNan(int nIdx)
{
    return(isnan(m_pData[nIdx])!=0);
}

char* CFmVector::GetCommaString(const char* szFormat)
{
    if (!m_pszBuf)
    {
        m_pszBuf = (char*)Calloc(1024+128, char);
        memset( m_pszBuf, 0, 1024+128);
    }
    else
        memset( m_pszBuf, 0, 1024+128);

    char szTempBuf[1024]={0};
    for (int i=0; i<m_nActLen; i++)
    {
        snprintf( szTempBuf, 512, szFormat, m_pData[i] );
        strcat(m_pszBuf, ", ");
        strncat(m_pszBuf, szTempBuf, 1024);

        if ( strlen(m_pszBuf) > 1024 )
        {
            strcat(m_pszBuf, "...");
            break;
        }
    }

    return (m_pszBuf+2);
}


void CFmVector::SetNames( CFmVectorStr* pNames)
{
    if ( m_pNames ) destroy( m_pNames );
    CFmNewTemp fmRef;
    m_pNames = new (fmRef) CFmVectorStr(pNames);
}

CFmVectorStr* CFmVector::GetNames()
{
    return m_pNames;
}

char* CFmVector::GetName(int idx )
{
    if (idx<0 ||idx>=m_nActLen)
        return NULL;

    if (m_pNames)
        return m_pNames->Get(idx);
    else
        return(NULL);
}


SEXP GetSEXP(CFmVector* pVct)
{
    SEXP sep;
    PROTECT( sep = allocVector(REALSXP, pVct->GetLength() ) );
    double* sep_r = REAL(sep);

    for (int i=0; i<pVct->GetLength(); i++)
        sep_r[i] = pVct->GetData()[i];

    if ( pVct->GetNames() )
    {
        SEXP exp_name;
        PROTECT(exp_name = allocVector(STRSXP, pVct->GetLength()));
        for (int i=0; i<pVct->GetLength(); i++)
            if (pVct->GetName(i))
                SET_STRING_ELT(exp_name, i, mkChar(pVct->GetName(i)));

        setAttrib(sep, R_NamesSymbol, exp_name);
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return sep;
}

int GetVector(SEXP pExp, CFmVector* pVct)
{
    pVct->Resize(0);
    double* sep_r = REAL(pExp);

    for (int i=0; i< GET_LENGTH(pExp); i++)
        pVct->Put(sep_r[i]);

    return(0);
}

CFmMatrix& CFmVector::GetRow()
{
	if (_DV) Rprintf(" GetRow.\n");

	// first check for a valid multiplication operation
    CFmMatrix* pTmp = CFmMatrix::FindReuseMatrix( 1, m_nActLen );

	double	 value ;
	for (int i = 0 ; i < m_nActLen ; ++i)
		pTmp->Set(0, i, m_pData[i] ) ;

	if ( this->IsReusable()) this->Release();

	return *pTmp ;
}

CFmMatrix& CFmVector::GetCol()
{
	if (_DV) Rprintf(" GetRow.\n");

	// first check for a valid multiplication operation
    CFmMatrix* pTmp = CFmMatrix::FindReuseMatrix( m_nActLen,1 );

	double	 value ;
	for (int i = 0 ; i < m_nActLen ; ++i)
		pTmp->Set(i, 0, m_pData[i] ) ;

	if ( this->IsReusable()) this->Release();

	return *pTmp ;
}

void CFmVector::StatCache(int* pnTotal, int* pnUsed)
{
	if (g_pReused==NULL)
	{
		*pnTotal = 0;
		*pnUsed = 0;
		return;
	}

	for(int i=0; i<REUSE_VECTOR_COUNT; i++)
	{
        CFmVector* pTmp = (CFmVector* )(g_pReused[i]);
		if ( pTmp==NULL )
			continue;

		*pnTotal = *pnTotal + 1;
		if( pTmp->GetRefer() > 0)
			*pnUsed = *pnUsed + 1;
	}

	return;
}

void destroy(CFmVector* p)
{
	CFmNewTemp fmRef;
	p->~CFmVector();
	operator delete(p, fmRef);
}
