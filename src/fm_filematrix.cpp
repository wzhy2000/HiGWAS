#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "fm_filematrix.h"
#include "fm_vector.h"
#include "fm_matrix.h"
#include "fm_rlogger.h"
#include "fm_linux.h"
#include "fm_err.h"
#include "fm_new.h"


CFmFileMatrix::CFmFileMatrix(char* szFile, bool bAppend, bool bAutodel=TRUE)
{
	m_bAutodel = bAutodel;
    m_nNumCols = 0;
    m_nNumRows = 0;
    m_pFile = NULL;
    m_pszFile = Strdup(szFile);

    m_nCacheCol0 = -1;
    m_nCacheCol1 = -1;

    CFmNewTemp refNew;
	m_pmatCache = new (refNew) CFmMatrix(0,100);

    OpenFile(bAppend);
}

CFmFileMatrix::~CFmFileMatrix()
{
    _log_info( _HI_, "CFmFileMatrix[%d,%d], file=%s\n",
            m_nNumRows, m_nNumCols, m_pszFile);

    if (m_pFile)
    {
	    RewriteHeader();
	    fflush( m_pFile);
		fclose(m_pFile);
        if (m_bAutodel) unlink(m_pszFile);
    }

    if (m_pszFile)  Free(m_pszFile);
    if (m_pmatCache) destroy( m_pmatCache );
}

char* CFmFileMatrix::GetFileName()
{
	return(m_pszFile);
}

int CFmFileMatrix::OpenFile(bool bAppend)
{
    char szModel[10]={"r+b"};
    if (!bAppend)
        strcpy(szModel, "w+b");

    m_pFile = fopen( m_pszFile, szModel);
    if (m_pFile==NULL)
    {
        if (bAppend)
            m_pFile = fopen( m_pszFile, "wb");

        if (m_pFile==NULL)
            return (ERR_OPEN_FILE);
    }

    if (bAppend)
    {
        FMMAT_FILE_FMT fmt;
        size_t cnt = fread( &fmt, sizeof(FMMAT_FILE_FMT), 1, m_pFile );
        if (cnt<1)
            return ERR_FILE_DATA;

        m_nNumCols = fmt.nNumCols;
        m_nNumRows = fmt.nNumRows;
    }
    else
    {
        RewriteHeader();
    }

    return(0);
}

int CFmFileMatrix::RewriteHeader()
{
    FMMAT_FILE_FMT fmt;
    memset(&fmt, 0, sizeof(FMMAT_FILE_FMT));
    strcpy(fmt.szHeader, "FUNMAP FILE MATRIX.\n");
    fmt.nNumCols = m_nNumCols;
    fmt.nNumRows = m_nNumRows;

    if ( fseek( m_pFile, 0, SEEK_SET) )
        return ERR_FILE_OPERATE;

    size_t cnt = fwrite( &fmt, sizeof(FMMAT_FILE_FMT), 1, m_pFile );
    if (cnt<1)
        return ERR_FILE_DATA;

    return(0);
}


int CFmFileMatrix::GetNumCols()
{
    return(m_nNumCols);
}

int CFmFileMatrix::GetNumRows()
{
    return(m_nNumRows);
}

int CFmFileMatrix::Rbind(CFmVector& vct)
{
    if (m_nNumCols==0)
    {
        m_nNumCols = vct.GetLength();
        RewriteHeader();
    }

    if ( vct.GetLength() != m_nNumCols )
        return(ERR_DATA_MISPATH);

    if ( fseek( m_pFile, 0, SEEK_END) )
        return ERR_FILE_OPERATE;

    int nRet = fwrite ( vct.GetData(), 1, vct.GetLength()*sizeof(double), m_pFile );
    if (nRet<=0)
        return ERR_FILE_OPERATE;

    m_nNumRows++;
    if (m_nNumRows%10==0)
        fflush( m_pFile);

    return(0);
}

int CFmFileMatrix::GetRow(int nRow, CFmVector& vct)
{
    long nPos = sizeof(FMMAT_FILE_FMT) + m_nNumCols*sizeof(double)*nRow;
    if ( fseek(m_pFile, nPos, SEEK_SET) )
        return ERR_FILE_OPERATE;

    double* pDouble = new double[m_nNumCols];
    size_t rlen = fread( pDouble, sizeof(double), m_nNumCols, m_pFile );
    if(rlen <= 0)
     	return ERR_FILE_OPERATE;

    vct.Resize(0);
    for (int i=0;i<m_nNumCols;i++)
        vct.Put(pDouble[i]);

    delete [] pDouble;
    return(0);
}

int CFmFileMatrix::GetCol(int nCol, CFmVector& vct)
{
    vct.Resize(0);
    for(int i=0; i<m_nNumRows; i++)
    {
        long nPos = sizeof(FMMAT_FILE_FMT) + (m_nNumCols*i + nCol)*sizeof(double);
        if ( fseek(m_pFile, nPos, SEEK_SET) )
            return ERR_FILE_OPERATE;

        double fv = 0;
        size_t cnt = fread( &fv, sizeof(double), 1, m_pFile );
        if (cnt<=0)
            return ERR_FILE_OPERATE;

        vct.Put(fv);
    }

    return(0);
}

int CFmFileMatrix::GetCacheCol(int nCol, CFmVector& vct)
{
	if( nCol>=m_nNumCols ) return(ERR_NULL_DATA);

    vct.Resize(0);
    if( nCol>=m_nCacheCol0 && nCol<m_nCacheCol1)
    {
		vct = m_pmatCache->GetCol(nCol-m_nCacheCol0);
		return(0);
	}

    fflush( m_pFile);

	m_pmatCache->Resize(m_nNumRows, 100);
	m_nCacheCol0 = (int)floor(nCol/100.0)*100;
	m_nCacheCol1 = (int)floor(nCol/100.0)*100+100;
	if( m_nCacheCol1> m_nNumCols ) m_nCacheCol1=m_nNumCols;

    for(int i=0; i<m_nNumRows; i++)
    {
        long nPos = sizeof(FMMAT_FILE_FMT) + (m_nNumCols*i + m_nCacheCol0)*sizeof(double);
        if ( fseek(m_pFile, nPos, SEEK_SET) )
            return ERR_FILE_OPERATE;

        double fv[100]={0.0};
        size_t cnt = fread( &fv, sizeof(double)*(m_nCacheCol1-m_nCacheCol0), 1, m_pFile );
        if (cnt<=0)
            return ERR_FILE_OPERATE;

		m_pmatCache->SetRow(i, fv );
    }

	vct = m_pmatCache->GetCol( nCol-m_nCacheCol0 );

    return(0);
}


double CFmFileMatrix::GetAt(int nRow, int nCol)
{
    long nPos = sizeof(FMMAT_FILE_FMT) + (m_nNumCols*nRow + nCol)*sizeof(double);
    if ( fseek(m_pFile, nPos, SEEK_SET) )
        return ERR_FILE_OPERATE;

    double fv = 0;
    size_t cnt = fread( &fv, sizeof(double), 1, m_pFile );
    if (cnt<=0)
        return ERR_FILE_OPERATE;

    return(fv);
}

int CFmFileMatrix::SetAt(int nRow, int nCol, double fval)
{
    long nPos = sizeof(FMMAT_FILE_FMT) + (m_nNumCols*nRow + nCol)*sizeof(double);
    if ( fseek(m_pFile, nPos, SEEK_SET))
        return ERR_FILE_OPERATE;

    size_t cnt = fwrite( &fval, sizeof(double), 1, m_pFile );
    if (cnt<=0)
        return ERR_FILE_OPERATE;

    return(0);
}

CFmMatrix* LoadFileMatrix(char* szFile)
{
	FILE* pFile = fopen( szFile, "rb");
    if (pFile==NULL)
    {
        _log_info(_HI_, "LoadFileMatrix: can not open the file of FileMatrix(%s).", szFile );
        return( NULL );
    }

    FMMAT_FILE_FMT fmt;
    size_t cnt = fread( &fmt, sizeof(FMMAT_FILE_FMT), 1, pFile );
    if (cnt<1)
    {
        _log_info(_HI_, "LoadFileMatrix: no data in the header of Filematrix.");
        return( NULL );
    }

    double* pBuf = Calloc( fmt.nNumCols, double );
	memset( pBuf, 0, sizeof(double)*fmt.nNumCols );

	CFmNewTemp refNew;
    CFmMatrix* pMat = new (refNew) CFmMatrix(0, 0);
	for (long int j=0; j<fmt.nNumRows; j++)
	{
		size_t nSize = fread(pBuf, sizeof(double)*fmt.nNumCols, 1, pFile);
		if (nSize != 1)
		{
			Rprintf("nSize=%d nPos=%d\n", nSize,  sizeof(double)*fmt.nNumCols*j );
			return( NULL );
		}

		pMat->Cbind( pBuf, fmt.nNumCols );
	}

    fclose(pFile);
	Free(pBuf);

	return(pMat);
}

void destroy(CFmFileMatrix* p)
{
	CFmNewTemp  fmRef;
	p->~CFmFileMatrix();
	operator delete(p, fmRef);
}
