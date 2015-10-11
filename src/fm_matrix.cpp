/* matrix.cpp  -	CFmMatrix Class
 *
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>

#include "fm_matrix.h"
#include "fm_rlogger.h"
#include "fm_err.h"
#include "fm_new.h"

#define MAX_MARGIN_COL  20
#define MAX_MARGIN_ROW  20

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
//#define new DEBUG_NEW
#endif

#define Strdup(X) strcpy(Calloc(strlen(X)+1, char), X) ;
#define MI(rows, cols, nrow, ncol) (nrow+ncol*rows)

//std::list<CFmMatrix*> CFmMatrix::g_matList;
CFmMatrix** CFmMatrix::g_pReused =NULL;
int CFmMatrix::g_nObjCount=0;

int CFmMatrix::_DM=0;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction CFmMatrix
//////////////////////////////////////////////////////////////////////

// default constructor, create a 1 * 1 array
CFmMatrix::CFmMatrix( bool bReused )
{
	if (bReused)
    {
        m_nObjId = g_nObjCount++;
    }
	else
		m_nObjId = -1;
#ifdef _DEBUG
    TRACE("Creating CFmMatrix object %1d - default constructor\n", m_nObjId) ;
#endif

    m_pId = NULL;
    m_nRefer = 0;
	m_bReuse = bReused;
	m_nNumCols = 1 ;
	m_nNumRows = 1 ;
    m_nMaxRows = MAX_MARGIN_ROW ;
    m_nMaxCols = MAX_MARGIN_COL ;

    m_pData = AllocateDouble(m_nMaxRows, m_nMaxCols) ;

    m_pRowNames = NULL;
    m_pColNames = NULL;

    AddRef();
}

CFmMatrix::CFmMatrix(const CFmMatrix& pOther)
{
    Rprintf("Copy & construct matrix.............\n");
    throw ("Copy constructor is not good for CFmMatrix");
}

// copy other's data
CFmMatrix::CFmMatrix(const CFmMatrix *pOther, bool bReused )
{
	if (bReused)
    {
        m_nObjId = g_nObjCount++;
    }
    else
		m_nObjId = -1;
#ifdef _DEBUG
    TRACE("Creating CFmMatrix object %1d - default constructor\n", m_nObjId) ;
#endif


    m_pId = NULL;
    m_nRefer = 0;
	m_bReuse = bReused;
    m_nNumCols = pOther->m_nNumCols;
    m_nNumRows = pOther->m_nNumRows;
    m_nMaxRows = pOther->m_nMaxRows;
    m_nMaxCols = pOther->m_nMaxCols;

    m_pData = AllocateDouble(m_nMaxRows, m_nMaxCols) ;
    memcpy(m_pData, pOther->m_pData, sizeof(double)*m_nNumCols*m_nNumRows );

    m_pRowNames = NULL;
    m_pColNames = NULL;

    CFmNewTemp fmRef;
    if( pOther->m_pRowNames) m_pRowNames = new (fmRef) CFmVectorStr( pOther->m_pRowNames );
    if( pOther->m_pColNames) m_pColNames = new (fmRef) CFmVectorStr( pOther->m_pColNames );

    AddRef() ;
}

CFmMatrix::CFmMatrix(int nRows, int nCols, bool bReused )
{
	if (bReused)
    {
        m_nObjId = g_nObjCount++;
    }
    else
		m_nObjId = -1;
#ifdef _DEBUG
    TRACE("Creating CFmMatrix object %1d - default constructor\n", m_nObjId) ;
#endif

    m_pId = NULL;
    m_nRefer = 0;
	m_bReuse = bReused;
	m_nNumCols = nCols ;
	m_nNumRows = nRows ;
    m_nMaxCols = nCols>=MAX_MARGIN_COL?nCols:MAX_MARGIN_COL;
    m_nMaxRows = nRows>=MAX_MARGIN_ROW?nRows:MAX_MARGIN_ROW;

    m_pData = AllocateDouble(m_nMaxRows, m_nMaxCols) ;

    m_pRowNames = NULL;
    m_pColNames = NULL;

	AddRef() ;
}

// construct a square matrix with 1.0's on the diagonal if required
CFmMatrix::CFmMatrix(int size, bool set_diagonal, double init_val, bool bReused )
{
	if (bReused)
    {
        m_nObjId = g_nObjCount++;
    }
    else
		m_nObjId = -1;
#ifdef _DEBUG
    TRACE("Creating CFmMatrix object %1d - default constructor\n", m_nObjId) ;
#endif

    m_pId = NULL;
    m_nRefer = 0;
	m_bReuse = bReused;
	m_nNumCols = size ;
	m_nNumRows = size ;
    m_nMaxCols = size>=MAX_MARGIN_COL?size:MAX_MARGIN_COL;
    m_nMaxRows = size>=MAX_MARGIN_ROW?size:MAX_MARGIN_ROW;

    m_pData = AllocateDouble(m_nMaxRows, m_nMaxCols) ;

    m_pRowNames = NULL;
    m_pColNames = NULL;

	AddRef();

	// set the dialognal if required
	if (set_diagonal)
		for (int i = 0 ; i < size ; ++i)
			Set(i, i, init_val) ;
}

CFmMatrix::CFmMatrix(int nRows, int nCols, int nMaxRows, int nMaxCols, bool bReused )
{
	if (bReused)
    {
        m_nObjId = g_nObjCount++;
    }
    else
		m_nObjId = -1;

#ifdef _DEBUG
    TRACE("Creating CFmMatrix object %1d - default constructor\n", m_nObjId) ;
#endif

    m_pId = NULL;
    m_nRefer = 0;
	m_bReuse = bReused;
	m_nNumCols = nRows ;
	m_nNumRows = nCols ;
	m_nMaxRows = nMaxRows;
	m_nMaxCols = nMaxCols;

    m_pData = AllocateDouble(m_nMaxRows, m_nMaxCols) ;

    m_pRowNames = NULL;
    m_pColNames = NULL;

	AddRef();
}

CFmMatrix::~CFmMatrix()
{
#ifdef _DEBUG
    TRACE("Destroying CFmMatrix object %1d\n", m_nObjId) ;
#endif
//    Rprintf("Destroying CFmMatrix object %X\n", this ) ;

	Release() ;
}

void CFmMatrix::SetId(const char* szId)
{
    if (m_pId) Free(m_pId);
    m_pId = Strdup(szId);
}

void CFmMatrix::FreeMemory()
{
    if (m_pId) Free(m_pId);
    m_pId = NULL;

    FreeDouble(m_pData);
    m_pData = NULL ;

    if (m_pRowNames) destroy( m_pRowNames);
    if (m_pColNames) destroy( m_pColNames);
}

void CFmMatrix::FreeDouble( double* pData )
{
//Rprintf("%X DOUBLEs Freed!\n", pData);
    Free(pData);
}

#ifdef _DEBUG
void CFmMatrix::Dump(CDumpContext& dc) const
{
	UNUSED_PARAMETER(dc) ;
    TRACE("CFmMatrix object    #%1d\n", m_nObjId) ;
	TRACE("Num columns     : %1d\n", m_nNumCols) ;
	TRACE("Num rows        : %1d\n", m_nNumRows) ;
	TRACE("Data pointer    : %lx\n", m_pData) ;
	TRACE("Reference count : %1d\n", GetReferenceCount()) ;
	for (int i = 0 ; i < m_nNumRows ; ++i)
	{
		TRACE("Row %2d,", i) ;
		for (int j = 0 ; j < m_nNumCols ; ++j)
		{
			TRACE("%e,", Get(j, i)) ;
		}

		TRACE("\n") ;
		// this is to allow all element data to be traced for very large matrices!
		Sleep(m_nNumCols *2) ;
	}
}
#endif

void CFmMatrix::Show(const char* szName)
{
    Rprintf("\n\n\n%s: Resuable=%d Refer=%d\n",szName,m_bReuse, m_nRefer );
	if( m_pColNames !=NULL)
	{
		Rprintf("\nColNames:\t");
		for(int i=0; i<m_pColNames->GetLength(); i++)
			Rprintf("%s\t",m_pColNames->Get(i) );
		Rprintf("\n");
	}
	if( m_pRowNames !=NULL)
	{
		Rprintf("\nRowNames:\t");
		for(int i=0; i<m_pRowNames->GetLength(); i++)
			Rprintf("%s\t",m_pRowNames->Get(i) );
		Rprintf("\n");
	}
    show_matrix2( m_pData, m_nNumRows, m_nNumCols, szName );
}

double* CFmMatrix::AllocateDouble( int nRows, int nCols )
{
    double	*pData = NULL;
//*MEMTEST:Rprintf("ALLOC MEMORY(CFmMatrix): %d[%d,%d]\n", (nCols * nRows + 1)*sizeof(double), nRows, nCols );

    if( nCols * nRows * sizeof(double)>1024*1024)
        _log_debug( _HI_ , "MEMORY: Try to allocate big memory(%.3fM bytes) for %s", (nCols * nRows + 1)*sizeof(double)/1024.0/1024.0, m_pId);

    pData = Calloc( nCols * nRows + 1, double);
//Rprintf("%d DOUBLEs allocated(%X)!\n",  nCols * nRows + 1, pData);

	if (pData==NULL)
	{
        _log_fatal( _HI_ , "MEMORY: failed to allocate %d bytes to CFmMatrix[%d,%d] for %s.", (nCols * nRows + 1)*sizeof(double), nRows, nCols, m_pId );
	}

	memset(pData, 0, sizeof(double) * (nCols * nRows + 1)) ;
	return pData ;
}

//R: r = Me + Mo;
CFmMatrix& CFmMatrix::operator+(CFmMatrix &other)
{
	if (_DM) Rprintf("op:  r = Me + Mo;\n");

	// first check for a valid addition operation
	if (m_nNumCols !=0 && m_nNumRows!=0)
	{
		if (m_nNumCols != other.m_nNumCols)
			throw("Invalid operation in op(+)") ;

		if (m_nNumRows != other.m_nNumRows)
			throw("Invalid operation in op(+)") ;
	}
	else
	{
		if (other.m_nNumRows * other.m_nNumCols>m_nMaxRows*m_nMaxCols)
			throw("Insufficient memory.") ;

		m_nNumRows = other.m_nNumRows;
		m_nNumCols = other.m_nNumCols;
	}

    CFmMatrix* pTmp = FindReuseMatrix( m_nNumRows, m_nNumCols );

	for (int i = 0 ; i < m_nNumRows; ++i)
		for (int j = 0 ; j < m_nNumCols ; ++j)
			pTmp->Set(i, j, Get(i, j) + other.Get(i, j)) ;

	if ( other.IsReusable()) other.Release();
	if ( IsReusable()) Release();

	return *pTmp ;
}

//R: r = Me - Mo;
CFmMatrix& CFmMatrix::operator-( CFmMatrix &other)
{
	if (_DM) Rprintf("op:  r = Me - Mo;\n");

	// first check for a valid subtraction operation
	if (m_nNumCols != other.m_nNumCols)
		throw("Invalid operation in op(-)");
	if (m_nNumRows != other.m_nNumRows)
		throw("Invalid operation in op(-)") ;

	// construct the object we are going to return
    CFmMatrix* pTmp = FindReuseMatrix( m_nNumRows, m_nNumCols );

	// now subtract the other matrix
	for (int i = 0 ; i < m_nNumRows; ++i)
		for (int j = 0 ; j < m_nNumCols ; ++j)
			pTmp->Set(i, j, Get(i, j) - other.Get(i, j)) ;

	if ( other.IsReusable()) other.Release();
	if ( this->IsReusable()) this->Release();

	return *pTmp ;
}

//R: r = Me %*% Mo;
CFmMatrix& CFmMatrix::operator*( CFmMatrix &other)
{
	if (_DM) Rprintf("op:  r = Me * Mo;\n");

	// first check for a valid multiplication operation
	if (m_nNumCols != other.m_nNumRows)
	{
		throw("Matrices do not have common size in op(* matrix)");
	}

    CFmMatrix* pTmp = FindReuseMatrix( m_nNumRows, other.m_nNumCols );

	// e.g.
	// [A][B][C]   [G][H]     [A*G + B*I + C*K][A*H + B*J + C*L]
	// [D][E][F] * [I][J] =   [D*G + E*I + F*K][D*H + E*J + F*L]
	//             [K][L]
	//
    /*double	 value ;
	for (int i = 0 ; i < pTmp->m_nNumRows; i++)
		for (int j = 0 ; j < pTmp->m_nNumCols ; j++)
			{
				value = 0.0 ;
				for (int k = 0 ; k < m_nNumCols; k++)
					value += Get(i, k) * other.Get(k, j) ;

				pTmp->Set(i, j, value) ;
            }*/


    matprod( m_pData, GetNumRows(), GetNumCols(), 'N',
            other.GetData(), other.GetNumRows(), other.GetNumCols(), 'N', pTmp->GetData() );

	if ( other.IsReusable()) other.Release();
	if ( this->IsReusable()) this->Release();

	return *pTmp ;
}

//R: r = Me %*% Vo;
CFmMatrix& CFmMatrix::operator*(CFmVector &other)
{
	if (_DM) Rprintf("op:  r = Me * V;\n");

	// first check for a valid multiplication operation
	if (m_nNumCols != other.GetLength() )
    {
        Rprintf("[%d,%d]*%d", m_nNumRows, m_nNumCols, other.GetLength() );
		throw( "Matrices do not have common size in op(* vector)" );
    }

    CFmMatrix* pTmp = FindReuseMatrix( m_nNumRows, 1 );

    /*double	 value ;
	for (int i = 0 ; i < m_nNumRows ; ++i)
	{
		value = 0.0 ;
		for (int k = 0 ; k < m_nNumCols ; ++k)
			value += Get(i, k) * other.Get(k) ;

		pTmp->Set(i, 0, value) ;
    }*/

    matprod( m_pData, GetNumRows(), GetNumCols(), 'N',
            other.GetData(), other.GetLength(), 1, 'N', pTmp->GetData() );


	if ( other.IsReusable()) other.Release();
	if ( this->IsReusable()) this->Release();

	return *pTmp ;
}

//R: b = (Me == Mo);
bool CFmMatrix::operator==( CFmMatrix &other)
{
	if (_DM) Rprintf("op:  b = Me == Mo;\n");

	if (&other == this)
		return true ;

	if (m_pData == other.m_pData)
		return true ;

	// different dimensions
	if (m_nNumCols != other.m_nNumCols || m_nNumRows != other.m_nNumRows)
		return false ;

	// buffers are the same
	if (memcmp(m_pData, other.m_pData, sizeof(double) * m_nNumCols * m_nNumRows) == 0)
		return true ;

	if ( other.IsReusable()) other.Release();
	if ( IsReusable()) Release();

	return false ;
}

//R: r = Me + v;
CFmMatrix& CFmMatrix::operator+(double value)
{
	if (_DM) Rprintf("op:  r = Me + v;\n");

	// construct the object we are going to return
    CFmMatrix* pTmp = FindReuseMatrix( m_nNumRows, m_nNumCols );

	for (int i = 0 ; i < m_nNumRows ; i++)
		for (int j = 0 ; j < m_nNumCols ; j++)
			pTmp->Set(i, j, Get(i, j) + value ) ;

	if ( IsReusable()) Release();

	return *pTmp;
}

//R: r = Me - v;
CFmMatrix& CFmMatrix::operator-(double value)
{
	if (_DM) Rprintf("op:  r = Me - v;\n");

	// construct the object we are going to return
    CFmMatrix* pTmp = FindReuseMatrix( m_nNumRows, m_nNumCols );

	for (int i = 0 ; i < m_nNumRows ; i++)
		for (int j = 0 ; j < m_nNumCols ; j++)
			pTmp->Set(i, j, Get(i, j) - value ) ;

	if ( IsReusable()) Release();

	return *pTmp;

}

CFmMatrix& CFmMatrix::operator*(double value)
{
	if (_DM) Rprintf("op:  r = Me * v;\n");

	// construct the object we are going to return
    CFmMatrix* pTmp = FindReuseMatrix( m_nNumRows, m_nNumCols );

	for (int i = 0 ; i < m_nNumRows ; i++)
		for (int j = 0 ; j < m_nNumCols ; j++)
			pTmp->Set(i, j, Get(i, j) * value ) ;

	if ( this->IsReusable()) this->Release();

	return *pTmp;
}

CFmMatrix& CFmMatrix::operator/(double value)
{
	if (_DM) Rprintf("op:  r = Me / v;\n");

	// construct the object we are going to return
    CFmMatrix* pTmp = FindReuseMatrix( m_nNumRows, m_nNumCols );

	for (int i = 0 ; i < m_nNumRows ; i++)
		for (int j = 0 ; j < m_nNumCols ; j++)
			pTmp->Set(i, j, Get(i, j) / value ) ;

	if ( this->IsReusable()) this->Release();

	return *pTmp;
}

CFmMatrix& CFmMatrix::operator^(double value)
{
	if (_DM) Rprintf("op:  r = Me ^ v;\n");

	// construct the object we are going to return
    CFmMatrix* pTmp = FindReuseMatrix( m_nNumRows, m_nNumCols );

	for (int i = 0 ; i < m_nNumRows ; i++)
		for (int j = 0 ; j < m_nNumCols ; j++)
			pTmp->Set(i, j, pow(Get(i, j), value) ) ;

	if ( this->IsReusable()) this->Release();

	return *pTmp;
}


//R: Me = f;
void CFmMatrix::operator=( double value)
{
	if (_DM) Rprintf("op:  Me = v;\n");

	m_nNumCols = 0;
	m_nNumRows = 0;
	for (int i = 0 ; i < m_nMaxCols*m_nMaxRows; i++)
		m_pData[i] = value;
}


//R: Me = Mo;
void CFmMatrix::operator=( CFmMatrix &other)
{
	if (_DM) Rprintf("op:  Me = Mo;\n");

	// first check for a valid addition operation
	if (m_nMaxCols*m_nMaxRows < other.m_nNumCols*other.m_nNumRows)
    {
        double* pData = AllocateDouble(other.m_nMaxRows, other.m_nMaxCols);
        FreeDouble( m_pData );
        m_pData = pData;
        m_nMaxRows = other.m_nMaxRows;
        m_nMaxCols = other.m_nMaxCols;
    }

	m_nNumCols = other.m_nNumCols;
	m_nNumRows = other.m_nNumRows;
	for (int i = 0 ; i < m_nNumRows ; i++)
		for (int j = 0 ; j <m_nNumCols ; j++)
			Set(i, j, other.Get(i, j)) ;

	if (other.IsReusable() )
		other.Release();
}

//R: Me += Mo;
void CFmMatrix::operator+=( CFmMatrix &other)
{
	if (_DM) Rprintf("op:  Me += Mo;\n");

	// first check for a valid addition operation
    if (m_nNumCols==0 && m_nNumRows==0)
    {
        if (m_nMaxCols>=other.m_nNumCols &&
            m_nMaxRows>=other.m_nNumRows)
        {
            m_nNumCols=other.m_nNumCols;
            m_nNumRows=other.m_nNumRows;
        }
        else
            throw( "Insuffient matrix space for op(+=)" );
    }
    else
    {
	if (m_nNumCols != other.m_nNumCols)
		throw( "Invalid operation in op(+=)" );
	if (m_nNumRows != other.m_nNumRows)
		throw( "Invalid operation in op(+=)" );
    }

	for (int i = 0 ; i < m_nNumRows; ++i)
		for (int j = 0 ; j < m_nNumCols; ++j)
			Set(i, j, Get(i, j) + other.Get(i, j)) ;

	if (other.IsReusable() )
		other.Release();
}

//R: Me -= Mo;
void CFmMatrix::operator-=( CFmMatrix &other)
{
	if (_DM) Rprintf("op:  Me -= Mo;\n");

	// first check for a valid subtraction operation
	if (m_nNumCols != other.m_nNumCols)
		throw( "Invalid operation in op(-=)" );
	if (m_nNumRows != other.m_nNumRows)
		throw( "Invalid operation in op(-=)" );

	for (int i = 0 ; i < m_nNumRows ; ++i)
		for (int j = 0 ; j < m_nNumCols ; ++j)
			Set(i, j, Get(i, j) - other.Get(i, j) ) ;

	if (other.IsReusable() )
		other.Release();
}

//R: Me *= f;
void CFmMatrix::operator*=( CFmMatrix &other)
{
	if (_DM) Rprintf("op:  Me *= Mo;\n");

	// first check for a valid multiplication operation
	if (m_nNumRows != other.m_nNumCols)
		throw( "Matrices do not have common size in op(*= vector)" );

    CFmMatrix* pTmp = FindReuseMatrix( m_nNumRows, other.m_nNumCols );

	double	 value ;
	for (int i = 0 ; i < pTmp->m_nNumRows; i++)
		for (int j = 0 ; j < pTmp->m_nNumCols ; j++)
			{
				value = 0.0 ;
				for (int k = 0 ; k < m_nNumCols; k++)
					value += Get(i, k) * other.Get(k, j) ;

				pTmp->Set(i, j, value) ;
			}

	m_nNumRows = pTmp->m_nNumRows;
	m_nNumCols = pTmp->m_nNumCols;
	for (int i = 0 ; i < pTmp->m_nNumRows; i++)
		for (int j = 0 ; j < pTmp->m_nNumCols ; j++)
			Set(i,j, pTmp->Get(i,j));

	if ( other.IsReusable()) other.Release();
	if ( pTmp->IsReusable()) pTmp->Release();
}

//R: Me += f;
void CFmMatrix::operator+=(double value)
{
	if (_DM) Rprintf("op:  Me += f;\n");

	if (m_nNumCols==0 || m_nNumRows==0)
	{
		for(int i=0; i<m_nMaxRows*m_nMaxCols; i++)
			m_pData[i] += value;

		return;
	}

	// just multiply the elements by the value
	for (int i = 0 ; i < m_nNumRows ; ++i)
		for (int j = 0 ; j < m_nNumCols ; ++j)
			Set(i, j, Get(i, j) + value) ;
}

//R: Me -= f;
void CFmMatrix::operator-=(double value)
{
	if (_DM) Rprintf("op:  Me -= f;\n");

	if (m_nNumCols==0 || m_nNumRows==0)
	{
		for(int i=0; i<m_nMaxRows*m_nMaxCols; i++)
			m_pData[i] -= value;

		return;
	}

	// just multiply the elements by the value
	for (int i = 0 ; i < m_nNumRows ; ++i)
		for (int j = 0 ; j < m_nNumCols ; ++j)
			Set(i, j, Get(i, j) - value) ;
}

//R: Me *= f;
void CFmMatrix::operator*=(double value)
{
	if (_DM) Rprintf("op:  Me *= f;\n");

	if (m_nNumCols==0 || m_nNumRows==0)
	{
		for(int i=0; i<m_nMaxRows*m_nMaxCols; i++)
			m_pData[i] *= value;

		return;
	}

	// just multiply the elements by the value
	for (int i = 0 ; i < m_nNumRows ; ++i)
		for (int j = 0 ; j < m_nNumCols ; ++j)
			Set(i, j, Get(i, j) / value) ;
}

//R: Me /= f;
void CFmMatrix::operator/=(double value)
{
	if (_DM) Rprintf("op:  Me /= f;\n");

	if (m_nNumCols==0 || m_nNumRows==0)
	{
		for(int i=0; i<m_nMaxRows*m_nMaxCols; i++)
			m_pData[i] /= value;

		return;
	}

	// just multiply the elements by the value
	for (int i = 0 ; i < m_nNumRows ; ++i)
		for (int j = 0 ; j < m_nNumCols ; ++j)
			Set(i, j, Get(i, j) / value) ;
}

//R: Me ^= f;
void CFmMatrix::operator^=(double value)
{
	if (_DM) Rprintf("op:  Me ^= f;\n");

	if (m_nNumCols==0 || m_nNumRows==0)
	{
		for(int i=0; i<m_nMaxRows*m_nMaxCols; i++)
			m_pData[i] = pow( m_pData[i], value) ;

		return;
	}

	// just multiply the elements by the value
	for (int i = 0 ; i < m_nNumRows ; ++i)
		for (int j = 0 ; j < m_nNumCols ; ++j)
			Set(i, j, pow( Get(i, j),  value ) ) ;
}

void CFmMatrix::ValueMultiple( CFmMatrix &other)
{
    if (_DM) Rprintf("op:  Me *= Mo;\n");

    if (m_nNumRows != other.m_nNumRows ||
        m_nNumCols != other.m_nNumCols)
        throw( "Matrices do not have same size in Me *= Mo" );

    for (int i = 0 ; i < m_nNumRows ; ++i)
        for (int j = 0 ; j < m_nNumCols ; ++j)
            Set(i, j, Get(i, j) * other.Get(i,j) );

    if ( other.IsReusable()) other.Release();
}


#ifdef _DEBUG
// release version is in-line
double CFmMatrix::Get(int nRow, int nCol) const
{
	return m_pData[MI(m_nNumRows, m_nNumCols, nRow, nCol)] ;
}
#endif

//
// To implement the reuse of matrix
//
void CFmMatrix::AddRef()
{
	m_nRefer ++;
    //if (m_nRefer==1 && !m_bReuse)
    //    g_matList.push_back( this );

    //if (!m_bReuse)
    //    Rprintf("CFmMatrix:[REUSE:%s]%d, %s,%d[%d,%d]\n", m_bReuse?"YES":"NO", g_matList.size(), m_pId, m_nObjId, m_nNumCols, m_nNumRows );
}

void CFmMatrix::GlobalDump()
{
    //std::list<CFmMatrix*>::iterator it;
    Rprintf( "CFmMatrix::GlobalDump:\n" );
    //for ( it=g_matList.begin() ; it != g_matList.end(); it++ )
    //    Rprintf( "%s:%d[%d,%d]\n", (*it)->m_pId, (*it)->m_nObjId, (*it)->m_nNumRows, (*it)->m_nNumCols );
}

void CFmMatrix::Release()
{
	m_nRefer --;
	if ( m_nRefer == 0)
	{
		if (!m_bReuse)
		{
            //Rprintf("Before CFmMatrix:%d, %s,%d[%d,%d]\n", g_matList.size(), m_pId, m_nObjId, m_nNumCols, m_nNumRows );
            //g_matList.remove(this);
            //Rprintf("After CFmMatrix:%d\n", g_matList.size());

            FreeMemory();
		}
		else
		{
			if (_DM) Rprintf("***RECYCLE M(%d):\n", m_nObjId);
		}

	}
}

CFmMatrix& CFmMatrix::GetAbs()
{
	// construct the object we are going to return
    CFmMatrix* pTmp = FindReuseMatrix( m_nNumRows, m_nNumCols  );

	for (int i = 0 ; i < m_nNumRows ; i++)
			for (int j = 0 ; j < m_nNumCols ; j++)
				pTmp->Set(i, j, abs( (int)Get(i, j) ) ) ;

	if (IsReusable() ) Release();

	return (*pTmp);
}

double CFmMatrix::RowProd(int nRow1, int nRow2)
{
	// construct the object we are going to return
	double sum = 0;
	for (int i = 0 ; i < m_nNumCols ; i++)
		sum += Get(nRow1, i)*Get(nRow2,i) ;

	if (IsReusable() ) Release();

	return (sum);
}

double CFmMatrix::ColProd(int nCol1, int nCol2)
{
    // construct the object we are going to return
    double sum = 0;
    for (int i = 0 ; i < m_nNumRows ; i++)
        sum += Get(i, nCol1)*Get(i, nCol2) ;

    if (IsReusable() ) Release();

    return (sum);
}

CFmMatrix& CFmMatrix::GetTransposed()
{
	// construct the object we are going to return
    CFmMatrix* pTmp = FindReuseMatrix( m_nNumCols, m_nNumRows );

	for (int i = 0 ; i < m_nNumRows ; i++)
			for (int j = 0 ; j < m_nNumCols ; j++)
				pTmp->Set(j, i, Get(i, j)) ;

	if (IsReusable() ) Release();

	return (*pTmp);
}

void CFmMatrix::Transpose()
{
	// construct the object we are going to return
    CFmMatrix& pTmp = GetTransposed();

    m_nNumCols = pTmp.m_nNumCols ;
    m_nNumRows = pTmp.m_nNumRows ;

    // copy across the transposed data
	for (int i = 0 ; i < m_nNumRows ; ++i)
			for (int j = 0 ; j <  m_nNumCols; ++j)
				Set(i, j, pTmp.Get(i, j)) ;

    CFmVectorStr* vct = m_pRowNames;
    m_pRowNames = m_pColNames;
    m_pColNames = vct;

	if (pTmp.IsReusable() )
		pTmp.Release();
}

double CFmMatrix::GetSum()
{
    double fVal = 0;
    for(int i=0;i<m_nNumRows;i++)
        for(int j=0;j<m_nNumCols;j++)
            fVal += Get(i, j);

    return(fVal);
}

// matrix inversion will only work on square matrices
CFmMatrix& CFmMatrix::GetInverted( bool bCheck )
{
	if (m_nNumCols != m_nNumRows)
		throw( "Matrix must be square." );

    CFmMatrix* pDiag = FindReuseMatrix( m_nNumCols, m_nNumRows );
	pDiag->Square(m_nNumCols, true, 1.0);

	int ret = solve(m_pData, m_nNumCols, pDiag->m_pData, m_nNumCols, m_nNumCols);
    if (ret!=0)
        throw("error in solve method.");

    if (bCheck)
    {
        pDiag->AddRef();
        this->AddRef();
        static CFmMatrix fmId(0,0);
        fmId = (*pDiag) * (*this);
        if ( fabs( fmId.GetSum()- 1.0*fmId.GetNumRows() ) >1E-4 )
        {
            fmId.Show("fmId");
            Rprintf("Wrong inverted value:%f\n", fmId.GetSum());
        }
    }

	if (IsReusable() ) Release();

	return(*pDiag);
}

// matrix inversion will only work on square matrices
int CFmMatrix::Invert(  bool bCheck )
{
	if (m_nNumCols != m_nNumRows)
		throw( "Matrix must be square." );

    CFmMatrix* pDiag = FindReuseMatrix( m_nNumCols, m_nNumRows );
	pDiag->Square(m_nNumCols, true, 1.0);

	int ret = solve(m_pData, m_nNumCols, pDiag->m_pData, m_nNumCols, m_nNumCols);
    if (ret!=0 )
    {
        if (bCheck)
        {
            pDiag->Release();
            return(-1);
        }
        else
            throw("error in solve method.");
    }

    if (bCheck)
    {
        pDiag->AddRef();
        this->AddRef();
        static CFmMatrix fmId(0,0);
        fmId = (*pDiag) * (*this);
        if ( fabs( fmId.GetSum()- 1.0*fmId.GetNumRows() ) >1E-4 )
        {
            //fmId.Show("fmId");
            Rprintf("Wrong inverted value:%f\n", fmId.GetSum());
            pDiag->Release();
            return(-1);
        }
    }

	memcpy(m_pData, pDiag->m_pData, m_nNumCols*m_nNumCols*sizeof(double) );
	pDiag->Release();

    return(0);
}

CFmMatrix& CFmMatrix::SubMatrix(int row_start, int row_size, int col_start, int col_size)
{
	// make sure the requested sub matrix is in the current matrix
	if (col_start + col_size > m_nNumCols)
		throw( "Sub matrix is not contained in source");

	if (row_start + row_size > m_nNumRows)
		throw( "Sub matrix is not contained in source") ;

    CFmMatrix* pTmp = FindReuseMatrix( row_size, col_size );

	for (int i = 0; i < row_size; i++)
		for (int j = 0; j <col_size; j++)
			pTmp->Set(i, j, Get(row_start + i,  col_start+ j)) ;


	if (IsReusable() ) Release();

	return *pTmp ;
}

void CFmMatrix::SetSubMatrix(CFmMatrix &other, int row_start, int col_start )
{
	if (col_start + other.m_nNumCols> m_nNumCols)
		throw( "Sub matrix is not contained in source" );

	if (row_start + other.m_nNumRows> m_nNumRows)
		throw( "Sub matrix is not contained in source" );

	for (int i = 0 ; i < other.m_nNumRows ; i++)
		for (int j = 0 ; j < other.m_nNumCols ; j++)
			Set( row_start+ i, col_start + j, other.Get(i, j)) ;

	if ( other.IsReusable()) other.Release();
}

void CFmMatrix::CopyMatrix( CFmMatrix* other, int row_start, int row_size, int col_start, int col_size)
{
	if ( col_start!=-1)
	{
		if ( col_size> m_nMaxCols)
			throw( "over column size" );

		if (col_start + col_size> other->m_nNumCols)
			throw( "Sub matrix is not contained in source" );
	}
	else
	{
		col_start = 0;
		col_size = other->m_nNumCols;
	}

	if ( row_start!=-1)
	{
		if (row_size> m_nMaxRows)
			throw( "over row size" );

		if (row_start + row_size> other->m_nNumRows)
			throw( "Sub matrix is not contained in source" );
	}
	else
	{
		row_start = 0;
		row_size = other->m_nNumRows;
	}


	m_nNumRows = row_size;
	m_nNumCols = col_size;

	for (int i = 0 ; i < row_size ; i++)
		for (int j = 0 ; j < col_size ; j++)
		{
			if ( !Set( i, j, other->Get(row_start+i, col_start+j)) )
				Rprintf("Failed to set %d, %d\n", i, j);
		}

	if ( other->IsReusable()) other->Release();
}

void CFmMatrix::ResetRow( int nRow, CFmVector& vct)
{
    if ( nRow >= m_nNumRows && nRow<0)
        throw("over maximum rows' size." );

    if (m_nNumCols != vct.GetLength() )
        throw( "Cannot ResetRow(), not same size" );

    // now add the other matrix
    for (int j = 0 ; j <m_nNumCols; ++j)
        Set(nRow, j, vct.Get(j)) ;

    if ( vct.IsReusable()) vct.Release();
}

void CFmMatrix::SetRow( int nRow, CFmVector& vct, char* szRowName )
{
	if ( nRow >= m_nNumRows)
		throw("over maximum rows' size." );

	if (m_nNumCols != vct.GetLength() )
        throw( "Cannot SetRow(), not same size" );

	// now add the other matrix
	for (int j = 0 ; j <m_nNumCols; ++j)
		Set(nRow, j, vct.Get(j)) ;

    SetRowName(nRow, szRowName);

	if ( vct.IsReusable()) vct.Release();
}

void CFmMatrix::SetRow( int nRow, double* fBuf , char* szRowName)
{
	if ( nRow >= m_nNumRows)
		throw("over maximum rows' size." );

	// now add the other matrix
	for (int j = 0 ; j <m_nNumCols; ++j)
		Set(nRow, j, fBuf[j] ) ;

    SetRowName(nRow, szRowName);
}

void CFmMatrix::SetRow( int nRow, double fVal , char* szRowName)
{
    if ( nRow >= m_nNumRows)
        throw("over maximum rows' size." );

    // now add the other matrix
    for (int j = 0 ; j <m_nNumCols; ++j)
        Set(nRow, j, fVal ) ;

    SetRowName(nRow, szRowName);
}

void CFmMatrix::SetCol( int nCol, CFmVector& vct , char* szColName)
{
	if ( nCol >= m_nNumCols)
		throw("over maximum columns' size." );

	if (m_nNumRows != vct.GetLength() )
        throw( "Cannot SetCol(), not same size" );

	// now add the other matrix
	for (int i = 0 ; i <m_nNumRows; i++)
		Set(i, nCol, vct.Get(i)) ;

    SetColName(nCol, szColName);

	if ( vct.IsReusable()) vct.Release();
}

void CFmMatrix::SetCol( int nCol, double* fBuf, char* szColName )
{
	if ( nCol >= m_nNumCols)
		throw("over maximum columns' size." );

	// now add the other matrix
	for (int i = 0 ; i <m_nNumRows; i++)
		Set(i, nCol, fBuf[i]) ;

    SetColName(nCol, szColName);
}

// concatinate the other matrix to ourselves
void CFmMatrix::Cbind(CFmMatrix &other)
{
	if (m_nNumRows!=0 && m_nNumCols!=0 )
	{
		if (m_nNumRows != other.m_nNumRows)
            throw( "Cannot concatenate matrices(Cbind), not same size" );
	}
	else
    {
        m_nNumRows = other.m_nNumRows;
        if (m_nMaxRows < m_nNumRows)
        {
            m_nMaxRows = m_nNumRows;
            double* pData= AllocateDouble( m_nMaxRows, m_nMaxCols);
            FreeDouble( m_pData );
            m_pData = pData;
        }
    }

    if ( (m_nNumCols + other.m_nNumCols) > m_nMaxCols )
        EnlargeCols( m_nNumCols + other.m_nNumCols + MAX_MARGIN_COL );

    int oldCols = m_nNumCols;
    m_nNumCols += other.m_nNumCols;

    // now add the other matrix
    for (int i = 0 ; i < other.m_nNumRows ; i++)
        for (int j = 0 ; j <other.m_nNumCols; j++)
        Set(i, j + oldCols, other.Get(i, j)) ;

    for(int i=m_nNumCols - other.m_nNumCols; i<m_nNumCols; i++ )
        SetColName( i, other.GetColName( i - oldCols ) );


	if ( other.IsReusable()) other.Release();
}


/*void CFmMatrix::Rbind(CFmMatrix &other)
{
	if (m_nNumCols !=0 && m_nNumRows !=0)
	{
		if (m_nNumCols != other.m_nNumCols)
            throw( "Cannot concatenate matrices(Rbind), not same size" );
	}
	else
		m_nNumCols = other.m_nNumCols;

	if (m_nNumRows + other.m_nNumRows> m_nMaxRows)
		throw( "over maximum rows' size." );

	// now add the other matrix
	for (int i = 0 ; i < other.m_nNumRows ; ++i)
		for (int j = 0 ; j <m_nNumCols; ++j)
			Set(i+ m_nNumRows, j , other.Get(i, j)) ;

	m_nNumRows += other.m_nNumRows;

    for(int i=m_nNumRows - other.m_nNumRows; i<m_nNumRows; i++ )
        SetRowName( i, other.GetRowName( i - (m_nNumRows - other.m_nNumRows) ) );

	if ( other.IsReusable()) other.Release();
}*/


// concatinate the other matrix to ourselves
void CFmMatrix::Cbind(CFmVector &other, char* szColName)
{
	if (m_nNumRows!=0 && m_nNumCols!=0 )
	{
		if (m_nNumRows != other.GetLength())
        {
            Rprintf("%d!=%d", m_nNumRows, other.GetLength() );
            other.Show("Vector");
            throw("Cannot concatenate matrices(Cbind), not same size" );
        }
	}
	else
    {
        m_nNumRows = other.GetLength();
        if (m_nMaxRows < m_nNumRows)
        {
            m_nMaxRows = m_nNumRows;
            double* pData= AllocateDouble( m_nMaxRows, m_nMaxCols);
            FreeDouble( m_pData );
            m_pData = pData;
        }
    }

	if (m_nNumCols + 1 > m_nMaxCols)
        EnlargeCols(m_nNumCols*12/10 );

    m_nNumCols += 1;
    for (int i = 0 ; i < m_nNumRows ; i++)
        Set(i, m_nNumCols-1, other.Get(i)) ;

    SetColName( m_nNumCols-1, szColName);
	if ( other.IsReusable())
		other.Release();
}


/*void CFmMatrix::Rbind(CFmVector &other, char* szRowName)
{
	if (m_nNumCols !=0 && m_nNumRows !=0)
	{
		if (m_nNumCols != other.GetLength())
            throw( "Cannot concatenate matrices(Rbind), not same size" );
	}
	else
		m_nNumCols = other.GetLength();

	if (m_nNumRows + 1 > m_nMaxRows)
		throw( "over maximum rows' size." );

	// now add the other matrix
	for (int i = 0 ; i < m_nNumCols ; ++i)
			Set(m_nNumRows+0, i, other.Get(i)) ;

    m_nNumRows += 1;
    SetColName(m_nNumRows-1, szRowName);


	if ( other.IsReusable())
		other.Release();
}*/


void CFmMatrix::Cbind(const double *pData, int nLen, char* szColName)
{
	if (m_nNumCols !=0 && m_nNumRows !=0)
	{
		if (m_nNumRows != nLen )
            throw( "Cannot concatenate matrices(Cbind), not same size" );
	}
	else
    {
        m_nNumRows = nLen;
        if (m_nMaxRows < m_nNumRows)
        {
            m_nMaxRows = m_nNumRows;
            double* pData= AllocateDouble( m_nMaxRows, m_nMaxCols);
            FreeDouble( m_pData );
            m_pData = pData;
        }
    }

	if (m_nNumCols + 1 > m_nMaxCols)
        EnlargeCols( m_nNumCols *12/10 );

    m_nNumCols += 1;
    // now add the new row
	for (int i = 0 ; i < m_nNumRows ; ++i)
        Set(i, m_nNumCols-1, pData[i]) ;

    SetColName(m_nNumCols-1, szColName);
}

/*void CFmMatrix::Rbind(const double *pData, int nLen, char* szRowName)
{
	if (m_nNumCols !=0 && m_nNumRows !=0)
	{
		if (m_nNumCols != nLen )
            throw( "Cannot concatenate matrices(Rbind), not same size" );
	}
	else
		m_nNumCols = nLen;

	if (m_nNumRows + 1 > m_nMaxRows)
		throw( "over maximum columns' size." );

	// now add the new row
	for (int i = 0 ; i < m_nNumCols ; ++i)
		Set(m_nNumRows, i,  pData[i]) ;

    m_nNumRows += 1;
    SetRowName(m_nNumRows-1, szRowName);
}*/

void CFmMatrix::Cbind(const double val, int nLen, char* szColName)
{
	if (m_nNumCols !=0 && m_nNumRows !=0)
	{
		if (m_nNumRows != nLen )
            throw( "Cannot concatenate matrices(Cbind), not same size" );
	}
	else
    {
        m_nNumRows = nLen;
        if (m_nMaxRows < m_nNumRows)
        {
            m_nMaxRows = m_nNumRows;
            double* pData= AllocateDouble( m_nMaxRows, m_nMaxCols);
            FreeDouble( m_pData);
            m_pData = pData;
        }
    }

	if (m_nNumCols + 1 > m_nMaxCols)
    {
        EnlargeCols( m_nNumCols *12/10 );
    }

    m_nNumCols += 1;
    // now add the new row
	for (int i = 0 ; i < m_nNumRows ; ++i)
        Set(i, m_nNumCols-1, val) ;

    SetColName(m_nNumCols-1, szColName);
}

/*void CFmMatrix::Rbind(const double val, int nLen, char* szRowName)
{
	if (m_nNumCols !=0 && m_nNumRows !=0)
	{
		if (m_nNumCols != nLen )
            throw( "Cannot concatenate matrices(Rbind), not same size" );
	}
	else
		m_nNumCols = nLen;

	if (m_nNumRows + 1 > m_nMaxRows)
		throw( "over maximum columns' size." );

	// now add the new row
	for (int i = 0 ; i < m_nNumCols ; ++i)
		Set(m_nNumRows, i,  val) ;

    m_nNumRows -= 1;
    SetRowName(m_nNumRows, szRowName);
}
*/

bool CFmMatrix::SetRowName(int nRow, char* szRowName)
{
    if (nRow<0 || nRow>= m_nNumRows)
    {
        _log_error(_HI_, "MATRIX: failed to set rowname for %d row", nRow);
        return(false);
    }

    if (szRowName)
    {
		CFmNewTemp fmRef;
        if (m_pRowNames==NULL)
            m_pRowNames = new (fmRef) CFmVectorStr( m_nNumRows, m_nMaxRows );

        if (m_nNumRows > m_pRowNames->GetLength())
            m_pRowNames->Resize(m_nNumRows);

        m_pRowNames->Set( nRow, szRowName );
    }
    return(true);
}

bool CFmMatrix::SetColName(int nCol, char* szColName)
{
    if (nCol<0 || nCol>= m_nNumCols)
    {
        _log_error(_HI_, "MATRIX: failed to set colname for %d col", nCol);
        return(false);
    }

    if (szColName)
    {
		CFmNewTemp fmRef;
        if (m_pColNames==NULL)
            m_pColNames = new (fmRef) CFmVectorStr( m_nNumCols, m_nMaxCols );

        if (m_nNumCols > m_pColNames->GetLength())
            m_pColNames->Resize(m_nNumCols);

        m_pColNames->Set( nCol, szColName);
    }

    return(true);
}

char* CFmMatrix::GetRowName(int nRow)
{
    if (nRow<0 || nRow>= m_nNumRows)
        return(NULL);
    else
    {
        if (m_pRowNames==NULL)
            return NULL;
        return(m_pRowNames->Get(nRow));
    }
}

char* CFmMatrix::GetColName(int nCol)
{
    if (nCol<0 || nCol>= m_nNumCols)
        return(NULL);
    else
    {
        if (m_pColNames==NULL)
            return NULL;
        return(m_pColNames->Get(nCol));
    }
}

void CFmMatrix::Square(int nSize, bool diagonal, double int_value)
{
    if (m_nMaxRows<nSize || m_nMaxCols < nSize)
	{
        double* pNew = AllocateDouble(nSize, nSize);
		m_nMaxRows = nSize;
		m_nMaxCols = nSize;

        FreeMemory();

        m_pData = pNew;
	}

    m_nNumRows = nSize;
	m_nNumCols = nSize;

    if (m_pRowNames) m_pRowNames->Reset(m_nNumRows);
    if (m_pColNames) m_pColNames->Reset(m_nNumCols);

	for (int i=0; i<nSize; i++)
		for (int j=0; j<nSize; j++)
		{
			Set(i, j, 0.0);
			if (diagonal && i==j )
				Set(i,j, int_value);
		}

	return;
}

CFmVector& CFmMatrix::GetRowNonzero( int nRow )
{
	if (nRow<0 || nRow>= m_nNumRows)
        throw( "invalid row value." );

    CFmVector* pVct = CFmVector::FindReuseVector( m_nNumCols );

	int jj=0;
	for (int j=0; j<m_nNumCols; j++)
		if (Get(nRow,j)!=0)
		{
			pVct->Set(jj, Get(nRow,j) );
			jj++;
		}

	pVct->SetLength( jj );

	if (IsReusable() ) Release();

	return(*pVct);
}

CFmVector& CFmMatrix::GetRowNonNan( int nRow )
{
        if (nRow<0 || nRow>= m_nNumRows)
        throw( "invalid row value." );

    CFmVector* pVct = CFmVector::FindReuseVector( m_nNumCols );

        int jj=0;
        for (int j=0; j<m_nNumCols; j++)
                if (!isnan( Get(nRow,j)) )
                {
                        pVct->Set(jj, Get(nRow,j) );
                        jj++;
                }

        pVct->SetLength( jj );

        if (IsReusable() ) Release();

        return(*pVct);
}

void CFmMatrix::ResetNanvalue(double fNan)
{
    for (int i=0; i<m_nNumCols*m_nNumRows; i++)
        if (m_pData[i]== fNan )
            m_pData[i] = NA_REAL;

}

int CFmMatrix::GetNanCount()
{
	int k=0;
	for (int i=0; i<m_nNumCols*m_nNumRows; i++)
        if (isnan(m_pData[i]))
			k++;

	if (IsReusable() ) Release();

	return k;
}

double CFmMatrix::GetDet()
{
    if (m_nNumRows != m_nNumCols)
        throw("'A' must be a square matrix");

	return( GetDeterminant( m_pData, m_nNumRows ) );
}

double CFmMatrix::GetMean()
{
    double fSum = 0;
    int nSum = 0;
    for(int i=0;i <m_nNumCols*m_nNumRows; i++)
    {
        if ( !ISNAN(m_pData[i]) )
        {
            nSum ++;
            fSum += m_pData[i];
        }
    }

    fSum /= nSum;
    return (fSum);
}

double CFmMatrix::GetSd()
{
    double mu = GetMean();

    double fSum = 0;
    int nSum = 0;
    for(int i=0;i <m_nNumCols*m_nNumRows; i++)
    {
        if ( !ISNAN(m_pData[i]) )
        {
            nSum ++;
            fSum += pow( (m_pData[i]-mu), 2 );
        }
    }

    double ret = sqrt(fSum/nSum);
    return (ret);
}


double CFmMatrix::SumColumn(int column)
{
	//ASSERT(column >= 0) ;					// bad column
	//ASSERT(column < m_nNumCols) ;				// bad column
	double	sum = 0.0 ;

    for (int i = 0 ; i < m_nNumRows ; ++i)
        sum += Get(i, column) ;

	if (IsReusable() ) Release();

	return sum ;
}

double CFmMatrix::SumRow(int row)
{
	//ASSERT(row >= 0) ;						// bad row
	//ASSERT(row < m_nNumRows) ;					// bad row
	double	sum = 0.0 ;

    for (int i = 0 ; i < m_nNumCols ; ++i)
        sum += Get( row, i ) ;

	if (IsReusable() ) Release();

	return sum ;
}

// returns the minimum value in a row of the matrix
double CFmMatrix::GetRowMin(int row)
{
	//ASSERT(row >= 0) ;
	//ASSERT(row < m_nNumRows) ;
	double	value = Get(0, row) ;
	for (int i = 1 ; i < m_nNumCols ; ++i)
		{
		if (Get(i, row) < value)
			value = Get(i, row) ;
		}

	if (IsReusable() ) Release();

	return value ;
}

// returns the maximum value in a row of the matrix
double CFmMatrix::GetRowMax(int row)
{
	//ASSERT(row >= 0) ;
	//ASSERT(row < m_nNumRows) ;
	double	value = Get(0, row) ;
	for (int i = 1 ; i < m_nNumCols ; ++i)
		{
		if (Get(i, row) > value)
			value = Get(i, row) ;
		}

	if (IsReusable() ) Release();

	return value ;
}

// returns the minimum value in a column of the matrix
double CFmMatrix::GetColumnMin(int column)
{
	//ASSERT(column >= 0) ;
	//ASSERT(column < m_nNumCols) ;
	double	value = Get(column, 0) ;
	for (int i = 1 ; i < m_nNumRows ; ++i)
		{
		if (Get(column, i) < value)
			value = Get(column, i) ;
		}

	if (IsReusable() ) Release();

	return value ;
}

// returns the maximum value in a column of the matrix
double CFmMatrix::GetColumnMax(int column)
{
	//ASSERT(column >= 0) ;
	//ASSERT(column < m_nNumCols) ;
	double	value = Get(column, 0) ;
	for (int i = 1 ; i < m_nNumRows ; ++i)
    {
		if (Get(column, i) > value)
			value = Get(column, i) ;
    }

	if (IsReusable() ) Release();

	return value ;
}

int CFmMatrix::WriteAsCSVFile(const char* filename, bool bAppend, const char* szTag, const char* szFloatFmt) const
{
    char szfmt[256]={","};
    if (szFloatFmt==NULL)
        strcpy(szfmt, ",%.8f");
    else
        strcat(szfmt, szFloatFmt);

    FILE* fp = NULL;
    if (!bAppend)
        fp = fopen( filename, "wt");
    else
        fp = fopen( filename, "a+");

    if (fp==NULL)
    {
        _log_error(_HI_, "The file can not be created");
        return(ERR_CREATE_FILE);
    }

    if (szTag!=NULL)
    {
        fprintf(fp, "\n>>-----------------%s--------------------<<\n", szTag);
    }

    if (m_pColNames)
    {
        fprintf(fp, "\\");
        for (int j = 0; j < m_nNumCols ; j++)
        {
            if (m_pColNames->GetLength()>j)
                fprintf(fp, ",%s", m_pColNames->Get(j) );
            else
                fprintf(fp, ", ");
        }
    }
    else
    {
        for (int j = 0; j < m_nNumCols ; j++)
        {
            fprintf(fp, ",Col%d", j);
        }
    }

    fprintf(fp, "\n");

    for (int i = 0; i < m_nNumRows ; i++)
    {
        if (m_pRowNames && m_pRowNames->GetLength()>i )
           fprintf(fp, "%s", m_pRowNames->Get(i));
        else
           fprintf(fp, "%d", i+1);

        for (int j = 0; j < m_nNumCols ; j++)
            if ( ISNAN( Get(i,j)) )
                fprintf(fp, ",NA");
            else
                fprintf(fp, szfmt, Get(i,j)  );

        fprintf(fp, "\n");
    }

    fclose(fp) ;

    return(0);
}

#define MAX_LINE_WIDTH 256*256
// this function reads a file for "," separated values and creates a matrix object from it
int CFmMatrix::ReadFromCSVFile(const char* filename, bool bColName, bool bRowName)
{
    FILE* fp = fopen( filename, "rt");
    if (fp==NULL)
    {
        _log_error(_HI_, "The file can not be opened(%s)", filename);
        return(ERR_OPEN_FILE);
    }

    char aLine[MAX_LINE_WIDTH+1]={0};
    if( fgets( aLine, MAX_LINE_WIDTH, fp ) == NULL)
    {
        _log_error( _HI_, "Failed to read a line (%s)", filename);
        return(ERR_READ_FILE);
    }

    CFmVectorStr vctRowNames(0);
    CFmVectorStr vctColNames(0);
    int col_size = 0;

    // count how many elements in the first lines
    char seps1[] = ",\t";
    char* token = strtok( aLine, seps1 );
    while(token)
    {
        vctColNames.Put(token);

        token = strtok( NULL, seps1 );
        col_size++;
    }

    CFmVector tmpRow(col_size, 0.0);
    if (bRowName && !bColName)
        tmpRow.Resize(col_size-1);

    CFmMatrix t_mat(tmpRow.GetLength(), 0);
    if (!bColName)
    {
        if (!bRowName)
        {
            for (int i=0; i<col_size; i++)
                tmpRow[i] = atof( vctColNames[i]);
            t_mat.Cbind(tmpRow);
        }
        else
        {
            for (int i=0; i<col_size-1; i++)
                tmpRow[i] = atof( vctColNames[i+1]);
            t_mat.Cbind(tmpRow);

            vctRowNames.Put( vctColNames[0]);
        }
    }

    CFmVectorStr vctTmpRow(0);
    while(!feof(fp))
    {
        if( fgets( aLine, MAX_LINE_WIDTH, fp ) == NULL)
            break;

        vctTmpRow.Reset(0);
        char* token = strtok( aLine, seps1 );
        int nNonEmty=0;
        while(token)
        {
            if (strlen(token)!=0 &&
                strncmp(token, "\n", 1)!=0 &&
                strncmp(token, "\r", 1)!=0 &&
                strncmp(token, "\t", 1)!=0)
                nNonEmty ++;

            vctTmpRow.Put(token);
            token = strtok( NULL, seps1 );
        }

        if (nNonEmty==0)
            break;

        if (!bRowName)
        {
            for (int i=0; i<col_size && i<vctTmpRow.GetLength(); i++)
                tmpRow[i] = atof( vctTmpRow[i] );
            t_mat.Cbind(tmpRow);
        }
        else
        {
            for (int i=1; i<col_size && i<vctTmpRow.GetLength(); i++)
            {
                tmpRow[i-1] = atof( vctTmpRow[i]);
            }
            t_mat.Cbind(tmpRow);
            vctRowNames.Put( vctTmpRow[0]);
        }
    }

    fclose(fp);

    m_nNumCols = t_mat.GetNumRows();
    m_nNumRows = t_mat.GetNumCols();
    m_nMaxRows = m_nNumRows;
    m_nMaxCols = m_nNumCols;

    if (m_pData) FreeDouble( m_pData );
    m_pData = AllocateDouble( m_nNumRows, m_nNumCols );

    for(int i=0; i<m_nNumRows; i++)
        for(int j=0; j<m_nNumCols; j++)
            m_pData[MI(m_nNumRows, m_nNumCols, i, j)] = t_mat.Get(j, i);

   	if(bRowName)
    	for(int i=0; i<m_nNumRows; i++)
    	    SetRowName(i, vctRowNames.Get(i));

	if(bColName)
	    for(int i=0; i<m_nNumCols; i++)
	        SetColName(i, vctColNames.Get(i));

    return(0);
}

CFmVector& CFmMatrix::GetRow(int nRow)
{
    CFmVector* pTmp = CFmVector::FindReuseVector( m_nNumCols );

	for(int i=0; i<m_nNumCols; i++)
		pTmp->Set(i, Get(nRow, i));

	if (IsReusable() ) Release();

	return *pTmp;
}

CFmVector& CFmMatrix::GetCol(int nCol)
{
    CFmVector* pTmp = CFmVector::FindReuseVector( m_nNumRows );

	for(int i=0; i<m_nNumRows; i++)
		pTmp->Set(i, Get(i, nCol));

	if (IsReusable() ) Release();

	return *pTmp;
}

bool CFmMatrix::IsFitted(int nRows, int nCols)
{
    return( m_nMaxCols>= nCols && m_nMaxRows >= nRows );
}

void CFmMatrix::Resize(int nRows, int nCols, bool reset)
{
    if (m_nMaxCols*m_nMaxRows  < nCols*nRows )
    {
        double *pNew = AllocateDouble(nRows, nCols );
        memcpy( pNew, m_pData, m_nNumRows*m_nNumCols*sizeof(double) );

        FreeDouble( m_pData );
        m_pData = pNew;

        m_nMaxRows = nRows;
        m_nMaxCols = nCols;
    }

    m_nNumRows = nRows;
    m_nNumCols = nCols;
    if (m_pColNames) m_pColNames->Reset(m_nNumCols);
    if (m_pRowNames) m_pRowNames->Reset(m_nNumRows);

    if (reset)
        memset( m_pData, 0, sizeof(double)*m_nNumRows*m_nNumCols);
}

#define REUSE_MATRIX_COUNT	1000
CFmMatrix* CFmMatrix::FindReuseMatrix(int nRows, int nCols)
{
	if (g_pReused==NULL)
	{
        g_pReused = (CFmMatrix**)Calloc( REUSE_MATRIX_COUNT, CFmMatrix* );
        memset(g_pReused, 0, sizeof(CFmMatrix*)*REUSE_MATRIX_COUNT );

//**MEMTEST: Rprintf("ALLOC MEMORY: %d\n", sizeof(CFmMatrix*)*REUSE_MATRIX_COUNT);
	}

	for(int i=0; i<REUSE_MATRIX_COUNT; i++)
	{
        CFmMatrix* pTmp = (CFmMatrix* )(g_pReused[i]);
		if ( pTmp!=NULL &&
			 pTmp->GetRefer()==0 &&
			 pTmp->IsFitted(nRows, nCols))
		{
			pTmp->Resize(nRows, nCols, true);
			pTmp->AddRef();

			if (_DM) Rprintf("***REUSE M(%d)\n", pTmp->m_nObjId);
			return(pTmp);
		}
        if (_DM && pTmp) Rprintf("***------TMP M(%d, %d)\n", pTmp->m_nObjId, pTmp->GetRefer());
    }

	CFmNewTemp fmRef;
    CFmMatrix* pMat = new (fmRef)CFmMatrix(nRows, nCols, true);
    if (_DM)
        Rprintf("***NEW M(%d): Len=%d\n", pMat->m_nObjId, nRows*nCols);
	if ( pMat->m_nObjId > 5000 )
		_DM=1;
	if ( pMat->m_nObjId > 5000*2 )
        throw("Resue Matrics have big problem\n...");

	for(int i=0; i<REUSE_MATRIX_COUNT; i++)
	{
        CFmMatrix* pTmp = (CFmMatrix* )(g_pReused[i]);
		if ( pTmp==NULL)
		{
			g_pReused[i] = pMat;
			break;
		}
	}

	return(pMat);
}

bool CFmMatrix::RemoveRows(CFmVector& nRows)
{
    for (int i=0; i<nRows.GetLength(); i++)
    {
        if ( nRows.Get(i)<0 || nRows.Get(i)>=m_nNumRows )
            return(false);
    }

    int newNumRow = m_nNumRows - nRows.GetLength();
    double* pDstData = m_pData;

    for(int j=0; j<m_nNumCols;j++)
    {
        for(int i=0, k=0; i<m_nNumRows; i++)
        {
            if (nRows.Find(i)>=0)
            {
            }
            else
            {
                pDstData[MI(newNumRow, m_nNumCols, k, j)] = m_pData[ MI(m_nNumRows, m_nNumCols, i, j)];
                k++;
            }
        }
    }

    //delete[] m_pData;
    //m_pData = pData;
    m_nNumRows = newNumRow;

    if(m_pRowNames)
        m_pRowNames->RemoveElements(nRows);

    return(true);
}

void CFmMatrix::EnlargeCols( int nMaxCols )
{
    double* pData = AllocateDouble(m_nMaxRows, nMaxCols);
    memcpy(pData, m_pData, sizeof(double)*m_nNumRows*m_nNumCols );
    FreeDouble( m_pData );

    m_nMaxCols = nMaxCols;
    m_pData = pData;
}

CFmVectorStr*   CFmMatrix::GetColNames()
{
    return (m_pColNames);
}

CFmVectorStr*  CFmMatrix:: GetRowNames()
{
    return (m_pRowNames);
}

bool CFmMatrix::SetColNames(CFmVectorStr* pNames)
{
    if (pNames->GetLength()!=m_nNumCols)
        return (false);

    if (m_pColNames)
        destroy( m_pColNames);

	CFmNewTemp fmRef;
    m_pColNames = new (fmRef) CFmVectorStr( pNames );
    return(true);
}

bool CFmMatrix::SetRowNames(CFmVectorStr* pNames)
{
    if (pNames->GetLength()!=m_nNumRows)
        return (false);

    if (m_pRowNames)
        destroy( m_pRowNames);

	CFmNewTemp fmRef;
    m_pRowNames = new (fmRef)CFmVectorStr( pNames );
    return(true);
}

int CFmMatrix::FindColumn(const char* szColumn)
{
	if (m_pColNames==NULL)
		return(-1);

    for (int i=0;i<m_nNumCols; i++)
    {
        if ( strcasecmp ( m_pColNames->Get(i), szColumn)==0 )
            return(i);
    }

    return( -1 );
}

SEXP GetSEXP(CFmMatrix* pMat)
{
    SEXP sep;
    PROTECT( sep = allocMatrix(REALSXP, pMat->GetNumRows(), pMat->GetNumCols() ) );
    double* sep_r = REAL(sep);

    double* pData = pMat->GetData();
    for (int i=0;i<pMat->GetNumRows()*pMat->GetNumCols(); i++)
        sep_r[i] = pData[i];

    if (pMat->GetColNames() || pMat->GetRowNames())
    {
        SEXP dim_name, col_name, row_name;
        int nProtected=0;
        PROTECT(dim_name = allocVector(VECSXP, 2));

        if (pMat->GetColNames())
        {
            PROTECT(col_name = allocVector( STRSXP, pMat->GetNumCols()));
            nProtected++;
            for (int i=0; i<pMat->GetNumCols(); i++)
                if (pMat->GetColName(i))
                    SET_STRING_ELT(col_name, i, mkChar(pMat->GetColName(i)));
        }
        else
            col_name = R_NilValue;

        if (pMat->GetRowNames())
        {
            PROTECT(row_name = allocVector( STRSXP, pMat->GetNumRows()));
            nProtected++;
            for (int i=0; i<pMat->GetNumRows(); i++)
            {
                if (pMat->GetRowName(i))
                    SET_STRING_ELT(row_name, i, mkChar(pMat->GetRowName(i)));
            }
        }
        else
            row_name = R_NilValue;

        SET_VECTOR_ELT(dim_name, 1, col_name);
        SET_VECTOR_ELT(dim_name, 0, row_name);

        setAttrib(sep, R_DimNamesSymbol, dim_name);
        UNPROTECT(nProtected+1);
    }

	UNPROTECT(1);
    return sep;
}

int GetMatrix(SEXP pExp, CFmMatrix* pMat)
{
    SEXP dims = getAttrib(pExp, R_DimSymbol);
    int nRow=0, nCol=0;

    if (length(dims) == 2)
    {
        nRow = INTEGER(dims)[0];
        nCol = INTEGER(dims)[1];
    }
    else if(length(dims) == 1)
    {
        nRow = 1;
        nCol = INTEGER(dims)[0];
    }
    else
        return (-1);

    pMat->Resize(nRow, nCol);
    double* pExp_r = REAL(pExp);

    for(int i=0; i<nRow; i++)
        for(int j=0; j<nCol; j++)
    {
        pMat->Set(i, j, pExp_r[j*nRow+i] );
    }


	SEXP dimnames = getAttrib( pExp, R_DimNamesSymbol );
	SEXP colnames = VECTOR_ELT(dimnames, 1);

	for (int j = 0; j < length(colnames); j++)
		pMat->SetColName( j, (char*)CHAR( STRING_ELT( colnames, j) ) );

	SEXP rownames = VECTOR_ELT(dimnames, 0);
	for (int i = 0; i < length(rownames); i++)
		pMat->SetRowName( i, (char*)CHAR( STRING_ELT( rownames, i) ) );

    return(0);
}

void CFmMatrix::StatCache(int* pnTotal, int* pnUsed)
{
	if (g_pReused==NULL)
	{
		*pnTotal = 0;
		*pnUsed = 0;
		return;
	}

	for(int i=0; i<REUSE_MATRIX_COUNT; i++)
	{
        CFmMatrix* pTmp = (CFmMatrix* )(g_pReused[i]);
		if ( pTmp==NULL )
			continue;

		*pnTotal = *pnTotal + 1;
		if( pTmp->GetRefer()> 0)
			*pnUsed = *pnUsed + 1;
	}

	return;
}

void destroy(CFmMatrix* p)
{
	CFmNewTemp  fmRef;
	p->~CFmMatrix();
	operator delete(p, fmRef);
}
