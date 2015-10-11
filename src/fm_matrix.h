// fm_matrix.h: interface for the CFmMatrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_MATRIX_H_)
#define _FM_MATRIX_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//
// Class : CFmMatrix
//
// This is an encapsulation of a matrix class that allows all standard operations to be
// performed on it. It also includes some odd procedures for importing and exporting matrices.
//
// Data is stored for the matrix in a flat allocation with 1 additional entry at the end which is
// used as a reference counter for how many objects are using that data section. When one matrix
// differs from another when referencing the same memory, a copy is made at that point, and they
// go their separate ways. This is very useful when returing CFmMatrix objects which then need to be
// assigned as it avoid additional memory allocation/deallocation calls.
//
#define UNUSED_PARAMETER(x) x

#include <stdio.h>
#include <string.h>
#include <list>
#include "fm_vector.h"
#include "fm_vector_str.h"

class CFmMatrix
{
public:
    friend class	CFmVector ;		// used for operator[][]

    // construction and destruction
    // default constructor
    explicit CFmMatrix(const CFmMatrix& pOther);
    explicit CFmMatrix(bool bReused = false );
    explicit CFmMatrix(const CFmMatrix *pOther, bool bReused = false);
    // constructs an empty matrix of this size
    explicit CFmMatrix(int nRows, int nCols, bool bReused = false);
    // creates a square matrix
    explicit CFmMatrix(int nSize, bool set_diagonal = true, double init_value=0, bool bReused = false) ;
    explicit CFmMatrix(int nRows, int nCols, int nMaxRows, int nMaxCols, bool bReused = false) ;

    virtual	~CFmMatrix();

    inline int      GetRefer()	const { return m_nRefer; };
    inline bool     IsReusable()	const { return m_bReuse; };
    inline int      GetNumCols() const { return m_nNumCols ; };
    inline int      GetNumRows() const { return m_nNumRows ; };
    inline int      GetMaxCols() const { return m_nMaxCols ; };
    inline int      GetMaxRows() const { return m_nMaxRows ; };
    inline double*  GetData()const { return m_pData ; };
    inline unsigned int GetBytes()const { return m_nNumCols*m_nNumRows*sizeof(double); };
    inline bool     Set(int nRow, int nCol, double value)
    {
        if ( nCol<m_nNumCols && nRow < m_nNumRows )
        {
            m_pData[ nCol*m_nNumRows + nRow] = value ;
            return true ;
        }
        else
        {
            Rprintf("Error in matrix.set");
            throw("Error in matrix.set");
            return false;
        }
    };

    void            Show(const char* szName=NULL);

    // matrix mathematical operations
    CFmMatrix&		operator+( CFmMatrix &other) ;
    CFmMatrix&		operator-( CFmMatrix &other) ;
    CFmMatrix&		operator*( CFmMatrix &other) ;
    CFmMatrix&		operator/( CFmMatrix &other) ;

    CFmMatrix&		operator*( CFmVector &other) ;

    CFmMatrix&		operator+(double value) ;
    CFmMatrix&		operator-(double value) ;
    CFmMatrix&		operator*(double value) ;
    CFmMatrix&		operator/(double value) ;
    CFmMatrix&		operator^(double value) ;
    CFmMatrix&		operator>(double value) ;

    void		operator+=(double value) ;
    void		operator-=(double value) ;
    void		operator*=(double value) ;
    void		operator/=(double value) ;
    void		operator^=(double value) ;

    void		operator=( double value) ;
    void		operator=( CFmMatrix &other) ;
    void		operator+=( CFmMatrix &other) ;
    void		operator-=( CFmMatrix &other) ;
    void		operator*=( CFmMatrix &other) ;
    void		operator/=( CFmMatrix &other) ;
    bool		operator==( CFmMatrix &other) ;
    void        ValueMultiple( CFmMatrix &other) ;

    double		RowProd(int nRow1, int nRow2);
    double		ColProd(int nRow1, int nRow2);
    CFmMatrix&		GetAbs();
    // matrix transposition
    CFmMatrix&		GetTransposed();
    // matrix inversion
    CFmMatrix&		GetInverted(bool bCheck = false) ;

    // element access
#ifdef _DEBUG
    double          Get(int nRow, int nCol) const ;
#else
    inline double   Get(int nRow, int nCol) const { return m_pData[ nCol*m_nNumRows + nRow ] ; } ;
#endif
    bool            SetRowName(int nRow, char* szRowname);
    bool            SetColName(int nCol, char* szRowname);
    char*           GetRowName(int nRow);
    char*           GetColName(int nCol);
    CFmVectorStr*   GetColNames();
    CFmVectorStr*   GetRowNames();
    bool            SetColNames(CFmVectorStr* pNames);
    bool            SetRowNames(CFmVectorStr* pNames);

    void            SetId(const char* szId);
    double          SumColumn(int col)  ;
    double          SumRow(int row)  ;
    double          SumColumnSquared(int col)  ;
    double          SumRowSquared(int row)  ;
    double          GetRowMin(int row)  ;
    double          GetRowMax(int row)  ;
    double          GetColumnMin(int col)  ;
    double          GetColumnMax(int col)  ;

    int             FindColumn(const char* szColName);

    CFmVector&      GetRow(int nRow);
    CFmVector&      GetCol(int nCol);
    bool            RemoveRows(CFmVector& nRows);

    void            Transpose() ;
    int             Invert(bool bCheck = false) ;

    void            SetRow( int nRow, CFmVector& vct, char* szRowName=NULL );
    void            SetRow( int nRow, double* fBuf, char* szRowName=NULL );
    void            SetRow( int nRow, double fBuf, char* szRowName=NULL );
    void            SetCol( int nCol, CFmVector& vct, char* szColName=NULL );
    void            SetCol( int nCol, double* fBuf, char* szColName=NULL );
    void            ResetRow( int nCol, CFmVector& vct);

    // matrix concatenation
    //void          rbind(CFmMatrix &other) ;
    //void          Rbind(const double *pData, int nLen, char* szRowName=NULL ) ;
    //void          Rbind(CFmVector& vct, char* szRowName=NULL );
    //void          Rbind(const double fVal, int nLen, char* szRowName=NULL );
    void            Cbind(CFmMatrix &other) ;
    void            Cbind(const double *pData, int nLen, char* szColName=NULL ) ;
    void            Cbind(CFmVector& vct, char* szColName=NULL );
    void            Cbind(const double fVal, int nLen, char* szColName=NULL );

    // sub matrix extraction, setting
    CFmMatrix&      SubMatrix(int row_start, int row_size, int col_start, int col_size);
    void            SetSubMatrix(CFmMatrix &other, int row_start, int col_start);
    void            CopyMatrix( CFmMatrix* other, int row_start, int row_size, int col_start, int col_size);

    void            ResetNanvalue(double fNan);
    CFmVector&      GetRowNonzero( int nRow );
    CFmVector&      GetRowNonNan( int nRow );
    int             GetNanCount();
    double          GetDet();
    double          GetMean();
    double          GetSd();
    double          GetSum();

    void            Square( int nSize, bool diagonal = true, double int_value=0);
    void            Resize(int nRows, int nCols, bool reset=false);
    bool            IsFitted(int nRows, int nCols);

    // export functions/import
    int             WriteAsCSVFile(const char* filename, bool bAppend=true, const char* szTag=NULL, const char* szFloatFmt=NULL) const ;
    int             ReadFromCSVFile(const char* filename, bool bColName, bool bRowName);

    static void     GlobalDump();
public:
    static          int	_DM;
    static          int g_nObjCount;
    static void     StatCache(int* pnTotal, int* pnUsed);

protected:
    static          CFmMatrix** g_pReused;
    static          CFmMatrix* FindReuseMatrix( int nRow, int nCols );

private:
	// internal variables
    char*           m_pId;
    int             m_nRefer;
    int             m_nObjId;
    bool            m_bReuse;
    int             m_nNumCols;
    int             m_nNumRows;
    int             m_nMaxRows;
    int             m_nMaxCols;
    double*         m_pData ;
    CFmVectorStr*   m_pRowNames;
    CFmVectorStr*   m_pColNames;
    void            Release() ;	// decrements the m_pData reference count


private:
	// private internal functions
    double*         AllocateDouble(int nRows, int nCols) ;
    void            FreeMemory();
    void            EnlargeCols( int nMaxCols );
    void            FreeDouble( double* pData );

    // reference counting functions
    void            AddRef() ;	// increments the m_pData reference count

    // helper functions
    char*           GetRowAsText(int row) const ;
    static char*    ReadLine(FILE *file) ;		// reads a \r\n delimited line of text from a file
    static int      GetStringToken(char* source, char* destination, int start, char ch) ;
    static std::list<CFmMatrix*> g_matList;
};

SEXP GetSEXP(CFmMatrix* pMat);
int GetMatrix(SEXP pExp, CFmMatrix* pMat);
void destroy( CFmMatrix* p);

#endif // !defined(_FM_MATRIX_H_)
