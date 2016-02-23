// fm_vector.h: interface for the CFmVector class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_VECTOR_H_)
#define _FM_VECTOR_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//
// Class : CFmVector
//

#include "Rmethods.h"

class CFmMatrix;
class CFmVectorStr;

// this class is used to help the operator[][] on a matrix object work correctly
// it only provides the operator[]
class CFmVector
{
public:
    ///friend class	CFmMatrix ;

    explicit CFmVector(bool bReused = false );
    explicit CFmVector(CFmMatrix*, int nRow_or_Cols, bool bRow=true, bool bReused = false );
    explicit CFmVector(int nSize, double fInit, int nMaxSize= -1, bool bReused = false );
    explicit CFmVector(CFmVector* pVct, bool bReused = false );
    explicit CFmVector(const CFmVector& pOther);

public:
    virtual ~CFmVector();

    inline int	GetRefer()	const { return m_nRefer; } ;
    inline bool	IsReusable()	const { return m_bReuse; } ;
    inline double* GetData(){return m_pData;};
    inline int	GetLength(){return m_nActLen;};
    inline int	GetBytes(){return m_nActLen*sizeof(double);};
    inline double Get(int i){return m_pData[i];};
    inline double GetLast(){
        if (m_nActLen>0)
            return (m_pData[m_nActLen-1]);
        else
            throw("NULL vector");
    }

    inline void	Set(int i, double f){m_pData[i] = f;};

    inline void Show(const char* szName=NULL)
    {
        show_vector2(m_pData, m_nActLen, szName);
    }

    double operator[](int row) const
    {
        if (row>=m_nActLen || row<0) throw("wrong index of vector");
            return m_pData[row] ;
    };
    double& operator[](int row)
    {
        if (row>=m_nActLen || row<0) throw("wrong index of vector");
            return m_pData[row];
    };

    bool		operator==(CFmVector &other) ;

    CFmVector&		operator+(double other);
    CFmVector&		operator-(double other);
    CFmVector&		operator*(double other);
    CFmVector&		operator/(double other);
    CFmVector&		operator^(double other);

    CFmVector&		operator+(CFmVector& other);
    CFmVector&		operator-(CFmVector& other);
    CFmVector&		operator*(CFmVector& other);
    CFmVector&		operator/(CFmVector& other);

    CFmVector&		abs();
    CFmMatrix&		operator*(CFmMatrix& other);
    CFmMatrix&		GetRow();
    CFmMatrix&		GetCol();

    void		operator=(double other);
    void 		operator=(double* other);
    void		operator=(CFmVector& other);
    void		Set(CFmVector& other);

    CFmVector& 		GetTransposed();
    CFmVector& 		GetReciprocal();

    bool 		SetLength(int nLen);
    double		Sum();
    double          GetMean();
    double          GetMedian();
    double          GetVar();
    double          GetMax();
    double          GetMin();
    double          Prod(CFmVector& other);
    void            UniquePut(double fValue);
    void            Put(double fValue);
    void            Join(double fvalue);
    void            Remove(int nPos);
    void            RemoveNan();
    int             GetLengthNonvalue(double fvalue);
    int             GetLengthNonNan();
    bool            IsFitted(int nLen);
    bool            IsNan(int nIdx);
    void            Resize(int nLen, bool reset=true);
    int             Find(double fVal);
    void            Append(CFmVector& other);
    void            Sort( bool bDecreasing = false);
    bool            Order(CFmVector& order);
    bool            Rearrange(CFmVector& order);
    bool            RemoveElements(CFmVector& nRows);
    int             WriteAsCSVFile(const char* szCsvFile, bool bAppend=true, const char* szTag=NULL);
    int             LoadValueFromRData( const char* szRDataFile, const char* szVarname);
    char*           GetCommaString(const char* szFormat);

	// reference counting functions
    void            AddRef() ;	// increments the m_pData reference count
    void            Release() ;	// decrements the m_pData reference count
    void            SetNames( CFmVectorStr* pNames);
    CFmVectorStr*   GetNames();
    char*           GetName(int idx );

public:
    static          int _DV;
    static          int g_nObjCount;
    static void     StatCache(int* pnTotal, int* pnUsed);

    static          CFmVector** g_pReused;
    static          CFmVector* FindReuseVector( int nLen );

protected:
    double*         m_pData;
    int             m_nMaxLen;
    int             m_nActLen;
    int             m_nRefer;
    int             m_nObjId;
    bool            m_bReuse;

private:
    double*         AllocateMemory( int nLen );
    void            FreeMemory( double* p );
    char*           m_pszBuf;
    CFmVectorStr*   m_pNames;
} ;

SEXP GetSEXP(CFmVector* pVct);
int GetVector(SEXP pExp, CFmVector* pVct);
void destroy(CFmVector* p);


#endif // !defined(_FM_VECTOR_H_)
