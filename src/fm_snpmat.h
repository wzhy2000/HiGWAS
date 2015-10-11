// snpmat.h: interface for the CFmSnpMat class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_SNPMAT_H_)
#define _FM_SNPMAT_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <math.h>
//
// Class : CFmSnpMat
//

#define UNUSED_PARAMETER(x) x

#include <stdio.h>
#include <string.h>
#include "fm_vector.h"

class CFmVectorStr;
//class CFmVector;
class CFmMatrix;

class CFmSnpMat
{
public:
//    friend class	CFmVector ;

    CFmSnpMat(int nRows, int nCols, int nMaxRows=-1, int nMaxCols=-1) ;
    virtual	~CFmSnpMat();

    inline int		GetNumSubjs() const { return m_nNumSubjs ; };
    inline int		GetNumSnps() const { return m_nNumSnps ; };
    inline int		GetMaxSubjs() const { return m_nMaxSubjs ; };
    inline int		GetMaxSnps() const { return m_nMaxSnps ; };
    inline char* 	GetData()const { return m_pData ; };
    inline int      GetLength()const { return m_nNumSubjs*m_nNumSnps ; };
    inline unsigned int GetBytes()const { return m_nNumSubjs*m_nNumSnps*sizeof(char); };
    inline bool     Set(int nSnp, int nSubj, double value)
	{
        if ( nSubj<m_nNumSubjs && nSnp < m_nNumSnps )
		{
            m_pData[ nSubj*m_nBaseSnps + nSnp] = (char)value ;
            return true ;
		}
		else
        {
            throw("Error in matrix.set");
            return false;
        }
	};

    inline double Get_a(int nSnp, int nSubj) const
    {
		if (m_bSD)
        	return (double)(m_pData[ nSubj*m_nBaseSnps + nSnp])/(((CFmVector&)(*m_pSD))[nSnp]);
        else
        	return (double)(m_pData[ nSubj*m_nBaseSnps + nSnp]);
    };

    inline double Get_d(int nSnp, int nSubj) const
    {
	    double d = (double)(m_pData[ nSubj*m_nBaseSnps + nSnp]) ;
	    d = (fabs(d)*(-1.0)+1.0);
		if (m_bSD)
	        d = d /(((CFmVector&)(*m_pSD))[nSnp]);
        return d;
    };

    void			operator=( CFmSnpMat &other) ;
    bool            SetSnpInfo(int nSnp, char* szSnpName, char* szChr, char* szPos);
    bool            SetSubjName(int nSubj, char* szSubjName);
    bool            GetSnpInfo(int nSnp, char* szSnpName, char* szChr, char* szPos);
    bool            GetSnpInfo(CFmVectorStr* pVctSnp, CFmVector* pVctChr, CFmVector* pVctPos );
    char*           GetSubjName(int nSubj);

    CFmVector&      GetSubj( int nSubj );
    CFmVector&      GetSnp( int nSnp );
    void            SetSubj( int nSubj, double* fData, int nLen, char* szSubjName=NULL);
    void            SetSubj( int nSubj, CFmVector& vct, char* szSubjName=NULL);
    void            SetSnp( int nSnp, double* fData, int nLen, char* szSnpName, char* szChr=NULL, char* szPos=NULL);
    void            SetSnp( int nSnp, CFmVector& vct, char* szSnpName, char* szChr=NULL, char* szPos=NULL);

    CFmMatrix*      ToMatrix();
    CFmMatrix*      GetAdd();
    CFmMatrix*      GetDom();
    void            FromMatrix(CFmMatrix* mat);

	// export functions/import
    int             WriteAsCSVFile(const char* filename);
    int             ReadFromCSVFile(const char* filename, bool bColName, bool bRowName);

	void			sd();

private:
	// internal variables
    int				m_nNumSubjs;
    int				m_nNumSnps;
    int 		    m_nMaxSnps;
    int 		    m_nMaxSubjs;
    int             m_nBaseSnps;
    char*			m_pData ;
    CFmVectorStr*   m_pSnpInfos;
    CFmVectorStr*   m_pSubjInfos;

	bool		    m_bSD;
	CFmVector* 	    m_pSD;
private:
	// private internal functions
    char*			AllocateChar(int nSnps, int nSubjs) ;
    void			FreeMemory();
};

void destroy(CFmSnpMat* p);

#endif // !defined(_FM_SNPMAT_H_)
