// fm_vector_str.h: interface for the CFmVectorStr class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_VECTOR_STR_H_)
#define _FM_VECTOR_STR_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "fm_vector.h"

// this class is used to help the operator[][] on a matrix object work correctly
// it only provides the operator[]
class CFmVectorStr
{
public:
    CFmVectorStr(CFmVectorStr* pOther);
    CFmVectorStr(int nSize, int nMaxSize= -1);

public:
    virtual ~CFmVectorStr();

    char**  GetData();
    void    Reset(int size);
    void    Resize(int size);
    char*   Get(int idx);
    void	Set(int idx, char* str);
    char*   operator[](int idx) const;

    int     GetMaxStrlen();
    int     GetLength();
    int     GetBytes();
    void    Put(char* str);
    void    Append(CFmVectorStr* another);
    bool    Remove(int idx);
    bool    RemoveElements(CFmVector& nRows);
    void    Show(const char* szName=NULL);
    int     Find(const char* szName );
    void    SetNames( CFmVectorStr* pNames);
    CFmVectorStr* GetNames();
    char*  GetName(int idx );

    int     WriteAsCSVFile(const char* szCsvFile, bool bAppend=false, const char* szTag=NULL);

protected:
    char**   m_pData;
    int     m_nMaxLen;
    int     m_nActLen;
    unsigned int     m_nMaxStrlen;
    CFmVectorStr* m_pNames;
private:
    char** AllocateMemory( int nLen );

} ;

SEXP GetSEXP(CFmVectorStr* pVct);
void destroy(CFmVectorStr* p);

#endif // !defined(_FM_VECTOR_STR_H_)
