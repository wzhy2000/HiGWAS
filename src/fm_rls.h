// FM_res.h: interface for the result class
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_RLS_H_)
#define _FM_RLS_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <map>

#define VAR_TYPE_STRING 2
#define VAR_TYPE_DOUBLE 1
#define VAR_TYPE_INT    0

#define MAX_CONT_ITEM   32
#define DWORD unsigned int

class CFmVectorStr;
class CFmMatrix;
class CFmVector;

typedef struct res_file_var
{
    char varName[64];
    DWORD varType;
    DWORD dimRow;
    DWORD dimCol;
    DWORD nBytes;
}VARFMT;

typedef struct res_file_fmt
{
    char szHeader[256];
    VARFMT varCont[MAX_CONT_ITEM];
}RESFMT;

class CFmRls
{
public:
    CFmRls( const char* szRlsFile );
    virtual ~CFmRls();

    int SaveResData();
    int LoadResFile();
    int SaveRData( const char* szRdataFile );

    CFmMatrix* GetMatrix(char* szVar);
    CFmVector* GetVector(char* szVar);
    CFmVectorStr* GetVectorStr(char* szVar);

    void SetData(char* szVar, CFmMatrix* pMat );
    void SetData(char* szVar, CFmVector* pVec );
    void SetData(char* szVar, CFmVectorStr* pVecStr );

private:
    std::map<char*, CFmMatrix*> m_MatMap;
    std::map<char*, CFmVector*> m_VecMap;
    std::map<char*, CFmVectorStr*> m_VecStrMap;

    char* m_szRlsFile;
};

#endif
