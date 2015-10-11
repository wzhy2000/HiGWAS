// fm_rtDump.h: interface for the result class
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_RTDUMP_H_)
#define _FM_RTDUMP_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <map>

#define VAR_TYPE_STRING 2
#define VAR_TYPE_DOUBLE 1
#define VAR_TYPE_INT    0

#define MAX_CONT_ITEM   256
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
    DWORD nPos;
    double fValue;
}RTVARFMT;

typedef struct rtdump_file_fmt
{
    char szHeader[256];
    RTVARFMT varCont[MAX_CONT_ITEM];
}RTDUMPFMT;

class CFmRtDump
{
public:
    CFmRtDump();
    virtual ~CFmRtDump();

    int OpenFile( char *szRtdumpFile );
    int SaveFile();
    bool IsEmpty();

    int Load(const char* szVar, double* pVar);
    int Load(const char* szVar, CFmMatrix* pVar);
    int Load(const char* szVar, CFmVector* pVar);
    int Load(const char* szVar, CFmVectorStr* pVar);

    void Save(const char* szVar, double pval );
    void Save(const char* szVar, CFmMatrix* pMat );
    void Save(const char* szVar, CFmVector* pVec );
    void Save(const char* szVar, CFmVectorStr* pVecStr );

private:
    RTVARFMT* FindHeader(const char* szVar);
    bool WriteStrlistToFile(CFmVectorStr* pList, FILE* fp);
    int fwrite_hugedata(void* pData, size_t nDatLen, FILE* fp);

    std::map<char*, double> m_DoubleMap;
    std::map<char*, CFmMatrix*> m_MatMap;
    std::map<char*, CFmVector*> m_VecMap;
    std::map<char*, CFmVectorStr*> m_VecStrMap;

    char* m_szRtFile;
    FILE* m_fp;
    bool m_bSaved;
    bool m_bEmpty;
    RTDUMPFMT m_fmt;
};

#endif
