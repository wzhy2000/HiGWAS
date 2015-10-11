// fm_sys.h: interface for the CFmSys class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_RSOURCE_H_)
#define _FM_RSOURCE_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CFmVector;
class CFmMatrix;
class CFmVectorStr;

class CFmRSource
{
public:
    CFmRSource();
    virtual ~CFmRSource();

    int Run(char* srcFile);
    int Run1(char* szCmd);
    int Run2(char* srcFile, char* szCmd);

    int SetGlobalVariable(const char* szVar, CFmVector* pVct);
    int SetGlobalVariable(const char* szVar, CFmMatrix* pMat);
    int SetGlobalVariable(const char* szVar, CFmVectorStr* pVctStr);

    int GetGlobalVector(const char* szVar, CFmVector* pVct);
    int GetGlobalMatrix(const char* szVar, CFmMatrix* pMat);
};

SEXP run_script(char* szSrc);

#endif
