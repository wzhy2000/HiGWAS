#ifndef _FM_FILEMATRIX_H_
#define _FM_FILEMATRIX_H_

#include <stdint.h>

typedef struct fmmat_file_fmt
{
    char szHeader[128];
    uint64_t nNumCols;
    uint64_t nNumRows;
}FMMAT_FILE_FMT;

class CFmMatrix;
class CFmVector;

class CFmFileMatrix
{
public:
    CFmFileMatrix(char* szFile, bool bAppend, bool bAutodel);
    virtual ~CFmFileMatrix();

    char* GetFileName();

    int Rbind(CFmVector& vct);
    int GetRow(int nRow, CFmVector& vct);
    int GetCol(int nCol, CFmVector& vct);
    int GetCacheCol(int nCol, CFmVector& vct);
    double GetAt(int nRow, int nCol);
    int SetAt(int nRow, int nCol, double fval);
    int GetNumCols();
    int GetNumRows();
    bool Update();

protected:
    int OpenFile(bool bAppend);
    int RewriteHeader();

private:
    int m_nCacheCol0;
    int m_nCacheCol1;
    CFmMatrix* m_pmatCache;
    int m_nNumCols;
    int m_nNumRows;
    FILE* m_pFile;
    char* m_pszFile;
    bool m_bAutodel;
};

CFmMatrix* LoadFileMatrix(char* szFile);

void destroy(CFmFileMatrix* p);

#endif // _FM_FILEMATRIX_H_
