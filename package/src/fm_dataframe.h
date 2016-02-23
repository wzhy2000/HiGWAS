#ifndef _FM_DATAFRAME_H_
#define _FM_DATAFRAME_H_

class CFmVectorStr;
class CFmVector;
class CFmMatrix;

class CFmDataFrame
{
public:
    CFmDataFrame();
    virtual ~CFmDataFrame();

    bool Set(int nRow, int nCol, char* value);
    char* Get(int nRow, int nCol );
    int GetNumRow();
    int GetNumCol();
    int GetMaxRow();
    int GetMaxCol();
    int FindColumn(const char* szColumn);
    char* GetRowName(int idx);
    char* GetColName(int idx);
    CFmVector& GetFloatCol(int idx);
    CFmVectorStr* GetStringCol(int idx);
    CFmVectorStr* RGetRowNames();
    CFmVectorStr* RGetColNames();

    int Load(const char* szFile, bool bRowName, bool bColName);

    CFmMatrix* GetMatrix(CFmVector* fmCols);

private:
    int AllocMemory( int nMaxRows, int nMaxCols );

    int m_nNumRows;
    int m_nNumCols;
    int m_nMaxRows;
    int m_nMaxCols;
    char** m_pData;
    CFmVectorStr* m_pRowNames;
    CFmVectorStr* m_pColNames;

};

void destroy(CFmDataFrame* p);

#endif // _FM_DATAFRAME_H_
