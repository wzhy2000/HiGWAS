// fm_sys.h: interface for the CFmSys class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_SYS_H_)
#define _FM_SYS_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CFmSys
{
public:
    CFmSys();
    virtual ~CFmSys();

    static int GetTempId(char* szGrpId=NULL);
    static int GetTempFile(char* szFile, const char*szExt, int nLen);
    static int GetSiblingFile(char* szFile, const char* szOrgFile, int nNum, bool bRBase=FALSE);

protected:
    static int g_nSysPid;
    static int g_nTaskId;
};

bool _change_fileext(char* filename, const char* newext, char* newfile);

#endif
