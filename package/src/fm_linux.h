// linux_diff.h: interface for the functions which are different
//                between LINUX AND WINDOWS
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_LINUX_H_)
#define _FM_LINUX_H_

#ifdef __cplusplus
extern "C" {
#endif

//#ifdef LINUX

#define MAX_PATH 4096

#define Strdup(X) strcpy(Calloc(strlen(X)+1, char), X) ;

int stricmp(const char* s1, const char* s2);

//#endif

#ifdef __cplusplus
}
#endif


#endif

