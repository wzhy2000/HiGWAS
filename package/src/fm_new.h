// fm_new.h: interface for new/delete
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_NEW_H_)
#define _FM_NEW_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <R.h>
#include <cstddef>

class CFmNewTemp {
public:
    void * allocate ( size_t size)
    {
        void* p = Calloc( size, char );
        return(p);
    }
    void deallocate ( void * p)
    {
        Free(p);
    }
};

void operator delete (void * p, CFmNewTemp  & arena);
void* operator new (std::size_t size, CFmNewTemp  & arena);

#ifdef __cplusplus
}
#endif

#endif


