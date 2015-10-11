#include <exception> // for std::bad_alloc
#include <new>
#include <R.h>
#include "fm_new.h"


void* operator new (size_t size)
{
 	void *p = Calloc(size, char);
 	if (p==0)
 	 	throw std::bad_alloc();
//Rprintf("new op=%X, %d\n", p, size);
 	return p;
}


void operator delete (void *p)
{
//Rprintf("delete op=%X\n", p);
 	Free(p);
}

void* operator new[] (size_t size)
{
 	void *p = Calloc( size, char );
 	if (p==0)
 	 	throw std::bad_alloc();
//Rprintf("new[] op=%X, %d\n", p, size);
 	return p;
}

void operator delete[] (void *p)
{
//Rprintf("delete[] op=%X\n", p);
 	Free(p);
}

void *operator new (std::size_t size, CFmNewTemp & arena)
{
    void* p = arena.allocate(size) ;
//Rprintf("new A=%X, %d\n", p, size);
	return(p);
}

void operator delete (void * p, CFmNewTemp  & arena)
{
//Rprintf("delete A=%X\n", p);
    arena.deallocate(p) ;
}

