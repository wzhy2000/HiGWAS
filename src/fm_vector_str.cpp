/* fm_vector_str.cpp  -	CFmVectorStr Class
 *
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>

#include "fm_vector_str.h"
#include "fm_rlogger.h"
#include "fm_err.h"
#include "fm_new.h"

#define Strdup(X) strcpy(Calloc(strlen(X)+1, char), X) ;

CFmVectorStr::CFmVectorStr(CFmVectorStr* pOther)
{
    m_pNames = NULL;
    m_nMaxLen = pOther->m_nMaxLen;
    m_nActLen = pOther->m_nActLen;
    m_nMaxStrlen = 0;

    m_pData = AllocateMemory( m_nMaxLen+1 );
    for (int i=0; i<pOther->GetLength(); i++)
    {
        m_pData[i] = Strdup(pOther->Get(i));
        if ( strlen(pOther->Get(i)) > m_nMaxStrlen )
            m_nMaxStrlen = strlen(pOther->Get(i));
    }
}

CFmVectorStr::CFmVectorStr(int nSize, int nMaxSize)
{
    m_pNames = NULL;
    m_nMaxStrlen = 0;
    m_nActLen= nSize;
	if (nMaxSize<=0)
        m_nMaxLen= m_nActLen>100?m_nActLen:100;
	else
		m_nMaxLen= nMaxSize;

    if (m_nActLen>m_nMaxLen)
        m_nMaxLen = m_nActLen;

    m_pData = AllocateMemory( m_nMaxLen+1 );
}

CFmVectorStr::~CFmVectorStr()
{
    for (int i=0; i<m_nActLen; i++)
    {
		if (m_pData[i])
            Free(m_pData[i]);
	}

    Free(m_pData);
}

void CFmVectorStr::Reset(int size)
{
    for (int i=0; i<m_nActLen; i++)
    {
        if (m_pData[i])
            Free(m_pData[i]);
        m_pData[i] = NULL;
        m_nMaxStrlen = 0;
    }

    if (m_nMaxLen<size)
    {
        Free(m_pData);
        m_nMaxLen = size;
        m_pData = AllocateMemory( m_nMaxLen+1 );
    }

    m_nActLen = size;
}

void CFmVectorStr::Resize(int size)
{
    if (m_nActLen==size)
    {
    }
    else
    if (m_nActLen<size)
    {
        if (m_nMaxLen>=size)
            m_nActLen = size;
        else
        {
            m_nMaxLen = size;
            char** ppData = AllocateMemory( m_nMaxLen+1 );
            memcpy(ppData, m_pData, sizeof(char*) * (m_nActLen) );
            Free( m_pData );
            m_pData = ppData;
            m_nActLen = size;
        }
    }
    else
    {
        for (int i=size; i<m_nActLen; i++)
        {
            if (m_pData[i])
                Free(m_pData[i]);
            m_pData[i] = NULL;
        }

        m_nActLen = size;
    }
}


char** CFmVectorStr::GetData()
{
    return m_pData;
}

int	CFmVectorStr::GetLength()
{
    return m_nActLen;
}

int	CFmVectorStr::GetMaxStrlen()
{
    return m_nMaxStrlen;
}

int	CFmVectorStr::GetBytes()
{
    return m_nActLen*(m_nMaxStrlen+1);
}

char* CFmVectorStr::Get(int idx)
{
    if (idx>=m_nActLen || idx<0)
        throw( "wrong index of vector in CFmVectorStr::Get()");

    return (m_pData[idx]);
}

void CFmVectorStr::Set(int idx, char* str)
{
    if (idx>=m_nActLen || idx<0)
        throw(  "wrong index of vector in CFmVectorStr::Set()");

    if (m_pData[idx]) Free(m_pData[idx]);
    if(str ==NULL || strlen(str)==0 )
        return;

    if(strlen(str)>m_nMaxStrlen)
        m_nMaxStrlen = strlen(str);

    m_pData[idx] = Strdup(str);
}

char* CFmVectorStr::operator[](int idx) const
{
    if (idx>=m_nActLen || idx<0)
        throw( "wrong index of vector in CFmVectorStr[]");

    return (m_pData[idx]);
}

void CFmVectorStr::Put(char* str)
{
    if (m_nActLen>=m_nMaxLen)
    {
        int nInc = m_nMaxLen/5>20 ? m_nMaxLen/5 : 20;
        m_nMaxLen += nInc;
        char** pData = AllocateMemory( m_nMaxLen+1 );
        memcpy(pData, m_pData, sizeof(char*) * m_nActLen );
        Free( m_pData );
        m_pData = pData;
    }

    if (strlen(str)>m_nMaxStrlen)
        m_nMaxStrlen = strlen(str);

    m_pData[m_nActLen] = Strdup(str);
    m_nActLen++;
}

void CFmVectorStr::Append(CFmVectorStr* another)
{
    if (another==NULL)
    {
        _log_debug( _HI_, "Append(), another is NULL pointer.");
        return;
    }

    for(int i=0;i<another->GetLength(); i++)
        Put(another->Get(i));
}

bool CFmVectorStr::Remove(int idx)
{
    if (idx<0 || idx>=m_nActLen)
        return false;

    if ( m_pData[idx] ) Free(m_pData[idx]);

    if ( idx < (m_nActLen-1) )
        memmove( &m_pData[idx], &m_pData[idx+1], (m_nActLen-idx-1) * sizeof(char*) );

	m_pData[ m_nActLen-1 ] = NULL;
    m_nActLen--;

    return(true);
}

char** CFmVectorStr::AllocateMemory( int nLen )
{
    char** pData = Calloc( nLen+1, char*);
//**MEMTEST: Rprintf("ALLOC MEMORY: %d(CFmVectorStr)\n", (nLen + 1));

    if (pData==NULL)
    {
        _log_fatal( _HI_, "MEMORY: failed to allocate %d bytes to CFmVectorStr.", nLen+1);
    }

    memset(pData, 0, (nLen+1)*sizeof(char*) ) ;
    return pData ;
}

bool CFmVectorStr::RemoveElements(CFmVector& nRows)
{
    for (int i=0; i<nRows.GetLength(); i++)
    {
        if ( nRows.Get(i)<0 || nRows.Get(i)>=m_nActLen )
            return(false);
    }

    int k=0;
    for (int i=0; i<m_nActLen; i++)
    {
        if (nRows.Find(i)>=0)
        {
            //delete string data
            if (m_pData[i]) Free( m_pData[i] );
        }
        else
        {
            if (k!=i)
            {
                m_pData[k] = m_pData[i];
            }

            k++;
        }
    }

    for(;k<m_nActLen;k++)
        m_pData[k] = 0;

    m_nActLen -= nRows.GetLength();
    _log_debug(_HI_, "Length(%d).",  m_nActLen );

    return(true);
}

void CFmVectorStr::Show(const char* szName )
{
    Rprintf("CFmVectorStr[%d]: %s\n", m_nActLen, szName);
    for(int i=0;i<m_nActLen; i++)
    {
        if (m_pData[i])
            Rprintf("[%d]:(%d)%s\n", i, strlen(m_pData[i]), m_pData[i] );
        else
            Rprintf("[%d]:(0)\n");
    }

    Rprintf("\n");
}

int CFmVectorStr::Find(const char* szName )
{
    for(int i=0;i<m_nActLen; i++)
    {
        if (strcmp(m_pData[i], szName)==0)
            return(i);
    }

    return(-1);
}

void CFmVectorStr::SetNames( CFmVectorStr* pNames)
{
    if ( m_pNames ) destroy(m_pNames);
    m_pNames = pNames;
}

CFmVectorStr* CFmVectorStr::GetNames()
{
    return m_pNames;
}

char* CFmVectorStr::GetName(int idx )
{
    if (idx<0 ||idx>=m_nActLen)
        return NULL;

    if (m_pNames)
        return m_pNames->Get(idx);
    else
        return(NULL);
}

SEXP GetSEXP(CFmVectorStr* pVct)
{
    SEXP sep;
    PROTECT( sep = allocVector(STRSXP, pVct->GetLength() ) );
    for (int i=0; i<pVct->GetLength(); i++)
        SET_STRING_ELT( sep, i, mkChar( pVct->Get((i))) );

    if ( pVct->GetNames() )
    {
        SEXP exp_name;
        PROTECT(exp_name = allocVector(STRSXP, pVct->GetLength()));
        for (int i=0; i<pVct->GetLength(); i++)
            if (pVct->GetName(i))
                SET_STRING_ELT(exp_name, i, mkChar(pVct->GetName(i)));

        setAttrib(sep, R_NamesSymbol, exp_name);
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return sep;
}

int CFmVectorStr::WriteAsCSVFile(const char* szCsvFile, bool bAppend, const char* szTag)
{
    FILE* fp = NULL;
    if (!bAppend)
        fp = fopen( szCsvFile, "wt");
    else
        fp = fopen( szCsvFile, "a+");

    if (fp==NULL)
    {
        _log_error( _HI_ , "The file can not be created(%s)", szCsvFile);
        return(ERR_CREATE_FILE);
    }

    if (szTag!=NULL)
    {
        fprintf(fp, "\n>>-----------------%s--------------------<<\n", szTag);
    }

    for (int i = 0; i < m_nActLen ; i++)
    {
        fprintf(fp, "%s", m_pData[i] );
        fprintf(fp, "\n");
    }

    fclose(fp);
    return(0);
}

void destroy(CFmVectorStr* p)
{
	CFmNewTemp fmRef;
	p->~CFmVectorStr();
	operator delete(p, fmRef);
}
