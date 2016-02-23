#include <ctype.h>
#include <string.h>

#include "fm_linux.h"
#include "fm_dataframe.h"
#include "fm_matrix.h"
#include "fm_vector_str.h"
#include "fm_rlogger.h"
#include "fm_err.h"
#include "fm_new.h"


#define MAX_LINE_WIDTH 256*256

CFmDataFrame::CFmDataFrame()
{
    m_nNumRows = 0;
    m_nNumCols = 0;
    m_nMaxRows = 0;
    m_nMaxCols = 0;
    m_pData    = NULL;
    m_pRowNames = NULL;
    m_pColNames = NULL;
}

CFmDataFrame::~CFmDataFrame()
{
    m_nNumRows = 0;
    m_nNumCols = 0;

	if (m_pRowNames) destroy( m_pRowNames );
    if (m_pColNames) destroy( m_pColNames );

    for (int i=0; i<m_nNumRows; i++)
        for (int j=0; j<m_nNumCols; j++)
        {
			char* p = Get(i,j);
			if(p) Free( p );
		}

    Free( m_pData );
}

CFmMatrix* CFmDataFrame::GetMatrix(CFmVector* pVctCols=NULL)
{
	CFmVector fmCols(pVctCols);
	if(fmCols.GetLength()==0)
		for(int i=0;i<m_nNumCols;i++)
			fmCols.Put(i);

	CFmNewTemp refNew;
	CFmMatrix* pMat = new (refNew) CFmMatrix(m_nNumRows, fmCols.GetLength());

	CFmVectorStr fmNameCols(0);

	for(int i=0; i<fmCols.GetLength(); i++)
	{
		CFmVector& pVct = GetFloatCol( fmCols.Get(i) );
		pMat->SetCol(i, pVct);

		if(m_pColNames)
			fmNameCols.Put( m_pColNames->Get( fmCols.Get(i) ) );
	}

	if(m_pRowNames)	pMat->SetRowNames( m_pRowNames);
	if(m_pColNames)	pMat->SetColNames( &fmNameCols );

    return(pMat);
}

bool CFmDataFrame::Set(int nRow, int nCol, char* value)
{
    if ( nCol<m_nNumCols && nRow < m_nNumRows )
    {
        m_pData[  nRow*m_nNumCols + nCol] = Strdup(value) ;
        return true ;
    }
    else
    {
        Rprintf("Error in DataFrame.set");
        throw("Error in DataFrame.set");
        return false;
    }
}

char* CFmDataFrame::Get(int nRow, int nCol )
{
    if ( nCol<m_nNumCols && nRow < m_nNumRows )
    {
        return( m_pData[ nRow*m_nNumCols + nCol] ) ;
    }
    else
    {
        Rprintf("Error in DataFrame.set");
        throw("Error in DataFrame.set");
        return NULL;
    }
}

int CFmDataFrame::GetNumRow()
{
    return(m_nNumRows);
}

int CFmDataFrame::GetNumCol()
{
    return(m_nNumCols);
}

int CFmDataFrame::GetMaxRow()
{
    return(m_nMaxRows);
}

int CFmDataFrame::GetMaxCol()
{
    return(m_nMaxCols);
}

char* CFmDataFrame::GetRowName(int idx)
{
    return (m_pRowNames->Get( idx ) );
}

char* CFmDataFrame::GetColName(int idx)
{
    return (m_pColNames->Get( idx ) );
}

int CFmDataFrame::FindColumn(const char* szColumn)
{
	if(m_pColNames==NULL)
		return(-1);

    for (int i=0;i<m_nNumCols; i++)
    {
        if ( strcasecmp ( m_pColNames->Get(i), szColumn)==0 )
            return(i);
    }

    return( -1 );
}

CFmVector& CFmDataFrame::GetFloatCol(int idx)
{
	CFmNewTemp refNew;
    CFmVector* pNew= new (refNew) CFmVector( 0, 0.0 );
    for (int i=0;i<m_nNumRows; i++)
    {
        char* szVal = Get(i, idx) ;
        char* pEnd = NULL;
        double f = strtod ( szVal, &pEnd);
        if ( f==0.0 )
        {
            if (pEnd <= szVal || pEnd > szVal + strlen(szVal) )
                pNew->Put( R_NaN );
            else
                pNew->Put( f );
        }
        else
            pNew->Put( f );
    }

    return *pNew;
}

CFmVectorStr* CFmDataFrame::GetStringCol(int idx)
{
	CFmNewTemp refNew;
    CFmVectorStr* pNew= new (refNew) CFmVectorStr( 0, 100 );
    for (int i=0;i<m_nNumRows; i++)
        pNew->Put( Get(i, idx) );

    return(pNew);
}

CFmVectorStr* CFmDataFrame::RGetRowNames()
{
    return(m_pRowNames);
}

CFmVectorStr* CFmDataFrame::RGetColNames()
{
    return(m_pColNames);
}

int CFmDataFrame::AllocMemory( int nMaxRows, int nMaxCols )
{
    if( nMaxRows * nMaxCols * sizeof(char*)>1024*1024)
        _log_debug( _HI_ , "MEMORY: Try to allocate big memory(%.3fM bytes)", (nMaxRows * nMaxCols + 1)*sizeof(double)/1024.0/1024.0);

    char** pData = (char**) Calloc( nMaxRows * nMaxCols + 1, char* ) ;
    if (pData==NULL)
    {
        _log_fatal( _HI_ , "MEMORY: failed to allocate %d bytes to CFmMatrix[%d,%d].", (nMaxRows * nMaxCols + 1)*sizeof(double), nMaxRows, nMaxCols );
        return(-1);
    }

    memset(pData, 0, sizeof(char*) * (nMaxCols * nMaxRows + 1)) ;

    if (m_pData)
    {
        memcpy(pData, m_pData, sizeof(char*) * (m_nNumCols * m_nNumRows) );
        Free( m_pData );
    }

    m_pData = pData;
    m_nMaxRows = nMaxRows;
    m_nMaxCols = nMaxCols;

    return(0);
}

char* _strsep(char **stringp, const char *delim)
{
    char *s;
    const char *spanp;
    int c, sc;
    char *tok;

    if ((s = *stringp) == NULL)
        return (NULL);
    for (tok = s;;) {
        c = *s++;
        spanp = delim;
        do {
            if ((sc = *spanp++) == c) {
                if (c == 0)
                    s = NULL;
                else
                    s[-1] = 0;
                *stringp = s;
                return (tok);
            }
        } while (sc != 0);
    }
    /* NOTREACHED */
}


// this function reads a file for "," separated values and creates a matrix object from it
int CFmDataFrame::Load(const char* filename, bool bRowName, bool bColName)
{
    _log_debug(_HI_, "CFmDataFrame::Load:%s, bRowName=%d, bColName=%d", filename, bRowName, bColName);

    FILE* fp = fopen( filename, "rt");
    if (fp==NULL)
    {
        _log_error(_HI_, "The file can not be opened(%s)", filename);
        return(ERR_OPEN_FILE);
    }

    char aLine[MAX_LINE_WIDTH+1]={0};
    if( fgets( aLine, MAX_LINE_WIDTH, fp ) == NULL)
    {
        _log_error( _HI_, "Failed to read a line (%s)", filename);

        fclose(fp);
        return(ERR_READ_FILE);
    }

    CFmVectorStr vctRowNames(0);
    CFmVectorStr vctColNames(0);
    int col_size = 0;

    // count how many elements in the first lines
    const char seps1[] = ",\t\n\r";

    char* running = Strdup (aLine);
    char* runnptr = running;
    char* token = _strsep((char**)&running, seps1 );
    while(token)
    {
		if(strlen(token)!=0)
		{
			vctColNames.Put( token);
		    col_size++;
        }

    	token = _strsep( (char**)&running, seps1 );
    }

    Free(runnptr);

    _log_debug(_HI_, "CFmDataFrame::Load:  col_size: %d", col_size);

    CFmVector tmpRow(col_size, 0.0);
    if (bRowName && !bColName)
        tmpRow.Resize(col_size-1);

    int nMaxRows = 20;
    int nMaxCols = tmpRow.GetLength();
    m_nNumCols = tmpRow.GetLength();
    m_nNumRows = 0;

    AllocMemory( nMaxRows, nMaxCols );

    if (!bColName)
    {
        if (!bRowName)
        {
            m_nNumCols = col_size;
            for (int i=0; i<col_size; i++)
                Set(0, i, vctColNames[i]);
        }
        else
        {
            m_nNumCols = col_size-1;
            for (int i=1; i<col_size+1; i++)
                Set(0, i-1, vctColNames[i]);
            vctRowNames.Put( vctColNames[0] );
        }

        m_nNumRows++;
    }

    CFmVectorStr vctTmpRow(0);
    while(!feof(fp))
    {
        if( fgets( aLine, MAX_LINE_WIDTH, fp ) == NULL)
            break;

        char* running = Strdup (aLine);

        vctTmpRow.Reset(0);
        char* token = _strsep( &running, seps1 );
        int nNonEmty = 0;

        while(token)
        {
            if (strlen(token)!=0 &&
                strncmp(token, "\n", 1)!=0 &&
                strncmp(token, "\r", 1)!=0 &&
                strncmp(token, "\t", 1)!=0)
                nNonEmty ++;
            vctTmpRow.Put(token);
            token = _strsep( &running, seps1 );
        }

        Free(running);
        if (nNonEmty==0)
            break;

        m_nNumRows++;
        if (m_nNumRows >= m_nMaxRows )
            AllocMemory( m_nNumRows+20, m_nMaxCols );

        if (!bRowName)
        {
            for (int i=0; i<col_size && i<vctTmpRow.GetLength(); i++)
                Set( m_nNumRows-1, i,  vctTmpRow[i] );
        }
        else
        {
            for (int i=1; i<col_size+1 && i<vctTmpRow.GetLength(); i++)
            {
                Set( m_nNumRows-1, i-1,  vctTmpRow[i] );
            }

            vctRowNames.Put( vctTmpRow[0]);
        }
    }

    fclose(fp);

	CFmNewTemp refNew;
    if (bRowName)
        m_pRowNames = new (refNew) CFmVectorStr( &vctRowNames );

    if (bColName)
        m_pColNames = new (refNew) CFmVectorStr( &vctColNames );

    _log_debug(_HI_, "CFmDataFrame::Load:%d,%d", m_nNumRows, m_nNumCols);

    return(0);
}


char *ltrim(char *s)
{
    while(isspace(*s))
        s++;
    return s;
}

char *rtrim(char *s)
{
    char* back = s + strlen(s);
    while(isspace(*--back));
    *(back+1) = '\0';
    return s;
}

char *trim(char *s)
{
    return rtrim(ltrim(s));
}

void destroy(CFmDataFrame* p)
{
	CFmNewTemp  fmRef;
	p->~CFmDataFrame();
	operator delete(p, fmRef);
}
