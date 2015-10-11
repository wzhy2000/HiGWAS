/* BLS_res.cpp  -	BLS Result Object
 *
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>
#include <Rmath.h>

#include "fm_matrix.h"
#include "fm_vector.h"
#include "fm_vector_str.h"
#include "fm_rlogger.h"
#include "fm_err.h"
#include "fm_linux.h"
#include "fm_rtdump.h"


#define _t(x) ((x).GetTransposed())

CFmRtDump::CFmRtDump()
{
    m_szRtFile = NULL;
    m_fp = NULL;
    m_bSaved = FALSE;
    m_bEmpty = TRUE;
}

CFmRtDump::~CFmRtDump()
{
    if (m_bSaved)
        SaveFile();

    if (m_szRtFile)
        Free(m_szRtFile);
}


void CFmRtDump::Save(const char* szVar, double val  )
{
    m_DoubleMap.insert ( std::pair<char*,double>( (char*)szVar, val ) );
    m_bSaved = TRUE;
}

void CFmRtDump::Save(const char* szVar, CFmMatrix* pMat  )
{
    m_MatMap.insert ( std::pair<char*,CFmMatrix*>( (char*)szVar, pMat ) );
    m_bSaved = TRUE;
}

void CFmRtDump::Save(const char* szVar, CFmVector* pVec  )
{
    m_VecMap.insert ( std::pair<char*,CFmVector*>( (char*)szVar, pVec ) );
    m_bSaved = TRUE;
}

void CFmRtDump::Save(const char* szVar, CFmVectorStr* pVecStr )
{
    m_VecStrMap.insert ( std::pair<char*,CFmVectorStr*>( (char*)szVar, pVecStr ) );
    m_bSaved = TRUE;
}

void SetVarInfo(RTDUMPFMT* pFmt, int nIdx, CFmMatrix* pMat, const char* szVar,  int nPos, unsigned int nType )
{
    _log_debug(_HI_, "SaveVarInfo: Var=%s, Matrix=[%d,%d] nPos=%d", szVar, pMat->GetNumRows(), pMat->GetNumCols(), nPos );

    strcpy( pFmt->varCont[ nIdx ].varName, szVar );
    pFmt->varCont[ nIdx ].varType = nType;
    pFmt->varCont[ nIdx ].dimRow  = pMat->GetNumRows();
    pFmt->varCont[ nIdx ].dimCol  = pMat->GetNumCols();
    pFmt->varCont[ nIdx ].nBytes  = pMat->GetBytes();
    pFmt->varCont[ nIdx ].nPos    = nPos;
}

void SetVarInfo(RTDUMPFMT* pFmt, int nIdx, CFmVector* pVct, const char* szVar,  int nPos, unsigned int nType, double fVal  )
{
    if (pVct)
        _log_debug(_HI_, "SaveVarInfo: Var=%s, Vector=[%d], nPos=%d", szVar, pVct->GetLength(), nPos);
    else
        _log_debug(_HI_, "SaveVarInfo: Var=%s, val=%f", szVar, fVal );

    strcpy( pFmt->varCont[ nIdx ].varName, szVar );
    pFmt->varCont[ nIdx ].varType = nType;
    if (pVct)
    {
        pFmt->varCont[ nIdx ].dimRow  = pVct->GetLength();
        pFmt->varCont[ nIdx ].nBytes  = pVct->GetBytes();
    }
    else
    {
        pFmt->varCont[ nIdx ].dimRow  = 1;
        pFmt->varCont[ nIdx ].nBytes  = 0;
    }

    pFmt->varCont[ nIdx ].dimCol  = 1;
    pFmt->varCont[ nIdx ].fValue  = fVal;
    pFmt->varCont[ nIdx ].nPos    = nPos;
}

void SetVarInfo(RTDUMPFMT* pFmt, int nIdx, CFmVectorStr* pVct, const char* szVar, int nPos, unsigned int nType)
{
    _log_debug(_HI_, "SaveVarInfo: Var=%s, VectorStr=[%d, strlen=%d] nPos=%d", szVar, pVct->GetLength(), pVct->GetMaxStrlen(), nPos );

    strcpy( pFmt->varCont[ nIdx ].varName, szVar );
    pFmt->varCont[ nIdx ].varType = nType;
    pFmt->varCont[ nIdx ].dimRow  = pVct->GetLength();
    pFmt->varCont[ nIdx ].dimCol  = pVct->GetMaxStrlen()+1;
    pFmt->varCont[ nIdx ].nBytes  = pVct->GetBytes();
    pFmt->varCont[ nIdx ].nPos    = nPos;
}

bool CFmRtDump::WriteStrlistToFile(CFmVectorStr* pList, FILE* fp)
{
    int nMaxStrlen = pList->GetMaxStrlen();
    char* pBuf = new char[ nMaxStrlen + 1 ];

    for(int i=0; i<pList->GetLength(); i++)
    {
        memset(pBuf, 0, nMaxStrlen+1);
        if ( pList->Get(i) != NULL )
            strcpy(pBuf, pList->Get(i) );

        if( fwrite(pBuf,  nMaxStrlen+1, 1, fp) != 1)
        {
            _log_error( _HI_, "WriteStrlistToFile: Failed to write items of the CFmVectorStr");

            delete[] pBuf;
            return (false);
        }
    }

    delete[] pBuf;
    return(true);
}

int CFmRtDump::fwrite_hugedata(void* pData, size_t nDatLen, FILE* fp)
{
    size_t nBlock = 256*256; /* 65536 bytes*/
    size_t nSize = 0;
    while( nSize < nDatLen )
    {
        if ( nSize + nBlock > nDatLen)
            nBlock = nDatLen - nSize;

        size_t nRet = fwrite( (char*)pData+nSize, nBlock, 1, fp );
        if ( nRet!= 1 )
        {
            _log_info( _HI_, "Failed to write a huge data, Error code:%d", ferror(fp) );
            return(ERR_WRITE_FILE);
        }

        nSize += nBlock;
    }

    return(0);
}

int CFmRtDump::SaveFile( )
{
    if (m_fp)
    {
        fclose(m_fp);
        m_fp = NULL;
    }

    if (!m_bSaved) return(0);

    _log_info(_HI_, "SaveFile: Start to save the results to file(%s).", m_szRtFile );

    memset(&m_fmt, 0, sizeof(RTDUMPFMT));
    strcpy(m_fmt.szHeader, "fGWAS1.0,32\r\nURL: Http://statgen.psu.edu/\r\nCenter of Statistical Genetics, Penn State University.\n");

    int nIndex = 0;

    std::map<char*, double>::iterator doubleIt;
    DWORD nPos = sizeof(RTDUMPFMT);
    for ( doubleIt=m_DoubleMap.begin() ; doubleIt != m_DoubleMap.end(); doubleIt++ )
    {
        char* szVar = (*doubleIt).first;
        double val = (*doubleIt).second;

        SetVarInfo( &m_fmt, nIndex, (CFmVector*)NULL, szVar, nPos, VAR_TYPE_DOUBLE, val);
        nIndex++;
    }

    std::map<char*,CFmVectorStr*>::iterator vecStrIt;
    for ( vecStrIt=m_VecStrMap.begin() ; vecStrIt != m_VecStrMap.end(); vecStrIt++ )
    {
        char* szVar = (*vecStrIt).first;
        CFmVectorStr* pVec = (*vecStrIt).second;

        SetVarInfo( &m_fmt, nIndex, pVec, szVar, nPos, VAR_TYPE_STRING);
        nPos += m_fmt.varCont[nIndex].nBytes;
        nIndex++;
    }

    std::map<char*,CFmVector*>::iterator vecIt;
    for ( vecIt=m_VecMap.begin() ; vecIt != m_VecMap.end(); vecIt++ )
    {
        char* szVar = (*vecIt).first;
        CFmVector* pVec = (*vecIt).second;

        SetVarInfo( &m_fmt, nIndex, pVec, szVar, nPos, VAR_TYPE_DOUBLE, 0);
        nPos += m_fmt.varCont[nIndex].nBytes;
        nIndex++;
    }

    std::map<char*,CFmMatrix*>::iterator matIt;
    for ( matIt=m_MatMap.begin() ; matIt != m_MatMap.end(); matIt++ )
    {
        char* szVar = (*matIt).first;
        CFmMatrix* pMat = (*matIt).second;

        SetVarInfo( &m_fmt, nIndex, pMat, szVar, nPos, VAR_TYPE_DOUBLE);
        nPos += m_fmt.varCont[nIndex].nBytes;
        nIndex++;
    }

    _log_info(_HI_, "SaveFile: variable count: %d.", nIndex );

    FILE* fp;
    if ((fp=fopen( m_szRtFile, "wb"))==NULL)
    {
        _log_error(_HI_, "Failed to open a DUMP file, file name=%s", m_szRtFile);
        return( ERR_CREATE_FILE );
    }

    fwrite(&m_fmt, sizeof(RTDUMPFMT),1 , fp);

    for ( vecStrIt=m_VecStrMap.begin() ; vecStrIt != m_VecStrMap.end(); vecStrIt++ )
    {
        CFmVectorStr* pVecStr = (*vecStrIt).second;
        WriteStrlistToFile( pVecStr, fp );
    }

    for ( vecIt=m_VecMap.begin() ; vecIt != m_VecMap.end(); vecIt++ )
    {
        CFmVector* pVec = (*vecIt).second;
        fwrite_hugedata(pVec->GetData(), pVec->GetBytes(), fp);
    }

    for ( matIt=m_MatMap.begin() ; matIt != m_MatMap.end(); matIt++ )
    {
        CFmMatrix* pMat = (*matIt).second;
        fwrite_hugedata(pMat->GetData(), pMat->GetBytes(), fp);
    }

    fclose(fp);

    m_bSaved = FALSE;

    _log_info(_HI_, "SaveFile: Successful to save the results to file(%s).", m_szRtFile );

    return(0);
}

RTVARFMT* CFmRtDump::FindHeader(const char* szVar)
{
    for (int i=0;i<MAX_CONT_ITEM;i++)
        if ( strcmp( m_fmt.varCont[i].varName, szVar ) == 0 )
            return &(m_fmt.varCont[i]);

    return NULL;
}

int CFmRtDump::Load(const char* szVar, CFmVectorStr* pVct )
{
    RTVARFMT* pHeader = FindHeader(szVar);
    if (!pHeader)
        return(-1);

    if ( pHeader->varType != VAR_TYPE_STRING )
        return( ERR_FILE_DATA);

    _log_info(_HI_, "Load: Read var(%s) from the position %d ", pHeader->varName, pHeader->nPos);

    pVct->Reset( pHeader->dimRow );

    fseek( m_fp, pHeader->nPos, SEEK_SET );
    char* pTemp = new char[ pHeader->dimCol + 1];
    if (!pTemp)
    {
        _log_error(_HI_, "Load: Filed to allocate %d bytes to a charaters buffer", pHeader->dimCol + 1 );
        return( ERR_MEM_ALLOC );
    }

    for (unsigned int i=0; i<pHeader->dimRow; i++)
    {
        if ( fread(pTemp, 1 , pHeader->dimCol, m_fp)== pHeader->dimCol )
        {
            pVct->Set(i, pTemp);
        }
        else
        {
            _log_error(_HI_, "Load: Filed to read a result file, nPos=%d, nBytes=%d.", pHeader->nPos, pHeader->nBytes);
            delete[] pTemp;
            return( ERR_READ_FILE );
        }
    }

    delete[] pTemp;
    return(0);
}

int CFmRtDump::Load(const char* szVar, CFmVector* pVct)
{
    RTVARFMT* pHeader = FindHeader(szVar);
    if (!pHeader)
        return(-1);

    if (pHeader->varType != VAR_TYPE_DOUBLE )
        return( ERR_FILE_DATA);

    _log_info(_HI_, "Load: Read var(%s) from the position %d ", pHeader->varName, pHeader->nPos);

    pVct->Resize(pHeader->dimRow);

    fseek( m_fp, pHeader->nPos, SEEK_SET);
    unsigned int nBytes = pHeader->nBytes;
    if ( fread(pVct->GetData(), 1 , nBytes, m_fp) != nBytes )
    {
        _log_error(_HI_, "Load: Filed to read a result file, nPos=%d, nBytes=%d.", pHeader->nPos, nBytes);
        return( ERR_READ_FILE );
    }

    return(0);
}

int CFmRtDump::Load(const char* szVar, CFmMatrix* pMat)
{
    RTVARFMT* pHeader = FindHeader(szVar);
    if (!pHeader)
        return(-1);

    if (pHeader->varType != VAR_TYPE_DOUBLE )
        return( ERR_FILE_DATA);

    _log_info(_HI_, "Load: Read var(%s) from the position %d ", pHeader->varName, pHeader->nPos);

    pMat->Resize(pHeader->dimRow, pHeader->dimCol, true);

    fseek( m_fp, pHeader->nPos, SEEK_SET);
    unsigned int nBytes = pHeader->nBytes;
    if ( fread(pMat->GetData(), 1 , nBytes, m_fp)!=nBytes)
    {
        _log_error(_HI_, "Load: Filed to read a result file, nPos=%d, nBytes=%d.", pHeader->nPos, nBytes);
        return( ERR_READ_FILE );
    }

    return(0);
}

int CFmRtDump::Load(const char* szVar, double* pVal)
{
    RTVARFMT* pHeader = FindHeader(szVar);
    if (!pHeader)
        return(-1);

    if (pHeader->varType != VAR_TYPE_DOUBLE )
        return( ERR_FILE_DATA);

    _log_info(_HI_, "Load: Read var(%s): %f ", pHeader->varName, pHeader->fValue);

    if (pHeader->dimCol==1 && pHeader->dimRow==1)
    {
        *pVal = pHeader->fValue;
        return(0);
    }
    else
    {
        _log_error(_HI_, "Load: Filed to read a double value, Row,col=(%d,%d)",pHeader->dimRow, pHeader->dimCol );
        return( ERR_READ_FILE );
    }

}

int CFmRtDump::OpenFile( char *szRtdumpFile )
{
    _log_info(_HI_, "OpenFile: Open a result file , file name=%s", szRtdumpFile);

    m_szRtFile = Strdup(szRtdumpFile);

    if ((m_fp=fopen(m_szRtFile, "r+b"))==NULL)
    {
        _log_error(_HI_, "OpenFile: Failed to open a result file , file name=%s", m_szRtFile);
        m_bEmpty = TRUE;
        return( ERR_OPEN_FILE );
    }

    if ( fread(&m_fmt, 1, sizeof(m_fmt), m_fp)!= sizeof(m_fmt))
    {
        _log_error(_HI_, "OpenFile: Failed to read a result file , file name=%s", m_szRtFile);
        fclose(m_fp);
        m_fp = NULL;
        return( ERR_READ_FILE );
    }

    m_bEmpty = FALSE;

    _log_info(_HI_, "OpenFile: Open a result file successful");

    return(0);
}

bool CFmRtDump::IsEmpty()
{
    return m_bEmpty;
}


