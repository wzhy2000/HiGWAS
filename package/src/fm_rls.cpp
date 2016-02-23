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
#include "fm_linux.h"
#include "fm_vector.h"
#include "fm_vector_str.h"
#include "fm_rlogger.h"
#include "fm_err.h"
#include "fm_rls.h"
#include "fm_new.h"



#define _t(x) ((x).GetTransposed())

CFmRls::CFmRls(const char* szRlsFile)
{
    m_szRlsFile = Strdup(szRlsFile);
}

CFmRls::~CFmRls()
{
    Free(m_szRlsFile);
}

CFmMatrix* CFmRls::GetMatrix(char* szVar)
{
    std::map<char*,CFmMatrix*>::iterator matIt;
    for ( matIt=m_MatMap.begin() ; matIt != m_MatMap.end(); matIt++ )
    {
        char *szVarX = (*matIt).first;
        if ( strcmp(szVarX, szVar)==0 )
            return( (*matIt).second );
    }
    return(NULL);
}

CFmVector* CFmRls::GetVector(char* szVar)
{
    std::map<char*,CFmVector*>::iterator matIt;
    for ( matIt=m_VecMap.begin() ; matIt != m_VecMap.end(); matIt++ )
    {
        char *szVarX = (*matIt).first;
        if ( strcmp(szVarX, szVar)==0 )
            return( (*matIt).second );
    }
    return(NULL);
}

CFmVectorStr* CFmRls::GetVectorStr(char* szVar)
{
    std::map<char*,CFmVectorStr*>::iterator matIt;
    for ( matIt=m_VecStrMap.begin() ; matIt != m_VecStrMap.end(); matIt++ )
    {
        char *szVarX = (*matIt).first;
        if ( strcmp(szVarX, szVar)==0 )
            return( (*matIt).second );
    }

    return(NULL);
}

void CFmRls::SetData(char* szVar, CFmMatrix* pMat  )
{
    m_MatMap.insert ( std::pair<char*,CFmMatrix*>( szVar, pMat ) );
}

void CFmRls::SetData(char* szVar, CFmVector* pVec  )
{
    m_VecMap.insert ( std::pair<char*,CFmVector*>( szVar, pVec ) );
}

void CFmRls::SetData(char* szVar, CFmVectorStr* pVecStr )
{
    m_VecStrMap.insert ( std::pair<char*,CFmVectorStr*>( szVar, pVecStr ) );
}

void SetVarInfo(RESFMT* pFmt, int nIdx, CFmMatrix* pMat, const char* szVar, unsigned int nType )
{
    _log_debug(_HI_, "SaveVarInfo: Var=%s, Matrix=[%d,%d]", szVar, pMat->GetNumRows(), pMat->GetNumCols() );

    strcpy( pFmt->varCont[ nIdx ].varName, szVar );
    pFmt->varCont[ nIdx ].varType = nType;
    pFmt->varCont[ nIdx ].dimRow  = pMat->GetNumRows();
    pFmt->varCont[ nIdx ].dimCol  = pMat->GetNumCols();
    pFmt->varCont[ nIdx ].nBytes  = pMat->GetBytes();
}

void SetVarInfo(RESFMT* pFmt, int nIdx, CFmVector* pVct, const char* szVar, unsigned int nType )
{
    _log_debug(_HI_, "SaveVarInfo: Var=%s, Vector=[%d]", szVar, pVct->GetLength());

    strcpy( pFmt->varCont[ nIdx ].varName, szVar );
    pFmt->varCont[ nIdx ].varType = nType;
    pFmt->varCont[ nIdx ].dimRow  = pVct->GetLength();
    pFmt->varCont[ nIdx ].dimCol  = 1;
    pFmt->varCont[ nIdx ].nBytes  = pVct->GetBytes();
}

void SetVarInfo(RESFMT* pFmt, int nIdx, CFmVectorStr* pVct, const char* szVar, unsigned int nType )
{
    _log_debug(_HI_, "SaveVarInfo: Var=%s, VectorStr=[%d, strlen=%d]", szVar, pVct->GetLength(), pVct->GetMaxStrlen() );

    strcpy( pFmt->varCont[ nIdx ].varName, szVar );
    pFmt->varCont[ nIdx ].varType = nType;
    pFmt->varCont[ nIdx ].dimRow  = pVct->GetLength();
    pFmt->varCont[ nIdx ].dimCol  = pVct->GetMaxStrlen()+1;
    pFmt->varCont[ nIdx ].nBytes  = pVct->GetBytes();
}

bool WriteStrlistToFile(CFmVectorStr* pList, FILE* fp)
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

int fwrite_hugedata(void* pData, size_t nDatLen, FILE* fp)
{
    size_t nBlock = 64*64*64; /* 262144 bytes*/
    size_t nSize = 0;
    while( nSize < nDatLen )
    {
        if ( nSize + nBlock > nDatLen)
            nBlock = nDatLen - nSize;

        size_t nRet = fwrite( (char*)pData+nSize, 1, nBlock, fp );
        if ( nRet!= nBlock )
        {
            _log_error( _HI_, "Failed to write a huge data, Error code:%d", ferror(fp) );
            return(ERR_WRITE_FILE);
        }

        nSize += nRet;
    }

    return(0);
}

int CFmRls::SaveResData( )
{
    _log_debug(_HI_, "SaveResData: Start to save the results to file(%s).", m_szRlsFile );

    RESFMT* pFmt = new RESFMT;
    memset(pFmt, 0, sizeof(RESFMT));
    strcpy(pFmt->szHeader, "fGWAS1.0,32\r\nURL: Http://statgen.psu.edu/\r\nCenter of Statistical Genetics, Penn State University.\n");

    int nIndex = 0;
    std::map<char*,CFmVectorStr*>::iterator vecStrIt;
    for ( vecStrIt=m_VecStrMap.begin() ; vecStrIt != m_VecStrMap.end(); vecStrIt++ )
    {
        char* szVar = (*vecStrIt).first;
        CFmVectorStr* pVec = (*vecStrIt).second;
        SetVarInfo( pFmt, nIndex, pVec, szVar, VAR_TYPE_STRING );
        nIndex++;
    }

    std::map<char*,CFmVector*>::iterator vecIt;
    for ( vecIt=m_VecMap.begin() ; vecIt != m_VecMap.end(); vecIt++ )
    {
        char* szVar = (*vecIt).first;
        CFmVector* pVec = (*vecIt).second;
        SetVarInfo( pFmt, nIndex, pVec, szVar, VAR_TYPE_DOUBLE);
        nIndex++;
    }

    std::map<char*,CFmMatrix*>::iterator matIt;
    for ( matIt=m_MatMap.begin() ; matIt != m_MatMap.end(); matIt++ )
    {
        char* szVar = (*matIt).first;
        CFmMatrix* pMat = (*matIt).second;
        SetVarInfo( pFmt, nIndex, pMat, szVar, VAR_TYPE_DOUBLE);
        nIndex++;
    }

    FILE* fp;
    if ((fp=fopen( m_szRlsFile, "wb"))==NULL)
    {
        delete pFmt;
        _log_error(_HI_, "Failed to open a result file to save result, file name=%s", m_szRlsFile);
        return( ERR_CREATE_FILE );
    }

    fwrite(pFmt, sizeof(RESFMT),1 , fp);

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
    delete pFmt;

    _log_debug(_HI_, "SaveResData: Successful to save the results to file(%s).", m_szRlsFile );

    return(0);
}

int LoadVarInfo(FILE* fp, CFmVectorStr* pVct, int nPos, VARFMT* pHeader)
{
    if (pHeader->varType != VAR_TYPE_STRING )
        return( ERR_FILE_DATA);

    _log_debug(_HI_, "LoadVarInfo: VString, Read var(%s) from the position %d ", pHeader->varName, nPos);

    pVct->Reset( pHeader->dimRow );

    fseek( fp, nPos, SEEK_SET);
    char* pTemp = new char[ pHeader->dimCol + 1];
    if (!pTemp)
    {
        _log_error(_HI_, "MEMORY: Filed to allocate %d bytes to a charaters buffer", pHeader->dimCol + 1 );
        return( ERR_MEM_ALLOC );
    }

    for (unsigned int i=0; i<pHeader->dimRow; i++)
    {
        if ( fread(pTemp, 1 , pHeader->dimCol, fp)== pHeader->dimCol )
        {
            pVct->Set(i, pTemp);
        }
        else
        {
            _log_error(_HI_, "LoadVarInfo: Filed to read a result file, nPos=%d, nBytes=%d.", nPos, pHeader->nBytes);
            delete[] pTemp;
            return( ERR_READ_FILE );
        }
    }

    delete[] pTemp;
    return(0);
}

int LoadVarInfo(FILE* fp, CFmVector* pVct, int nPos, VARFMT* pHeader)
{
    if (pHeader->varType != VAR_TYPE_DOUBLE )
        return( ERR_FILE_DATA);

    _log_debug(_HI_, "LoadVarInfo: Vector, Read var(%s) from the position %d ", pHeader->varName, nPos);

    pVct->Resize(pHeader->dimRow);

    fseek( fp, nPos, SEEK_SET);
    unsigned int nBytes = pHeader->nBytes;
    if ( fread(pVct->GetData(), 1 , nBytes, fp) != nBytes )
    {
        _log_error(_HI_, "LoadVarInfo: Filed to read a result file, nPos=%d, nBytes=%d.", nPos, nBytes);
        return( ERR_READ_FILE );
    }

    return(0);
}

int LoadVarInfo(FILE* fp, CFmMatrix* pMat, int nPos, VARFMT* pHeader)
{
    if (pHeader->varType != VAR_TYPE_DOUBLE )
        return( ERR_FILE_DATA);

    _log_info(_HI_, "LoadVarInfo: Matrix, Read var(%s) from the position %d ", pHeader->varName, nPos);

    pMat->Resize(pHeader->dimRow, pHeader->dimCol, true);

    fseek( fp, nPos, SEEK_SET);
    unsigned int nBytes = pHeader->nBytes;
    if ( fread(pMat->GetData(), 1 , nBytes, fp)!=nBytes)
    {
        _log_error(_HI_, "LoadVarInfo: Filed to read a result file, nPos=%d, nBytes=%d.", nPos, nBytes);
        return( ERR_READ_FILE );
    }

    return(0);
}

int CFmRls::LoadResFile()
{
    _log_info(_HI_, "LoadResFile: Open a result file , file name=%s", m_szRlsFile);

    FILE* fp;

    if ((fp=fopen(m_szRlsFile, "rb"))==NULL)
    {
        _log_error(_HI_, "LoadResFile: Filed to open a result file , file name=%s", m_szRlsFile);
        return( ERR_OPEN_FILE );
    }

    RESFMT fmt;
    if ( fread(&fmt, 1, sizeof(fmt), fp)!= sizeof(fmt))
    {
        _log_error(_HI_, "LoadResFile: Filed to read a result file , file name=%s", m_szRlsFile);
        fclose(fp);
        return( ERR_READ_FILE );
    }

	CFmNewTemp refNew;

    int nRet = 0;
    int nPos = sizeof(fmt);
    for(int i=0;i<MAX_CONT_ITEM;i++)
    {
        nRet = 0;
        char* szOrgVar = fmt.varCont[i].varName;

        if (strlen(szOrgVar)==0)
            continue;

        char* szVar = Strdup(szOrgVar);
        if ( fmt.varCont[i].dimCol == 1 )
        {
            CFmVector* pVec = new (refNew) CFmVector( fmt.varCont[i].dimRow, 0 );
            nRet += LoadVarInfo(fp, pVec, nPos,  &(fmt.varCont[i]) );
            m_VecMap.insert( std::pair<char*,CFmVector*>( szVar, pVec ) );
        }
        else
        if ( fmt.varCont[i].varType == VAR_TYPE_STRING )
        {
            CFmVectorStr* pVecStr = new  (refNew) CFmVectorStr( fmt.varCont[i].dimRow );
            nRet += LoadVarInfo(fp, pVecStr, nPos, &(fmt.varCont[i]) );
            m_VecStrMap.insert( std::pair<char*,CFmVectorStr*>( szVar, pVecStr ) );
        }
        else if (fmt.varCont[i].dimRow>0 && fmt.varCont[i].dimCol>1)
        {
            CFmMatrix* pMat = new  (refNew) CFmMatrix( (int)(fmt.varCont[i].dimRow), (int)(fmt.varCont[i].dimCol) );
            nRet += LoadVarInfo(fp, pMat, nPos, &(fmt.varCont[i]) );
            m_MatMap.insert( std::pair<char*,CFmMatrix*>( szVar, pMat ) );
        }

		Free(szVar);

        if(nRet!=0)
            break;

        nPos += fmt.varCont[i].nBytes;
    }

    fclose(fp);

    _log_info(_HI_, "LoadResFile: Open a result file %s", nRet==0?"successful":"unsuccessful");

    return(nRet);
}


int CFmRls::SaveRData( const char* szRdataFile )
{
    FILE* fp;
    if ((fp=fopen(szRdataFile, "wb"))==NULL)
    {
        _log_error(_HI_, "Filed to create a RDATA to save result, file name=%s", szRdataFile);
        return( ERR_CREATE_FILE );
    }

    int nSize = 0;
    nSize += m_MatMap.size();
    nSize += m_VecMap.size();
    nSize += m_VecStrMap.size();

    _log_info(_HI_, "SaveRData: Start to save RDATA to file(%s).count=%d", szRdataFile, nSize);

    SEXP list, names;
    PROTECT(list=allocVector(VECSXP, nSize));
    PROTECT(names = allocVector(STRSXP, nSize));

    int nIndex=0;
    std::map<char*,CFmVectorStr*>::iterator vecStrIt;
    for ( vecStrIt=m_VecStrMap.begin() ; vecStrIt != m_VecStrMap.end(); vecStrIt++ )
    {
        SEXP exp_vec;
        CFmVectorStr* pVecStr = (*vecStrIt).second;
        PROTECT(exp_vec = GetSEXP( pVecStr ));
        SET_VECTOR_ELT(list, nIndex, exp_vec);
        SET_STRING_ELT(names, nIndex, mkChar( (*vecStrIt).first ));
        nIndex++;
    }

    std::map<char*,CFmVector*>::iterator vecIt;
    for ( vecIt=m_VecMap.begin() ; vecIt != m_VecMap.end(); vecIt++ )
    {
        SEXP exp_vec;
        CFmVector* pVec = (*vecIt).second;
        PROTECT(exp_vec = GetSEXP( pVec ));
        SET_VECTOR_ELT(list,  nIndex, exp_vec);
        SET_STRING_ELT(names, nIndex, mkChar( (*vecIt).first ));
        nIndex++;
    }

    std::map<char*,CFmMatrix*>::iterator matIt;
    for ( matIt=m_MatMap.begin() ; matIt != m_MatMap.end(); matIt++ )
    {
        SEXP exp_mat;
        CFmMatrix* pMat = (*matIt).second;
        PROTECT(exp_mat = GetSEXP(pMat));
        SET_VECTOR_ELT(list,  nIndex, exp_mat);
        SET_STRING_ELT(names, nIndex, mkChar( (*matIt).first ));
        nIndex++;
    }

    setAttrib(list, R_NamesSymbol, names);
    R_SaveToFileV( list, fp, true);
    fclose(fp);

    UNPROTECT( nSize + 2);

    _log_info(_HI_, "SaveRData: Successfule to save file(%s).", szRdataFile);

    return(0);
}

void destroy(CFmRls* p)
{
	CFmNewTemp  fmRef;
	p->~CFmRls();
	operator delete(p, fmRef);
}


