/* fm_mapmat.cpp  -	CFmMapMat Class
 *
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>

#include "fm_err.h"
#include "fm_rlogger.h"
#include "fm_vector_str.h"
#include "fm_vector.h"
#include "fm_matrix.h"
#include "fm_mapmat.h"
#include "fm_dataframe.h"
#include "fm_linux.h"
#include "fm_new.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
//#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction CFmSnpMat
//////////////////////////////////////////////////////////////////////

#ifdef _DEBUG
int CFmSnpMat::m_NextObjectNumber = 1 ;
#endif

CFmMapMat::CFmMapMat()
{
    m_pMap = NULL;
    m_pMarkers = NULL;
    m_nChr = 0;
}


CFmMapMat::~CFmMapMat()
{
    if(m_pMap) destroy( m_pMap );
    if(m_pMarkers) destroy( m_pMarkers );
}

int CFmMapMat::GetNumChrs()
{
    return(m_nChr);
}

int CFmMapMat::GetStartMarker(int nChr)
{
    for (int i=0; i<m_pMap->GetNumRows(); i++)
    {
        if ( m_pMap->Get(i,0) == nChr)
            return(i);
    }

    return(-1);

}

int CFmMapMat::GetNumMarkers(int nChr)
{
    int nCnt=0;
    for (int i=0; i<m_pMap->GetNumRows(); i++)
    {
        if ( m_pMap->Get(i,0) == nChr)
            nCnt++;
    }

    return(nCnt);
}

double CFmMapMat::GetChrLength(int nChr)
{
    double fLen=0;
    for (int i=0; i<m_pMap->GetNumRows(); i++)
    {
        if ( m_pMap->Get(i,0) == nChr && m_pMap->Get(i,1) > fLen)
            fLen = m_pMap->Get(i,1);
    }

    return(fLen);
}

double CFmMapMat::GetMapLength()
{
    return(0);
}

int CFmMapMat::GetNumMarkers()
{
    return(0);
}

double CFmMapMat::GetPosition(int nChr, int nIndex)
{
    int nCnt=0;
    for (int i=0; i<m_pMap->GetNumRows()-1; i++)
    {
        if (nCnt >= nIndex )
            return( m_pMap->Get(i,1) );

        if ( m_pMap->Get(i,0) == nChr)
           nCnt++;
    }

    return(-1);
}

double CFmMapMat::GetDistance(int nChr, int nLeft, int nRight)
{
    int nCnt=0;
    double fLeft=-1;
    double fRight=-1;
    for (int i=0; i<m_pMap->GetNumRows()-1; i++)
    {
        if (nCnt >= nLeft && fLeft<0)
            fLeft = m_pMap->Get(i,1);

        if ( m_pMap->Get(i,0) == nChr)
           nCnt++;

        if (nCnt >= nRight && fRight<0)
            fRight = m_pMap->Get(i+1,1);

        if (fRight>=0 && fLeft>=0)
            return( fRight-fLeft);
    }

    return(-1);
}

char* CFmMapMat::GetMarkerName(int nIdx)
{
    return m_pMarkers->Get(nIdx);
}

char* CFmMapMat::GetMarkerName(int nChr, int nMarker)
{
    int nCnt=0;
    for (int i=0; i<m_pMap->GetNumRows(); i++)
    {
        if ( m_pMap->Get(i,0) == nChr && nCnt == nMarker)
            return(m_pMarkers->Get(i));
        if ( m_pMap->Get(i,0) == nChr)
            nCnt++;
    }

    return(NULL);
}

int CFmMapMat::LoadCSV(const char* filename)
{
    _log_info( _HI_, "Phenotype: Phenotype file(%s) is being loaded", filename);

    CFmDataFrame df;
    int ret = df.Load(filename, true, true);
    if (ret)
        return ret;

	CFmNewTemp refNew;

    int nIdx = df.FindColumn("Marker");
    if (nIdx>=0)
        m_pMarkers = new (refNew) CFmVectorStr( df.GetStringCol(nIdx));

    int nIdx1 = df.FindColumn("grp_idx");
    int nIdx2 = df.FindColumn("Distance");
    if (nIdx1>=0 && nIdx2>=0)
    {
        m_pMap = new (refNew) CFmMatrix(0,0);
        CFmVector& Vct1 = df.GetFloatCol(nIdx1);
        CFmVector& Vct2 = df.GetFloatCol(nIdx2);
        m_nChr = (int)Vct1.GetMax();
        m_pMap->Cbind( Vct1 );
        m_pMap->Cbind( Vct2 );
    }

    _log_info( _HI_, "CFmMapMat::Load %d %d %d", nIdx, nIdx1, nIdx2);

    return(0);
}

int CFmMapMat::WriteCSV(const char* filename, bool bAppend, const char* szTag)
{
    FILE* fp = NULL;
    if (!bAppend)
        fp = fopen( filename, "wt");
    else
        fp = fopen( filename, "a+");

    if (fp==NULL)
    {
        _log_error(_HI_, "The file can not be created");
        return(ERR_CREATE_FILE);
    }

    if (szTag!=NULL)
    {
        fprintf(fp, "\n>>-----------------%s--------------------<<\n", szTag);
    }

    fprintf(fp, "marker, grp_idx, Distance\n" );
    for (int i = 0; i < m_pMarkers->GetLength() ; i++)
    {
        fprintf(fp, "%s", m_pMarkers->Get(i));
        fprintf(fp, ",%.4f", m_pMap->Get(i,0) );
        fprintf(fp, ",%.4f", m_pMap->Get(i,1)  );
        fprintf(fp, "\n");
    }

    fclose(fp) ;

    return(0);
}



CFmDatPheno::CFmDatPheno(char* szFile_pheno )
{
    m_pPhenoY = NULL;
    m_pCovarX = NULL;
    m_pTimes  = NULL;

    m_nSubjN = 0;
    m_nMesuQ = 0;

    m_szFile_pheno = Strdup(szFile_pheno);
}

CFmDatPheno::~CFmDatPheno()
{
	if (m_szFile_pheno) Free(m_szFile_pheno);
    if (m_pTimes) destroy( m_pTimes );
}

int CFmDatPheno::Load()
{
    _log_info( _HI_, "CFmDatPheno::Load Phenotype file(%s) is being loaded", m_szFile_pheno);

    CFmDataFrame df;
    int ret = df.Load(m_szFile_pheno, true, true);
    if (ret)
        return ret;

	CFmNewTemp refNew;

    m_pPhenoY = new (refNew) CFmMatrix(0,0);
    m_pTimes = new (refNew) CFmVector(0, 0.0);
    for (int i=0; i<df.GetNumCol(); i++ )
    {
        if ( strcmp(df.GetColName(i), "id")==0 ||
             strcmp(df.GetColName(i), "ID")==0 ||
             strcmp(df.GetColName(i), "Id")==0 ||
             strcmp(df.GetColName(i), "X")==0 ||
             strcmp(df.GetColName(i), "x")==0 )
            continue;

        CFmVector& Vct = df.GetFloatCol(i);
        m_pPhenoY->Cbind( Vct );
        m_pTimes->Put(atof(df.GetColName(i)));
    }

    for (int i=0; i<df.GetNumCol(); i++ )
    {
        if ( strcmp(df.GetColName(i), "X")==0 ||
             strcmp(df.GetColName(i), "x")==0 )
        {
            CFmVector& Vct = df.GetFloatCol(i);
            m_pCovarX = new (refNew) CFmVector( &Vct );
        }
    }

    *m_pTimes = (*m_pTimes)*(m_pTimes->GetLength());

    m_nMesuQ = m_pPhenoY->GetNumCols();
    m_nSubjN = m_pPhenoY->GetNumRows();

    _log_info( _HI_, "CFmDatPheno::Load m_nMesuQ=%d, m_nSubjN=%d", m_nMesuQ, m_nSubjN);

    return(0);
}

CFmDatGeno::CFmDatGeno(char* szFile_snp )
{
    m_pSnpMat = 0;
    m_pSnpNames = 0;
    m_szFile_snp = Strdup(szFile_snp);
}

CFmDatGeno::~CFmDatGeno()
{
	if(m_szFile_snp) Free(m_szFile_snp);
}

int CFmDatGeno::Load()
{
    _log_info(_HI_, "m_nMesuQ: Start to load SNP data into a matrix from %s", m_szFile_snp);

	CFmNewTemp refNew;
    m_pSnpMat = new (refNew) CFmMatrix(0, 0);
    int ret = m_pSnpMat->ReadFromCSVFile( m_szFile_snp, true, true);
    if (ret!=0)
    {
        _log_info(_HI_, "m_nMesuQ: Failed to load the SNP data into a matrix from %s",m_szFile_snp);
        return(ret);
    }

    _log_info(_HI_, "m_nMesuQ: The SNP data is loaded into the matrix[%d,%d]",
                    m_pSnpMat->GetNumRows(), m_pSnpMat->GetNumCols());

    return(0);
}

void destroy(CFmMapMat* p)
{
	CFmNewTemp  fmRef;
	p->~CFmMapMat();
	operator delete(p, fmRef);
}

void destroy(CFmDatPheno* p)
{
	CFmNewTemp  fmRef;
	p->~CFmDatPheno();
	operator delete(p, fmRef);
}

void destroy(CFmDatGeno* p)
{
	CFmNewTemp  fmRef;
	p->~CFmDatGeno();
	operator delete(p, fmRef);
}
