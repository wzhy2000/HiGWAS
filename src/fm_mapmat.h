// fm_mapmat.h: interface for the CFmMapMat class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_MAPMAT_H_)
#define _FM_MAPMAT_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//
// Class : CFmSnpMat
//

#define UNUSED_PARAMETER(x) x

#include <stdio.h>
#include <string.h>

class CFmVectorStr;
class CFmVector;
class CFmMatrix;

class CFmMapMat
{
public:
    CFmMapMat() ;
    virtual	~CFmMapMat();

    int  GetNumChrs();
    int  GetNumMarkers();
    int  GetNumMarkers(int nChr);
    int  GetStartMarker(int nChr);

    double GetMapLength();
    double GetChrLength(int nChr);
    double GetDistance(int nChr, int nLeft, int nRight);
    double GetPosition(int nChr, int nIndex);

    char* GetMarkerName(int nChr, int nMarker);
    char* GetMarkerName(int nTotalIdx);

    int  LoadCSV(const char* filename);
    int  WriteCSV(const char* filename, bool bAppend=false, const char* szTag=NULL);

private:
	// private internal functions
    char*			AllocateChar(int nRows, int nCols) ;
    void			FreeMemory();

    CFmMatrix * m_pMap;
    CFmVectorStr* m_pMarkers;
    int m_nChr;
};

class CFmDatPheno
{
public:
    CFmDatPheno(char* szFile_pheno );
    virtual ~CFmDatPheno();

    int Load();

    CFmMatrix* m_pPhenoY;     //Matrix[N,Q];
    CFmVector* m_pCovarX;     //VECTOR[N];
    CFmVector* m_pTimes;

    int m_nSubjN;
    int m_nMesuQ;

private:
    char* m_szFile_pheno;
};

class CFmDatGeno
{
public:
    CFmDatGeno(char* szFile_snp );
    virtual ~CFmDatGeno();

    int Load();

    CFmMatrix* m_pSnpMat;
    CFmVectorStr* m_pSnpNames;

private:
    char* m_szFile_snp;
};


void destroy(CFmMapMat* p);
void destroy(CFmDatPheno* p);
void destroy(CFmDatGeno* p);

#endif // !defined(_FM_MAPMAT_H_)
