/* CFmDat.cpp  -	GLS Data Object
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <stdio.h>
#include <math.h>

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>
#include <Rmath.h>

#include "fm_linux.h"
#include "fm_rlogger.h"
#include "fm_err.h"
#include "fm_dataframe.h"
#include "fm_dat.h"
#include "fm_new.h"

CFmDat_Plink::CFmDat_Plink(char* szFile_tped, char* szFile_tfam )
{
    m_sTped_file = Strdup(szFile_tped);
    m_sTfam_file = Strdup(szFile_tfam);
    m_nSnpP = 0;
	m_pSubIds = NULL;
	m_pSnpNames = NULL;
}

CFmDat_Plink::~CFmDat_Plink()
{
    Free(m_sTped_file);
    Free(m_sTfam_file);

	if(m_pSnpNames) destroy( m_pSnpNames );
	if(m_pSubIds) destroy( m_pSubIds );

    _log_debug(_HI_, "CFmDat_Plink is released successfully.");
}

int CFmDat_Plink::Load(char* szPresigFile)
{
    _log_info( _HI_, "LoadPlink: Start to load TPed/Tfam,tped=%s", m_sTped_file);
    int nRet = LoadTped( m_sTped_file, szPresigFile );
    if (nRet)
        return nRet;

    nRet = LoadTfam( m_sTfam_file );
    if (nRet)
        return nRet;

    m_nSnpP = m_PackedSNP.GetNumSnps();
    _log_info( _HI_, "LoadPlink: Start to remove rare SNPs,snp=%d", m_nSnpP);
    nRet = m_PackedSNP.RemoveRareSNP( 0.01 );
    if (nRet!=0)
    {
        _log_info(_HI_, "LoadSnpFile: Failed to remove rare SNPs");
        return(nRet);
    }

    m_nSnpP = m_PackedSNP.GetNumSnps();
    _log_info( _HI_, "LoadPlink: non-rare SNP:%d", m_nSnpP );

	CFmNewTemp refNew;
    m_pSnpNames = new (refNew) CFmVectorStr(m_nSnpP);
    for (int i=0;i<m_nSnpP;i++)
    {
        char szBuf[64]={0};
        SNPINFO* pInfo;
        if (m_PackedSNP.GetSnpInfo(i, &pInfo))
        {
            snprintf(szBuf, 63, "%s", pInfo->szSnpId);
            m_pSnpNames->Set(i, szBuf);
        }
    }

    _log_info( _HI_, "LoadPlink: successful");
    return nRet;
}

/* TPED Column:
1 chr
2 snp id
3 genetic distance
4 base-pair position
5/6 allel1
7/8 allel2 */

bool CFmDat_Plink::Extract_snpinfo( char* szLine, CFmVectorStr* pSigSnp )
{
    //char gen_type[] = {'C','A','T','G'};
    char szFullId [256]={0};
    char szChrId [256]={0};
    char szSnpId [256]={0};
    char szGenDist[256]={0};
    char szBasePair[256]={0};
    char snpPairs[MAX_SUBJ_INDATA*2+1]={0};

    char seps1[] = " ,\t";
    char* token = strtok( szLine, seps1 );
    int nCol = 0;
    while( token != NULL && nCol< MAX_SUBJ_INDATA*2 + 4 )
    {
       /* While there are tokens in "string" */
       if (nCol<4)
       {
            switch(nCol)
            {
                case 0: strcpy( szChrId, token);
                    break;
                case 1: strcpy( szSnpId, token);
                    break;
                case 2: strcpy( szGenDist, token);
                    break;
                case 3: strcpy( szBasePair, token);
                    break;
            }
       }
       else
           snpPairs[nCol-4] = token[0];

       /* Get next token: */
       token = strtok( NULL, seps1 );
       nCol++;
    }

    if ( pSigSnp!=NULL && pSigSnp->GetLength()>0 )
    {
        sprintf(szFullId, "%s/%s", szChrId, szSnpId);
        /* the SNP is not selected by the significant subset*/
        if( pSigSnp->Find(szFullId)<0)
            return(true);
    }

    int nSubj = strlen(snpPairs)/2;
    int nSnpIdx = -1;
    if ( m_PackedSNP.GetNumSnps() == 0 )
    {
        m_PackedSNP.InitData(nSubj, 1);
        nSnpIdx =0;
    }
    else
    {
        if ( m_PackedSNP.GetNumSubjs() != nSubj)
        {
            _log_error(_HI_ , "The TPED file does not have same subjects with SNP data.");
            return(false);
        }
        else
        {
            nSnpIdx = m_PackedSNP.AppendSnp();
            if (nSnpIdx<0)
            {
                _log_error(_HI_ , "Memory is not enough to append SNP(%d).", nSnpIdx);
                return(false);
            }
        }
    }

    char BigA ='0';
    for (int j=0; j<nSubj*2; j++)
    {
        BigA = snpPairs[j];
        if (BigA!='0') break;
    }

    //Debug
    //Rprintf("%d, %d, %s, %s, %s, %s, %c\n", nSnpIdx, nSubj, szChrId, szSnpId, szGenDist, szBasePair, BigA);

    for(int i=0; i<nSubj; i++)
    {
        int snp_d = 1;
        if (snpPairs[2*i]=='0' && snpPairs[2*i+1]=='0')
            snp_d = 3;
        else if (snpPairs[2*i]==BigA && snpPairs[2*i+1]==BigA)
            snp_d = 0;
        else if (snpPairs[2*i]!=BigA && snpPairs[2*i+1]!=BigA)
            snp_d = 2;
        else
            snp_d = 1;

        m_PackedSNP.Set( i, nSnpIdx, snp_d);
        if ( m_PackedSNP.Get( i, nSnpIdx ) != snp_d )
        {
            _log_fatal(_HI_, "Extract_snpinfo: Failed to check snp %d(%d != %d)",  i, snp_d, m_PackedSNP.Get( i, nSnpIdx ) );
        }
    }

    m_PackedSNP.SetSnpInfo(nSnpIdx, (char*)szSnpId, (char*)szChrId, (char*)szGenDist, (char*)szBasePair, BigA );
    return(true);
}

int CFmDat_Plink::LoadTped( char* szFile_tped, char* szPresig )
{
    _log_info( _HI_, "TPED file is being loaded. File=%s PreSig=%s", szFile_tped, szPresig);

    CFmVectorStr fmPreSig(0);
    if (szPresig!=NULL && strlen(szPresig)>0)
    {
        _log_info( _HI_, "Load preselect significant SNP. File=%s", szPresig);
        CFmDataFrame df;
        df.Load( szPresig, TRUE, TRUE );
        fmPreSig.Append( df.RGetRowNames() );

        _log_info( _HI_, "%d SNP in PreSig File", fmPreSig.GetLength());
    }

    FILE* fp = fopen( szFile_tped, "rt");
    if ( fp ==NULL )
    {
        _log_error(_HI_, "Failed to open the TPED file(%s)", szFile_tped);

        return( ERR_OPEN_FILE );
    }

    char aLine[MAX_SUBJ_INDATA*32]={0};
    while( !feof( fp ) )
    {
        if(  fgets( aLine, MAX_SUBJ_INDATA*32, fp ) != NULL)
            if ( !Extract_snpinfo( aLine, &fmPreSig ) )
                return( ERR_FILE_DATA );
    }

    fclose(fp);

    _log_info( _HI_, "TPED file is loaded successful. Subjects:%d, SNPs:%d",
                    m_PackedSNP.GetNumSubjs(), m_PackedSNP.GetNumSnps() );


    return(0);
}

/* TFAM Column:
1 Family ID
2 Individual ID
3 Paternal ID
4 Maternal ID
5 Sex (1=male; 2=female..)
6 Phenotype */

int CFmDat_Plink::LoadTfam( char* szFile_tfam )
{
    _log_info( _HI_, "TFAM file is being loaded. File=%s", szFile_tfam);

    FILE* fp = fopen( szFile_tfam, "rt");
    if (fp == NULL)
    {
        _log_error(_HI_, "Failed to open the TFAM file(%s).", szFile_tfam);
        return( ERR_OPEN_FILE );
    }

	CFmNewTemp refNew;
	m_pSubIds = new (refNew) CFmVectorStr(0, 1000);
    int nSubj=0;
    char aLine[256]={0};
    char seps1[] = " ,\t";
    while( !feof( fp ) )
    {
        if( fgets( aLine, 256, fp ) != NULL)
        {
            char* token = strtok( aLine, seps1 );
            int i=1;
            while( token != NULL && i<=6)
            {
               //just get subId
               if (i==2)
               {
                    m_pSubIds->Put( token );
               }

               token = strtok( NULL, seps1 );
               i++;
           }

			nSubj++;
        }
    }

    _log_info( _HI_, "TFAM file is loaded into memory.(Subj:%d)", nSubj);
    fclose(fp);

    return(0);
}


CFmDat_Simple::CFmDat_Simple(char* szFile_snp )
{
    m_szFile_snp = Strdup(szFile_snp);
	m_pSubIds = NULL;
    m_pSnpNames= NULL;

}

CFmDat_Simple::~CFmDat_Simple()
{
    Free(m_szFile_snp);
	if(m_pSubIds) destroy( m_pSubIds );
	if(m_pSnpNames) destroy( m_pSnpNames );

    _log_debug(_HI_, "CFmDat_Simple is released successfully.");
}

int CFmDat_Simple::Load(char* szPresigFile)
{
    _log_info(_HI_, "LoadSnpFile: Start to load SNP data into a matrix from %s", m_szFile_snp);

    CFmMatrix mat(0,0);
    int ret = mat.ReadFromCSVFile( m_szFile_snp, true, true);
    if (ret!=0)
    {
        _log_error(_HI_, "LoadSnpFile: Failed to load the SNP data into a matrix from %s",m_szFile_snp);
        return(ret);
    }

    _log_info(_HI_, "LoadSnpFile: The SNP data is loaded into the matrix[%d,%d]",
                mat.GetNumRows(), mat.GetNumCols()-2);
    if ( m_PackedSNP.GetNumSnps() == 0 )
        m_PackedSNP.InitData(mat.GetNumCols()-2, mat.GetNumRows());

    for(int iSnp=0; iSnp<mat.GetNumRows(); iSnp++)
    {
        for(int jSub=0; jSub<mat.GetNumCols()-2; jSub++)
        {
            int snp_d = (int)(mat.Get( iSnp, jSub+2));
            if (snp_d != 0 && snp_d != 1 && snp_d!= 2)
                snp_d=3;

            m_PackedSNP.Set( jSub, iSnp, snp_d);
            if ( m_PackedSNP.Get( jSub, iSnp ) != snp_d )
            {
                _log_fatal(_HI_, "LoadSnpFile: Failed to check snp %d(%d != %d)", iSnp, snp_d, m_PackedSNP.Get(jSub, iSnp ) );
            }
        }

        char szChr[64] = {"0"};
        char szPos[64] = {"0"};
        sprintf(szChr, "%d", (int)mat.Get( iSnp, 0) );
        sprintf(szPos, "%d", (int)mat.Get( iSnp, 1) );

        m_PackedSNP.SetSnpInfo(iSnp, mat.GetRowName(iSnp), (char*)szChr, (char*)szPos, (char*)szPos, 'C' );
    }

    _log_info(_HI_, "LoadSnpFile: SNP data is load into Data object successfully. SNP=%d, subj=%d",
                m_PackedSNP.GetNumSnps(),m_PackedSNP.GetNumSubjs() );


    int nSnpP = m_PackedSNP.GetNumSnps();
    _log_info( _HI_, "LoadSimple: Start to remove rare SNPs,snp=%d", nSnpP);
    ret = m_PackedSNP.RemoveRareSNP( 0.01 );
    if (ret!=0)
    {
        _log_info(_HI_, "LoadSnpFile: Failed to remove rare SNPs");
        return(ret);
    }

    nSnpP = m_PackedSNP.GetNumSnps();

    _log_info( _HI_, "LoadSimple: non-rare SNP:%d", nSnpP );

	CFmNewTemp refNew;
    m_pSnpNames = new (refNew) CFmVectorStr(nSnpP);
    for (int i=0;i<nSnpP;i++)
    {
        SNPINFO* pInfo;
        if (m_PackedSNP.GetSnpInfo(i, &pInfo))
        {
            m_pSnpNames->Set(i, pInfo->szSnpId);
        }
    }

    return(0);
}

CFmDat_Pheno::CFmDat_Pheno(char* szFile_pheno, bool bZnorm, char* szYname, char* szZname, char* szXname )
{
    m_szFile_pheno = Strdup(szFile_pheno);
    m_bZnorm = bZnorm;
    m_pAttachedPhe = NULL;

	m_pszXname = NULL;
	m_pszZname = NULL;
	m_pszYname = Strdup(szYname);

	if(szXname) m_pszXname = Strdup(szXname);
	if(szZname) m_pszZname = Strdup(szZname);

	Init();
}

void CFmDat_Pheno::Init()
{
    m_pPhenoY = NULL ;
    m_pPhenoZ = NULL ;
    m_pPhenoZ0 = NULL ;
    m_pCovars = NULL ;
	m_pSubjs = NULL;

    m_nSubjN = 0;
    m_nMesuQ = 1;

	CFmNewTemp refNew;
	m_pXnames = new (refNew) CFmVectorStr(0);
    char* pch = strtok (m_pszXname,",");
    while (pch != NULL)
    {
        m_pXnames->Put(pch);
        pch = strtok (NULL, ",");
    }

    //m_pXnames->Show("m_pXnames");
}

CFmDat_Pheno::CFmDat_Pheno(CFmMatrix* pFmPhe, bool bZnorm, char* szYname, char* szZname, char* szXname )
{
    m_szFile_pheno = NULL;
    m_pAttachedPhe = pFmPhe;
    m_bZnorm = bZnorm;

	m_pszXname = NULL;
	m_pszZname = NULL;
	m_pszYname = Strdup(szYname);

	if(szXname) m_pszXname = Strdup(szXname);
	if(szZname) m_pszZname = Strdup(szZname);

	Init();
}

CFmDat_Pheno::~CFmDat_Pheno()
{
    if (m_szFile_pheno) Free( m_szFile_pheno );
    if (m_pszXname) Free( m_pszXname);
    if (m_pszYname) Free( m_pszYname);
    if (m_pszZname) Free( m_pszZname);

    if (m_pXnames) destroy( m_pXnames );
    if (m_pSubjs) destroy( m_pSubjs );
    if (m_pPhenoY) destroy( m_pPhenoY );
    if (m_pPhenoZ) destroy( m_pPhenoZ );
    if (m_pPhenoZ0) destroy( m_pPhenoZ0 );
    if (m_pCovars) destroy( m_pCovars );

    _log_debug(_HI_, "CFmDat_Pheno is released successfully.");
}

int CFmDat_Pheno::LoadLongdt( CFmPackedSNP* pPackedSNP , CFmVectorStr* pFamSubjs)
{
    //include response value(y);
    if ( m_pszYname ==NULL || m_pszZname == NULL )
    {
        _log_error(_HI_, "No response value or  covariate in the model parameter");
        return( ERR_PARAM_VALUE );
    }

	CFmMatrix* pFmPhe = NULL;
	if( m_szFile_pheno != NULL)
	{
	    CFmDataFrame df;
		int ret = df.Load(m_szFile_pheno, false, true);
		if (ret!=0)
		{
			_log_error(_HI_, "Failed to open the Phenotype file(%s)", m_szFile_pheno);
			return( ERR_OPEN_FILE );
		}

		CFmVector vctCol(0, 0.0);
		for(int i=1;i<df.GetNumCol();i++)	vctCol.Put(i);
		pFmPhe = df.GetMatrix( &vctCol );
		pFmPhe->SetRowNames( df.GetStringCol(0) );
	}
	else
		pFmPhe = m_pAttachedPhe;

    m_nSubjN = pPackedSNP->GetNumSubjs();
    if ( pFmPhe->GetNumRows() != m_nSubjN )
    {
        _log_error(_HI_, "The id's count is greater than the count in SNP data, %d!=%d", pFmPhe->GetNumRows(), m_nSubjN);
        return( ERR_UNKNOWN );
    }

    CFmNewTemp refNew;
    m_pSubjs = new (refNew) CFmVectorStr( pFmPhe->GetRowNames() );

    // Data from PLINK command
	CFmVector subjOrd(0, 0.0);
    if (pFamSubjs)
    {
        for(int i=0;i<m_pSubjs->GetLength(); i++)
        {
            int nPos = pFamSubjs->Find ( m_pSubjs->Get(i) );
            if ( nPos>=0  )
                subjOrd.Put( nPos+1 );
            else
                _log_error(_HI_, "Failed to find %s in phenotype", m_pSubjs->Get(i));
        }
    }
    // Data from SIMPLE command
    else
    {
        for(int i=0;i<m_pSubjs->GetLength(); i++)
            subjOrd.Put( i+1 );
    }

    char* szVarY = m_pszYname;
    CFmVector fmVctPosY(0, 0.0);
    for(int i=1; i<= pFmPhe->GetNumCols(); i++)
    {
        char szLongY[256] = {0};
        sprintf(szLongY, "%s_%d", szVarY, i);
        int nIdx = pFmPhe->FindColumn( szLongY );
        if (nIdx<0)
            continue;

        fmVctPosY.Put(i);
    }

    m_nMesuQ = fmVctPosY.GetLength();
    if (m_nMesuQ<1)
    {
        _log_error(_HI_, "Failed to find the column(%s)", szVarY);
        return( ERR_PARAM_VALUE );
    }

    // Y matrix
    if (m_pPhenoY) destroy( m_pPhenoY );

    m_pPhenoY = new (refNew) CFmMatrix( m_nSubjN, m_nMesuQ );

    for(int i=0; i<fmVctPosY.GetLength(); i++)
    {
        char szLongY[256] = {0};
        sprintf(szLongY, "%s_%d", szVarY, (int)(fmVctPosY.Get(i)) );
        int nIdx = pFmPhe->FindColumn( szLongY );
        if (nIdx<0)
            _log_fatal( _HI_, "Why?(DataFrame.nIdx=%d)", nIdx );

        CFmVector& VctY = pFmPhe->GetCol(nIdx);
        VctY.Rearrange( subjOrd );

        m_pPhenoY->SetCol(i, VctY );
    }


    _log_prompt( _HI_, "Hypothesis Z=%s", m_pszZname);
    if (m_pszZname == NULL)
        return( ERR_PARAM_VALUE );

    if (m_pPhenoZ) destroy( m_pPhenoZ );
    if (m_pPhenoZ0) destroy( m_pPhenoZ0 );

    for(int i=0; i<fmVctPosY.GetLength(); i++)
    {
        char szLongZ[256] = {0};
        sprintf(szLongZ, "%s_%d", m_pszZname, (int)(fmVctPosY.Get(i)) );
        int nIdx = pFmPhe->FindColumn( szLongZ );
        if (nIdx<0)
        {
            _log_error( _HI_, "Hypothesis Z doesn't match Y(%d!=%d).", fmVctPosY.GetLength(), m_nMesuQ );
            return( ERR_PARAM_VALUE );
        }
    }

    m_pPhenoZ = new (refNew) CFmMatrix( m_nSubjN, m_nMesuQ );
    for( int i=0; i<fmVctPosY.GetLength(); i++)
    {
        char szLongZ[256] = {0};
        sprintf(szLongZ, "%s_%d", m_pszZname, (int)(fmVctPosY.Get(i)) );
        int nIdx = pFmPhe->FindColumn( szLongZ );
        if (nIdx<0)
            _log_fatal( _HI_, "Why?(DataFrame.nIdx=%d)", nIdx );

        CFmVector& VctZ = pFmPhe->GetCol(nIdx);
        VctZ.Rearrange( subjOrd );
        m_pPhenoZ->SetCol(i, VctZ );
    }

    CFmVector phe_null(0, 0.0);
    CFmVector vct(0, 0.0);
    CFmVector vct_ord(0, 0.0);

    for(int i=0; i<m_pPhenoZ->GetNumRows(); i++)
    {
        vct_ord.Resize(0, true);
        for(int j=0; j<m_pPhenoZ->GetNumCols(); j++)
            if ( !isnan( m_pPhenoZ->Get(i,j) ) && !isnan( m_pPhenoY->Get(i,j) ) )
            	vct_ord.Put(j);

         vct.Resize(0, true);
         for( int j=0; j<m_pPhenoZ->GetNumCols(); j++)
         	vct.Put(R_NaN);

         for( int j=0; j<vct_ord.GetLength(); j++)
            vct.Set( j, m_pPhenoZ->Get(i, (int)(vct_ord[j])) );
         m_pPhenoZ->SetRow(i, vct);

         vct.Resize(0, true);
         for( int j=0; j<m_pPhenoY->GetNumCols(); j++)
         	vct.Put(R_NaN);

         for( int j=0; j<vct_ord.GetLength(); j++)
            vct.Set( j, m_pPhenoY->Get(i, (int)(vct_ord[j])) );
         m_pPhenoY->SetRow(i, vct);
    }

    // Covars matrix
    m_pCovars = new (refNew) CFmMatrix( m_nSubjN, m_pXnames->GetLength() + 1);
	m_pCovars->SetColName( 0, (char*)"Mu" );
    for(int i=0; i<m_nSubjN; i++)
        m_pCovars->Set( i, 0, 1.0 );

    for(int i=1; i<m_pXnames->GetLength()+1;i++)
    {
		char* szVarX = m_pXnames->Get( i-1 );
        int nIdx = pFmPhe->FindColumn( szVarX );
        if (nIdx<0)
        {
            _log_error(_HI_, "Failed to find the column(%s)", szVarX);
            return( ERR_PARAM_VALUE );
        }
        else
        {
            _log_debug(_HI_, "Covariate(%d), nCol=%d, %s", i, nIdx, szVarX);
            CFmVector& VctX = pFmPhe->GetCol(nIdx);
            VctX.Rearrange( subjOrd );
            m_pCovars->SetCol(i, VctX, szVarX );
        }
    }

    RemoveMissing(pPackedSNP);

    m_pPhenoZ0 = new (refNew) CFmMatrix(m_pPhenoZ);

    if (!m_bZnorm)
    {
        // -- Find minimum and maximum value in the age data
        double Z_max=0;
        double Z_min=100;
        for (int i=0; i<m_pPhenoZ->GetNumRows(); i++)
        for (int j=0; j<m_pPhenoZ->GetNumCols(); j++)
        {
            double tmp_z = m_pPhenoZ->Get(i,j);
            if ( tmp_z != 0 )
            {
                if (tmp_z>Z_max) Z_max = tmp_z;
                if (tmp_z<Z_min) Z_min = tmp_z;
            }
        }

        _log_debug( _HI_, "Phenotype: bZNorm=%d Range=(%.2f, %.2f)",
                    m_bZnorm?1:0, Z_min, Z_max);

        for (int i=0; i<m_pPhenoZ->GetNumRows(); i++)
        for (int j=0; j<m_pPhenoZ->GetNumCols(); j++)
        {
            double tmp_z = m_pPhenoZ->Get(i,j);
            if (isnan(tmp_z))
                continue;

            double tp = (tmp_z - Z_min)*2/(Z_max-Z_min)-1;
            m_pPhenoZ->Set( i, j, tp );
        }
    }

	if( m_szFile_pheno != NULL)
		destroy( pFmPhe );

    _log_info( _HI_, "Phenotype: %d subjects are read from phenotype file", m_pPhenoY->GetNumRows() );

    return(0);
}

int CFmDat_Pheno::LoadNonlongdt( CFmPackedSNP* pPackedSNP, CFmVectorStr* pFamSubjs )
{
    //include response value(y);
    if ( m_pszYname == NULL )
    {
        _log_error(_HI_, "No response value or  covariate in the model parameter");
        return( ERR_OPEN_FILE );
    }

	CFmMatrix* pFmPhe=NULL;
	if( m_szFile_pheno!=NULL)
	{
	    CFmDataFrame df;
		int ret = df.Load(m_szFile_pheno, false, true);
		if (ret!=0)
		{
			_log_error(_HI_, "Failed to open the Phenotype file(%s)", m_szFile_pheno);
			return( ERR_OPEN_FILE );
		}

		CFmVector vctCol(0, 0.0);
		for(int i=1;i<df.GetNumCol();i++)	vctCol.Put(i);
		pFmPhe = df.GetMatrix( &vctCol );
		pFmPhe->SetRowNames( df.GetStringCol(0) );
	}
	else
		pFmPhe = m_pAttachedPhe;

	CFmNewTemp refNew;
	m_pSubjs = new (refNew) CFmVectorStr( pFmPhe->GetRowNames() );

	CFmVector subjOrd(0, 0.0);
	if (pFamSubjs)
	{
		for(int i=0;i<m_pSubjs->GetLength(); i++)
		{
			int nPos = pFamSubjs->Find ( m_pSubjs->Get(i) );
			if ( nPos>=0  )
				subjOrd.Put( nPos+1 );
			else
				_log_error(_HI_, "Failed to find %s in phenotype", m_pSubjs->Get(i));
		}
	}

    m_nMesuQ = 1;
    m_nSubjN = pPackedSNP->GetNumSubjs();
    if ( pFmPhe->GetNumRows() != m_nSubjN )
    {
        _log_error(_HI_, "The id's count is greater than the count in SNP data(%d:%d)", pFmPhe->GetNumRows(), m_nSubjN);
        return(-1);
    }

    if (m_pPhenoY) destroy( m_pPhenoY );

    // Y matrix
    m_pPhenoY = new (refNew) CFmMatrix( m_nSubjN, 1 );
    int nIdx = pFmPhe->FindColumn( m_pszYname );
    if (nIdx<0)
    {
        _log_error(_HI_, "Failed to find the column(%s)", m_pszYname );
        return( ERR_OPEN_FILE );
    }
    else
    {
        CFmVector& VctY = pFmPhe->GetCol(nIdx);

		if (subjOrd.GetLength()>0) VctY.Rearrange( subjOrd );
        m_pPhenoY->SetCol(0, VctY );
    }

    // Covar matrix
    m_pCovars = new (refNew) CFmMatrix( m_nSubjN, m_pXnames->GetLength() + 1  );
	m_pCovars->SetColName( 0, (char*)"Mu" );
    for(int i=0; i<m_nSubjN; i++)
		m_pCovars->Set( i, 0, 1.0);

    for(int i=1; i< m_pXnames->GetLength() + 1;i++)
    {
        char* szVarX = m_pXnames->Get(i-1 );

        int nIdx = pFmPhe->FindColumn( szVarX );
        if (nIdx<0)
        {
            _log_error(_HI_, "Failed to find the column(%s)", szVarX);
            return( ERR_OPEN_FILE );
        }
        else
        {
            CFmVector& VctX = pFmPhe->GetCol(nIdx);
			if (subjOrd.GetLength()>0)
                VctX.Rearrange( subjOrd );
            m_pCovars->SetCol(i, VctX, szVarX );
        }
    }

	RemoveMissing(pPackedSNP);

	if( m_szFile_pheno != NULL)
		destroy( pFmPhe );

    _log_info( _HI_, "Phenotype: %d subjects are read from phenotype file", m_pPhenoY->GetNumRows() );

    return(0);
}


int CFmDat_Pheno::RemoveMissing( CFmPackedSNP* pPackedSNP )
{
    CFmVector phe_null(0, 0.0);

    // Y Matrix
	CFmVector vctY(0, 0.0);

	if (m_pPhenoY->GetNumCols()==1)
    {
		for(int i=0; i<m_pPhenoY->GetNumCols(); i++)
		{
			vctY = m_pPhenoY->GetCol(i );
			for(int i=0;i<vctY.GetLength(); i++)
				if ( isnan( vctY.Get(i) ) )
					phe_null.UniquePut(i);
		}
	}
	else
	{
		CFmVector vct_ord(0, 0.0);
		// Z Matrix
		for(int i=0; i<m_pPhenoZ->GetNumRows(); i++)
		{
			vct_ord.Resize(0, true);
			for(int j=0; j<m_pPhenoZ->GetNumCols(); j++)
				//if ( m_pPhenoZ->Get(i,j) != 0.0 &&
				//	 m_pPhenoY->Get(i,j) != 0.0 )
				if ( !isnan( m_pPhenoZ->Get(i,j)) &&
					 !isnan( m_pPhenoY->Get(i,j)) )
				vct_ord.Put(j);

			if (vct_ord.GetLength()==0)
					phe_null.UniquePut(i);
		}
	}

    // Covar Matrix
    for(int i=0; i<m_pCovars->GetNumCols(); i++)
    {
        vctY = m_pCovars->GetCol(i );
		for(int j=0;j<vctY.GetLength(); j++)
			if ( isnan( vctY.Get(j) ) )
				phe_null.UniquePut(j);
    }

    _log_info( _HI_, "RemoveMissing: %d subjects will be removed.", phe_null.GetLength() );
    if (phe_null.GetLength()>0)
    {
        if (m_pPhenoZ) m_pPhenoZ->RemoveRows( phe_null );

        m_pCovars->RemoveRows( phe_null );
        m_pPhenoY->RemoveRows( phe_null );
        pPackedSNP->RemoveSubjs( phe_null );
	    m_nSubjN = pPackedSNP->GetNumSubjs();
    }

    _log_info( _HI_, "RemoveMissing: %d subjects are left after removing empty rows", m_pPhenoY->GetNumRows() );

    return(0);
}

void destroy(CFmDat_Plink* p)
{
	CFmNewTemp  fmRef;
	p->~CFmDat_Plink();
	operator delete(p, fmRef);
}

void destroy(CFmDat_Simple* p)
{
	CFmNewTemp  fmRef;
	p->~CFmDat_Simple();
	operator delete(p, fmRef);
}

void destroy(CFmDat_Pheno* p)
{
	CFmNewTemp  fmRef;
	p->~CFmDat_Pheno();
	operator delete(p, fmRef);
}
