/* gls_dat.cpp  -	GLS Data Object
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

#include "fm_rlogger.h"
#include "fm_err.h"
#include "fm_dat.h"
#include "fm_simulate.h"
#include "fm_snpmat.h"
#include "fm_matrix.h"
#include "fm_vector.h"
#include "fm_packedsnp.h"
#include "fm_dat.h"
#include "fm_new.h"

#include "gls_cfg.h"
#include "gls_dat.h"
#include "gls_par.h"
#include "gls_test.h"

GLS_dat::GLS_dat(CMDOPTIONS* pCmd)
{
    m_pSnpNames = NULL;
    m_pPhenoY = NULL;
    m_pPhenoZ = NULL;
    m_pPhenoZ0= NULL;
    m_pCovars = NULL;
    m_pSimuSnps = NULL;
    m_pSnpNames = NULL;

    m_pPlink = NULL;
    m_pSimple = NULL;
    m_pSimulate = NULL;

    m_bSimulate = false;
    m_nSubjN  = 0;
    m_nSnpP   = 0;
    m_nMesuTime = 0;
    m_nRound = 1;

    m_sPheno_file = NULL;
    m_sGeno_file = NULL;

    m_pZRange = NULL;
    m_pCmd = pCmd;
}

GLS_dat::~GLS_dat()
{
    ResetAllObjects();
    _log_debug(_HI_, "GLS_dat is released successfully.");
}

void GLS_dat::ResetAllObjects()
{
    if (m_pPlink) destroy( m_pPlink );
    if (m_pSimple) destroy( m_pSimple );
    if (m_pSimulate) destroy( m_pSimulate );
    if (m_pPhenoY) destroy( m_pPhenoY );
    if (m_pPhenoZ) destroy( m_pPhenoZ );
    if (m_pPhenoZ0) destroy( m_pPhenoZ0 );
    if (m_pCovars) destroy( m_pCovars );

    m_pSnpNames = NULL;
    m_pPhenoY = NULL;
    m_pPhenoZ = NULL;
    m_pPhenoZ0 = NULL;
    m_pCovars = NULL;
    m_pSimuSnps = NULL;
    m_pSnpNames = NULL;

    m_bSimulate = false;
    m_nSubjN  = 0;
    m_nSnpP   = 0;
    m_nMesuTime = 0;
    m_nRound = 1;

    if (m_sPheno_file) Free( m_sPheno_file );
    m_sPheno_file = NULL;
    if (m_sGeno_file) Free( m_sGeno_file );
    m_sGeno_file = NULL;
}

char* GLS_dat::GetSnpName(int idx)
{
    if (idx<0 || idx>=m_pSnpNames->GetLength())
        return(NULL);

    return (m_pSnpNames->Get(idx));
}

int GLS_dat::LoadPlink( char* szFile_tped, char* szFile_tfam, char* szFile_pheno, bool bZnorm  )
{
    ResetAllObjects();

    m_sPheno_file = Strdup(szFile_pheno);
    m_sGeno_file = Strdup(szFile_tped);

    _log_info( _HI_, "LoadPlink: Start to read the TPED file(%s)", szFile_tped);

	CFmNewTemp refNew;
    m_pPlink = new (refNew) CFmDat_Plink( szFile_tped, szFile_tfam);
    int ret = m_pPlink->Load();
    if (ret!=0)
    {
        _log_info(_HI_, "LoadPlink: Failed to load %s", szFile_tped);
        return(ret);
    }

    m_pPackedSNP = &(m_pPlink->m_PackedSNP);
    m_pSnpNames = m_pPlink->m_pSnpNames;
    m_nSubjN = m_pPackedSNP->GetNumSubjs();
    m_nSnpP = m_pPackedSNP->GetNumSnps();
    m_nRound =1;
    m_bSimulate = false;

    _log_info( _HI_, "LoadPlink: Start to read the phenotype file(%s, model=%s,%s,%s Znorm=%d)", szFile_pheno, m_pCmd->szYname, m_pCmd->szZname, m_pCmd->szXname, bZnorm?1:0 );

    CFmDat_Pheno* pPheno = new (refNew) CFmDat_Pheno( szFile_pheno, bZnorm, m_pCmd->szYname, m_pCmd->szZname, m_pCmd->szXname );
    ret = pPheno->LoadLongdt( m_pPackedSNP, m_pPlink->m_pSubIds  );
    if (ret!=0)
    {
        _log_info(_HI_, "LoadSnpFile: Failed to load %s", szFile_pheno);
        return(ret);
    }

    m_nSubjN  = pPheno->m_nSubjN;
    m_pPhenoY = new (refNew) CFmMatrix( pPheno->m_pPhenoY, false );
    m_pPhenoZ = new (refNew) CFmMatrix( pPheno->m_pPhenoZ, false );
    m_pPhenoZ0 = new (refNew) CFmMatrix( pPheno->m_pPhenoZ0, false );
    m_pCovars = new (refNew) CFmMatrix( pPheno->m_pCovars );
    m_nMesuTime = pPheno->m_nMesuQ;

    //set 0.0 as NaN
    m_pPhenoY->ResetNanvalue( 0.0 );
    m_pPhenoZ->ResetNanvalue( 0.0 );
    m_pPhenoZ0->ResetNanvalue( 0.0 );

    destroy( pPheno );

    _log_info( _HI_, "LoadPlink: Successful, nSubjN=%d", m_nSubjN  );
    return(0);
}

int GLS_dat::LoadSimple( char* szFile_snp, char* szFile_pheno, bool bZnorm  )
{
    ResetAllObjects();

    m_sPheno_file = Strdup(szFile_pheno);
    m_sGeno_file = Strdup(szFile_snp);

    _log_info( _HI_, "LoadSimple: Start to read the SNP file(%s)", szFile_snp);

	CFmNewTemp refNew;
    m_pSimple = new (refNew) CFmDat_Simple( szFile_snp );
    int ret = m_pSimple->Load();
    if (ret!=0)
    {
        _log_info(_HI_, "LoadSimple: Failed to load %s", szFile_snp);
        return(ret);
    }

    m_pPackedSNP = &(m_pSimple->m_PackedSNP);
    m_pSnpNames = m_pSimple->m_pSnpNames;
    m_nSubjN = m_pPackedSNP->GetNumSubjs();
    m_nSnpP = m_pPackedSNP->GetNumSnps();
    m_nRound =1;
    m_bSimulate = false;

    _log_info( _HI_, "LoadSimple: Start to read the phenotype file(%s, Znorm=%d)", szFile_pheno, bZnorm?1:0 );

    CFmDat_Pheno* pPheno = new (refNew) CFmDat_Pheno( szFile_pheno, bZnorm, m_pCmd->szYname, m_pCmd->szZname, m_pCmd->szXname );
    ret = pPheno->LoadLongdt( m_pPackedSNP, m_pSimple->m_pSubIds );
    if (ret!=0)
    {
        _log_info(_HI_, "LoadSimple: Failed to load %s", szFile_pheno);
        return(ret);
    }

    m_nSubjN   = pPheno->m_nSubjN;
    m_pPhenoY  = new (refNew) CFmMatrix( pPheno->m_pPhenoY, false );
    m_pPhenoZ  = new (refNew) CFmMatrix( pPheno->m_pPhenoZ, false );
    m_pPhenoZ0 = new (refNew) CFmMatrix( pPheno->m_pPhenoZ0, false );
    m_pCovars  = new (refNew) CFmMatrix( pPheno->m_pCovars, false );
    m_nMesuTime= pPheno->m_nMesuQ;

    destroy( pPheno );

    _log_info( _HI_, "LoadSimple: Successful, nSubjN=%d", m_nSubjN  );
    return(0);
}

int GLS_dat::AttachSnpmat( CFmMatrix* pFmPhe, CFmMatrix* pFmSnp, bool bZnorm  )
{
    ResetAllObjects();

    m_sPheno_file = Strdup("SNPMAT_dataset.csv");
    m_sGeno_file = Strdup("SNPMAT_dataset.csv");

    _log_info( _HI_, "AttachSnpmat: Start to read the SNPMAT" );

	CFmNewTemp refNew;
    m_pPackedSNP = new (refNew) CFmPackedSNP();
    m_pPackedSNP->InitData( pFmSnp->GetNumCols()-2, pFmSnp->GetNumRows() );

    for(int iSnp=0; iSnp<pFmSnp->GetNumRows(); iSnp++)
    {
        for(int jSub=0; jSub<pFmSnp->GetNumCols()-2; jSub++)
        {
            int snp_d = (int)(pFmSnp->Get( iSnp, jSub+2));
            if (snp_d != 0 && snp_d != 1 && snp_d!= 2)
                snp_d=3;

            m_pPackedSNP->Set( jSub, iSnp, snp_d);
            if ( m_pPackedSNP->Get( jSub, iSnp ) != snp_d )
            {
                _log_fatal(_HI_, "CFmPackedSNP: Failed to check snp %d(%d != %d)", iSnp, snp_d, m_pPackedSNP->Get(jSub, iSnp ) );
            }
        }

        char szChr[64] = {"0"};
        char szPos[64] = {"0"};
        sprintf(szChr, "%d", (int)pFmSnp->Get( iSnp, 0) );
        sprintf(szPos, "%d", (int)pFmSnp->Get( iSnp, 1) );

        m_pPackedSNP->SetSnpInfo(iSnp, pFmSnp->GetRowName(iSnp), (char*)szChr, (char*)szPos, (char*)szPos, 'C' );
    }

    int nSnpP = m_pPackedSNP->GetNumSnps();
    int ret = m_pPackedSNP->RemoveRareSNP( 0.01 );
    if (ret!=0)
    {
		_log_error( _HI_, "CFmPackedSNP: Failed to remove rare SNPs");
		return( ERR_SNPFILE_LONG );
	}
	else
	    _log_prompt( _HI_, "Total SNPs: %d, Rare SNPs: %d, Rare MAF= %.3f", nSnpP, nSnpP - m_pPackedSNP->GetNumSnps(), 0.01 );

    m_pSnpNames = new (refNew) CFmVectorStr( m_pPackedSNP->GetNumSnps() );
    for (int i=0;i<m_pPackedSNP->GetNumSnps();i++)
    {
        SNPINFO* pInfo;
        if (m_pPackedSNP->GetSnpInfo(i, &pInfo))
        {
            m_pSnpNames->Set(i, pInfo->szSnpId);
        }
    }

    m_nSubjN = m_pPackedSNP->GetNumSubjs();
    m_nSnpP  = m_pPackedSNP->GetNumSnps();
    m_nRound = 1;
    m_bSimulate = false;

    _log_info( _HI_, "AttachSnpmat: Start to read the phenotype file(%s, Znorm=%d)", m_sPheno_file, bZnorm?1:0 );

    CFmVectorStr* pSubids = new (refNew) CFmVectorStr( pFmSnp->GetColNames() );
    pSubids->Remove(1);
    pSubids->Remove(0);

    CFmDat_Pheno* pPheno = new (refNew) CFmDat_Pheno( pFmPhe, bZnorm, m_pCmd->szYname, m_pCmd->szZname, m_pCmd->szXname );
    ret = pPheno->LoadLongdt( m_pPackedSNP,  pSubids );
    if (ret!=0)
    {
        _log_info(_HI_, "AttachSnpmat: Failed to load SNPMAT");
        return(ret);
    }

    m_nSubjN    = pPheno->m_nSubjN;
    m_pPhenoY   = new (refNew) CFmMatrix( pPheno->m_pPhenoY, false );
    m_pPhenoZ   = new (refNew) CFmMatrix( pPheno->m_pPhenoZ, false );
    m_pPhenoZ0  = new (refNew) CFmMatrix( pPheno->m_pPhenoZ0, false );
    m_pCovars   = new (refNew) CFmMatrix( pPheno->m_pCovars, false );
    m_nMesuTime = pPheno->m_nMesuQ;

    destroy( pPheno );
    destroy( pSubids );

    _log_info( _HI_, "AttachSnpmat: Successful, nSubjN=%d", m_nSubjN  );

    return(0);
}

int GLS_dat::LoadPhenoOnly(char* szFile_pheno, bool bZnorm )
{
    _log_info( _HI_, "LoadPhenoOnly: Start to read the phenotype file(%s, Znorm=%d)", szFile_pheno, bZnorm?1:0 );

    CFmVectorStr* pSubjs = NULL;
    if (m_pPlink)
        pSubjs = m_pPlink->m_pSubIds;
    if (m_pPlink)
        pSubjs = m_pSimple->m_pSubIds;

	CFmNewTemp refNew;
    CFmDat_Pheno* pPheno = new (refNew) CFmDat_Pheno( szFile_pheno, bZnorm, m_pCmd->szYname, m_pCmd->szZname, m_pCmd->szXname );
    int ret = pPheno->LoadNonlongdt( m_pPackedSNP, pSubjs );
    if (ret!=0)
    {
        _log_info(_HI_, "LoadPhenoOnly: Failed to load %s", szFile_pheno);
        return(ret);
    }

    m_nSubjN  = pPheno->m_nSubjN;
    m_pPhenoY = new (refNew) CFmMatrix( pPheno->m_pPhenoY );
    m_pCovars = new (refNew) CFmMatrix( pPheno->m_pCovars );
    return(0);
}

int GLS_dat::Simulate_test(GLS_par* par, CMDOPTIONS* pCmd )
{
    ResetAllObjects();
    m_sPheno_file = Strdup( "simu.pheno.LS");
    m_sGeno_file = Strdup( "simu.geno.LS");

    _log_info( _HI_, "Simulation_test: Initialize data object. nTimes:%d", m_nMesuTime );

    //--generate snp data
    _log_info( _HI_, "Simulation_test: Generate genotypical data.");
    proc_gls_simu_geno_test( m_pSimuSnps );

    //--generate phenotype data
    _log_info( _HI_, "Simulation_test: Generate phenotypical data.");
    proc_gls_simu_pheno_test( m_pPhenoY, m_pPhenoZ, m_pCovars, m_pZRange );

	CFmNewTemp refNew;
    m_pSnpNames = new (refNew) CFmVectorStr(m_nSnpP);
    for (int i=0;i<m_nSnpP;i++)
    {
        char szBuf[64]={0};
        snprintf(szBuf, 63, "SNP%d", i+1);
        m_pSnpNames->Set(i, szBuf);
    }

    if (strlen(pCmd->szPheoutFile))
    {
    }

    if (strlen(pCmd->szSnpoutFile))
    {
    }

    return(0);
}

int GLS_dat::Simulate( GLS_par* par, CMDOPTIONS* pCmd )
{
    ResetAllObjects();

    m_sPheno_file = Strdup( "simu.pheno.LS");
    m_sGeno_file = Strdup( "simu.geno.LS");

    _log_info( _HI_, "Simulate: Start to Simulate");

    if ( par->simu_n * par->simu_p* 1  == 0)
        return (ERR_PARAM_VALUE);

    SIMU_PARAM_LONGDT simu_par;
    simu_par.simu_grps      = par->simu_grps;
    simu_par.simu_n         = par->simu_n;
    simu_par.simu_p         = par->simu_p;
    simu_par.simu_snp_rho   = par->simu_snp_rho;
    simu_par.simu_snp_miss  = par->simu_snp_missing;
    simu_par.simu_rho       = par->simu_rho;
    simu_par.simu_sigma2    = par->simu_sigma2;
    simu_par.simu_a_len     = par->simu_a_len;
    simu_par.simu_d_len     = par->simu_d_len;
    simu_par.simu_covar_len = par->simu_covar_len;
    simu_par.sig_p          = par->sig_p;

    memcpy( simu_par.simu_mu          ,par->simu_mu           ,sizeof(par->simu_mu));
    memcpy( simu_par.simu_covar_effect,par->simu_covar_effect ,sizeof(par->simu_covar_effect));
    memcpy( simu_par.simu_covar_range ,par->simu_covar_range  ,sizeof(par->simu_covar_range));
    memcpy( simu_par.simu_a_pos       ,par->simu_a_pos        ,sizeof(par->simu_a_pos));
    memcpy( simu_par.simu_a_effect    ,par->simu_a_effect     ,sizeof(par->simu_a_effect));
    memcpy( simu_par.simu_d_pos       ,par->simu_d_pos        ,sizeof(par->simu_d_pos));
    memcpy( simu_par.simu_d_effect    ,par->simu_d_effect     ,sizeof(par->simu_d_effect));
    memcpy( simu_par.simu_z_range     ,par->simu_z_range      ,sizeof(par->simu_z_range));
    memcpy( simu_par.simu_z_count     ,par->simu_z_count      ,sizeof(par->simu_z_count));

	CFmNewTemp refNew;
	CFmSimulate* m_pSimulate = new (refNew) CFmSimulate( &simu_par );
    int nRet = m_pSimulate->Simulate( pCmd->szSnpoutFile, pCmd->szPheoutFile, (char*)"0" );
    if (nRet)
        return(nRet);

    m_pSimuSnps = m_pSimulate->m_pSimuSnps;
    m_pSnpNames = m_pSimulate->m_pSnpNames;
    m_nSubjN = m_pSimulate->m_nSubjN;
    m_nSnpP = m_pSimulate->m_nSnpP;
    m_nRound =1;
    m_bSimulate = true;

    _log_info( _HI_, "Simulate: successful" );

    return(nRet);
}

int GLS_dat::GetPartialSNP( CFmSnpMat* pMat, int nSnpStart, int nSnpCnt )
{
    _log_debug( _HI_, "GetPartialSNP(%d,%d), m_bSimulate=%d",nSnpStart, nSnpCnt, m_bSimulate?1:0 );
    if (m_bSimulate)
    {
        char szChr[16]={"1"};
        char szPos[16]={0};

        for (int i=nSnpStart;i<nSnpStart+nSnpCnt; i++)
        {
            CFmVector& vct = m_pSimuSnps->GetRow(i);
            sprintf( szPos, "%d", i);
            pMat->SetSnp(i, vct, m_pSimuSnps->GetRowName(i), szChr, szPos);
        }

        return(0);
    }

    if (pMat->GetMaxSubjs()!=m_nSubjN || pMat->GetMaxSnps() < nSnpCnt )
    {
        _log_debug( _HI_, "CFmMatrix does not have enough space(%d,%d)",nSnpCnt, m_nSubjN );
        return( ERR_MATRIX_DAT);
    }

    double* fSNPs = new double[m_nSubjN];
    for (int i=0; i<nSnpCnt;i++)
    {
        m_pPackedSNP->GetSnpRow( nSnpStart+i, fSNPs);

        SNPINFO* pInfo=NULL;
        char* pSnpId = NULL;
        char* pChr = NULL;
        char* pPos = NULL;

        if ( m_pPackedSNP->GetSnpInfo(nSnpStart+i, &pInfo) && pInfo)
        {
            pSnpId = pInfo->szSnpId;
            pChr   = pInfo->szChrId;
            pPos   = pInfo->szBasePair;
        }

        for(int k=0;k<m_nSubjN;k++)
            fSNPs[k] = fSNPs[k]-1;

        pMat->SetSnp(i, fSNPs, m_nSubjN, pSnpId, pChr, pPos );
    }

    return(0);
}

int GLS_dat::SelectRefitGenos(CFmVector* pSelSnp, CFmSnpMat* pMat, bool bGenNorm)
{
    _log_debug(_HI_, "SelectRefitGenos:%d bSimulate=%d bGenNorm=%d", pSelSnp->GetLength(),
               m_bSimulate?1:0, bGenNorm?1:0);


    _log_debug(_HI_, "SelectRefitGenos:pMat[%d,%d]", pMat->GetNumSnps(),pMat->GetNumSubjs() );

    double* fSNPs = new double[m_nSubjN];
    for(int i=0; i<pSelSnp->GetLength(); i++)
    {
        if (m_bSimulate)
        {
            CFmVector vct( m_pSimuSnps->GetNumCols(), 0 );
            vct= m_pSimuSnps->GetRow( (int)pSelSnp->Get(i));
            char *szRowName= m_pSimuSnps->GetRowName( (int)pSelSnp->Get(i));
            char szPos[16]={0};
            sprintf(szPos, "%d", i);
            pMat->SetSnp(i, vct, szRowName, (char*)"1", szPos);
        }
        else
        {
            if ( m_pPackedSNP->GetSnpRow( (int)pSelSnp->Get(i), fSNPs) )
            {
                if (bGenNorm)
                    for(int k=0;k<m_nSubjN;k++)
                        fSNPs[k] = fSNPs[k]-1;

                SNPINFO *pInfo = NULL;
                char* pSnpId = NULL;
                char* pChr   = NULL;
                char* pPos   = NULL;

                if ( m_pPackedSNP->GetSnpInfo( (int)pSelSnp->Get(i), &pInfo) && pInfo)
                {
                    pSnpId = pInfo->szSnpId;
                    pChr = pInfo->szChrId;
                    pPos   = pInfo->szBasePair;
                }

                _log_debug(_HI_, "SelectRefitGenos:pMat[%d,%s, Chr:%s, Pos:%s]", i, pSnpId, pChr, pPos);
                pMat->SetSnp(i, fSNPs, m_nSubjN, pSnpId, pChr, pPos);
            }
            else
            {
                delete[] fSNPs;
                return( ERR_NULL_DATA );
            }
        }
    }

    delete[] fSNPs;

    return(0);
}

int GLS_dat::AppendSnp(CFmSnpMat* pMat)
{
    _log_debug(_HI_, "AppendSnp: pMat[%d,%d]", pMat->GetNumSnps(), pMat->GetNumSubjs() );

    int nSubj = pMat->GetNumSnps();
    for(int k=0; k<pMat->GetNumSubjs(); k++ )
    {
        int nSnpIdx = -1;
        if ( m_pPackedSNP->GetNumSnps() == 0 )
        {
            m_pPackedSNP->InitData(nSubj, 1);
            nSnpIdx =0;
        }
        else
        {
            if ( m_pPackedSNP->GetNumSubjs() != nSubj)
            {
                _log_error(_HI_ , "AppendSnp: The TPED file does not have same subjects with SNP data.");
                return( ERR_MATRIX_DAT );
            }
            else
            {
                nSnpIdx = m_pPackedSNP->AppendSnp();
                if (nSnpIdx<0)
                {
                    _log_error(_HI_ , "Memory is not enough to append SNP(%d).", nSnpIdx);
                    return(nSnpIdx);
                }
            }
        }

        for(int i=0; i<nSubj; i++)
        {
            int snpd = pMat->Get_a(k, i);
            m_pPackedSNP->Set( i, nSnpIdx, snpd);
            if ( m_pPackedSNP->Get( i, nSnpIdx ) != snpd )
            {
                _log_fatal(_HI_, "AppendSnp: Failed to check snp %d(%d != %d)",  i, snpd, m_pPackedSNP->Get( i, nSnpIdx ) );
            }
        }

        char szSnpName[256]={0};
        char szChr[256]={0};
        char szPos[256]={0};
        char szEmpty[]="";

        pMat->GetSnpInfo( k, szSnpName, szChr, szPos );
        m_pPackedSNP->SetSnpInfo(nSnpIdx, szChr, szSnpName, szEmpty, szPos, '0' );
    }

    m_nSubjN = m_pPackedSNP->GetNumSubjs();
    m_nSnpP = m_pPackedSNP->GetNumSnps();

    _log_info(_HI_ , "AppendSnp: m_PackedSNP[ SnpP=%d, SubjN=%d ]", m_nSnpP, m_nSubjN );

    return( 0 );
}

int GLS_dat::Summary(char* szOutFile)
{
    char szTemp[MAX_PATH*256] = {0};

    sprintf(szTemp, "%s \nSUMMARY(Data)\n",     szTemp);
    sprintf(szTemp, "%s -------------------------------------------\n", szTemp );
    sprintf(szTemp, "%s DATA Simulate    = %s\n", szTemp, m_bSimulate?"Yes":"No");
    sprintf(szTemp, "%s DATA Subjects    = %d\n", szTemp, m_nSubjN);
    sprintf(szTemp, "%s DATA Snps        = %d\n", szTemp, m_nSnpP);
    sprintf(szTemp, "%s DATA Times       = %d\n", szTemp, m_nMesuTime);
    sprintf(szTemp, "%s DATA PhenFile    = %s\n", szTemp, m_sPheno_file);
    sprintf(szTemp, "%s DATA GenoFile    = %s\n", szTemp, m_sGeno_file);
    sprintf(szTemp, "%s -------------------------------------------\n", szTemp );

    _log_info(_HI_, szTemp);

    if(szOutFile!=NULL)
    {
        FILE* fp = fopen(szOutFile, "at");
        if (fp==NULL)
        {
            _log_error(_HI_, "Failed to open file(%s) to append the summary information", szOutFile);
        }
        else
        {
            fwrite(szTemp, strlen(szTemp),1, fp);
            fclose(fp);
        }
    }

    return(0);
}

void destroy(GLS_dat* p)
{
	CFmNewTemp  fmRef;
	p->~GLS_dat();
	operator delete(p, fmRef);
}
