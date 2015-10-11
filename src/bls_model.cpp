/* BLS_model.cpp  -	LS2 Computational model
 *
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */


#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>

#include "Rmethods.h"
#include "fm_rlogger.h"
#include "fm_pcf.h"
#include "fm_sys.h"
#include "fm_rtdump.h"
#include "fm_rsource.h"
#include "fm_err.h"
#include "fm_new.h"

#include "bls_model.h"
#include "bls_dat.h"
#include "bls_cfg.h"

#define _t(x) ((x).GetTransposed())

BLS::BLS()
{
    m_nDataType = RUNMODE_UNK;
    m_bRefit = false;

    m_pDat = NULL;
    m_pRes = NULL;
    m_pCfg = NULL;
}

BLS::~BLS()
{
    if (m_pDat) destroy( m_pDat );
	if (m_pRes) destroy( m_pRes );
    //dont delete this point because its creator is not this class.
    //if (m_pCfg) destroy( m_pCfg );

    _log_debug(_HI_, "BLS is released successfully.");
}

int BLS::LoadSimulate( CMDOPTIONS *pCmd, BLS_par* pPar )
{
	int ret=0;

    m_pCmd = pCmd;
    m_nDataType = RUNMODE_SIM;

    _log_prompt(_HI_, "Load the parameter file:%s", pCmd->szParFile);

    CFmPcf pcf;
    pcf.UpdatePcfFile(PCF_PAR_READ );

    pPar->Summary();

    _log_prompt(_HI_, "Start the data generation.");

	CFmNewTemp  refNew;
    m_pDat = new (refNew) BLS_dat( pCmd );

    pcf.UpdatePcfFile(PCF_DAT_LOAD );
    //Simulate_test??
    if ( (ret = m_pDat->Simulate( pPar, pCmd)) !=0)
    {
        _log_error(_HI_, "Failed to do data generation.");
        return(ret);
    }

    _log_prompt(_HI_, "Data generation is successful.");
    return(0);
}

int BLS::LoadPlink( CMDOPTIONS *pCmd )
{
	//!!!TEST for R
	R_set_seed( 100 );

    m_pCmd = pCmd;
    m_nDataType = RUNMODE_PLINK;

    _log_info(_HI_, "Check TPED/TFAM/phenotype parameter.");

    if (strlen(pCmd->szTpedFile)>=MAX_PATH)
        return( ERR_TPEDFILE_LONG );

    if (strlen(pCmd->szTfamFile)>=MAX_PATH)
        return( ERR_TFAMFILE_LONG );

    if (strlen(pCmd->szPheFile)>=MAX_PATH)
        return( ERR_TPHEFILE_LONG );

    _log_prompt(_HI_, "Load the data sets from PLINK files.\nTPED=%s\nTFAM=%s\nPHE=%s\n",
                pCmd->szTpedFile, pCmd->szTfamFile, pCmd->szPheFile );

    CFmPcf pcf;
    pcf.UpdatePcfFile(PCF_DAT_LOAD );

	CFmNewTemp  refNew;
    m_pDat = new (refNew) BLS_dat(pCmd);
    int ret = m_pDat->LoadPlink( pCmd->szTpedFile,
                                 pCmd->szTfamFile,
                                 pCmd->szPheFile,
                                 pCmd->bZNormalized,
                                 pCmd->szPreSigFile);
    if (ret!=0)
    {
        _log_error(_HI_, "Failed to load the PLINK files.");
        return(ret);
    }


    _log_prompt(_HI_, "PLINK data are loaded successfully.");
    m_pDat->Summary();

    return(0);
}

int BLS::LoadSimple( CMDOPTIONS *pCmd )
{
    m_pCmd = pCmd;
    m_nDataType = RUNMODE_SIMPLE;

    _log_info(_HI_, "Check SNP/phenotype parameter.");

    if (strlen(pCmd->szSnpFile)>=MAX_PATH)
        return( ERR_TPEDFILE_LONG );

    if (strlen(pCmd->szPheFile)>=MAX_PATH)
        return( ERR_TPHEFILE_LONG );

    _log_prompt(_HI_, "Load the data sets from simple SNP files");

    CFmPcf pcf;
    pcf.UpdatePcfFile(PCF_DAT_LOAD );

	CFmNewTemp  refNew;
    m_pDat = new (refNew) BLS_dat(pCmd);
    int ret = m_pDat->LoadSimple(pCmd->szSnpFile,
                                 pCmd->szPheFile,
                                 pCmd->bZNormalized,
                                 pCmd->szPreSigFile);
    if (ret!=0)
    {
        _log_error(_HI_, "Failed to load the PLINK files.");
        return(ret);
    }

    _log_prompt(_HI_, "Simple data are loaded successfully.");

    m_pDat->Summary();
    return(0);
}

int BLS::LoadSnpmat( CFmMatrix* pFmPhe, CFmMatrix* pFmSnp, CMDOPTIONS *pCmd )
{
    m_pCmd = pCmd;
    m_nDataType = RUNMODE_SIMPLE;

    _log_info(_HI_, "Check SNP/phenotype parameter.");

    if ( pFmPhe->GetNumRows() != pFmSnp->GetNumCols()-2 )
        return( ERR_TPEDFILE_LONG );

    CFmPcf pcf;
    pcf.UpdatePcfFile(PCF_DAT_LOAD );

	CFmNewTemp  refNew;
    m_pDat = new (refNew) BLS_dat(pCmd);
    int ret = m_pDat->AttachSnpmat(pFmPhe, pFmSnp, pCmd->bZNormalized);
    if (ret!=0)
    {
        _log_error(_HI_, "Failed to load the SNPMAT data.");
        return(ret);
    }

    _log_prompt(_HI_, "SNPMAT data are loaded successfully.");

    m_pDat->Summary();
    return(0);
}

int BLS::Varsel(BLS_cfg* pCfg)
{
    // Create Result object
    if (m_pRes!=NULL) destroy( m_pRes );

	m_pCfg = pCfg;

    _log_info( _HI_, "VARSEL: PREC(3)" );

    int ret = proc_prec( 3 );
    if (ret != 0 )
       _log_error( _HI_, "VARSEL: Failed to PREC:(%d)", ret );

    int P = m_pDat->m_nSnpP;
    int N = m_pDat->m_nSubjN;
    int nMcmcSnps = m_pCfg->m_nMcmcSnps;
    int nMcmcBlocks = (int)ceil(P*1.0/nMcmcSnps);

    _log_prompt(_HI_, "VARSEL: SnpP:%d, SubjN:%d, nMcmcSnps: %d nMcmcBlocks: %d", P, N, nMcmcSnps, nMcmcBlocks );

	CFmNewTemp  refNew;

    m_bRefit = false;
    m_pRes = new (refNew) BLS_res( m_pCmd, m_pCfg );

    m_pRes->InitVarsel( m_pCmd->nSimuRound, P, N, nMcmcBlocks, nMcmcSnps, m_pCfg->m_nMcmcIter, m_pDat->m_pCovars->GetColNames() );

    CFmPcf pcf;
    pcf.UpdatePcfFile( PCF_VARSEL, 0, nMcmcBlocks, -1, -1);

    // MCMC produre are called by the partial SNPs;
    _log_info( _HI_, "VARSEL: MCMC procedure START (SNP:%d)", P );
    pcf.UpdatePcfFile( PCF_VARSEL, 0, 1, -1, -1);

    CFmFileMatrix *r_pFileMat = NULL;

    if (strlen(m_pCmd->szRtdump)>0)
    {
        sprintf(m_szFileMat, "%s%s",m_pCmd->szRtdump, "bls.fmat");
        r_pFileMat = new (refNew) CFmFileMatrix( m_szFileMat, FALSE, TRUE);
    }
    else
    {
		CFmSys::GetTempFile(m_pCmd->szTmpFmat, "", MAX_PATH);
        sprintf(m_szFileMat, "%s%s",m_pCmd->szTmpFmat, ".bls.fmat");
		r_pFileMat = new (refNew) CFmFileMatrix( m_szFileMat, FALSE, TRUE);
    }

    CFmVector vctSnpChr(0, 0.0);
    CFmVector vctSnpPos(0, 0.0);
    CFmVectorStr vctSnpName(0);

    try
    {
        CFmSnpMat gen(m_pDat->m_nSnpP, m_pDat->m_nSubjN);   /* [P, N]*/
        m_pDat->GetPartialSNP( &gen, 0, m_pDat->m_nSnpP );
        gen.GetSnpInfo(&vctSnpName, &vctSnpChr, &vctSnpPos);

        int ret = proc_mcmc( *(m_pDat->m_pPhenoY), *(m_pDat->m_pCovars), gen, r_pFileMat );
        if (ret!=0)
        {
            _log_error( _HI_, "VARSEL: MCMC procedure return error code(%d).", ret);
            return(ret);
        }

		CFmMatrix* gen_pA = gen.GetAdd();//[P, N]
		gen_pA->Transpose(); //[P, N]->[N, P]
		m_pRes->SetResult( false, gen_pA, m_pDat->m_pPhenoY->GetVar(), m_pDat->m_pCovars->GetNumCols(),
							&vctSnpName, &vctSnpChr, &vctSnpPos, r_pFileMat );

		destroy( gen_pA );

    }
    catch(const char* str)
    {
        pcf.UpdatePcfFile( PCF_EXCEPT, 0, nMcmcBlocks,  -1, -1);

        //-- keep the trace tag into the log file.
        _log_error( _HI_, "VARSEL: Exception is caught in the procedure of MCMC.Tracetag=%s, Exception=%s", m_szTraceTag, str);
            return( ERR_EXCEPTION );
    }


    _log_debug( _HI_, "VARSEL: pFileMat: FileMatrix[%d,%d] ", r_pFileMat->GetNumRows(), r_pFileMat->GetNumCols());


    if (strlen( m_pCmd->szRdataFile) )
    {
        _log_prompt(_HI_, "VARSEL: Save the results into RDATA file(rdata.out=%s)", m_pCmd->szRdataFile );
        m_pRes->SaveRData( m_pCmd->szRdataFile);
    }

    if (strlen( m_pCmd->szFigOutFile)>0 )
    {
        _log_prompt(_HI_, "VARSEL: Save the figures into PDF file(fig.out=%s)", m_pCmd->szFigOutFile );
        ExportFig(m_pCmd->szFigOutFile, false );
    }

    if (strlen( m_pCmd->szSigOutFile)>0 )
    {
        _log_prompt(_HI_, "VARSEL: Save the significant SNP into CSV file(sig.out=%s)", m_pCmd->szSigOutFile);

		char szSig1[MAX_PATH]={0};
		CFmSys::GetSiblingFile( szSig1, m_pCmd->szSigOutFile, 1 );
        m_pRes->SaveSigFile( szSig1, false );
    }

    destroy( r_pFileMat );

    _log_prompt(_HI_, "VARSEL: End" );

    return 0;
}

int BLS::Refit(BLS_cfg* pCfg)
{
    if (m_pRes==NULL)
    {
        _log_error( _HI_, "REFIT: No result object for the refit");
        return( ERR_NULL_DATA );
    }

	m_pCfg = pCfg;

    _log_info( _HI_, "REFIT: Start" );

    int ret = m_pRes->InitRefit( m_pCfg->m_nMcmcIter, m_pDat->m_pCovars->GetColNames());
    if ( ret!=0 )
    {
        _log_error( _HI_, "REFIT: Failed to select siginificant snp." );
        return( ret );
    }

    CFmVector SnpList(0, 0.0);
    ret = m_pRes->GetRefitSnp( &SnpList );
    if (ret)
    {
        _log_error( _HI_, "REFIT: Failed to get refitting SNPs" );
        return(ret);
    }

    if (SnpList.GetLength()==0 || m_pDat->m_nSubjN==0)
    {
        _log_info( _HI_, "REFIT: No snps or subjects are selected to refit.");
        return( ERR_NO_ITEMS );
    }

    CFmSnpMat gen(SnpList.GetLength(), m_pDat->m_nSubjN );
    ret = m_pDat->SelectRefitGenos(&SnpList, &gen, true);
    if ( ret !=0 )
    {
        _log_info( _HI_, "REFIT: Failed to extract SNP data info from the data object");
        return( ret );
    }

    /* Now gen is [N,P] */
    int N = m_pDat->m_nSubjN;
    int P = SnpList.GetLength();
    int nMcmcSnps = SnpList.GetLength();
    int nMcmcBlocks = 1;

	CFmNewTemp  refNew;

    _log_prompt(_HI_, "REFIT: SnpP:%d, SubjN:%d, nMcmcSnps: %d nMcmcBlocks: %d", P, N, nMcmcSnps, nMcmcBlocks );
    m_bRefit = true;

    CFmPcf pcf;
    pcf.UpdatePcfFile( PCF_REFIT, 0, 1, -1, -1);

    CFmFileMatrix *r_pFileMat = NULL;

    if (strlen(m_pCmd->szRtdump)>0)
    {
        sprintf(m_szFileMat, "%s%s",m_pCmd->szRtdump, "bls.fmat");
        r_pFileMat = new (refNew) CFmFileMatrix( m_szFileMat, FALSE, TRUE);
    }
    else
    {
		CFmSys::GetTempFile(m_pCmd->szTmpFmat, "", MAX_PATH);
        sprintf(m_szFileMat, "%s%s",m_pCmd->szTmpFmat, "bls.fmat");
		r_pFileMat = new (refNew) CFmFileMatrix( m_szFileMat, FALSE, TRUE);
    }

    try
    {
        int ret = proc_mcmc( *m_pDat->m_pPhenoY, *(m_pDat->m_pCovars), gen, r_pFileMat );

        if (ret!=0)
        {
            _log_info( _HI_, "REFIT: Failed to do the MCMC procedure.");

            destroy( r_pFileMat );
            return(ret);
        }

		CFmVector vctSnpChr(0, 0.0);
		CFmVector vctSnpPos(0, 0.0);
		CFmVectorStr vctSnpName(0);
		gen.GetSnpInfo(&vctSnpName, &vctSnpChr, &vctSnpPos);

		CFmMatrix* gen_pA = gen.GetAdd();//[P, N]
		gen_pA->Transpose(); //[P, N]->[N, P]
		m_pRes->SetResult( true, gen_pA, m_pDat->m_pPhenoY->GetVar(), m_pDat->m_pCovars->GetNumCols(),
						 &vctSnpName, &vctSnpChr, &vctSnpPos, r_pFileMat );
		destroy( gen_pA );
    }
    catch(const char* str)
    {
        pcf.UpdatePcfFile( PCF_EXCEPT, 1, 1,  -1, -1);

        //-- keep the trace tag into the log file.
        _log_error( _HI_, "REFIT: An exception is caught in the MCMC procedure. Tracetag=%s, Exception=%s", m_szTraceTag, str);

        destroy( r_pFileMat );
        return( ERR_EXCEPTION );
    }

    _log_debug( _HI_, "REFIT: pFileMat: FileMatrix[%d,%d]", r_pFileMat->GetNumRows(), r_pFileMat->GetNumCols() );

    m_pRes->GetSigSNP();
    m_pRes->Summary();

    if (strlen( m_pCmd->szRdataFile) )
    {
        _log_prompt(_HI_, "REFIT: Save the results to RDATA file(rdata.out=%s).", m_pCmd->szRdataFile);
        m_pRes->SaveRData( m_pCmd->szRdataFile);
    }

    if (strlen( m_pCmd->szFigOutFile)>0 )
    {
        _log_prompt(_HI_, "REFIT: Save the figures into PDF file(fig.out=%s)", m_pCmd->szFigOutFile );
        ExportFig(m_pCmd->szFigOutFile, true );
    }

    if (strlen( m_pCmd->szSigOutFile)>0 )
    {
        _log_prompt(_HI_, "REFIT: Save the significant SNP into CSV file(sig.out=%s)", m_pCmd->szSigOutFile);

		char szSig1[MAX_PATH]={0};
		CFmSys::GetSiblingFile( szSig1, m_pCmd->szSigOutFile, 2 );
        m_pRes->SaveSigFile( szSig1, true );
    }

    destroy( r_pFileMat );

    _log_prompt(_HI_, "REFIT: End" );

    return 0;
}

int BLS::proc_mcmc( CFmVector& Y0, CFmMatrix& Covs, CFmSnpMat& gen, CFmFileMatrix* r_pFileMat)
{
	CFmNewTemp  refNew;

    //!!!TEST for R
    //int nSeed = 100;
    int nSeed = m_pCfg->GetTempId();
    R_set_seed( nSeed );

    _log_info( _HI_, "set.seed(%d)",nSeed);

    CFmPcf pcf;
    int N = Y0.GetLength();   //*[N]
    int P = gen.GetNumSnps(); //*[N,P]
    int R = m_pCfg->m_nMcmcIter;

//**SEQTEST
    _log_debug( _HI_, "proc_mcmc() MCMC: N=%d, P=%d, Q=%d, R=%d", N, P, 1, R);

    CFmVector Y(&Y0);

    double gammaTau = 1;
    double deltaTau = 1;
    double gammaTau_st = 1;
    double deltaTau_st = 1;
    double lambda2 = 1;
    double lambda_st2 = 1;

    CFmVector sigma_Coefs( Covs.GetNumCols(), 1 );
    CFmVector Coefs( Covs.GetNumCols(), 0 );

 	double mu = Y.GetMean();
 	double sigma2 = (Y - mu).GetVar();
    static CFmVector Y_eps(0, 0.0);
    Y_eps = Y;


    CFmVector a(P, 0.0);
    CFmVector tau2(P, 1);
    CFmVector d(P, 0.0);
    CFmVector tau_st2(P, 1);

    CFmVector tmp(  P, 0.0);
    CFmVector aVar( P, 0.0);
    CFmVector dVar( P, 0.0);
    CFmVector tmp3( N, 0.0);
    CFmVector spdY_by_d(N, 0.0);
    CFmVector spdY_by_a(N, 0.0);
    CFmVector spd_pD_col(N, 0.0);
    CFmVector spd_pA_col(N, 0.0);

//---------------------------------------------------
// PART X in the R code
//---------------------------------------------------
    CFmMatrix* gen_pA = gen.GetAdd();//[P, N]
    CFmMatrix* gen_pD = gen.GetDom();//[P, N]
    gen_pA->Transpose(); //[P, N]->[N, P]
    gen_pD->Transpose(); //[P, N]->[N, P]

    _log_debug( _HI_, "proc_mcmc() gen_pA:[%d, %d]", gen_pA->GetNumRows(), gen_pA->GetNumCols() );

    CFmVector tmp_genA( P, 0.0);
    for (int j=0; j<P; j++)
       tmp_genA[j] = gen_pA->ColProd(j, j);

    CFmVector tmp_genD( P, 0.0);
    for (int j=0; j<P; j++)
       tmp_genD[j] = gen_pD->ColProd(j, j);

    strcpy(m_szTraceTag, "X");

    //**SEQTEST
    //_log_debug( _HI_, "PART X.....McmcIter: %d", m_pCfg->m_nMaxIter);

    GetRNGstate();
    for (int round=0; round<m_pCfg->m_nMcmcIter; round++)
    {
        if ( round % m_pCfg->m_nMcmcHint ==0 || round<=0 )
        {
            PCFDAT* pDat = pcf.GetPcfBlock();
            char szBuf[256]={0};
            sprintf(szBuf, "Section:%d Round:%d Progress:%.3f%% Total:%.0fs Left:%.0fs\n",
                    m_pCfg->m_nSectId, round,
                    pDat->fProgress*100,
                    pDat->fELapsedSeconds + pDat->fEstimatedSeconds,
                    pDat->fEstimatedSeconds );
            //Rprintf( szBuf );
            _log_info( _HI_, szBuf );
        }

//---------------------------------------------------
// part A in the R code
// #--updating alpha
// if (length(X) > 1 )
// {
//    alphaVar <- 1/(t(X)%*%X/sigma2+1/sigmaAlpha);
//    alphaMu  <- alphaVar*( t(X) %*% (Y-Z*beta-Xi%*%t(a)-Zeta%*%t(d) )/sigma2);
//    alpha <- rnorm(1, mean=alphaMu, sd=sqrt(alphaVar));
//}
//---------------------------------------------------
//---------------------------------------------------
// part B in the R code
//  #--updating beta
//  if (length(Z) > 1 )
//  {
//      betaVar <- 1/(t(Z)%*%Z/sigma2+1/sigmaBeta);
//      betaMu  <- betaVar*( t(Z) %*% (Y-X*alpha-Xi%*%t(a)-Zeta%*%t(d)) /sigma2);
//      beta <- rnorm(1, mean=betaMu, sd=sqrt(betaVar));
//  }
//---------------------------------------------------
        strcpy(m_szTraceTag, "AB");

        //**SEQTEST
        //_log_debug( _HI_, "PART AB.....round=%d/%d ", round, R);

        static CFmVector X0(0, 0.0);
        for(int i=0; i<Covs.GetNumCols(); i++ )
        {
            X0 = Covs.GetCol(i);
            double alphaVar = 1/( X0.Prod(X0)/sigma2 + 1 /sigma_Coefs[i] );

            Y_eps = Y;
            for (int j=0; j<Covs.GetNumCols(); j++ )
				if (i!=j)
                    Y_eps = Y_eps - Covs.GetCol(j) * Coefs[j];

            Y_eps = Y_eps - ((*gen_pA)*a).GetCol(0) - ((*gen_pD)*d).GetCol(0);

            double alphaMu  = alphaVar*( X0.Prod(Y_eps)/sigma2);
            Coefs[i]  = rnorm( alphaMu, sqrt(alphaVar) );
       }

//---------------------------------------------------
// part C in the R code
//  #--Updating additive effects: a, for each SNP
//  tmp <- c();
//  for (j in 1:p)
//      tmp <- c( tmp, t(Xi[,j])%*%Xi[,j] );
//
//  if (!refit) tmp  <- (tmp + 1./tau2);
//
//  tmp3 <- Y - X*alpha - Z*beta - Zeta%*%t(d);
//  aVar <- sigma2/tmp;
//  for (j in 1:p)
//  {
//      aMu_j <- t(Xi[,j])%*%( tmp3 - Xi%*%t(a) + Xi[,j]*a[j] )/tmp[j];
//      a[j] <- rnorm( 1, mean=aMu_j, sd=sqrt( aVar[j] ) );
//  }
//---------------------------------------------------
        strcat(m_szTraceTag, "C");

        //**SEQTEST
        //_log_debug( _HI_, "PART C.....round=%d/%d", round, R );

        //for (int j=0; j<P; j++)
        //   tmp[j] = (*gen_pA).ColProd(j, j);

        tmp = tmp_genA;
        if (!m_bRefit)
            tmp = tmp + tau2.GetReciprocal();
        aVar = tmp.GetReciprocal() * sigma2;

        Y_eps = Y;
        for (int j=0; j<Covs.GetNumCols(); j++ )
            Y_eps = Y_eps - Covs.GetCol(j) * Coefs[j];
        tmp3 = Y_eps - ((*gen_pD) * d).GetCol(0);

        spdY_by_a = ((*gen_pA)*a).GetCol(0);

        if (m_pCmd->bAddUsed)
			for (int j=0; j<P; j++)
			{
				spd_pA_col = gen_pA->GetCol(j);

				double aMu_j = spd_pA_col.Prod( tmp3 - spdY_by_a + spd_pA_col*a[j] )/tmp[j];
				double old_aj = a[j];
				double new_a = rnorm( aMu_j, sqrt( aVar[j] ) );

				//isnana(new_a) <==> new_a!=new_a
				if ( new_a != new_a )
				{
					_log_debug( _HI_, "PART C_test.....j=%d, aMu_j=%f, var_a=%f", j, aMu_j, aVar[j]);

					aVar.WriteAsCSVFile("aVar.csv");
					a.WriteAsCSVFile("a.csv");
					d.WriteAsCSVFile("d.csv");
					tmp_genA.WriteAsCSVFile("tmp_genA.csv");
					tmp.WriteAsCSVFile("tmp.csv");
					tau2.WriteAsCSVFile("tau2.csv");
					tmp3.WriteAsCSVFile("tmp3.csv");
					spdY_by_a.WriteAsCSVFile("spdY_by_a.csv");
					spd_pA_col.WriteAsCSVFile("spd_pA_col.csv");

					_log_fatal( _HI_, "FAILED by NAN Problem");
				}

				a[j] = new_a;

				{
					int N_AA = 0;
					int N_aa = 0;
					for (int k=0;k<N;k++)
					{
						if (spd_pA_col[k]>0) N_AA++;
						if (spd_pA_col[k]<0) N_aa++;
					}

					if (N_AA==0 || N_aa==0) a[j] = 0.0;
				}


				spdY_by_a = spdY_by_a + spd_pA_col*(a[j]-old_aj);
			}

//---------------------------------------------------
// part D in the R code
//  #--Updating tau2 underlying a in non-refit status
//  if (!refit)
//  {
//      for (j in 1:p)
//      {
//          InvTau2_1 <- sqrt(lambda2 * sigma2)/abs(a[j]);
//          tau2[j] <- 1/func_invGau(InvTau2_1, lambda2);
//      }
//    #--Updating lambda2 underlying tau2
//    lambda2 <- rgamma( 1, p+gammaTau, scale=sum(tau2)/2 + deltaTau );
//  }
//---------------------------------------------------
        strcat(m_szTraceTag, "D");

        //**SEQTEST
        //_log_debug( _HI_, "PART D.....round=%d/%d", round, R);

        if (!m_bRefit)
        {
            for (int j=0; j<P; j++)
            {
				if (a[j] == 0.0 )
				{
					tau2[j] = 0.0;
					continue;
				}

                double InvTau2_1 = sqrt( lambda2 * sigma2)/fabs(a[j]);
                double _tau2 = 1/func_invGau(InvTau2_1, lambda2);
                if ( _tau2 <= 0 || _tau2 != _tau2 )
                {
                    _log_debug( _HI_, "PART D_test.....j=%d, tau2=%f, lambda2=%f, sigma2=%f, a[j]=%f, mu=%f",
                                     j, tau2[j], lambda2, sigma2, a[j], InvTau2_1);
                }
                else
                    tau2[j] = _tau2;
            }

            lambda2 = rgamma( P+gammaTau, tau2.Sum()/2 + deltaTau );
            if ( lambda2 <= 0 || isnan(lambda2) )
            {
                _log_debug( _HI_, "PART D_test.....tau2*=%f, P=%f, gammaTau_st=%f, deltaTau_st=%f",
                             tau2.Sum(), P, gammaTau_st, deltaTau_st);
                _log_fatal( _HI_, "FAILED by NAN Problem");
            }
        }

//---------------------------------------------------
// part E in the R code
//  #--Updating dominant effects: d, for each SNP
//  tmp <- c();
//  for (j in 1:p)
//    tmp <- c(tmp, t(Zeta[,j])%*%Zeta[,j] );
//
//  if (!refit) tmp  <- (tmp + 1./tau_st2);
//
//  tmp3 <- Y - X*alpha - Z*beta - Xi%*%t(a);
//  dVar <- sigma2/tmp;
//  for (j in 1:p)
//  {
//      dMu_j <- t(Zeta[,j])%*%(tmp3-Zeta%*%t(d)+Zeta[,j]*d[j])/tmp[j];
//      d[j]  <- rnorm( 1, mean=dMu_j, sd=sqrt(dVar[j]));
//  }
//---------------------------------------------------
        strcat(m_szTraceTag, "E");
        //**SEQTEST
        //_log_debug( _HI_, "PART E.....round=%d/%d", round, R);

        //for (int j=0; j<P; j++)
        //    tmp[j] = (*gen_pD).ColProd(j, j);

        tmp = tmp_genD;
        if (!m_bRefit)
             tmp = tmp + tau_st2.GetReciprocal();
        dVar = tmp.GetReciprocal() * sigma2;

        Y_eps = Y;
        for (int j=0; j<Covs.GetNumCols(); j++ )
            Y_eps = Y_eps - Covs.GetCol(j) * Coefs[j];
        tmp3 = Y_eps - ((*gen_pA) * a).GetCol(0);

        spdY_by_d = ((*gen_pD)*d).GetCol(0);

        if (m_pCmd->bDomUsed)
			for (int j=0; j<P; j++)
			{
				spd_pD_col = gen_pD->GetCol(j);
				double dMu_j = spd_pD_col.Prod( tmp3 - spdY_by_d + spd_pD_col*d[j] )/tmp[j];
				double old_dj = d[j];
				double new_d = rnorm( dMu_j, sqrt(dVar[j]));

				if ( new_d != new_d )
				{
					_log_debug( _HI_, "PART E_test.....j=%d, dMu_j=%f, var_d=%f", j, dMu_j, dVar[j]);

					dVar.WriteAsCSVFile("dVar.csv");
					a.WriteAsCSVFile("a.csv");
					d.WriteAsCSVFile("d.csv");
					tmp_genD.WriteAsCSVFile("tmp_genD.csv");
					tmp.WriteAsCSVFile("tmp.csv");
					tau_st2.WriteAsCSVFile("tau_st2.csv");
					tmp3.WriteAsCSVFile("tmp3.csv");
					spdY_by_d.WriteAsCSVFile("spdY_by_d.csv");
					spd_pD_col.WriteAsCSVFile("spd_pD_col.csv");

					_log_fatal( _HI_, "FAILED by NAN Problem");
				}


				d[j]  = new_d;
				spdY_by_d = spdY_by_d + spd_pD_col*(d[j]-old_dj);
			}

//---------------------------------------------------
// part F in the R code
//  #--Updating tau_st2 underlying d in non-refit status
//  if (!refit)
//  {
//      for (j in 1:p)
//      {
//          InvTau2_1 = sqrt(lambda_st2*sigma2)/abs(d[j]);
//          tau_st2[j] = 1/func_invGau(InvTau2_1, lambda_st2);
//      }
//
//    #--Updating lambda_star2 underlying tau_st2
//    lambda_st2 <- rgamma( 1, p + gammaTau_st, scale=sum(tau_st2)/2+deltaTau_st );
//}
//---------------------------------------------------
        strcat(m_szTraceTag, "F");
        //**SEQTEST
        //_log_debug( _HI_, "PART F.....round=%d/%d", round, R);

        if (!m_bRefit)
        {
            for (int j=0; j<P; j++)
            {
                double InvTau2_1 = sqrt(lambda_st2*sigma2)/fabs(d[j]);
                double _tau_st2 = 1/func_invGau(InvTau2_1, lambda_st2);

                if ( _tau_st2 <= 0 || _tau_st2 != _tau_st2 )
                {
                    _log_debug( _HI_, "PART F_test.....j=%d, tau2*=%f, lambda2*=%f, sigma2=%f, d[j]=%f, mu=%f",
                                 j, tau_st2[j], lambda_st2, sigma2, d[j], InvTau2_1);
                }
                else
                    tau_st2[j] = _tau_st2;
            }

            lambda_st2 = rgamma( P + gammaTau_st, tau_st2.Sum()/2+deltaTau_st );
            if ( lambda_st2 <= 0 || isnan(lambda_st2) )
            {
                _log_debug( _HI_, "PART F_test.....tau2*=%f, P=%f, gammaTau_st=%f, deltaTau_st=%f",
                             tau_st2.Sum(), P, gammaTau_st, deltaTau_st);
                _log_fatal( _HI_, "FAILED by NAN Problem");
            }
        }

//---------------------------------------------------
// part G in the R code
//  #--Updating residual variance: sigma2
//  tmp0   <- Y - X*alpha - Z*beta - Xi%*%t(a) - Zeta%*%t(d);
//  # scale parameter
//  tmp    <- t(tmp0)%*%tmp0/N;
//  sigma2 <- 1/rchisq(1,N)*N*tmp;
//---------------------------------------------------
        strcat(m_szTraceTag, "G");
        //**SEQTEST
        //_log_debug( _HI_, "PART G.....round=%d/%d", round, R);

        Y_eps = Y;
        for (int j=0; j<Covs.GetNumCols(); j++ )
            Y_eps = Y_eps - Covs.GetCol(j) * Coefs[j];
        tmp3 = Y_eps- ((*gen_pA)*a).GetCol(0) - ((*gen_pD)*d).GetCol(0);

        //**SEQTEST
        //_log_debug( _HI_, "PART G_test.....Y_eps=%f, a=%f, d=%f", Y_eps.Sum(), a.Sum(), d.Sum());

        sigma2 = 1/rchisq(N) * tmp3.Prod(tmp3);

//---------------------------------------------------
// part H in the R code
//
//  r_sigma2  <- rbind( r_sigma2, sigma2);
//  r_a       <- rbind( r_a, a);
//  r_d       <- rbind( r_d, d);
//  r_alpha   <- rbind( r_alpha, alpha);
//  r_beta    <- rbind( r_beta, beta);
//  r_tau2    <- rbind( r_tau2, tau2);
//  r_lambda2 <- rbind( r_lambda2, lambda2);
//  r_tau_st2 <- rbind( r_tau_st2, tau_st2);
//  r_lambda_st2 <- rbind( r_lambda_st2, r_lambda_st2);
//---------------------------------------------------
        strcat(m_szTraceTag, "H");
        //**SEQTEST
        //_log_debug( _HI_, "PART H.....round=%d/%d", round, R);

        //-- save results of current round
        CFmVector fmTmp(0, 0.0);
        if ( round >= m_pCfg->GetBurnInRound() )
        {
			fmTmp.Resize(0);
            fmTmp.Append( Coefs);
            fmTmp.Append( a );
            fmTmp.Append( d );
            int ret = r_pFileMat->Rbind( fmTmp );
            if (ret)
            {
                _log_error( _HI_, "PART H: Failed to CFmFileMatrix::Rbind(ret=%d)", ret);
                return(ret);
            }
        }

		//**SEQTEST
		//_log_debug(_HI_, "Round=%d, mu=%.4f, alpha=%.4f, lambda=%.4f, lambda*=%.4f, sigma2=%f\n", round, Coefs[0], Coefs[1], lambda2, lambda_st2, sigma2 );

		//-- save the trace tag into the log file.
        if (round % 100 == 0)
            pcf.UpdatePcfFile( PCF_GOING, -1, -1, round+1, R);
    }

    PutRNGstate();

    destroy( gen_pA );
    destroy( gen_pD );

    return(0);
}


/*Sample code in Java:
public double inverseGaussian(double mu, double lambda) {
       Random rand = new Random();
       double v = rand.nextGaussian();   // sample from a normal distribution with a mean of 0 and 1 standard deviation
       double y = v*v;
       double x = mu + (mu*mu*y)/(2*lambda) - (mu/(2*lambda)) * Math.sqrt(4*mu*lambda*y + mu*mu*y*y);
       double test = rand.nextDouble();  // sample from a uniform distribution between 0 and 1
       if (test <= (mu)/(mu + x))
              return x;
       else
              return (mu*mu)/x;
}*/

double BLS::func_invGau(double mu, double lambda)
{
    //squared normal, i.e., chi-square with df=1
    double y = pow(rnorm(0, 1), 2);
    double x    = mu + 0.5*mu/lambda * ( mu*y - sqrt(4*mu*lambda*y + mu*mu*y*y) );
    if (runif(0, 1) <= (mu/(mu+x)))
        return(x);
    else
        return (mu*mu)/x;
}

int GetCoeff( CFmVector* pY, CFmVector* pX, CFmVector* pCoeff)
{
    CFmVector b1(0, 0.0);
    b1 = ( (*pY) - pY->GetMean() )*( (*pX) - pX->GetMean() );
    CFmVector b0(0, 0.0);
    b0 = ( (*pX) - pX->GetMean() )*( (*pX) - pX->GetMean() );
    float b = 0;
    if (b0.Sum()!=0 )
        b = b1.Sum()/b0.Sum();

    pCoeff->Resize(0);
    pCoeff->Put( pY->GetMean() - b*pX->GetMean() );
    pCoeff->Put( b );

    return(0);
}

CFmMatrix* BLS::GetPrePcaAd()
{
    CFmRSource srcR;

    CFmVector vct(0, 0.0);
    CFmVector vct_ca(0, 0.0);
    CFmVector vct_cd(0, 0.0);
    CFmVector vct_ret(0, 0.0);

    for(int i=0; i<m_pDat->m_nSnpP;i++)
    {
        vct.Resize(0);
        m_pDat->m_pPackedSNP->GetSnpRow( i, &vct );

        vct = vct - 1;
        GetCoeff( m_pDat->m_pPhenoY, &vct, &vct_ret);
        vct_ca.Put( vct_ret[1] );

        vct = vct.abs()*(-1)+1;
        GetCoeff( m_pDat->m_pPhenoY, &vct, &vct_ret );
        vct_cd.Put( vct_ret[1] );
    }


    CFmVector vct_sa(0, 0.0);
    vct_sa = vct_ca.abs();
    vct_sa.Sort(true);
    int nMarker = (int)floor( 0.03*m_pDat->m_nSnpP )-1;
    if (nMarker<0) nMarker = 0;
    double fMarker_a = vct_sa[nMarker];

    CFmVector vct_sd(0, 0.0);
    vct_sd = vct_cd.abs();
    vct_sd.Sort(true);
    double fMarker_d = vct_sd[nMarker];

    _log_info( _HI_, "PREC2(ma=%.4f, md=%.4f)", fMarker_a, fMarker_d );

	CFmNewTemp refNew;
    CFmMatrix* pSnp = new (refNew) CFmMatrix(0, 0);
    for( int i=0; i<m_pDat->m_nSnpP; i++ )
    {
        if(fabs(vct_ca[i]) >= fMarker_a)
        {
            vct.Resize(0);
            m_pDat->m_pPackedSNP->GetSnpRow( i, &vct );
            vct = vct - 1;
            pSnp->Cbind( vct );
        }
    }

    for( int i=0; i<m_pDat->m_nSnpP; i++ )
    {

        if(fabs(vct_cd[i]) >= fMarker_d)
        {
            vct.Resize(0);
            m_pDat->m_pPackedSNP->GetSnpRow( i, &vct );
            vct = vct - 1;
            vct = vct.abs()*(-1)+1;
            pSnp->Cbind( vct );
        }
    }

    //Row:SNP Col:Individual
    //pSnp->Transpose();

    return pSnp;
}

int BLS::proc_prec( int nOrder )
{
    CFmMatrix* pPcaAd = GetPrePcaAd();
    if (pPcaAd==NULL)
        return(-1);

    char szCmd[ MAX_PATH ]={0};
    strcpy(szCmd, "proc_prec<-function( m ) 				\
    {\n														\
		if (m>length(pca[1,]))	m<- length(pca[1,]);\n 		\
		reduced_pred_pca <- prcomp(pca, scale = TRUE);\n 	\
		reduced_pca1 <- predict(reduced_pred_pca)[,c(1:m)];\n \
		prec_y <<- lm(py~reduced_pca1)$fitted.values;\n 	\
	}\n" );

    sprintf(szCmd, "%s proc_prec(%d)\n", szCmd, nOrder );

    CFmRSource srcR;
    srcR.SetGlobalVariable( "pca", pPcaAd );
    srcR.SetGlobalVariable( "py",  m_pDat->m_pPhenoY );

    int ret = srcR.Run1(szCmd );
    if ( ret !=0 )
        _log_prompt(_HI_, "!!proc_prec, Failed to save result to %s", szCmd );

    CFmVector y_vect(0, 0.0);
    ret = srcR.GetGlobalVector("prec_y", &y_vect );
    if ( ret==0 )
        *(m_pDat->m_pPhenoY) =  y_vect;
    else
        _log_prompt(_HI_, "!!proc_prec, Failed to load result from %s", szCmd );

    destroy( pPcaAd );

    return(ret);
}

int BLS::ExportFig(char* szFigFile, bool bRefit )
{
    char szRFile[ MAX_PATH ]={0};
    strcpy(szRFile,m_pCmd->szAppName);
    char* pos = strrchr( szRFile,  '/' );
    if (pos==NULL)
        pos = strrchr( szRFile,  '\\' );
    if (pos==NULL)
    {
        strcpy(szRFile, "hgwas.r");
        pos = szRFile;
    }
    else
        pos++;

    strcpy(pos, "hgwas.r");

    _log_debug( _HI_, "ExportFig, R File=%s", szRFile );

    // export manhattan figure 1
    char szFig1[MAX_PATH]={0};
    CFmSys::GetSiblingFile(szFig1, szFigFile, bRefit?2:1, TRUE );

    char szTmpFile[MAX_PATH]={0};
    sprintf( szTmpFile, "%s.%s", szFig1, "csv" );
    m_pRes->SaveCsv4Fig( szTmpFile, bRefit );

    char szCmd[MAX_PATH ]={0};
    sprintf(szCmd, "plot_adh2(\"%s\",\"%s\")", szTmpFile, szFig1 );

    CFmRSource src2;
    int ret = src2.Run2(szRFile, szCmd );
    if ( ret ==0 )
        _log_prompt(_HI_, "ExportFig,  Save the Add/Dom/H2 figure(%s)", szFig1 );
    else
        _log_prompt(_HI_, "!!ExportFig, Failed to Export the Add/Dom/H2 figure(%s)", szFig1 );

    //unlink(szTmpFile);
    return(ret);
}



SEXP BLS::GetRObj()
{
	if (m_pRes)
		return(m_pRes->GetRObj());
	else
		return(R_NilValue);
}

void destroy(BLS* p)
{
	CFmNewTemp  fmRef;
	p->~BLS();
	operator delete(p, fmRef);
}
