/* gls_model.cpp  -	LS2 Computational model
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
#include "fm_rsource.h"
#include "fm_filematrix.h"
#include "fm_err.h"
#include "fm_pcf.h"
#include "fm_new.h"

#include "gls_model.h"
#include "gls_dat.h"
#include "gls_R.h"
#include "gls_par.h"
#include "gls_res.h"
#include "gls_cfg.h"

#define _t(x) ((x).GetTransposed())

GLS::GLS()
{
    m_nDataType = RUNMODE_UNK;
    m_bRefit = false;

	m_pDat = NULL;
	m_pRes = NULL;
    m_pCfg = NULL;
}

GLS::~GLS()
{
    if (m_pDat) destroy( m_pDat );
    if (m_pRes) destroy( m_pRes );
    //if (m_pCfg) destroy( m_pCfg );

    _log_info(_HI_, "GLS is released successfully.");
}

int GLS::LoadSimulate( CMDOPTIONS *pCmd, GLS_par* pPar  )
{
	int ret=0;

    m_pCmd = pCmd;
    m_nDataType = RUNMODE_SIM;

    pPar->Summary();

    CFmPcf pcf;
    pcf.UpdatePcfFile(PCF_PAR_READ );

    _log_prompt(_HI_, "Start the data generation.");

	CFmNewTemp refNew;
    m_pDat = new (refNew) GLS_dat(pCmd);

    pcf.UpdatePcfFile(PCF_DAT_LOAD );
    //Simulate_test??
    if ( (ret = m_pDat->Simulate(pPar, pCmd)) !=0)
    {
        _log_error(_HI_, "Failed to do data generation.");
        return(ret);
    }

    _log_prompt(_HI_, "Data generation is successful.");
    return(0);
}

int GLS::LoadPlink( CMDOPTIONS *pCmd )
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

    _log_prompt(_HI_, "Load the data sets from PLINK files.");

    CFmPcf pcf;
    pcf.UpdatePcfFile(PCF_DAT_LOAD );

	CFmNewTemp refNew;
    m_pDat = new (refNew) GLS_dat( pCmd );
    int ret = m_pDat->LoadPlink( pCmd->szTpedFile,
                                 pCmd->szTfamFile,
                                 pCmd->szPheFile,
                                 pCmd->bZNormalized);
    if (ret!=0)
    {
        _log_error(_HI_, "Failed to load the PLINK files.");
        return(ret);
    }


    _log_prompt(_HI_, "PLINK data are loaded successfully.");
    m_pDat->Summary();

    return(0);
}

int GLS::LoadSimple( CMDOPTIONS *pCmd )
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

	CFmNewTemp refNew;
    m_pDat = new (refNew) GLS_dat( pCmd );
    int ret = m_pDat->LoadSimple(pCmd->szSnpFile,
                                 pCmd->szPheFile,
                                 pCmd->bZNormalized);
    if (ret!=0)
    {
        _log_error(_HI_, "Failed to load the SIMPLE files.");
        return(ret);
    }

    _log_prompt(_HI_, "Simple data are loaded successfully.");

    m_pDat->Summary();
    return(0);
}

int GLS::LoadSnpmat( CFmMatrix* pFmPhe, CFmMatrix* pFmSnp, CMDOPTIONS *pCmd )
{
    m_pCmd = pCmd;
    m_nDataType = RUNMODE_SIMPLE;

    _log_info(_HI_, "Check SNP/phenotype parameter.");

    if ( pFmPhe->GetNumRows() != pFmSnp->GetNumCols()-2 )
        return( ERR_TPEDFILE_LONG );

    CFmPcf pcf;
    pcf.UpdatePcfFile(PCF_DAT_LOAD );

	CFmNewTemp refNew;
    m_pDat = new (refNew) GLS_dat( pCmd );
    int ret = m_pDat->AttachSnpmat( pFmPhe, pFmSnp, pCmd->bZNormalized);
    if (ret!=0)
    {
        _log_error(_HI_, "Failed to attach the SNP dataset.");
        return(ret);
    }

    _log_prompt(_HI_, "SNP data are loaded successfully.");

    m_pDat->Summary();
    return(0);
}

int GLS::Varsel(GLS_cfg *pCfg)
{
    // Create Result object
    if (m_pRes!=NULL)
		destroy( m_pRes );

	m_pCfg = pCfg;

    int P = m_pDat->m_nSnpP;
    int N = m_pDat->m_nSubjN;

    _log_prompt( _HI_, "VARSEL: SnpP:%d, SubjN:%d, nMcmcIter:%d", P, N, m_pCfg->m_nMcmcIter );

    m_bRefit = false;

	CFmNewTemp refNew;
    m_pRes = new (refNew) GLS_res(m_pCmd, m_pCfg);
    m_pRes->InitVarsel( m_pCmd->nSimuRound, P, N, m_pCfg->m_nMcmcIter);

    CFmPcf pcf;
    pcf.UpdatePcfFile( PCF_VARSEL, 0, 1, -1, -1);

	CFmFileMatrix* pMatRet=NULL;
	if (m_pCmd->szMatadFile)
		pMatRet = new (refNew) CFmFileMatrix( m_pCmd->szMatadFile, FALSE, FALSE);

    try
    {
        CFmSnpMat gen( P , m_pDat->m_nSubjN);   /* [P, N]*/
        m_pDat->GetPartialSNP( &gen, 0, P );

        int ret = proc_mcmc( *(m_pDat->m_pPhenoY), *(m_pDat->m_pPhenoZ), *(m_pDat->m_pPhenoZ0), *(m_pDat->m_pCovars), gen, pMatRet );
        if (ret!=0)
        {
        	_log_error( _HI_, "VARSEL: MCMC procedure return error code(%d).", ret);
        	return(ret);
        }

        CFmVectorStr vctSnpName(0);
        CFmVector vctSnpChr(0, 0.0);
        CFmVector vctSnpPos(0, 0.0);
        gen.GetSnpInfo(&vctSnpName, &vctSnpChr, &vctSnpPos);

	    m_pRes->SetMcmcResults( false, &vctSnpName, &vctSnpChr, &vctSnpPos, pMatRet, m_pDat->m_pCovars->GetNumCols() );

    }
    catch(const char* str)
    {
        pcf.UpdatePcfFile( PCF_EXCEPT, 1, 1,  -1, -1);
        if (pMatRet) destroy( pMatRet );

        //-- keep the trace tag into the log file.
        _log_error( _HI_, "VARSEL: Exception is caught in the procedure of MCMC.Tracetag=%s, Exception=%s", m_szTraceTag, str );
        return( ERR_EXCEPTION );
    }

    if (pMatRet) destroy( pMatRet );

    if (strlen( m_pCmd->szRdataFile) )
    {
        _log_prompt( _HI_, "VARSEL: Save the results into RDATA file(%s)", m_pCmd->szRdataFile );
        m_pRes->SaveRData( m_pCmd->szRdataFile);
    }

    if (strlen( m_pCmd->szFigOutFile)>0 )
    {
        _log_prompt( _HI_, "VARSEL: Save the figures into PDF file(%s)", m_pCmd->szFigOutFile );
        ExportFig(m_pCmd->szFigOutFile, false);
    }

    if (strlen( m_pCmd->szSigOutFile)>0 )
    {
        _log_prompt( _HI_, "VARSEL: Save the significant SNP into CSV file(sig.out=%s)", m_pCmd->szSigOutFile);

        char szSig1[MAX_PATH]={0};
        CFmSys::GetSiblingFile( szSig1, m_pCmd->szSigOutFile, 1 );
        m_pRes->SaveSigFile( szSig1, false );
    }

    _log_prompt( _HI_, "VARSEL: End");

    return 0;
}

int GLS::Refit(GLS_cfg *pCfg)
{
    if (m_pRes==NULL)
    {
        _log_error( _HI_, "REFIT: No result object for the refit");
        return( ERR_NULL_DATA );
    }

	m_pCfg = pCfg;

    int ret = m_pRes->InitRefit( m_pCfg->m_nMcmcIter );
    if ( ret!=0 )
    {
        _log_error( _HI_, "REFIT: Failed to select siginificant snp." );
        return( ret );
    }

    CFmVector* pSnpList = NULL;
    ret = m_pRes->GetRefitSnp( &pSnpList );
    if (ret)
    {
        _log_error( _HI_, "REFIT: Failed to get refitting SNPs" );
        return(ret);
    }

    if (pSnpList->GetLength()==0 || m_pDat->m_nSubjN==0)
    {
        _log_info( _HI_, "REFIT: No snps or subjects are selected to refit.");
        return( ERR_NO_ITEMS );
    }

    CFmSnpMat gen(pSnpList->GetLength(), m_pDat->m_nSubjN );
    ret = m_pDat->SelectRefitGenos(pSnpList, &gen, true);
    if ( ret !=0 )
    {
        _log_info( _HI_, "REFIT: Failed to extract SNP data info from the data object");
        return( ret );
    }

    /* Now gen is [N,P] */
    int P = pSnpList->GetLength();
    int N = m_pDat->m_nSubjN;

    _log_prompt( _HI_, "REFIT: SnpP:%d, SubjN:%d ", P, N );
    m_bRefit = true;

    CFmPcf pcf;
    pcf.UpdatePcfFile( PCF_REFIT, 0, 1, -1, -1);

	CFmNewTemp refNew;
	CFmFileMatrix* pMatRet=NULL;
	if (m_pCmd->szMatadFile)
		pMatRet = new (refNew) CFmFileMatrix( m_pCmd->szMatadFile, FALSE, FALSE);

    try
    {
        int ret = proc_mcmc( *(m_pDat->m_pPhenoY), *(m_pDat->m_pPhenoZ), *(m_pDat->m_pPhenoZ0), *(m_pDat->m_pCovars), gen, pMatRet );

        if (ret!=0)
        {
            _log_info( _HI_, "REFIT: Failed to do the MCMC procedure.");
            return(ret);
        }

        CFmVector vctSnpChr(0, 0.0);
        CFmVector vctSnpPos(0, 0.0);
        CFmVectorStr vctSnpName(0);
        gen.GetSnpInfo(&vctSnpName, &vctSnpChr, &vctSnpPos);

        m_pRes->SetMcmcResults( true, &vctSnpName, &vctSnpChr, &vctSnpPos, pMatRet, m_pDat->m_pCovars->GetNumCols());
    }
    catch(const char* str)
    {
        pcf.UpdatePcfFile( PCF_EXCEPT, 1, 1,  -1, -1);
        if (pMatRet) destroy( pMatRet );

        //-- keep the trace tag into the log file.
        _log_error( _HI_, "REFIT: An exception is caught in the MCMC procedure. Tracetag=%s, Exception=%s", m_szTraceTag, str);
        return( ERR_EXCEPTION );
    }

    if (pMatRet) destroy( pMatRet );

    _log_info( _HI_, "REFIT: Select Sig. SNPs.");

    m_pRes->GetSigSNP();
    m_pRes->Summary();

    if (strlen( m_pCmd->szRdataFile) )
    {
        _log_prompt( _HI_, "REFIT: Save the results to RDATA file(%s).", m_pCmd->szRdataFile);
        m_pRes->SaveRData( m_pCmd->szRdataFile);
    }

    if (strlen( m_pCmd->szFigOutFile)>0 )
    {
        _log_prompt( _HI_, "REFIT: Save the figures into PDF file(%s)", m_pCmd->szFigOutFile );
        ExportFig(m_pCmd->szFigOutFile, true);
    }

    if (strlen( m_pCmd->szSigOutFile)>0 )
    {
        _log_prompt( _HI_, "REFIT: Save the significant SNP into CSV file(sig.out=%s)", m_pCmd->szSigOutFile);

        char szSig1[MAX_PATH]={0};
        CFmSys::GetSiblingFile( szSig1, m_pCmd->szSigOutFile, 2 );
        m_pRes->SaveSigFile( szSig1, true );
    }

    _log_prompt( _HI_, "REFIT: End");

    return 0;
}

//**
int GLS::proc_mcmc( CFmMatrix& Y,  CFmMatrix& Z,  CFmMatrix& Z0, CFmMatrix& X, CFmSnpMat& gen, CFmFileMatrix* pMatRet)
{
    //!!!TEST for R
    int nSeed = m_pCfg->GetTempId();
    //int nSeed = 100;
    R_set_seed( nSeed );
    _log_prompt( _HI_, "set.seed(%d)",nSeed);

    CFmPcf pf;

    int N = Y.GetNumRows();   //*[N,Q]
    int P = gen.GetNumSnps(); //*[N,P]
    int Q = Y.GetNumCols();
    int R = m_pCfg->m_nMcmcIter;
    int nC = X.GetNumCols()-1;

//**SEQTEST Rprintf("DEBUG:N=%d, P=%d, Q=%d, R=%d\n ", N, P, Q, R);

    int sum_Ji = Y.GetNumRows()*Y.GetNumCols()-Y.GetNanCount();

    double gammaTau = 1;
    double deltaTau = 1;
    double gammaTau_st = 1;
    double deltaTau_st = 1;
    double sigma2 = 8;
    double rho=0.4;

    double lambda2 = 1;
    double lambda2_x = 1;
    double lambda_st2 = 1;
    double lambda_st2_x = 1;
    CFmVector tau2(P, 1);
    CFmVector tau2_x(P, 1);
    CFmVector tau_st2(P, 1);
    CFmVector tau_st2_x(P, 1);

    CFmMatrix sigmaMu( LG, true, 1);
    CFmMatrix sigmaAlpha(LG, true, 1);
    CFmMatrix alpha( nC+1, LG);
    CFmMatrix a(P,LG);
    CFmMatrix d(P,LG);

    CFmMatrix r_sigma2(0, 0, 1, 1);
    CFmMatrix r_tau2(0, 0, P, 1);
    CFmMatrix r_lambda2(0, 0, 1, 1);
    CFmMatrix r_tau_st2(0, 0, P, 1);
    CFmMatrix r_lambda_st2(0, 0, 1, 1);

	CFmNewTemp refNew;
    CFmMatrix** all_corMat = Calloc(N, CFmMatrix* );
    for (int i=0; i<N; i++)
        all_corMat[i] = new (refNew) CFmMatrix(0, 0, Q, Q);

    CFmMatrix** all_corMat_MH = Calloc(N, CFmMatrix* );
    for (int i=0; i<N; i++)
        all_corMat_MH[i] = new (refNew)  CFmMatrix(0, 0, Q, Q);

    CFmVector** all_yi = Calloc(N, CFmVector* );
    for (int i=0; i<N; i++)
        all_yi[i] = new  (refNew) CFmVector(Q, 0);

    CFmVector** all_rd = Calloc(N, CFmVector* );
    for (int i=0; i<N; i++)
        all_rd[i] = new  (refNew) CFmVector(Q, 0);

    CFmMatrix** all_ui = Calloc(N, CFmMatrix* );
    for (int i=0; i<N; i++)
        all_ui[i] = new  (refNew) CFmMatrix( 0, 0, Q, LG);

//---------------------------------------------------
// PART A in the R code
//---------------------------------------------------
strcpy(m_szTraceTag, "A");
//**SEQTEST _log_debug( _HI_, "PART A.....%d", R);

    _log_debug( _HI_, "PART A.....Z: bZNorm=%d", m_pCmd->bZNormalized?1:0);

    // -- create Ui matrix, conditioning on standarized age
    CFmVector Z_tmp(1, 0, Q );
    CFmVector tp(1,0,Q );
    for( int i=0; i<N; i++)
    {
        tp = Z.GetRow(i);
        tp.RemoveNan();

        CFmVector tp1( tp.GetLength(), 1);
        (all_ui[i])->Cbind( tp1 );
        (all_ui[i])->Cbind( tp );
        (all_ui[i])->Cbind( ((tp^2)*3.0 - 1.0)/2.0 );
		if (LG>=4) (all_ui[i])->Cbind( ((tp^3)*5.0 - tp*3.0)/2.0 );
		if (LG>=5) (all_ui[i])->Cbind( ((tp^4)*35.0 - (tp^2)*30.0 +3 )/8.0 );
    }

    // -- initialize the residual Y
    CFmVector yi0(1,1);
    for (int i=0; i<N; i++)
    {
        yi0 = Y.GetRow(i);
        yi0.RemoveNan();
        (*all_yi[i]) = yi0;
        (*all_rd[i]) = yi0;
    }

    CFmMatrix sigma0( 0, 0, Q, Q );
    CFmMatrix sigma00( 0, 0, Q, Q );

    CFmMatrix Mu_Var_j(0, 0, Q, Q);
    CFmMatrix Mu_Mu_j(0, 0, Q, Q);
    CFmMatrix a_old( &a );
    CFmMatrix d_old( &d );
    CFmVector mu(LG,0);
    CFmMatrix aVar_j(0, 0, Q, Q);
    CFmMatrix aMu_j(0, 0, Q, Q);
    CFmMatrix dVar_j(0, 0, Q, Q);
    CFmMatrix dMu_j(0, 0, Q, Q);
    CFmMatrix alpha_Var_j(0, 0, Q, Q);
    CFmMatrix alpha_Mu_j(0, 0, Q, Q);

    CFmVector allRd_i(Q, 0.0);
    CFmMatrix sigma_old(Q,Q);
    CFmMatrix sigma_new(Q,Q);
    CFmMatrix exp_mat(Q,Q);
    CFmVector tmpv_sigma2(Q,Q);
    CFmMatrix tmpv_get2(Q,Q);
    CFmMatrix tmp( 0, 0, Q, Q);
    CFmVector mean_effect(0, 0.0, Q );
    CFmMatrix diag( LG, true, 1);

    CFmMatrix tmp2(0, 0, Q, Q);
    CFmMatrix tmp3(0, 0, Q, Q);
    CFmMatrix tmp4(0, 0, Q, Q);


    GetRNGstate();
    for (int round=0; round< R ; round++)
    {
        if ( round % m_pCfg->m_nMcmcHint ==0 || round<=0 )
        {
            PCFDAT* pcf = pf.GetPcfBlock();
            char szBuf[256]={0};
            sprintf(szBuf, "Section:%d Round:%d Progress:%.3f%% Total:%.0fs Left:%.0fs\n",
                    m_pCfg->m_nSectId, round,
                    pcf->fProgress*100,
                    pcf->fELapsedSeconds+pcf->fEstimatedSeconds,
                    pcf->fEstimatedSeconds );
            _log_info( _HI_, szBuf );
        }

//---------------------------------------------------
// part B in the R code
//---------------------------------------------------
strcpy(m_szTraceTag, "B");
//**SEQTEST _log_debug( _HI_, "PART B.....round=%d/%d", round, R);

        //-- new proposal of rho
        double tmp5 = runif( fmax2(-1, rho - m_pCfg->m_fRhoTuning), fmin2(1, rho + m_pCfg->m_fRhoTuning) );

        //-- calculating correlation matrix
        for(int i=0; i<N; i++)
        {
            int m = Z.GetRow(i).GetLengthNonNan();
            sigma0.Resize(m, m, true);
            sigma00.Resize(m, m, true);
            for (int ii=0; ii<m; ii++)
                for (int jj=ii; jj<m; jj++)
                {
                    if (ii==jj)
                    {
                        sigma0.Set(ii, jj, 1);
                        sigma00.Set(ii, jj, 1);
                    }
                    else
                    {
                        sigma0.Set(ii,jj,    pow(rho, Z0.Get(i,jj) - Z0.Get(i,ii) ) );
                        sigma0.Set(jj, ii,   sigma0.Get(ii, jj) );
                        sigma00.Set(ii,jj,   pow(tmp5, Z0.Get(i,jj)- Z0.Get(i,ii) ) );
                        sigma00.Set(jj, ii,  sigma00.Get(ii, jj) );
                    }
                }

            (*all_corMat[i]) =  sigma0;
            (*all_corMat_MH[i]) = sigma00;
        }

//---------------------------------------------------
// part C in the R code
//---------------------------------------------------
strcat(m_szTraceTag, "C");
//**SEQTEST _log_debug( _HI_, "PART C.....round=%d/%d", round, R);

        // ---- update mu
        tmp2 = 0;
        tmp3 = 0;

        for (int i=0; i<N; i++)
        {
            tmp =  ((*all_corMat[i]) * sigma2).GetInverted() * (*all_ui[i]);

			tmp4.Resize(1, all_rd[i]->GetLength(), TRUE );
			if(nC>0)
			for( int nX=0; nX < nC; nX++ )
				tmp4 = tmp4 + alpha.GetRow(nX) * _t( (*all_ui[i]) ) * X.Get(i,nX+1);

           	tmp3 = tmp3 + ( (*all_rd[i]) - tmp4.GetRow(0) ) * tmp;
            tmp2 = tmp2 + _t( (*all_ui[i]) ) * tmp;
        }

        Mu_Var_j = ((sigmaMu.GetInverted()) + tmp2).GetInverted();
        Mu_Mu_j = Mu_Var_j * _t( tmp3 );

        double fBuf[LG];
        int ret = rmultnorm( 1, Mu_Mu_j.GetData(), Mu_Var_j.GetData(), LG, fBuf);
        mu = (double*)fBuf;

//---------------------------------------------------
// part E in the R code
//---------------------------------------------------
strcat(m_szTraceTag, "E");
//**SEQTEST _log_debug( _HI_, "PART E.....round=%d/%d", round, R);

        // ---- Calculating effect other than genetic effects: yi-\mu_i after updating a,d

        for (int i=0;i<N;i++)
        {
            mean_effect.Resize((all_yi[i])->GetLength(), true);
            mean_effect = 0.0;

            for (int jj=0; jj<P; jj++)
            {

                if (gen.Get_a(jj,i) != 0 )
                    mean_effect = mean_effect + ( (*all_ui[i])*a.GetRow(jj) ).GetCol(0)*gen.Get_a(jj,i);

                if (gen.Get_d(jj,i) != 0 )
                    mean_effect = mean_effect + ( (*all_ui[i])*d.GetRow(jj) ).GetCol(0)*gen.Get_d(jj,i);

            }

            (*all_rd[i]) = (*all_yi[i]) - mean_effect;
        }

        a_old = a;
        d_old = d;

//---------------------------------------------------
// part F in the R code
//---------------------------------------------------
strcat(m_szTraceTag, "F");
//**SEQTEST _log_debug( _HI_, "PART F.....round=%d/%d", round, R);

		CFmVector vctP(P, 0.0);

        // -- updating additive effects: a for each SNP
        //CFmMatrix diag(4, true, 1);
        if(m_pCmd->bAddUsed)
        for(int j=0; j<P; j++)
        {
            tmp2 = 0;
            tmp3 = 0;

            int N0 = 0;
            int N2 = 0;
            int N1 = 0;

            for(int i=0; i<N; i++)
            {
                double a0 = gen.Get_a(j,i);
                if (a0!=0)
                {
					tmp4.Resize(1, all_rd[i]->GetLength(), TRUE );
					if(nC>0)
					for( int nX=0; nX < nC; nX++ )
						tmp4 = tmp4 + (alpha.GetRow(nX) * _t( (*all_ui[i]) ) ) * X.Get(i,nX+1);

                    tmp = ( ((*all_corMat[i]) * sigma2).GetInverted() ) * a0 * (*all_ui[i]);
                    tmp3 = tmp3 + ( (*all_rd[i]) - (mu * _t( (*all_ui[i]) ) ).GetRow(0)
                                  - tmp4.GetRow(0)
                                  + ( a.GetRow(j) * _t( (*all_ui[i]) ) ).GetRow(0) * a0 ) * tmp;
                    tmp2 = tmp2 + _t( (*all_ui[i]) )* tmp * a0;
                }


                if( a0 >0 ) N0++;
                if( a0 <0 ) N2++;
				if( a0 == 0 ) N1++;
            }

            if (N0==0 || N2==0)
            {
                double fBuf[LG] = { 0.0 };
                a.SetRow( j,  (double*)fBuf);
                vctP[j] = 0;
            }
            else
            {
                if (N1==0) vctP[j] = 1; else vctP[j] = 2;

                if (!m_bRefit)
                {
                    if (N1>0)
						diag.Square(LG, true, 1/tau2[j]);
                    else
						diag.Square(LG, true, 1/tau2_x[j]);

                    aVar_j =  (diag/sigma2 + tmp2).GetInverted();
                }
                else
                {
                    diag.Square(LG, true, 1);
	                //aVar_j =  (diag/sigma2 + tmp2).GetInverted();
	                aVar_j =  (tmp2).GetInverted();
                }

                aMu_j = aVar_j* _t( tmp3 );
                rmultnorm(1, aMu_j.GetData(), aVar_j.GetData(), LG, fBuf);
                a.SetRow( j,  (double*)fBuf);

                // -- update all_rd[i] for the newly updated genetic effect.
                for (int i=0; i<N; i++)
                    if (gen.Get_a(j,i)!=0)
                        (*all_rd[i]) = (*all_rd[i]) + _t( (*all_ui[i]) * a_old.GetRow(j) - (*all_ui[i]) * a.GetRow(j) ).GetRow(0) * gen.Get_a(j,i);
           }

        }

//---------------------------------------------------
// part G in the R code
//---------------------------------------------------
//CFmMatrix::_DM=0;
strcat(m_szTraceTag, "G");
//**SEQTEST _log_debug( _HI_, "PART G.....round=%d/%d", round, R);
    if(!m_bRefit && m_pCmd->bAddUsed )
    {
		int a_P1 =0 ;
		int a_P2 =0 ;
        // --- Updating tau2 underlying a
        for (int j=0; j<P; j++)
        {
			//if (a.RowProd(j,j )!=0)
			if (vctP[j]==2)
			{
				double InvTau2_1 = sqrt( LG * lambda2 * sigma2/a.RowProd(j,j ) );
				tau2[j] = 1/func_invGau( InvTau2_1, LG*lambda2);
				tau2_x[j] = 0;
				a_P2 ++;
			}
			else if (vctP[j]==1)
			{
				double InvTau2_1 = sqrt( LG * lambda2_x * sigma2/a.RowProd(j,j ) );
				tau2_x[j] = 1/func_invGau( InvTau2_1, LG*lambda2_x);
				tau2[j] = 0;
				a_P1 ++;
			}
			else
			{
				tau2[j] = 0;
				tau2_x[j] = 0;
			}
        }

        //-- Updating lambda2 underlying tau2
        lambda2   = rgamma( (LG*a_P2 + a_P2)/2 + gammaTau, (LG/2.0 * tau2.Sum() + deltaTau) );
        lambda2_x = rgamma( (LG*a_P1 + a_P1)/2 + gammaTau, (LG/2.0 * tau2_x.Sum() + deltaTau) );
    }
    else
    {
		lambda2 = 0;
		lambda2_x = 0;
	}

//---------------------------------------------------
// part H in the R code
//---------------------------------------------------
strcat(m_szTraceTag, "H");
//**SEQTEST _log_debug( _HI_, "PART H.....round=%d/%d", round, R);

		int d_P = 0;
        //-- Updating dominant effects: d, for each SNP
        if(m_pCmd->bDomUsed)
        for(int j=0; j<P; j++)
        {
            tmp2 = 0;
            tmp3 = 0;

            int N0 = 0;
            int N2 = 0;
            int N1 = 0;

            for(int i=0; i<N; i++)
            {
                double d0 = gen.Get_d(j,i);
                double a0 = gen.Get_a(j,i);
                if (d0!=0)
                {
					tmp4.Resize(1, all_rd[i]->GetLength(), TRUE );
					if(nC>0)
					for( int nX=0; nX < nC; nX++ )
						tmp4 = tmp4 + ( alpha.GetRow(nX) * _t( (*all_ui[i]) ) ) * X.Get(i,nX+1);

                    //calculate mean and variance
                    tmp  = ( (*all_corMat[i]) * sigma2).GetInverted() * (*all_ui[i])*d0;
                    tmp3 = tmp3 + ( (*all_rd[i]) - (mu * _t( (*all_ui[i]) ) ).GetRow(0)
                                - tmp4.GetRow(0)
                                + ( d.GetRow(j) *_t( (*all_ui[i]) ) ).GetRow(0)*d0 ) * tmp;
                    tmp2 = tmp2 + _t( (*all_ui[i]) ) * tmp * d0;

                }

                if( a0 >0 ) N0++;
                if( a0 <0 ) N2++;
                if( a0 ==0 ) N1++;
            }

            if ( N1==0 )
            {
                double fBuf[LG] = { 0.0 };
                d.SetRow( j,  (double*)fBuf);
                vctP[j] = 0;
            }
            else
            {
                if (N0>0 && N2>0) vctP[j] = 2; else vctP[j] = 1;
                if (!m_bRefit)
                {
                    if (vctP[j] == 2)
	                {
						diag.Square(LG, true, 1/tau_st2[j] );
	                	dVar_j = ( tmp2 + diag /sigma2 ).GetInverted();
	                }
	                else
	                {
						diag.Square(LG, true, 1/tau_st2_x[j] );
	                	dVar_j = ( tmp2 + diag /sigma2 ).GetInverted();
	                }
                }
                else
                {
	                diag.Square(LG, true, 1 );
	                //dVar_j = ( diag/sigma2 + tmp2 ).GetInverted();
	                dVar_j = ( tmp2 ).GetInverted();
                }

                dMu_j  = dVar_j * _t(tmp3);
                ret = rmultnorm( 1, dMu_j.GetData(), dVar_j.GetData(), LG, fBuf);
                d.SetRow(j, fBuf);

                // Update all_rd{i} for the newly updated genetic effect
                for (int i=0; i<N; i++)
                    if (gen.Get_d(j,i)!=0)
                        (*all_rd[i]) = (*all_rd[i]) + _t( (*all_ui[i]) * d_old.GetRow(j) - (*all_ui[i]) * d.GetRow(j) ).GetRow(0)* gen.Get_d(j,i);
            }
       }

//---------------------------------------------------
// part I in the R code
//---------------------------------------------------
strcat(m_szTraceTag, "I");
//**SEQTEST _log_debug( _HI_, "PART I.....round=%d/%d", round, R);
    if(!m_bRefit && m_pCmd->bDomUsed)
    {
		int d_P1 =0 ;
		int d_P2 =0 ;
        // -- Updating tau_st2 underlying d
        for (int j=0; j<P; j++)
        {
			//if (d.RowProd(j,j )!=0)
			if (vctP[j]==2)
			{
				double InvTau2_1 = sqrt( LG * lambda_st2 * sigma2/d.RowProd(j,j ) );
				tau_st2[j] = 1/func_invGau(InvTau2_1, LG*lambda_st2 );
				tau_st2_x[j] = 0;
				d_P2 ++;
			}
			else if (vctP[j]==1)
			{
				double InvTau2_1 = sqrt( LG * lambda_st2_x * sigma2/d.RowProd(j,j ) );
				tau_st2_x[j] = 1/func_invGau(InvTau2_1, LG*lambda_st2_x );
				tau_st2[j] = 0;
				d_P1 ++;
			}
			else
			{
				tau_st2[j] = 0;
				tau_st2_x[j] = 0;
			}
        }


        // -- Updating lambda_st2 underlying tau_star2
        lambda_st2 = rgamma( (LG * d_P2 + d_P2)/2+gammaTau_st, (LG/2.0*tau_st2.Sum() + deltaTau_st) );
        lambda_st2_x = rgamma( (LG * d_P1 + d_P1)/2+gammaTau_st, (LG/2.0*tau_st2_x.Sum() + deltaTau_st) );
    }
    else
    {
		lambda_st2 = 0;
        lambda_st2_x = 0;
	}

//---------------------------------------------------
// part J in the R code
//---------------------------------------------------
strcat(m_szTraceTag, "J");
//**SEQTEST _log_debug( _HI_, "PART J.....round=%d/%d", round, R);

        // -- Updating single covariates: alpha
        for(int nX=0; nX<nC; nX ++)
        {
			tmp2 = 0;
			tmp3 = 0;
			for (int i=0; i<N; i++)
			{
				tmp4.Resize(1, all_rd[i]->GetLength(), TRUE );
				if ( nC > 0)
				for( int nX2=0; nX2 < nC; nX2++ )
					if ( nX2 != nX )
						tmp4 = tmp4 + ( alpha.GetRow(nX2) * _t( (*all_ui[i]) ) ) * X.Get(i,nX2+1);

				// calculate mean and variance
				tmp  = ( (*all_corMat[i]) * sigma2 ).GetInverted() * X.Get(i, nX + 1) * (*all_ui[i]);
				tmp3 = tmp3 + ( (*all_rd[i]) - (mu * _t( (*all_ui[i]) ) ).GetRow(0) - tmp4.GetRow(0)) * tmp;
				tmp2 = tmp2 + _t( (*all_ui[i]) ) * tmp * X.Get(i,nX+1);
			}

			alpha_Var_j = ( tmp2 + sigmaAlpha.GetInverted()).GetInverted();
			alpha_Mu_j  = alpha_Var_j * _t(tmp3);
			rmultnorm( 1, alpha_Mu_j.GetData(), alpha_Var_j.GetData(), LG, fBuf );
			alpha.SetRow(nX, fBuf );
		}

//---------------------------------------------------
// part K in the R code
//---------------------------------------------------
strcat(m_szTraceTag, "K");
//**SEQTEST _log_debug( _HI_, "PART K.....round=%d/%d", round, R);

        //-- Updating residual variance: sigma2
        double sigma2_scale = 0;
        tmpv_sigma2.Resize(Q, true);
        tmpv_get2.Resize( Q, Q, true);

        for(int i=0; i<N; i++)
        {
			tmp4.Resize(1, all_rd[i]->GetLength(), TRUE );
			if (nC > 0)
			for( int nX=0; nX < nC; nX++ )
				tmp4 = tmp4 + ( alpha.GetRow( nX ) * _t( (*all_ui[i]) ) ) * X.Get(i,nX+1);

            // calculate scale parameter
            tmpv_sigma2 = (*all_rd[i]) - ( mu* _t( (*all_ui[i]) ) ).GetRow(0) - tmp4.GetRow(0);
            tmpv_get2 = (tmpv_sigma2 * ( (*all_corMat[i]).GetInverted() ) * _t( tmpv_sigma2) );
            sigma2_scale = sigma2_scale + tmpv_get2.Get(0,0) ;
        }

        sigma2_scale = sigma2_scale/(sum_Ji*1.0);
        double f_rchisq=rchisq( sum_Ji*1.0 );
        sigma2 = 1/f_rchisq*sum_Ji*1.0*sigma2_scale;

//---------------------------------------------------
// part L in the R code
//---------------------------------------------------
strcat(m_szTraceTag, "L");
//**SEQTEST _log_debug( _HI_, "PART L.....round=%d/%d", round, R);

        //-- Updating rho, proposal distribution Q
        double Qnew = dunif( tmp5, fmax2(-1, rho  - m_pCfg->m_fRhoTuning), fmin2(1, rho  + m_pCfg->m_fRhoTuning), FALSE );
        double Qold = dunif( rho,  fmax2(-1, tmp5 - m_pCfg->m_fRhoTuning), fmin2(1, tmp5 + m_pCfg->m_fRhoTuning), FALSE );

        // target density with new sample
        double det_diff = 1;
        double exp_diff = 0;
        allRd_i.Resize(Q, true);
        sigma_old.Resize( Q, Q, true);
        sigma_new.Resize( Q, Q, true);
        exp_mat.Resize( Q, Q, true);

        for (int i=0; i<N; i++)
        {
            sigma_old = (*all_corMat[i])*sigma2;
            sigma_new = (*all_corMat_MH[i])*sigma2;

			tmp4.Resize(1, all_rd[i]->GetLength(), TRUE );
			if ( nC > 0 )
			for( int nX=0; nX < nC; nX++ )
				tmp4 = tmp4 + ( alpha.GetRow( nX) * _t( (*all_ui[i]) ) ) * X.Get(i,nX+1);

            // calculate scale parameter
            allRd_i  = (*all_rd[i]) - (mu * _t((*all_ui[i]))).GetRow(0) - tmp4.GetRow(0);
            exp_mat  = allRd_i * sigma_new.GetInverted() * _t( allRd_i );
            exp_diff = exp_diff + exp_mat.Get(0,0);
            exp_mat  = allRd_i * sigma_old.GetInverted() * _t( allRd_i );
            exp_diff = exp_diff - exp_mat.Get(0,0);
            det_diff = det_diff * sqrt( sigma_old.GetDet()/sigma_new.GetDet() );
        }

//---------------------------------------------------
// part M in the R code
//---------------------------------------------------
strcat(m_szTraceTag, "M");
//**SEQTEST _log_debug( _HI_, "PART M.....round=%d/%d", round, R);

        double Pnew_Pold = det_diff * exp( -exp_diff/2.0 );
        double accpp     = fmin2( 1 , Pnew_Pold * Qold/ Qnew );
        double tmp6      = runif(0, 1);
        rho       		 = tmp5*(tmp6<accpp) + rho*(tmp6>=accpp);

        //-- save results of current round
        if ( round >= m_pCfg->GetBurnInRound() )
        {
			CFmVector fmVct(0, 0.0);
			CFmVector fmTmp0(0, 0.0);

			fmVct.Append(mu);

			for(int nX=0; nX<nC; nX++)
			{
				fmTmp0 = alpha.GetRow(nX);
				fmVct.Append(fmTmp0);
			}

			for(int i=0; i<a.GetNumRows(); i++)
			{
				fmTmp0 = a.GetRow(i);
				fmVct.Append( fmTmp0 );

				fmTmp0 = d.GetRow(i);
				fmVct.Append( fmTmp0 );
			}

			pMatRet->Rbind(fmVct);
        }

        _log_debug(_HI_, "Round=%d sigma2=%f rho=%.4f tmp5=%f lambda2=%f,%f,%f,%f mu=%.4f, %.4f, %.4f ",  round, sigma2, rho, tmp5, lambda2, lambda2_x, lambda_st2, lambda_st2_x, mu[0], mu[1], mu[2] );

        //-- save the trace tag into the log file.
        if ( round%50 == 0 )
            pf.UpdatePcfFile( PCF_GOING, -1, -1, round+1, R);
    }

    for (int i=0; i<N; i++) destroy( all_corMat[i]);   Free( all_corMat );
    for (int i=0; i<N; i++) destroy( all_corMat_MH[i]); Free( all_corMat_MH );
    for (int i=0; i<N; i++) destroy( all_yi[i]);  Free( all_yi );
    for (int i=0; i<N; i++) destroy( all_rd[i]);  Free( all_rd );
    for (int i=0; i<N; i++) destroy( all_ui[i]);  Free( all_ui );

    PutRNGstate();

    return(0);
}

double GLS::func_invGau(double theta, double chi)
{
    //squared normal, i.e., chi-square with df=1
    double chisq1 = pow(rnorm(0, 1), 2);

	double y1    = theta + 0.5*theta/chi * ( theta*chisq1 - sqrt(4*theta*chi*chisq1 + theta*theta*chisq1*chisq1) );
	double y2    = theta*theta/y1;
	double out_1 = runif(0, 1) < (theta/(theta+y1));
	double value = out_1*y1+(1-out_1)*y2;

    return(value);
}

int GLS::ExportFig(char* szFigFile, bool bRefit )
{
    char szRFile[ MAX_PATH ]={0};
    strcpy( szRFile,m_pCmd->szAppName );
    char* pos = strrchr( szRFile,  '/' );
    if (pos==NULL)
        pos = strrchr( szRFile,  '\\' );
    if (pos==NULL)
    {
        strcpy(szRFile, "lgwas.r");
        pos = szRFile;
    }
    else
        pos++;

    strcpy(pos, "lgwas.r");

    _log_debug( _HI_, "ExportFig, R File=%s", szRFile );

    // export manhattan figure 1
    char szFig1[MAX_PATH]={0};
    CFmSys::GetSiblingFile(szFig1, szFigFile, bRefit?2:1, TRUE );

    char szTmpFile[MAX_PATH]={0};
    sprintf( szTmpFile, "%s.%s", szFig1, "csv" );
    m_pRes->SaveCsv4Fig( szTmpFile, bRefit );

    char szCmd[MAX_PATH ]={0};
    if (!bRefit)
        sprintf(szCmd, "plot_adh2(\"%s\",\"%s\")", szTmpFile, szFig1 );
    else
        sprintf(szCmd, "plot_sig_curve( \"%s\",\"%s\", %d )", szTmpFile, szFig1, LG );

    CFmRSource src2;
    int ret = src2.Run2(szRFile, szCmd );
    if ( ret ==0 )
        _log_prompt( _HI_, "ExportFig,  Save the Add/Dom figure(%s)", szFig1 );
    else
        _log_prompt( _HI_, "!!ExportFig, Failed to Export the Add/Dom figure(%s)", szFig1 );

    //unlink(szTmpFile);
    return(ret);
}

SEXP GLS::GetRObj()
{
	if (m_pRes)
		return(m_pRes->GetRObj());
	else
		return(R_NilValue);
}
