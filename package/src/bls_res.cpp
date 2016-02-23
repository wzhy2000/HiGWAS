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
#include "fm_rls.h"
#include "fm_filematrix.h"
#include "fm_new.h"
#include "fm_err.h"

#include "bls_cfg.h"
#include "bls_res.h"

#define _t(x) ((x).GetTransposed())

BLS_res::BLS_res(CMDOPTIONS *pCmd, BLS_cfg* pCfg)
{
    m_pCfg = pCfg;
	m_pCmd = pCmd;

    m_nSimuRound = 1;
    m_nSnpP 	 = 0;
    m_nSubjN 	 = 0;
    m_nMcmcSect  = 0;
    m_nMcmcSnp   = 0;
    m_nMcmcIter  = 0;
    m_nRefitSnp  = 0;

    m_pVarsel_SnpName = NULL;
    m_pVarsel_SnpChr  = NULL;
    m_pVarsel_SnpPos  = NULL;
    m_pVarsel_Coefs   = NULL;
    m_pVarsel_Ra      = NULL;
    m_pVarsel_Rd      = NULL;
    m_pVarsel_Rh2     = NULL;
    m_pVarsel_QBest   = NULL;

    m_pRefit_SnpName  = NULL;
    m_pRefit_SnpChr   = NULL;
    m_pRefit_SnpPos   = NULL;
    m_pRefit_Coefs    = NULL;
    m_pRefit_Ra       = NULL;
    m_pRefit_Rd       = NULL;
    m_pRefit_Rh2      = NULL;
    m_pRefit_QBest    = NULL;

	m_pRefitSnps  = NULL;
	m_pSigSnps    = NULL;
    m_pSigAddSnps = NULL;
    m_pSigDomSnps = NULL;
    m_pCovNames    = NULL;
}

BLS_res::~BLS_res()
{
    if(m_pVarsel_SnpName) destroy( m_pVarsel_SnpName );
    if(m_pVarsel_SnpChr) destroy( m_pVarsel_SnpChr );
    if(m_pVarsel_SnpPos) destroy( m_pVarsel_SnpPos );
    if(m_pVarsel_Coefs) destroy( m_pVarsel_Coefs );
    if(m_pVarsel_Ra)   destroy( m_pVarsel_Ra );
    if(m_pVarsel_Rd)   destroy( m_pVarsel_Rd );
	if(m_pVarsel_Rh2)   destroy( m_pVarsel_Rh2 );
	if(m_pVarsel_QBest) destroy( m_pVarsel_QBest );

    if(m_pRefit_SnpName) destroy( m_pRefit_SnpName );
    if(m_pRefit_SnpChr) destroy( m_pRefit_SnpChr );
    if(m_pRefit_SnpPos) destroy( m_pRefit_SnpPos );
    if(m_pRefit_Coefs) destroy( m_pRefit_Coefs );
    if(m_pRefit_Ra) destroy( m_pRefit_Ra );
    if(m_pRefit_Rd) destroy( m_pRefit_Rd );
    if(m_pRefit_Rh2) destroy( m_pRefit_Rh2 );
	if(m_pRefit_QBest) destroy( m_pRefit_QBest );

	if(m_pRefitSnps) destroy( m_pRefitSnps );
	if(m_pSigSnps) destroy( m_pSigSnps );
    if(m_pSigAddSnps) destroy( m_pSigAddSnps );
    if(m_pSigDomSnps) destroy( m_pSigDomSnps );
    if(m_pCovNames) destroy( m_pCovNames );

    _log_debug(_HI_, "BLS_res is released successfully.");
}


// mRound:   Simulation Round,
// nSnpP,    SNP count,
// nSujN,    Individual count,
// mcmcSect  MCMC section count,
// mcmcSnps  SNP count for every MCMC procedure,
// mcmcIter  MCMC maximum iteration
int BLS_res::InitVarsel(int mRound, int nSnpP, int SubjN, int mcmcSect, int mcmcSnps, int mcmcIter, CFmVectorStr* pCovNames)
{
    m_nSimuRound = mRound;
    m_nSnpP      = nSnpP;
    m_nSubjN     = SubjN;
    m_nMcmcSect  = mcmcSect;
    m_nMcmcSnp   = mcmcSnps;
    m_nMcmcIter  = mcmcIter;

	CFmNewTemp refNew;
	m_pCovNames  = new (refNew) CFmVectorStr(pCovNames);

    return(0);
}

char* BLS_res::GetVarselSnpName(int idx)
{
    if (idx<0 || idx>=m_pVarsel_SnpName->GetLength())
        return(NULL);

    return (m_pVarsel_SnpName->Get(idx));
}

int BLS_res::GetRefitSnp(CFmVector* pVct)
{
	if (m_pVarsel_Ra==NULL)
        return(ERR_NULL_DATA);

	for (long int i=0; i< m_nSnpP; i++)
    {
	    if ( m_pVarsel_Ra->Get(i,0) || m_pVarsel_Rd->Get(i,0) )
			pVct->Put(i);
    }

    return(0);
}


int BLS_res::GetMcmcResult(CFmFileMatrix* pFileMat, long int idx, double fQval, CFmVector& fmVctRet)
{
    int nMcmc= m_pCfg->m_nMcmcIter - m_pCfg->GetBurnInRound();
    int q1 = (int)floor( fQval * nMcmc);
    int q2 = (int)floor( ( 1- fQval)* nMcmc -1 );

	CFmVector fmVctTmp(0, 0.0);
	int ret = pFileMat->GetCacheCol(idx, fmVctTmp);
	if(ret!=0) return(ret);

	int bSig = 1;
    fmVctTmp.Sort();
	if (fmVctTmp[q1]<=0 && fmVctTmp[q2]>=0)
	    bSig = 0;

	CFmVector fmTmp(0, 0.0);
	for(long int q=q1; q<= q2; q++)
		fmTmp.Put( fmVctTmp[q] );

	fmVctRet.Resize(0);
	fmVctRet.Put( bSig * 1.0 );
	fmVctRet.Put( fmTmp.GetMedian() );
	fmVctRet.Put( fmTmp.GetMin() );
	fmVctRet.Put( fmTmp.GetMax() );

    return(0);
}

int BLS_res::SortoutResult( CFmFileMatrix* pFileMat, CFmMatrix* gen_pA, double y_var, int nSnp, int nCov,
						     CFmMatrix* pmatCoef, CFmMatrix* pmatA, CFmMatrix* pmatD, CFmVector* pvctH2 )
{
    _log_debug(_HI_, "SortoutResult(), NCOL(FileMat)=%d, SNP=%d nCov=%d", pFileMat->GetNumCols(), nSnp, nCov );

    CFmVector fmVctRet(0, 0.0);

	for(long int i=0; i< nCov; i++)
	{
		GetMcmcResult(pFileMat, i, m_pCfg->m_fQval_add, fmVctRet);
		pmatCoef->Cbind( fmVctRet );
	}

	for (long int i=0; i< nSnp; i++)
    {
		GetMcmcResult(pFileMat, nCov + i , m_pCfg->m_fQval_add, fmVctRet);
		pmatA->Cbind( fmVctRet );
		if ( fmVctRet[0] >0 )
			_log_debug(_HI_, "significant SNP [%d], ADD = %.4f, %.4f, %.4f, %.4f", i+1, fmVctRet[0], fmVctRet[1], fmVctRet[2], fmVctRet[3] );
	}

	for (long int i=0; i< nSnp; i++)
    {
		GetMcmcResult(pFileMat, nCov + nSnp + i , m_pCfg->m_fQval_dom, fmVctRet);
		pmatD->Cbind( fmVctRet );
		if ( fmVctRet[0] >0 )
			_log_debug(_HI_, "significant SNP [%d], DOM = %.4f, %.4f, %.4f, %.4f", i+1, fmVctRet[0], fmVctRet[1], fmVctRet[2], fmVctRet[3] );
    }

	pmatA->Transpose();
	pmatD->Transpose();
	pmatCoef->Transpose();

	CFmVector fmAdd(0, 0.0);
	fmAdd= pmatA->GetCol(1);
	CFmVector fmDom(0, 0.0);
	fmDom = pmatD->GetCol(1);
    GetH2( gen_pA, fmAdd, fmDom, pvctH2, y_var );

    return(0);
}


int BLS_res::SetResult(bool bRefit, CFmMatrix* gen_pA, double y_var, int nCov,
						CFmVectorStr* pVctSnp, CFmVector* pVctChr, CFmVector* pVctPos, CFmFileMatrix* pFileMat )
{
	int nSnp = pVctSnp->GetLength();
    _log_debug(_HI_, "SetResult: bRefit=%d, nSNP=%d, nCov=%d", bRefit?1:0, nSnp, nCov );

	CFmNewTemp refNew;

    if (bRefit)
    {
        m_pRefit_SnpName = new (refNew) CFmVectorStr(pVctSnp);
        m_pRefit_SnpChr = new (refNew) CFmVector(pVctChr);
        m_pRefit_SnpPos = new (refNew) CFmVector(pVctPos);
        m_pRefit_Coefs = new (refNew) CFmMatrix(0, 0);
        m_pRefit_Ra = new (refNew) CFmMatrix(0, 0);
        m_pRefit_Rd = new (refNew) CFmMatrix(0, 0);
        m_pRefit_Rh2 = new (refNew) CFmVector(0, 0.0);

 		SortoutResult( pFileMat, gen_pA, y_var, nSnp, nCov, m_pRefit_Coefs, m_pRefit_Ra, m_pRefit_Rd, m_pRefit_Rh2 );

 		m_pRefit_QBest = new (refNew) CFmMatrix(nSnp, 8);
 		SortoutQBest( pFileMat, m_pRefit_QBest, nSnp, nCov);
 		m_pRefit_QBest->SetRowNames(m_pRefit_SnpName);

    }
    else
    {
        m_pVarsel_SnpName = new (refNew) CFmVectorStr(pVctSnp);
        m_pVarsel_SnpChr = new (refNew) CFmVector(pVctChr);
        m_pVarsel_SnpPos = new (refNew) CFmVector(pVctPos);
        m_pVarsel_Coefs = new (refNew) CFmMatrix(0, 0);
        m_pVarsel_Ra = new (refNew) CFmMatrix(0, 0);
        m_pVarsel_Rd = new (refNew) CFmMatrix(0, 0);
        m_pVarsel_Rh2 = new (refNew) CFmVector(0, 0.0);

		SortoutResult( pFileMat, gen_pA, y_var, nSnp, nCov, m_pVarsel_Coefs, m_pVarsel_Ra, m_pVarsel_Rd, m_pVarsel_Rh2 );

 		m_pVarsel_QBest = new (refNew) CFmMatrix(nSnp, 8);
 		SortoutQBest( pFileMat, m_pVarsel_QBest, nSnp, nCov);
 		m_pVarsel_QBest->SetRowNames(m_pVarsel_SnpName);
    }

    return(0);
}

int BLS_res::InitRefit(int nMcmcIter, CFmVectorStr* pCovNames)
{
	CFmNewTemp refNew;

    if (m_pRefitSnps) destroy( m_pRefitSnps );
    m_pRefitSnps = new (refNew) CFmVector( 0, 0.0 );

	if (m_pCovNames==NULL)
		m_pCovNames = new (refNew) CFmVectorStr(pCovNames);

    _log_debug(_HI_, "InitRefit: VARSEL[%d,%d]", m_pVarsel_Ra->GetNumRows(), m_pVarsel_Ra->GetNumCols() );

	for (long int i=0; i< m_nSnpP; i++)
    {
		if ( m_pVarsel_Ra->Get(i,0) >0 || m_pVarsel_Rd->Get(i,0)>0 )
            m_pRefitSnps->Put(i);
    }

    if (m_pRefitSnps->GetLength()<=0)
    {
        _log_prompt( _HI_, "REFIT: No SNP is selected for the refit procedure");
        return(0);
    }

    _log_debug(_HI_, "InitRefit: %d SNPs are selected to refit.", m_pRefitSnps->GetLength() );

    m_nRefitSnp = m_pRefitSnps->GetLength();
    m_nMcmcSnp = m_pRefitSnps->GetLength();
    m_nMcmcSect = 1;
    m_nMcmcIter = nMcmcIter;

    return(0);
}

int BLS_res::GetSigSNP()
{
	CFmNewTemp refNew;

	int nSnp = m_nRefitSnp;

	if (m_pSigSnps) destroy( m_pSigSnps );
	if (m_pSigDomSnps) destroy( m_pSigDomSnps );
	if (m_pSigAddSnps) destroy( m_pSigAddSnps );

	m_pSigSnps = new (refNew) CFmVector(0, 0.0);
	m_pSigDomSnps = new (refNew) CFmVector(0, 0.0);
	m_pSigAddSnps = new (refNew) CFmVector(0, 0.0);

	_log_debug(_HI_, "GetSigSNP: NSnp=%d", nSnp );

	for (long int i=0; i< nSnp; i++)
    {
		if ( m_pRefit_Ra->Get(i,0) >0 || m_pRefit_Rd->Get(i,0)>0 )
            m_pSigSnps->Put(i);

		if ( m_pRefit_Rd->Get(i,0) >0 )
            m_pSigDomSnps->Put(i);

		if (  m_pRefit_Ra->Get(i,0) >0 )
            m_pSigAddSnps->Put(i);
    }

    _log_debug(_HI_, "SelectSigSNP: %d SNPs are selected to refit.", m_pSigSnps->GetLength() );

    return( m_pSigSnps->GetLength() );
}

int BLS_res::Summary(char* szOutFile)
{
    char szTemp[MAX_PATH*256] = {0};

    sprintf(szTemp, "%s \nSUMMARY(Result)\n", szTemp );
    sprintf(szTemp, "%s -------------------------------------------\n", szTemp );
    sprintf(szTemp, "%s RESULT Subjects    = %d\n",  szTemp, m_nSubjN);
    sprintf(szTemp, "%s RESULT Snps        = %d\n",  szTemp, m_nSnpP);
    sprintf(szTemp, "%s RESULT Q.add       = %.3f\n",szTemp, m_pCfg->m_fQval_add);
    sprintf(szTemp, "%s RESULT Q.dom       = %.3f\n",szTemp, m_pCfg->m_fQval_dom);
    sprintf(szTemp, "%s RESULT Refit SNP   = %d\n",  szTemp, m_nRefitSnp);
    for(int i=0; i<m_nRefitSnp && i<20 ; i++)
    {
        int idx = (int) ( m_pRefitSnps->Get(i) );
        sprintf(szTemp, "%s RESULT (VARSEL) [%02d]= %s, %.2f, %.2f\n",
                szTemp, idx+1, m_pVarsel_SnpName->Get(idx),
                m_pVarsel_Ra->Get(idx, 1), m_pVarsel_Rd->Get(idx, 1) );
    }

    CFmVector* pSigSnpA=m_pSigAddSnps;
    CFmVector* pSigSnpD=m_pSigDomSnps;

    sprintf(szTemp, "%s RESULT (REFIT) Sig. Additive SNP   = %d\n", szTemp, m_nRefitSnp);

    for(int i=0; i<pSigSnpA->GetLength() && i<20 ; i++)
    {
        int idx = (int)( pSigSnpA->Get(i) );
        sprintf(szTemp, "%s RESULT A[%02d]= %s, %.2f, %.2f\n",
                szTemp, idx+1, m_pRefit_SnpName->Get(idx),
                m_pRefit_Ra->Get(idx, 1), m_pRefit_Rh2->Get(idx) );
    }

    sprintf(szTemp, "%s RESULT (REFIT) Sig. Dominant SNP   = %d\n", szTemp, m_nRefitSnp);

    for(int i=0; i<pSigSnpD->GetLength() && i<20 ; i++)
    {
        int idx = (int)( pSigSnpD->Get(i) );
        sprintf(szTemp, "%s RESULT D[%02d]= %s, %.2f, %.2f\n",
                szTemp, i+1, m_pRefit_SnpName->Get(idx),
                m_pRefit_Rd->Get(idx, 1), m_pRefit_Rh2->Get(idx) );
    }
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

int BLS_res::SaveRData( char* szRdataFile )
{
    _log_info(_HI_, "SaveResData: Start to save the RDATA to file(%s).", szRdataFile);

    CFmRls rls(szRdataFile);

    // #NO1: adh2
    if (m_pVarsel_SnpName)
    {
		CFmMatrix matVs1(0,0);

        matVs1.Cbind(*m_pVarsel_SnpChr, (char*)("chr"));
        matVs1.Cbind(*m_pVarsel_SnpPos, (char*)("pos"));
        matVs1.Cbind(*m_pVarsel_Ra);
        matVs1.Cbind(*m_pVarsel_Rd);
        matVs1.Cbind(*m_pVarsel_Rh2, (char*)("h2"));

		matVs1.SetRowNames(m_pVarsel_SnpName);

		rls.SetData( (char*)"varsel", &matVs1);
	}

    if (m_pRefit_SnpName)
    {
		CFmMatrix matVs1(0,0);

        matVs1.Cbind(*m_pRefit_SnpChr, (char*)("chr"));
        matVs1.Cbind(*m_pRefit_SnpPos, (char*)("pos"));
        matVs1.Cbind(*m_pRefit_Ra );
        matVs1.Cbind(*m_pRefit_Rd );
        matVs1.Cbind(*m_pRefit_Rh2, (char*)("h2"));

		matVs1.SetRowNames(m_pRefit_SnpName);

		rls.SetData( (char*)"refit", &matVs1);
	}

    if(m_pVarsel_Coefs)
    {
        rls.SetData( (char*)"varsel.cov", m_pVarsel_Coefs);
    }

    if(m_pRefit_Coefs)
    {
        rls.SetData( (char*)"refit.cov", m_pRefit_Coefs);
    }

    if(m_pVarsel_QBest)
    {
        rls.SetData( (char*)"varsel.Qbest", m_pVarsel_QBest);
    }

    if(m_pRefit_QBest)
    {
        rls.SetData( (char*)"refit.Qbest", m_pRefit_QBest);
    }

    int nRet = rls.SaveRData( szRdataFile );
    if(nRet==0)
        _log_info(_HI_, "SaveResData: Successful to save the RDATA to file(%s).", szRdataFile);

    return(nRet);
}

int BLS_res::SaveCsv4Fig( char* szCsvFile, bool bRefit )
{
	CFmMatrix matVs1(0,0);

    if (!bRefit)
    {
        matVs1.Cbind(*m_pVarsel_SnpChr, (char*)("chr"));
        matVs1.Cbind(*m_pVarsel_SnpPos, (char*)("pos"));
        matVs1.Cbind(*m_pVarsel_Ra );
        matVs1.Cbind(*m_pVarsel_Rd );
        matVs1.Cbind(*m_pVarsel_Rh2, (char*)("h2"));

		matVs1.SetRowNames(m_pVarsel_SnpName);
	}

    if (bRefit)
    {
        matVs1.Cbind(*m_pRefit_SnpChr, (char*)("chr"));
        matVs1.Cbind(*m_pRefit_SnpPos, (char*)("pos"));
        matVs1.Cbind(*m_pRefit_Ra );
        matVs1.Cbind(*m_pRefit_Rd );
        matVs1.Cbind(*m_pRefit_Rh2, (char*)("h2"));

		matVs1.SetRowNames(m_pRefit_SnpName);
	}


    matVs1.WriteAsCSVFile(szCsvFile, false);
    return(0);

}

int BLS_res::SaveSigFile( char* szCsvFile, bool bRefit )
{
	_log_info(_HI_, "SaveSigFile: CSV File=%s bRefit=%d", szCsvFile,  bRefit?1:0 );

	CFmMatrix matVs1(0,0);

    if (!bRefit)
    {
        matVs1.Cbind(*m_pVarsel_SnpChr, (char*)("chr"));
        matVs1.Cbind(*m_pVarsel_SnpPos, (char*)("pos"));
        matVs1.Cbind(*m_pVarsel_Ra );
        matVs1.Cbind(*m_pVarsel_Rd );
        matVs1.Cbind(*m_pVarsel_Rh2, (char*)("h2"));

		matVs1.SetRowNames(m_pVarsel_SnpName);
	}

    if (bRefit)
    {
        matVs1.Cbind(*m_pRefit_SnpChr, (char*)("chr"));
        matVs1.Cbind(*m_pRefit_SnpPos, (char*)("pos"));
        matVs1.Cbind(*m_pRefit_Ra );
        matVs1.Cbind(*m_pRefit_Rd );
        matVs1.Cbind(*m_pRefit_Rh2, (char*)("h2"));

		matVs1.SetRowNames(m_pRefit_SnpName);
	}

	CFmMatrix matVs2(0,0);
	for(int i=0;i<matVs1.GetNumRows();i++)
	{
		if ( m_pVarsel_Ra->Get(i,0)>0 || m_pVarsel_Rd->Get(i,1)>0 )
		{
			matVs2.Cbind(matVs1.GetRow(i));
		}
	}

	matVs2.Transpose();
    matVs2.WriteAsCSVFile(szCsvFile, false);
    return(0);

}

SEXP BLS_res::GetRObj()
{
    _log_info(_HI_, "GetRObj: Start to save the matrix to R.");

	SEXP sRet, t;
	int nList = (m_pVarsel_SnpName != NULL?1:0) +
			  	(m_pRefit_SnpName  != NULL?1:0) +
				(m_pVarsel_Coefs   != NULL?1:0) +
				(m_pRefit_Coefs    != NULL?1:0) +
				(m_pVarsel_QBest   != NULL?1:0) +
				(m_pRefit_QBest    != NULL?1:0);

    PROTECT(sRet = t = allocList(nList));

    // #NO1: ret.vs
    if (m_pVarsel_SnpName)
    {
		CFmMatrix matVs1(0,0);

        matVs1.Cbind(*m_pVarsel_SnpChr, (char*)("chr"));
        matVs1.Cbind(*m_pVarsel_SnpPos, (char*)("pos"));
        matVs1.Cbind(*m_pVarsel_Ra );
        matVs1.Cbind(*m_pVarsel_Rd );
        matVs1.Cbind(*m_pVarsel_Rh2, (char*)("h2"));

		matVs1.SetRowNames(m_pVarsel_SnpName);
		SEXP expVS1 = GetSEXP(&matVs1);
		SETCAR( t, expVS1 );
		SET_TAG(t, install("varsel") );
		t = CDR(t);
	}

    if (m_pRefit_SnpName)
    {
		CFmMatrix matVs1(0,0);

        matVs1.Cbind(*m_pRefit_SnpChr, (char*)("chr"));
        matVs1.Cbind(*m_pRefit_SnpPos, (char*)("pos"));
        matVs1.Cbind(*m_pRefit_Ra );
        matVs1.Cbind(*m_pRefit_Rd );
        matVs1.Cbind(*m_pRefit_Rh2, (char*)("h2"));

		matVs1.SetRowNames(m_pRefit_SnpName);
		SEXP expVS1 = GetSEXP(&matVs1);
		SETCAR( t, expVS1 );
		SET_TAG(t, install("refit") );
		t = CDR(t);
     }

	if(m_pVarsel_Coefs)
	{
		SEXP expVS = GetSEXP(m_pVarsel_Coefs);
		SETCAR( t, expVS );
		SET_TAG(t, install("varsel_cov") );
		t = CDR(t);
	}

	if(m_pRefit_Coefs)
	{
		SEXP expVS = GetSEXP(m_pRefit_Coefs);
		SETCAR( t, expVS );
		SET_TAG(t, install("refit_cov") );
		t = CDR(t);
	}

	if(m_pVarsel_QBest)
	{
		SEXP expVS = GetSEXP(m_pVarsel_QBest);
		SETCAR( t, expVS );
		SET_TAG(t, install("varsel_Qbest") );
		t = CDR(t);
	}

	if(m_pRefit_QBest)
	{
		SEXP expVS = GetSEXP(m_pRefit_QBest);
		SETCAR( t, expVS );
		SET_TAG(t, install("refit_Qbest") );
		t = CDR(t);
	}

	UNPROTECT(1);

    _log_debug(_HI_, "GetRObj: Successful to save the matrix to R.");

    return( sRet );
}

//#-------------computing h2--------------
//Y  <- dat$prec_y;
//N  <- length(Y);
//p1 <- c();
//for (j in 1:p )
//{
//       tmp <- dat$snps[j,];
//        p1  <- c(p1, (2*length(which(tmp==1))+length(which(tmp==0)))/2/N );
//}
//
//p0 <- 1 - p1;
//var_a <- 2*p1*p0*(me_a + (p1-p0)*me_d)^2;
//var_d <- 4*(p1^2)*(p0^2)*(me_d^2);
//var_d <- var_d*3;
//r_h2  <- (var_a + var_d)/var(Y)*100;

int BLS_res::GetH2(CFmMatrix* gen_pA, CFmVector& vctA, CFmVector& vctD, CFmVector* pVctH2, double var_y )
{
    CFmVector tmp1(0, 0.0);
    CFmVector tmp_1(0, 0.0);
    pVctH2->Resize(0);

    int N = gen_pA->GetNumRows();

    for(int i=0; i<gen_pA->GetNumCols(); i++)
    {
        tmp1 = gen_pA->GetCol(i);
        tmp1 = (tmp1*tmp1 + tmp1)/2;

        tmp_1 = gen_pA->GetCol(i);
        tmp_1 = (tmp_1*tmp_1 - tmp_1)/2;

        double f0 = ( N - tmp1.Sum() - tmp_1.Sum() );
        double p1 = (2*tmp1.Sum() + f0)/(2.0*N);
        double p0 = 1-p1;

        double var_a = 2*p1*p0*pow(vctA[i] + (p1-p0)*vctD[i],2);
        double var_d = 4*(p1*p1)*(p0*p0)*(vctD[i]*vctD[i]);
        //double var_d = var_d*3;
        double r_h2  = (var_a + var_d)/var_y*100;
        pVctH2->Put(r_h2);
    }

    return(0);
}

int BLS_res::SortoutQBest(CFmFileMatrix* pMatRet, CFmMatrix* pQBest, int nSnp, int nCov)
{
	_log_debug(_HI_, "SortoutQBest: FileMatrix=%s NROW(pFileMat)=%d nSNP=%d nCov=%d ",
			pMatRet->GetFileName(), pMatRet->GetNumRows(), nSnp, nCov );

	CFmVector fmQBest(4, 0.0);

	for (long int i=0; i<nCov; i++)
		GetQBestInfo( pMatRet, i, fmQBest);

	CFmVector fmTemp(0, 0.0);
	for (long int i=0; i< nSnp; i++)
    {
		fmTemp.Resize(0);

		GetQBestInfo( pMatRet, nCov+ i, fmQBest);
		fmTemp.Append( fmQBest);

		GetQBestInfo( pMatRet, nCov+ nSnp +i, fmQBest);
		fmTemp.Append( fmQBest);

		pQBest->SetRow(i, fmTemp );
    }

    _log_debug(_HI_, "SortoutQBest: Successful to sort out.");
    return(0);
}

int BLS_res::GetQBestInfo(CFmFileMatrix* pMatRet, int idx, CFmVector& fmQBest)
{
    int nMcmc= m_pCfg->m_nMcmcIter - m_pCfg->GetBurnInRound();

	fmQBest.Resize( 4 );
	fmQBest[0] = R_NaN;
	fmQBest[1] = R_NaN;
	fmQBest[2] = R_NaN;
	fmQBest[3] = R_NaN;

	CFmVector fmVct(nMcmc, 0.0);
    int ret = pMatRet->GetCacheCol( idx, fmVct );
    if(ret!=0) return( ret );
    fmVct.Sort();

	CFmVector fmTmp(0, 0.0);
	for(int qx=0 ; qx < nMcmc/2; qx++)
		if(( fmVct[qx]>0 && fmVct[nMcmc-1-qx]>0 ) || ( fmVct[qx]<0 && fmVct[nMcmc-1-qx]<0 ) )
		{
			fmTmp.Resize(0);
			for(int qxi = qx; qxi <= nMcmc-1-qx; qxi++)
				fmTmp.Put(fmVct[qxi]);

			fmQBest[0] = qx/(nMcmc*1.0 );
			fmQBest[1] = fmTmp.GetMedian();
			fmQBest[2] = fmTmp.GetMin();
			fmQBest[3] = fmTmp.GetMax();

			break;
		}

	if(isnan(fmQBest[0]))
	{
		fmQBest[0] = 0.5;
		fmQBest[1] = fmVct[nMcmc/2];
		fmQBest[2] = fmVct[nMcmc/2];
		fmQBest[3] = fmVct[nMcmc/2];
	}

	return(0);
}

void destroy(BLS_res* p)
{
	CFmNewTemp  fmRef;
	p->~BLS_res();
	operator delete(p, fmRef);
}
