/* bls_RI.cpp  -	BLS application
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include "bls_R.h"

#define BOOLEAN_ELT(x,__i__)	LOGICAL(x)[__i__]
#define INTEGER_ELT(x,__i__)	INTEGER(x)[__i__]
#define NUMERIC_ELT(x,__i__)	REAL(x)[__i__]

void bls_simulate( char** pszPhe_out,
					 char** pszSnp_out,
					 int* pnSimu_grp,
					 int* pnSimu_n,
					 int* pnSimu_p,
					 double* pfSimu_snp_rho,
					 double* pfSimu_snp_missing,
					 double* pfSimu_rho,
					 double* pfSimu_sigma2,
					 double* pfSimu_mu,
					 int *pnSimu_covar_count,
					 double* pfSimu_covar_coeff,
					 int* pnSimu_sig_p,
					 int* pnSimu_a_len,
					 int* pnSimu_a_pos,
					 double* pfSimu_a_effect,
					 int* pnSimu_d_len,
					 int* pnSimu_d_pos,
					 double* pfSimu_d_effect,
					 double* pfSimu_covar_range,
					 double* pfSimu_t_range,
					 int* pnDebug,
					 int* err)
{
	int ret = blasso_simulate(  *pszPhe_out,
						 *pszSnp_out,
						 *pnSimu_grp,
						 *pnSimu_n,
						 *pnSimu_p,
						 *pfSimu_snp_rho,
						 *pfSimu_snp_missing,
						 *pfSimu_rho,
						 *pfSimu_sigma2,
						 *pfSimu_mu,
						 *pnSimu_covar_count,
						 pfSimu_covar_coeff,
						 *pnSimu_sig_p,
						 *pnSimu_a_len,
						 pnSimu_a_pos,
						 pfSimu_a_effect,
						 *pnSimu_d_len,
						 pnSimu_d_pos,
						 pfSimu_d_effect,
						 pfSimu_covar_range,
						 pfSimu_t_range,
					     *pnDebug);
	*err = ret;
}

SEXP bls_simple( SEXP szPhefile,
  		   	SEXP szSnpfile,
  		   	SEXP szYname,
  		   	SEXP szXname,
  		   	SEXP sbRefit,
  		   	SEXP sbAddUsed,
  		   	SEXP sbDomUsed,
  		   	SEXP snMcmcIter,
		   	SEXP sfBurnInRound,
		   	SEXP sfRhoTuning,
	        SEXP sfQval_add,
	        SEXP sfQval_dom,
			SEXP snDebug)
{
	const char* pszPhefile = CHAR(STRING_ELT(szPhefile,0));
	const char* pszSnpfile = CHAR(STRING_ELT(szSnpfile,0));
	const char* pszYname   = CHAR(STRING_ELT(szYname,0));
	const char* pszXname   = CHAR(STRING_ELT(szXname,0));
	bool bRefit            = BOOLEAN_ELT(sbRefit, 0);
	bool bAddUsed          = BOOLEAN_ELT(sbAddUsed, 0);
	bool bDomUsed          = BOOLEAN_ELT(sbDomUsed, 0);
	int  nMcmcIter         = round(NUMERIC_ELT( snMcmcIter, 0) );
	double fRhoTuning      = NUMERIC_ELT( sfRhoTuning, 0 );
	double fBurnInRound    = NUMERIC_ELT( sfBurnInRound, 0 );
	double fQval_add       = NUMERIC_ELT( sfQval_add, 0 );
	double fQval_dom       = NUMERIC_ELT( sfQval_dom, 0 );
	int    nDebug          = round(NUMERIC_ELT( snDebug, 0 ) );

	SEXP ret = blasso_simple( pszPhefile,
  		   			pszSnpfile,
  		   			pszYname,
  		   			pszXname,
  		   			bRefit,
  		   			bAddUsed,
  		   			bDomUsed,
  		   			nMcmcIter,
		   			fBurnInRound,
		   			fRhoTuning,
	        		fQval_add,
	        		fQval_dom,
	        		nDebug);
	return(ret);
}

SEXP bls_plink_tped( SEXP szPhefile,
  		   	SEXP szTpedfile,
  		   	SEXP szTfamfile,
  		   	SEXP szYname,
  		   	SEXP szXname,
  		   	SEXP sbRefit,
  		   	SEXP sbAddUsed,
  		   	SEXP sbDomUsed,
  		   	SEXP snMcmcIter,
		   	SEXP sfBurnInRound,
		   	SEXP sfRhoTuning,
	        SEXP sfQval_add,
	        SEXP sfQval_dom,
			SEXP snDebug )
{
	const char* pszPhefile = CHAR(STRING_ELT(szPhefile,0));
	const char* pzTpedfile = CHAR(STRING_ELT(szTpedfile,0));
	const char* pzTfamfile = CHAR(STRING_ELT(szTfamfile,0));
	const char* pszYname   = CHAR(STRING_ELT(szYname,0));
	const char* pszXname   = CHAR(STRING_ELT(szXname,0));
	bool bRefit            = BOOLEAN_ELT(sbRefit, 0);
	bool bAddUsed          = BOOLEAN_ELT(sbAddUsed, 0);
	bool bDomUsed          = BOOLEAN_ELT(sbDomUsed, 0);
	int  nMcmcIter         = round(NUMERIC_ELT( snMcmcIter, 0) );
	double fRhoTuning      = NUMERIC_ELT( sfRhoTuning, 0);
	double fBurnInRound    = NUMERIC_ELT( sfBurnInRound, 0);
	double fQval_add       = NUMERIC_ELT( sfQval_add, 0);
	double fQval_dom       = NUMERIC_ELT( sfQval_dom, 0);
	int   nDebug           = round(NUMERIC_ELT(snDebug, 0) );

	SEXP ret = blasso_plink_tped( pszPhefile,
  		   			pzTpedfile,
  		   			pzTfamfile,
  		   			pszYname,
  		   			pszXname,
  		   			bRefit,
  		   			bAddUsed,
  		   			bDomUsed,
  		   			nMcmcIter,
		   			fBurnInRound,
		   			fRhoTuning,
	        		fQval_add,
	        		fQval_dom,
	        		nDebug);
	return(ret);
}

SEXP bls_snpmat( SEXP szPheMat,
  		   	SEXP szSnpMat,
  		   	SEXP szYname,
  		   	SEXP szXname,
  		   	SEXP sbRefit,
  		   	SEXP sbAddUsed,
  		   	SEXP sbDomUsed,
  		   	SEXP snMcmcIter,
		   	SEXP sfBurnInRound,
		   	SEXP sfRhoTuning,
	        SEXP sfQval_add,
	        SEXP sfQval_dom,
			SEXP snDebug)
{
	const char* pszYname   = CHAR(STRING_ELT(szYname,0));
	const char* pszXname   = CHAR(STRING_ELT(szXname,0));
	bool bRefit            = BOOLEAN_ELT(sbRefit, 0);
	bool bAddUsed          = BOOLEAN_ELT(sbAddUsed, 0);
	bool bDomUsed          = BOOLEAN_ELT(sbDomUsed, 0);
	int  nMcmcIter         = round(NUMERIC_ELT( snMcmcIter, 0) );
	double fRhoTuning      = NUMERIC_ELT( sfRhoTuning, 0 );
	double fBurnInRound    = NUMERIC_ELT( sfBurnInRound, 0 );
	double fQval_add       = NUMERIC_ELT( sfQval_add, 0 );
	double fQval_dom       = NUMERIC_ELT( sfQval_dom, 0 );
	int   nDebug           = round(NUMERIC_ELT( snDebug, 0 ) );

	SEXP ret = blasso_snpmat( szPheMat,
  		   			szSnpMat,
  		   			pszYname,
  		   			pszXname,
  		   			bRefit,
  		   			bAddUsed,
  		   			bDomUsed,
  		   			nMcmcIter,
		   			fBurnInRound,
		   			fRhoTuning,
	        		fQval_add,
	        		fQval_dom,
	        		nDebug);
	return(ret);
}
