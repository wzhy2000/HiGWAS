/* gls_RI.cpp  -	GLS application
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include "gls_R.h"

#define BOOLEAN_ELT(x,__i__)	LOGICAL(x)[__i__]
#define INTEGER_ELT(x,__i__)	INTEGER(x)[__i__]
#define NUMERIC_ELT(x,__i__)	REAL(x)[__i__]

void gls_simulate( char** pszPhe_out,
					 char** pszSnp_out,
					 int* pnSimu_grp,
					 int* pnSimu_n,
					 int* pnSimu_p,
					 double* pfSimu_snp_rho,
					 double* pfSimu_snp_missing,
					 double* pfSimu_rho,
					 double* pfSimu_sigma2,
					 double* pfSimu_mu,
					 int* pnSimu_covar_len,
					 double* pfSimu_covar_range,
					 double* pfSimu_covar_effect,
					 int* pnSimu_sig_p,
					 int* pnSimu_add_len,
					 double* pfSimu_add_effect,
					 int* pnSimu_dom_len,
					 double* pfSimu_dom_effect,
					 double* pfSimu_z_range,
					 int* pnSimu_z_count,
					 int* pnDebug,
					 int* err)
{
	int ret = glasso_simulate(  *pszPhe_out,
						 *pszSnp_out,
						 *pnSimu_grp,
						 *pnSimu_n,
						 *pnSimu_p,
						 *pfSimu_snp_rho,
						 *pfSimu_snp_missing,
						 *pfSimu_rho,
						 *pfSimu_sigma2,
						 pfSimu_mu,
						 *pnSimu_covar_len,
						 pfSimu_covar_effect,
						 pfSimu_covar_range,
						 *pnSimu_sig_p,
						 *pnSimu_add_len,
						 pfSimu_add_effect,
						 *pnSimu_dom_len,
						 pfSimu_dom_effect,
						 pfSimu_z_range,
					     pnSimu_z_count,
					     *pnDebug);
	*err = ret;
}

SEXP gls_simple( SEXP szPhefile,
  		   	SEXP szSnpfile,
  		   	SEXP szYname,
  		   	SEXP szZname,
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
	const char* pszZname   = CHAR(STRING_ELT(szZname,0));
	const char* pszXname   = CHAR(STRING_ELT(szXname,0));

	bool bRefit    = BOOLEAN_ELT(sbRefit, 0);
	bool bAddUsed  = BOOLEAN_ELT(sbAddUsed, 0);
	bool bDomUsed  = BOOLEAN_ELT(sbDomUsed, 0);

	int  nMcmcIter      = round(NUMERIC_ELT( snMcmcIter, 0) );
	double fBurnInRound = NUMERIC_ELT( sfBurnInRound, 0);
	double fRhoTuning   = NUMERIC_ELT( sfRhoTuning,0);
	double fQval_add    = NUMERIC_ELT( sfQval_add, 0);
	double fQval_dom    = NUMERIC_ELT( sfQval_dom, 0);
	int   nDebug        = round(NUMERIC_ELT(snDebug, 0));

	SEXP ret = glasso_simple( pszPhefile,
  		   			pszSnpfile,
  		   			pszYname,
  		   			pszZname,
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

SEXP gls_plink_tped( SEXP szPhefile,
  		   	SEXP szTpedfile,
  		   	SEXP szTfamfile,
  		   	SEXP szYname,
  		   	SEXP szZname,
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

	bool bRefit   = BOOLEAN_ELT(sbRefit, 0);
	bool bAddUsed = BOOLEAN_ELT(sbAddUsed, 0);
	bool bDomUsed = BOOLEAN_ELT(sbDomUsed, 0);

	const char* pszYname = CHAR(STRING_ELT(szYname,0));
	const char* pszZname = CHAR(STRING_ELT(szZname,0));
	const char* pszXname = CHAR(STRING_ELT(szXname,0));

	int  nMcmcIter   = round( NUMERIC_ELT( snMcmcIter, 0) );
	double fBurnInRound = NUMERIC_ELT( sfBurnInRound, 0);
	double fRhoTuning = NUMERIC_ELT( sfRhoTuning,0);
	double fQval_add  = NUMERIC_ELT( sfQval_add, 0);
	double fQval_dom  = NUMERIC_ELT( sfQval_dom, 0);
	int    nDebug     = round( NUMERIC_ELT(snDebug, 0) );

	SEXP ret = glasso_plink_tped( pszPhefile,
  		   			pzTpedfile,
  		   			pzTfamfile,
  		   			pszYname,
  		   			pszZname,
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

SEXP gls_snpmat( SEXP smatPhe,
  		   	SEXP smatSnp,
  		   	SEXP szYname,
  		   	SEXP szZname,
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
	bool bRefit   = BOOLEAN_ELT(sbRefit, 0);
	bool bAddUsed = BOOLEAN_ELT(sbAddUsed, 0);
	bool bDomUsed = BOOLEAN_ELT(sbDomUsed, 0);

	const char* pszYname = CHAR(STRING_ELT(szYname,0));
	const char* pszZname = CHAR(STRING_ELT(szZname,0));
	const char* pszXname = CHAR(STRING_ELT(szXname,0));

	int  nMcmcIter    = round( NUMERIC_ELT( snMcmcIter, 0) );
	double fBurnInRound = NUMERIC_ELT( sfBurnInRound, 0);
	double fRhoTuning = NUMERIC_ELT( sfRhoTuning,0);
	double fQval_add  = NUMERIC_ELT( sfQval_add, 0);
	double fQval_dom  = NUMERIC_ELT( sfQval_dom, 0);
	int    nDebug     = round( NUMERIC_ELT( snDebug, 0) );

	SEXP ret = glasso_snpmat( smatPhe,
  		   			smatSnp,
  		   			pszYname,
  		   			pszZname,
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


