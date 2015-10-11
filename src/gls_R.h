// gls_app.h: interface for the glasso command.
//
//////////////////////////////////////////////////////////////////////

#if !defined(GLS_APP_H__INCLUDED_)
#define GLS_APP_H__INCLUDED_

#include <stdbool.h>
#include "fm_linux.h"
#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#ifdef __cplusplus
extern "C" {
#endif

#define _RUN_BOTH       0   //Step 1+2: both
#define _RUN_VARSEL     1   //Step 1: variable selection
#define _RUN_REFIT      2   //Step 2: Refit

#define MAX_RESIN_CNT      128

typedef struct cmd_param{
    bool bHelp;
    bool bVersion;
    bool bZNormalized;
    bool bAddUsed;
    bool bDomUsed;
    bool bRefit;
    int  nDebug;
    int  nRunmode;
    int  nSimuRound;
    char szAppName[MAX_PATH];
    char szParFile[MAX_PATH];
    char szTpedFile[MAX_PATH];
    char szTfamFile[MAX_PATH];
    char szVsretFile[MAX_PATH];
    char szLogFile[MAX_PATH];
    char szPheFile[MAX_PATH];
    char szSnpFile[MAX_PATH];
    char szRdataFile[MAX_PATH];
    char szFigOutFile[MAX_PATH];
    char szSigOutFile[MAX_PATH];
    char szPcfFile[MAX_PATH];
    char szCfgFile[MAX_PATH];
    char szResFile[MAX_PATH];
    char szYname[MAX_PATH];
    char szXname[MAX_PATH];
    char szZname[MAX_PATH];
    char szPheoutFile[MAX_PATH];
    char szSnpoutFile[MAX_PATH];
    char szMifFile[MAX_PATH];
    char szResinFile[MAX_PATH];
    char szMatadFile[MAX_PATH];
} CMDOPTIONS;


int glasso_simulate( const char* szPhe_out,
			 const char* szSnp_out,
			 int nSimu_grp,
			 int nSimu_n,
			 int nSimu_p,
			 double fSimu_snp_rho,
			 double fSimu_snp_missing,
			 double fSimu_rho,
			 double fSimu_sigma2,
			 double* pfSimu_mu,
			 int nSimu_covar_len,
			 double* pfSimu_covar_effect,
			 double* pfSimu_covar_range,
			 int nSimu_sig_p,
			 int nSimu_add_len,
			 double* pfSimu_add_effect,
			 int nSimu_dom_len,
			 double* pfSimu_dom_effect,
			 double* pfSimu_z_range,
			 int* pnSimu_longdt_points,
			 int nDebug);

SEXP glasso_simple( const char* pzPhefile,
  		   	const char*  pzSnpfile,
  		   	const char*  pzYname,
  		   	const char*  pzZname,
  		   	const char*  pzXname,
  		   	bool bRefit,
  		   	bool bAddUsed,
  		   	bool bDomUsed,
  		   	int nMcmcIter,
		   	double fBurnInRound,
			double fRhoTuning,
		   	double fQval_add,
		   	double fQval_dom,
		   	int nDebug);

SEXP glasso_plink_tped( const char* pzPhefile,
  		   	const char*  pzTpedfile,
  		   	const char*  pzTfamfile,
  		   	const char*  pzYname,
  		   	const char*  pzZname,
  		   	const char*  pzXname,
  		   	bool bRefit,
  		   	bool bAddUsed,
  		   	bool bDomUsed,
  		   	int nMcmcIter,
		   	double fBurnInRound,
			double fRhoTuning,
		   	double fQval_add,
		   	double fQval_dom,
		   	int nDebug);


SEXP glasso_snpmat( SEXP smatPhe,
  		   	SEXP smatSnp,
  		   	const char*  pzYname,
  		   	const char*  pzZname,
  		   	const char*  pzXname,
  		   	bool bRefit,
  		   	bool bAddUsed,
  		   	bool bDomUsed,
  		   	int nMcmcIter,
		   	double fBurnInRound,
			double fRhoTuning,
		   	double fQval_add,
		   	double fQval_dom,
		   	int nDebug);
#ifdef __cplusplus
}
#endif

#endif


