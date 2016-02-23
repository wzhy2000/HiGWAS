/* BLS_par.cpp  -	BLS Parameter Object for simulation
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

#include "Rmethods.h"
#include "fm_rlogger.h"
#include "fm_inifile.h"
#include "fm_err.h"
#include "fm_new.h"

#include "bls_dat.h"
#include "bls_cfg.h"
#include "bls_par.h"

BLS_par::BLS_par( CMDOPTIONS* pCmd )
{
    m_pCmd = pCmd;
    m_szParFile = NULL;

    LoadDefault();
}

BLS_par::~BLS_par()
{
    if(m_szParFile)
        Free(m_szParFile);

    _log_debug(_HI_, "BLS_par is released successfully.");
}

void BLS_par::LoadDefault()
{
    simu_grps    = 1;
    simu_n       = 100;
    simu_p       = 1000;
    simu_snp_rho = 0.1;
    simu_rho     = 0.4;
    simu_mu      = 26;
    simu_sigma2  = 3;

    simu_covar_len      = 2;
    simu_covar_coefs[0] = 0;
    simu_covar_coefs[1] = 0;

    simu_a_len   = 11;
    simu_a_pos[0]= 100;
    simu_a_pos[1]= 200;
    simu_a_pos[2]= 300;
    simu_a_pos[3]= 400;
    simu_a_pos[4]= 500;
    simu_a_pos[5]= 600;
    simu_a_pos[6]= 700;
    simu_a_pos[7]= 800;
    simu_a_pos[8]= 850;
    simu_a_pos[9]= 900;
    simu_a_pos[10]=950;
    simu_a_effect[0]= 1.2;
    simu_a_effect[1]= 1.2;
    simu_a_effect[2]= 1.2;
    simu_a_effect[3]= 0.8;
    simu_a_effect[4]= 0.8;
    simu_a_effect[5]= 0.8;
    simu_a_effect[6]= 0.4;
    simu_a_effect[7]= 0.4;
    simu_a_effect[8]= 1.2;
    simu_a_effect[9]= 0.8;
    simu_a_effect[10]=1.2;

    simu_d_len   = 11;
    simu_d_pos[0]=  50;
    simu_d_pos[1]= 150;
    simu_d_pos[2]= 250;
    simu_d_pos[3]= 350;
    simu_d_pos[4]= 450;
    simu_d_pos[5]= 550;
    simu_d_pos[6]= 650;
    simu_d_pos[7]= 750;
    simu_d_pos[8]= 850;
    simu_d_pos[9]= 900;
    simu_d_pos[10]=950;
	simu_d_effect[0]= 1.2;
    simu_d_effect[1]= 1.2;
    simu_d_effect[2]= 1.2;
    simu_d_effect[3]= 0.8;
    simu_d_effect[4]= 0.8;
    simu_d_effect[5]= 0.8;
    simu_d_effect[6]= 0.4;
    simu_d_effect[7]= 0.4;
    simu_d_effect[8]= 1.2;
    simu_d_effect[9]= 1.2;
    simu_d_effect[10]=0.8;

	simu_t_range[0] = -1;
	simu_t_range[1] = 1;
	simu_covar_range[0] = 0;
	simu_covar_range[1] = 1;

    sig_p = Get_SigP();
}

int BLS_par::Get_SigP()
{
    // simu_a_posi
    SEXP par_simu_a_posi, par_simu_d_posi;
    PROTECT( par_simu_a_posi = allocVector( INTSXP, simu_a_len ) );
    int* val = INTEGER(par_simu_a_posi);
    for (int i=0; i<simu_a_len; i++)
         val[i] = simu_a_pos[i];

    // simu_d_posi
    PROTECT( par_simu_d_posi = allocVector( INTSXP, simu_d_len ) );
    val = INTEGER(par_simu_d_posi);
    for (int i=0; i<simu_d_len; i++)
        val[i] = simu_d_pos[i];

    SEXP sig_pos0, e1;
    int errorOccurred;
    PROTECT(e1 = lang3(install("c"), par_simu_a_posi, par_simu_d_posi ) );
    PROTECT(sig_pos0 = R_tryEval(e1, R_GlobalEnv, &errorOccurred) );

    //Rf_duplicate doesnt work?
    //SEXP sigpos = Rf_duplicate( sig_pos0 );

    SEXP sig_pos, e2;
    PROTECT(e2 = lang2(install("unique"), sig_pos0 ) );
    PROTECT(sig_pos = R_tryEval(e2, R_GlobalEnv, &errorOccurred) );

    sig_p = LENGTH(sig_pos);

    UNPROTECT(6);

    return(sig_p);
}

int BLS_par::Load(char* szParFile)
{
    CFmIniFile iFile;
    if (!iFile.OpenIniFile(szParFile))
    {
        _log_error( _HI_, "Failed to open IniFile(%s).", szParFile);
        return( ERR_INF_FILE);
    }

    _log_info(_HI_, "LoadConfig: Start to read a parameter configuration file(%s).", szParFile);

    simu_grps  = iFile.ReadInt("blasso", "simu_grps", -1);
    simu_n     = iFile.ReadInt("blasso", "simu_n", -1);
    simu_p     = iFile.ReadInt("blasso", "simu_p", -1);
    simu_rho   = iFile.ReadDouble("blasso", "simu_rho", -1);
    simu_snp_rho= iFile.ReadDouble("blasso", "simu_snp_rho", -1);
    simu_sigma2= iFile.ReadDouble("blasso", "simu_sigma2", -1);
    simu_mu    = iFile.ReadDouble("blasso", "simu_mu", -1);

    simu_a_len = iFile.ReadInt("blasso", "simu_a_len", -1);
    simu_d_len = iFile.ReadInt("blasso", "simu_d_len", -1);

    double x_pos[100]={0};
    double x_eff[100]={0};

    const char* szaPos = iFile.ReadString("blasso", "simu_a_pos", "");
    iFile.SplitFloat(szaPos, x_pos);
    for (int i=0;i<simu_a_len;i++)
    {
        _log_debug(_HI_, "simu_a_pos:%.2f", x_pos[i]);
        simu_a_pos[i] = (int)round( x_pos[i] );
    }

    const char* szaEff = iFile.ReadString("blasso", "simu_a_effect", "");
    iFile.SplitFloat(szaEff, x_eff);
    for (int i=0;i<simu_a_len;i++)
    {
        _log_debug(_HI_, "simu_a_effect:%.2f", x_eff[i]);
        simu_a_effect[i] = x_eff[i];
    }

    const char* szdPos = iFile.ReadString("blasso", "simu_d_pos", "");
    iFile.SplitFloat(szdPos, x_pos);
    for (int i=0;i<simu_d_len;i++)
    {
        _log_debug(_HI_, "simu_d_pos:%.2f", x_pos[i]);
        simu_d_pos[i] = (int)round( x_pos[i] );
    }

    const char* szdEff = iFile.ReadString("blasso", "simu_d_effect", "");
    iFile.SplitFloat(szdEff, x_eff);
    for (int i=0;i<simu_d_len;i++)
    {
        _log_debug(_HI_, "simu_d_effect:%.2f", x_eff[i]);
        simu_d_effect[i] = x_eff[i] ;
    }

    simu_covar_len = iFile.ReadInt("blasso", "simu_covar_len", -1);
    const char* szdCoef = iFile.ReadString("blasso", "simu_covar_coeff", "");
    iFile.SplitFloat(szdCoef, x_eff);
    for (int i=0;i<simu_covar_len;i++)
    {
        _log_debug(_HI_, "simu_covar_coeff:%.2f", x_eff[i]);
        simu_covar_coefs[i] = x_eff[i] ;
    }

	//default T_range is (-1, 1);
    const char* sz_z_range = iFile.ReadString("blasso", "simu_t_range", "");
    iFile.SplitFloat(sz_z_range, x_pos);
    for (int i=0;i<2;i++)
    {
        _log_debug(_HI_, "simu_z_range:%.2f", x_pos[i]);
        simu_t_range[i] = x_pos[i];
    }

	//default cov_range is (-1, 1);
    const char* sz_cov_range = iFile.ReadString("blasso", "simu_cov_range", "");
    iFile.SplitFloat(sz_cov_range, x_pos);
    for (int i=0;i<2;i++)
    {
        _log_debug(_HI_, "simu_covar_range:%.2f", x_pos[i]);
        simu_covar_range[i] = x_pos[i];
    }

    sig_p = Get_SigP();
    m_szParFile = Strdup(szParFile);

    _log_info(_HI_, "parameter file is loaded successfully.(sig_p:%d)", sig_p);

	return(0);
}

int BLS_par::Summary( char* szOutFile )
{
    char szTemp[MAX_PATH*256] = {0};

    sprintf(szTemp, "%s \nSUMMARY(Parameter)\n", szTemp );
    sprintf(szTemp, "%s -------------------------------------------\n", szTemp );
    sprintf(szTemp, "%s PARAMETER parameter = %s\n", szTemp, m_szParFile );
    sprintf(szTemp, "%s PARAMETER simu_grps   = %d\n", szTemp, simu_grps );
    sprintf(szTemp, "%s PARAMETER simu_n      = %d\n", szTemp, simu_n );
    sprintf(szTemp, "%s PARAMETER simu_p      = %d\n", szTemp, simu_p );
    sprintf(szTemp, "%s PARAMETER simu_rho    = %.2f\n", szTemp, simu_rho );
    sprintf(szTemp, "%s PARAMETER simu_snp_rho    = %.2f\n", szTemp, simu_snp_rho );
    sprintf(szTemp, "%s PARAMETER simu_sigma2 = %.2f\n", szTemp, simu_sigma2 );
    sprintf(szTemp, "%s PARAMETER simu_mu     = %.2f\n", szTemp, simu_mu );
    sprintf(szTemp, "%s PARAMETER simu_t_range= %.2f, %.2f\n", szTemp, simu_t_range[0], simu_t_range[1] );
    sprintf(szTemp, "%s PARAMETER simu_cov_range= %.2f, %.2f\n", szTemp, simu_covar_range[0], simu_covar_range[1] );

    for(int i=0;i<simu_covar_len;i++)
        sprintf(szTemp, "%s PARAMETER covariate coeffient(%d) :%.2f\n",
                        szTemp, i, simu_covar_coefs[i] );

    for(int i=0;i<simu_a_len;i++)
        sprintf(szTemp, "%s PARAMETER additive    = Pos:%d, effect:%.2f\n",
                        szTemp, simu_a_pos[i], simu_a_effect[i] );

    for(int i=0;i<simu_d_len;i++)
        sprintf(szTemp, "%s PARAMETER dominant    = Pos:%d, effect:%.2f\n",
                        szTemp, simu_d_pos[i], simu_d_effect[i] );

    sprintf(szTemp, "%s PARAMETER sig_p       =%d\n", szTemp, sig_p);
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

void destroy(BLS_par* p)
{
	CFmNewTemp  fmRef;
	p->~BLS_par();
	operator delete(p, fmRef);
}
