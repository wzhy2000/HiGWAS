/* gls_par.cpp  -	GLS Parameter Object for simulation
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
#include "fm_err.h"
#include "fm_new.h"

#include "gls_dat.h"
#include "gls_cfg.h"
#include "gls_par.h"

GLS_par::GLS_par( CMDOPTIONS* pCmd )
{
    m_pCmd = pCmd;
    m_szParFile = NULL;
    LoadDefault();
}

GLS_par::~GLS_par()
{
    if(m_szParFile)
        Free(m_szParFile);

    _log_info(_HI_, "GLS_par is released successfully.");
}

void GLS_par::LoadDefault()
{
    simu_grps= 2;
    simu_n   = 100;
    simu_p   = 1000;
    simu_rho = 0.4;
    simu_snp_rho = 0.1;

    simu_mu[0]  = 26.79;
    simu_mu[1]  = -6.16;
    simu_mu[2]  = 3.75;
    simu_mu[3]  = -6.39;
    simu_sigma2 = 16;

    simu_covar_len = 1;
    simu_covar_effect[0][0] = 6;
    simu_covar_effect[0][1] = 0.29;
    simu_covar_effect[0][2] = -5.34;
    simu_covar_effect[0][3] = 6.49;

    simu_covar_range[0] = -1;
    simu_covar_range[1] = +1;

    simu_a_len = 3;
    simu_a_pos[0] = 1;
    simu_a_pos[1] = 2;
    simu_a_pos[2] = 3;
    simu_a_effect[0][0] = 2.08;
    simu_a_effect[0][1] = 1.77;
    simu_a_effect[0][2] = -4.11;
    simu_a_effect[0][3] = 1.09;
    simu_a_effect[1][0] = 2.34;
    simu_a_effect[1][1] = -0.4;
    simu_a_effect[1][2] = 1.48;
    simu_a_effect[1][3] = -9.43;
    simu_a_effect[2][0] = 2.80;
    simu_a_effect[2][1] = -4.5;
    simu_a_effect[2][2] = 2.00;
    simu_a_effect[2][3] = 0.00;

    simu_d_len = 3;
    simu_d_pos[0]   = 4;
    simu_d_pos[1]   = 5;
    simu_d_pos[2]   = 6;
    simu_d_effect[0][0] = 2.98;
    simu_d_effect[0][1] = -4.27;
    simu_d_effect[0][2] = 9.64;
    simu_d_effect[0][3] = 2.85;
    simu_d_effect[1][0] = 2.09;
    simu_d_effect[1][1] = 2.64;
    simu_d_effect[1][2] = 3.81;
    simu_d_effect[1][3] = 3.07;
    simu_d_effect[2][0] = 2.53;
    simu_d_effect[2][1] = -2.45;
    simu_d_effect[2][2] = 5.42;
    simu_d_effect[2][3] = -3.92;

    simu_z_range[0] = 30;
    simu_z_range[1] = 80;
    simu_z_count[0] = 5;
    simu_z_count[1] = 12;

    sig_p = Get_SigP();
}

int GLS_par::Get_SigP()
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

bool extract_value(char* optionv, const char* left, double* pValue, int nMaxLen, int* pActlen=NULL)
{
    if (strncmp(optionv, left, strlen(left))!=0)
        return false;

    char szTmp[256] ={0};
    strcpy( szTmp, optionv + strlen(left) );

    char* szRight =NULL;
    if (strncmp(szTmp, "=", 1)==0)
    {
        szRight = szTmp + 1;
    }
    else
    {
        char seps1[] = " \t";
        char* token = NULL;
        token = strtok( szTmp, seps1 );
        if ( token== NULL || strncmp(token, "=", 1)!=0 )
            return false;

        strcpy( szTmp, optionv + strlen(left) );
        szRight = strchr( szTmp, '=');
        szRight++;
    }

    int i=0;
    char seps2[] = " ,\t\n";
    char* token = strtok( szRight, seps2 );
    while( token != NULL && i<nMaxLen)
    {
       /* While there are tokens in "string" */
       pValue[i] = atof(token);
       /* Get next token: */
       token = strtok( NULL, seps2 );
       i++;
    }

    if (pActlen)
        *pActlen = i;
    return (true);
}

int GLS_par::Load(char* szParFile)
{
    FILE* fp = fopen(szParFile, "r");
    if (!fp)
    {
        _log_error(_HI_, "Failed to open the parameter file:%s", szParFile);
        return( ERR_OPEN_FILE );
    }

    _log_info(_HI_, "Start to read the parameter file:%s", szParFile);

    char aLine[256]={0};
    while( !feof( fp ) )
    {
        if( fgets( aLine, 255, fp ) == NULL)
        {
            if ( ferror(fp)!=0 )
            {
                _log_error(_HI_, "Failed to get a line from the file:%s", szParFile);
                return( ERR_OPEN_FILE );
            }
            else
                break;
        }

        _log_debug( _HI_, "Par Info: %s", aLine);

        double tmp;
        if ( extract_value(aLine, "simu_n", &tmp, 1) )
        {
            simu_n = (int)tmp;
            continue;
        }
        if ( extract_value(aLine, "simu_p", &tmp, 1) )
        {
            simu_p = (int)tmp;
            continue;
        }
        if ( extract_value(aLine, "simu_grps", &tmp, 1) )
        {
            simu_grps = (int)tmp;
            continue;
        }

        if ( extract_value(aLine, "simu_rho", &simu_rho, 1) ) continue;
        if ( extract_value(aLine, "simu_snp_rho", &simu_snp_rho, 1) ) continue;
        if ( extract_value(aLine, "simu_sigma2", &simu_sigma2, 1) ) continue;
        if ( extract_value(aLine, "simu_mu", simu_mu, 4) ) continue;

        if ( extract_value(aLine, "simu_a_len", &tmp, 1) )
        {
            simu_a_len = (int)tmp;
            continue;
        }
        if ( extract_value(aLine, "simu_d_len", &tmp, 1) )
        {
            simu_d_len = (int)tmp;
            continue;
        }
        if ( extract_value(aLine, "simu_covar_len", &tmp, 1) )
        {
            simu_covar_len = (int)tmp;
            continue;
        }

        //a
        double a_pos[50]={0};
        if ( extract_value(aLine, "simu_a_pos", a_pos, 50, &simu_a_len) )
        {
            for (int i=0;i<simu_a_len;i++)
            {
                _log_debug(_HI_, "a_pos:%.2f", a_pos[i]);
                simu_a_pos[i] = (int)round( a_pos[i] );
            }
             continue;
        }

        //d
        double d_pos[50]={0};
        if ( extract_value(aLine, "simu_d_pos", d_pos, 50, &simu_d_len) )
        {
            for (int i=0;i<simu_d_len;i++)
            {
                _log_debug(_HI_, "d_pos:%.2f", d_pos[i]);
                simu_d_pos[i] = (int)round( d_pos[i] );
            }
            continue;
        }

        double dTmp[2]={0};
        if ( extract_value(aLine, "simu_z_range", dTmp, 2 ) )
        {
            simu_z_range[0] = (int)round(dTmp[0]);
            simu_z_range[1] = (int)round(dTmp[1]);
            continue;
        }

        if ( extract_value(aLine, "simu_z_count", dTmp, 2 ) )
        {
            simu_z_count[0] = (int)round(dTmp[0]);
            simu_z_count[1] = (int)round(dTmp[1]);
            continue;
        }

        if ( extract_value(aLine, "simu_covar_range", dTmp, 2 ) )
        {
            simu_covar_range[0] = dTmp[0];
            simu_covar_range[1] = dTmp[1];
            continue;
        }

        bool bEffect_a = false;
        for(int i=0; i<simu_a_len; i++)
        {
            char szLeft[64]={0};
            sprintf(szLeft, "simu_a_effect[%d]", i);
            if ( extract_value(aLine, szLeft, (simu_a_effect[i]), 4 ) )
            {
                bEffect_a = true;
                break;
            }
        }

        if(bEffect_a) continue;

        bool bEffect_d = false;
        for(int i=0; i<simu_d_len; i++)
        {
            char szLeft[64]={0};
            sprintf(szLeft, "simu_d_effect[%d]", i);
            if ( extract_value(aLine, szLeft, (simu_d_effect[i]), 4 ) )
            {
                bEffect_d = true;
                break;
            }
        }

        if(bEffect_d) continue;

        bool bCovar_e = false;
        for(int i=0; i<simu_covar_len; i++)
        {
            char szLeft[64]={0};
            sprintf(szLeft, "simu_covar_effect[%d]", i);
            if ( extract_value(aLine, szLeft, (simu_covar_effect[i]), 4 ) )
            {
                bCovar_e = true;
                break;
            }
        }

        if(bCovar_e) continue;


    }
    fclose( fp );

    sig_p = Get_SigP();
    m_szParFile = Strdup(szParFile);

    _log_info(_HI_, "parameter file is loaded successfully.(sig_p:%d)", sig_p);

	return(0);
}

int GLS_par::Summary( char* szOutFile )
{
    char szTemp[MAX_PATH*256] = {0};

    sprintf(szTemp, "%s \nSUMMARY(Parameter)\n", szTemp );
    sprintf(szTemp, "%s -------------------------------------------\n", szTemp );
    sprintf(szTemp, "%s PARAMETER parameter = %s\n", szTemp, m_szParFile );
    sprintf(szTemp, "%s PARAMETER simu_grps   = %d\n", szTemp, simu_grps);
    sprintf(szTemp, "%s PARAMETER simu_n      = %d\n", szTemp, simu_n);
    sprintf(szTemp, "%s PARAMETER simu_p      = %d\n", szTemp, simu_p);
    sprintf(szTemp, "%s PARAMETER simu_snp_rho= %.2f\n", szTemp, simu_snp_rho);
    sprintf(szTemp, "%s PARAMETER simu_rho    = %.2f\n", szTemp, simu_rho);
    sprintf(szTemp, "%s PARAMETER simu_sigma2 = %.2f\n", szTemp, simu_sigma2);
    sprintf(szTemp, "%s PARAMETER simu_z_range= %d, %d\n",
                        szTemp, simu_z_range[0], simu_z_range[1]);
    sprintf(szTemp, "%s PARAMETER simu_z_count= %d, %d\n",
                        szTemp, simu_z_count[0], simu_z_count[1]);
    sprintf(szTemp, "%s PARAMETER simu_mu     = %.2f %.2f %.2f %.2f\n",
                        szTemp, simu_mu[0], simu_mu[1], simu_mu[2], simu_mu[3] );

	for(int i =0; i< simu_covar_len; i++)
    sprintf(szTemp, "%s PARAMETER simu_covar[%d] = %.2f %.2f %.2f %.2f\n",
                        szTemp, i+1, simu_covar_effect[i][0], simu_covar_effect[i][1], simu_covar_effect[i][2], simu_covar_effect[i][3] );

    for(int i=0;i<simu_a_len;i++)
    sprintf(szTemp, "%s PARAMETER additive    = Pos:%d, effect:%.2f %.2f %.2f %.2f\n",
                        szTemp, simu_a_pos[i],
                        simu_a_effect[i][0],
                        simu_a_effect[i][1],
                        simu_a_effect[i][2],
                        simu_a_effect[i][3] );

    for(int i=0;i<simu_d_len;i++)
    sprintf(szTemp, "%s PARAMETER dominant    = Pos:%d, effect:%.2f %.2f %.2f %.2f\n",
                        szTemp, simu_d_pos[i],
                        simu_d_effect[i][0],
                        simu_d_effect[i][1],
                        simu_d_effect[i][2],
                        simu_d_effect[i][3] );

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


void destroy(GLS_par* p)
{
	CFmNewTemp  fmRef;
	p->~GLS_par();
	operator delete(p, fmRef);
}
