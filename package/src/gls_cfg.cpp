/* GLS_cfg.cpp  -	GLS System Object
 *
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */


#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fm_rlogger.h"
#include "fm_err.h"
#include "fm_new.h"
#include "gls_cfg.h"

GLS_cfg::GLS_cfg():CFmSys()
{
	m_nSectId = 1;
	m_nTaskId = 11;
    m_nMcmcHint = 100;

    m_nMcmcIter = 2000;
    m_nMcmcSnps =100000;
    m_fRhoTuning = 0.095;
    m_fBurnInRound = 0.3; //30%
    m_fDetRefit = 0;
    m_fPerRefit = 0.10;
    m_fNumSig   = 0;
    m_fPerSig   = 0.02;
    m_fQval_add = 0.05;
    m_fQval_dom = 0.09;

    _log_debug( _HI_, "CFG File: NULL");
}


GLS_cfg::GLS_cfg(char* szCfgFile):CFmSys()
{
	m_nSectId = 1;
	m_nTaskId = 11;
    m_nMcmcHint = 100;

    m_nMcmcIter = 2000;
    m_nMcmcSnps =100000;
    m_fRhoTuning = 0.095;
    m_fBurnInRound = 0.3; //30%
    m_fDetRefit = 0;
    m_fPerRefit = 0.10;
    m_fNumSig   = 0;
    m_fPerSig   = 0.02;
    m_fQval_add = 0.05;
    m_fQval_dom = 0.09;

    _log_debug( _HI_, "CFG File: %s", szCfgFile);
    if(strlen(szCfgFile)!=0)
        Load(szCfgFile);
}

GLS_cfg::~GLS_cfg()
{
    _log_debug(_HI_, "GLS_cfg is released successfully.");
}

bool extract_value(char* optionv, const char* left, char* pValue)
{
    if (strncmp(optionv, left, strlen(left))!=0)
        return false;

    char szTmp[256] ={0};
    strcpy( szTmp, optionv + strlen(left) );

    char seps1[] = "=";
    char* token = strtok( szTmp, seps1 );
    if (token!=NULL)
    {
        strcpy(pValue, token);
        return (true);
    }

    return (false);
}

int GLS_cfg::Load(char* szCfgFile)
{
    FILE* fp = fopen(szCfgFile, "r");
    if (!fp)
    {
        _log_error(_HI_, "Failed to open the parameter file:%s", szCfgFile);
        return(-1);
    }

    char aLine[256]={0};
    while( !feof( fp ) )
    {
        if( fgets( aLine, 255, fp ) == NULL)
        {
            if ( ferror(fp)!=0 )
            {
                _log_error(_HI_, "Failed to get a line from the file:%s", szCfgFile);
                return(-1);
            }
            else
                break;
        }

        _log_debug( _HI_, "CFG Info: %s", aLine);

        char tmp[256];
        if ( extract_value(aLine, "max_iter", tmp ) )
        {
            m_nMcmcIter = atoi(tmp);
            continue;
        }
        if ( extract_value(aLine, "mcmc_snps", tmp ) )
        {
            m_nMcmcSnps = atoi(tmp);
            continue;
        }
        if ( extract_value(aLine, "rho_tuning", tmp ) )
        {
            m_fRhoTuning = atof(tmp);
            _log_debug( _HI_, "m_fRhoTuning: %.4f", m_fRhoTuning);
            if (m_fRhoTuning<0.01)
                m_fRhoTuning = 0.095;
            continue;
        }
        if ( extract_value(aLine, "refit_criterion", tmp ) )
        {
            m_fPerRefit = 0;
            m_fDetRefit = 0;

            if (strchr(tmp, '%') != NULL)
            {
                char* pos = strchr(tmp, '%');
                *pos = 0;
                m_fPerRefit = atof(tmp);
            }
            else
                m_fDetRefit = atof(tmp);

            continue;
        }
        if ( extract_value(aLine, "q_criterion", tmp ) )
        {
            m_fQval_add = atof(tmp);
            m_fQval_dom = atof(tmp)*1.8;
            continue;
        }

        if ( extract_value(aLine, "sig_criterion", tmp ) )
        {
            m_fPerSig = 0;
            if (strchr(tmp, '%') != NULL)
            {
                char* pos = strchr(tmp, '%');
                *pos = 0;
                m_fPerSig = atof(tmp)/100.0;
            }
            else
                m_fNumSig = atof(tmp);

            continue;
        }

        if ( extract_value(aLine, "burn_in_round", tmp ) )
        {
            m_fBurnInRound = atof(tmp);
            _log_debug( _HI_, "m_fBurnInRound: %f", m_fBurnInRound);
            continue;
        }
    }

    fclose( fp );
    return(0);
}

int GLS_cfg::GetBurnInRound()
{
    return((int)round(m_fBurnInRound*m_nMcmcIter));
}

void destroy(GLS_cfg* p)
{
	CFmNewTemp  fmRef;
	p->~GLS_cfg();
	operator delete(p, fmRef);
}
