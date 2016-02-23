// BLS_par.h: interface for the BLS_par class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(BLS_PAR_H__INCLUDED_)
#define BLS_PAR_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include "bls_R.h"

class BLS_par
{
public:
    BLS_par( CMDOPTIONS* pCmd );
    virtual ~BLS_par();

	int Load(char* szParFile);
    int Summary(char* szOutFile=NULL);

public:
    //how many chromosome(group) data
    int simu_grps;
    int simu_n;

    // the number of total predictors
    int simu_p;
    double simu_snp_rho;
    double simu_snp_missing;

    // pairwise correlation coefficient
    double simu_rho;
    double simu_mu;
    double simu_sigma2;
    int  simu_covar_len;
    double simu_covar_coefs[100];

    //a
    int simu_a_pos[100];
    int simu_a_len;
    double simu_a_effect[100];

    //d
    int simu_d_pos[100];
    int simu_d_len;
    double simu_d_effect[100];

    double simu_t_range[2];
    double simu_covar_range[2];
    int sig_p;

private:
    int Get_SigP();
    void LoadDefault();

    char* m_szParFile;
    CMDOPTIONS* m_pCmd;
};

void destroy( BLS_par* p );


#endif

