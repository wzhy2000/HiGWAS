// gls_par.h: interface for the GLS_par class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(GLS_PAR_H__INCLUDED_)
#define GLS_PAR_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "gls_R.h"

class GLS_par
{
public:
    GLS_par( CMDOPTIONS* pCmd );
    virtual ~GLS_par();

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
    double simu_mu[4];
    double simu_sigma2;

    int simu_covar_len;
    double simu_covar_range[2];
    double simu_covar_effect[100][4];

    //a
    int simu_a_pos[100];
    int simu_a_len;
    double simu_a_effect[100][4];

    //d
    int simu_d_pos[100];
    int simu_d_len;
    double simu_d_effect[100][4];

    //z
    int simu_z_range[2];
    int simu_z_count[2];

    int sig_p;

private:
    int Get_SigP();
    void LoadDefault();

    char* m_szParFile;
    CMDOPTIONS* m_pCmd;
};


void destroy(GLS_par* p);


#endif

