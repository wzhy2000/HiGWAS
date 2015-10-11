#ifndef _FM_SIMULATE_H_
#define _FM_SIMULATE_H_

#include <stdio.h>
#include <stdlib.h>

#include "fm_linux.h"

class CFmVectorStr;
class CFmMatrix;
class CFmVector;

typedef struct _simu_param0
{
    //how many chromosome(group) data
    int simu_grps;
    int simu_n;
    // the number of total predictors
    int simu_p;
    double simu_snp_rho;
    double simu_snp_miss;

    // pairwise correlation coefficient
    double simu_rho;
    double simu_mu;
    double simu_sigma2;

    int simu_covar_len;
    double* simu_covar_pcoefs;

    //a
    int simu_a_len;
    int* simu_a_ppos;
    double* simu_a_peffect;

    //d
    int simu_d_len;
    int* simu_d_ppos;
    double* simu_d_peffect;

    //z
    double* simu_t_prange;
    double* simu_covar_prange;
    int sig_p;
} SIMU_PARAM;

typedef struct _simu_param1
{
    //how many chromosome(group) data
    int simu_grps;
    int simu_n;
    // the number of total predictors
    int simu_p;
    double simu_snp_rho;
    double simu_snp_miss;

    // pairwise correlation coefficient
    double simu_rho;
    double simu_mu[4];
    double simu_sigma2;

    int simu_covar_len;
    double simu_covar_effect[100][4];
    double simu_covar_range[2];

    //a
    int simu_a_len;
    int simu_a_pos[100];
    double simu_a_effect[100][4];

    //d
    int simu_d_len;
    int simu_d_pos[100];
    double simu_d_effect[100][4];

    //z
    int simu_z_range[2];
    int simu_z_count[2];

    int sig_p;
} SIMU_PARAM_LONGDT;


class CFmSimulate
{
public:
    CFmSimulate( SIMU_PARAM* par );
    CFmSimulate( SIMU_PARAM_LONGDT* par );

    virtual ~CFmSimulate();

    CFmVectorStr* m_pSnpNames;
    CFmMatrix*  m_pSimuSnps;   //Matrix[P,N] for simulation data

    CFmMatrix* m_pCovar;      //Matrix[N,covar_len];
    CFmMatrix* m_pPhenoY;     //Matrix[N,Q];
    CFmMatrix* m_pCovarZ;     //Matrix[N,Q];
    CFmMatrix* m_pCovarX;     //VECTOR[N];
    CFmVector* m_pZRange;     //VECTOR[Q];

    int m_nSubjN;
    int m_nSnpP;
    int m_nMesuTime;

    int Simulate( char* szSnpoutFile, char* szPheoutFile, char* szGrpId );

private:
    int SavePhenoFile( char* szPheoutFile );
    int SaveSnpFile( char* szSnpoutFile );
    int Simu_geno( double* gen_r );
    int Simu_pheno( double* snp, CFmMatrix* pMat_cov, CFmVector* pVct_Y );
    int Simu_pheno_longdt( double* snp, double* ind_Y_r, double* ind_X_r, double* ind_Z_r, double* vct_Z_sample );

private:
    SIMU_PARAM* m_par;
    SIMU_PARAM_LONGDT* m_par_longdt;

    char m_sPheno_file[MAX_PATH];
    char m_sGeno_file[MAX_PATH];
};

void destroy( CFmSimulate* p);

#endif // _FM_SIMULATE_H_
