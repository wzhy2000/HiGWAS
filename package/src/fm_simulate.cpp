#include <stdio.h>
#include <math.h>

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>
#include <Rmath.h>

#include "fm_linux.h"
#include "fm_simulate.h"
#include "fm_matrix.h"
#include "fm_vector.h"
#include "fm_vector_str.h"
#include "fm_rlogger.h"
#include "fm_err.h"
#include "fm_sys.h"
#include "fm_new.h"


CFmSimulate::CFmSimulate( SIMU_PARAM* par  )
{
    m_par = par;
    m_par_longdt = NULL;

    m_nSubjN = par->simu_n;
    m_nSnpP  = par->simu_p;
    m_nMesuTime = 1;

	CFmNewTemp refNew;
    m_pPhenoY  = new (refNew) CFmMatrix( m_nSubjN, m_nMesuTime );
    m_pCovarX  = new (refNew) CFmMatrix( m_nSubjN, par->simu_covar_len );
    m_pSimuSnps= new (refNew) CFmMatrix( m_nSnpP, m_nSubjN);
    m_pCovarZ  = NULL;
    m_pZRange  = NULL;

    _log_info( _HI_, "%d, %d, %d, %d, %d", par->simu_p, par->simu_n, par->sig_p, par->simu_a_len, par->simu_d_len);

}

CFmSimulate::CFmSimulate( SIMU_PARAM_LONGDT* par  )
{
    m_par = NULL;
    m_par_longdt = par;

    m_nSubjN = par->simu_n;
    m_nSnpP  = par->simu_p;
    m_nMesuTime = par->simu_z_count[1];

	CFmNewTemp refNew;
    m_pPhenoY = new (refNew) CFmMatrix( m_nSubjN, m_nMesuTime );
    m_pCovarZ = new (refNew) CFmMatrix( m_nSubjN, m_nMesuTime );
    m_pCovarX = new (refNew) CFmMatrix( m_nSubjN, par->simu_covar_len );
    m_pSimuSnps= new (refNew) CFmMatrix( m_nSnpP, m_nSubjN);
    m_pZRange = new (refNew) CFmVector( m_nMesuTime, 0.0 );
}

CFmSimulate::~CFmSimulate()
{
	if(m_pPhenoY) destroy( m_pPhenoY );
	if(m_pCovarZ) destroy( m_pCovarZ );
	if(m_pCovarX) destroy( m_pCovarX );
	if(m_pSimuSnps) destroy( m_pSimuSnps );
	if(m_pZRange) destroy( m_pZRange );
}

int CFmSimulate::Simulate(char *szSnpoutFile, char* szPheoutFile, char* szGrpId)
{
    R_set_seed( CFmSys::GetTempId() );

    strcpy( m_sPheno_file, "simu.pheno.LS");
    strcpy( m_sGeno_file, "simu.geno.LS");

    //--generate snp data
    _log_info( _HI_, "Simulation: Generate genotypical data.");
    int ret = Simu_geno( m_pSimuSnps->GetData() );
    if(ret!=0)
        return(ret);

    //--generate phenotype data
    _log_info( _HI_, "Simulation: Generate phenotypical data.");
    if (m_par_longdt)
    {
        ret = Simu_pheno_longdt( m_pSimuSnps->GetData() ,
                        m_pCovarX->GetData(),
                        m_pPhenoY->GetData(),
                        m_pCovarZ->GetData(),
                        m_pZRange->GetData());
    }
    else
    {
		CFmVector vctY(m_nSubjN, 0.0);
        ret = Simu_pheno(m_pSimuSnps->GetData(), m_pCovarX, &vctY );
		m_pPhenoY->SetCol( 0, vctY, (char*)"y" );
    }

    if(ret!=0)
        return(ret);

    //--create a list for the snps' name.
	CFmNewTemp refNew;
    m_pSnpNames = new (refNew) CFmVectorStr(m_nSnpP);
    for (int i=0;i<m_nSnpP;i++)
    {
        char szBuf[64]={0};
        snprintf(szBuf, 63, "G%s-SNP%d", szGrpId, i+1);
        m_pSnpNames->Set(i, szBuf);
        m_pSimuSnps->SetRowName(i, szBuf);
    }

    //--save phenotype file based on -pheout parameter
    if (strlen( szPheoutFile ))
    {
        _log_debug( _HI_, "Start to save phenotype data to %s.", szPheoutFile);
        int ret = SavePhenoFile( szPheoutFile );
        if (!ret)
            _log_info( _HI_, "The phenotype data is saved into %s.", szPheoutFile);
        else
            _log_info( _HI_, "Failed to save phenotype data into %s.", szPheoutFile);

        strcpy( m_sPheno_file, szPheoutFile);
    }

    //--save phenotype file based on -snpout parameter
    if (strlen(szSnpoutFile))
    {
        _log_debug( _HI_, "Start to save genotype data to %s.", szSnpoutFile);
        int ret = SaveSnpFile( szSnpoutFile  );
        if (!ret)
            _log_info( _HI_, "The genotype data is saved into %s.", szSnpoutFile);
        else
            _log_info( _HI_, "Failed to save genotype data into %s.", szSnpoutFile);

        strcpy( m_sGeno_file, szSnpoutFile );

    }

    return(0);
}

int CFmSimulate::Simu_geno( double* gen_r )
{
    int simu_p, simu_n, sig_p, simu_a_len, simu_d_len;
    int* simu_a_pos, *simu_d_pos;
    double simu_snprho, simu_snpmiss;

    if (m_par)
    {
        simu_p     = m_par->simu_p;
        simu_n     = m_par->simu_n;
        sig_p      = m_par->sig_p;
        simu_a_len = m_par->simu_a_len;
        simu_d_len = m_par->simu_d_len;
        simu_a_pos = m_par->simu_a_ppos;
        simu_d_pos = m_par->simu_d_ppos;
        simu_snprho= m_par->simu_snp_rho;
        simu_snpmiss= m_par->simu_snp_miss;

        _log_info( _HI_, "%d, %d, %d, %d, %d", simu_p, simu_n, sig_p, simu_a_len, simu_d_len);
    }
    else
    {
        simu_p     = m_par_longdt->simu_p;
        simu_n     = m_par_longdt->simu_n;
        sig_p      = m_par_longdt->sig_p;
        simu_a_len = m_par_longdt->simu_a_len;
        simu_d_len = m_par_longdt->simu_d_len;
        simu_a_pos = m_par_longdt->simu_a_pos;
        simu_d_pos = m_par_longdt->simu_d_pos;
        simu_snprho= m_par_longdt->simu_snp_rho;
        simu_snpmiss= m_par_longdt->simu_snp_miss;
    }

    //R:-----------------------------------
    //R: c <- qnorm(0.25);
    double c = qnorm(0.25, 0, 1, TRUE, FALSE);

    //R: sigma_x <- matrix( seq(par$simu_rho, par$simu_rho, length=par$sig_p*par$sig_p), ncol=par$sig_p);
    //   sigma_x <- sigma_x + diag(seq(1-par$simu_rho, 1-par$simu_rho, length=par$sig_p ));
    int ncol = sig_p;
    SEXP sigma_x;
    PROTECT( sigma_x = allocMatrix( REALSXP, ncol, ncol ) );
    double* sigma_x_r = REAL( sigma_x);
    for (int i=0; i<ncol; i++)
    for (int j=0; j<ncol; j++)
    {
        if (i==j)
            sigma_x_r [MI(ncol, ncol, i, j)] = 1;
        else
            sigma_x_r [MI(ncol, ncol, i, j)] = simu_snprho;
    }

    //R:-----------------------------------
    //R: x <- rmvnorm ( n=par$simu_n, mean=seq(0,0,length=par$sig_p), sigma=sigma_x, method="chol");
    SEXP mean;
    PROTECT(mean=allocVector( REALSXP, ncol ) );
    double* mean_n = REAL( mean);
    for (int i=0; i<ncol; i++)
        mean_n[i] = 0;

    SEXP x;
    GetRNGstate();
    PROTECT( x = rmvnorm_chol( simu_n, mean, sigma_x ) );
    PutRNGstate();
    double* x_r = REAL(x);

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

    //R:-----------------------------------
    //R: sig_pos <- unique( c(par$simu_a_posi, par$simu_d_posi));
    SEXP a_pos, d_pos, sig_pos0, e1;
    PROTECT( a_pos = allocVector(INTSXP, simu_a_len));
    int* a_pos_r = INTEGER(a_pos);
    for (int i=0; i<simu_a_len; i++)
         a_pos_r[i] = simu_a_pos[i]-1;

    PROTECT( d_pos = allocVector(INTSXP, simu_d_len));
    int* d_pos_r = INTEGER(d_pos);
    for (int i=0; i<simu_d_len; i++)
         d_pos_r[i] = simu_d_pos[i]-1;

    int errorOccurred;
    PROTECT(e1 = lang3( install("c"), a_pos, d_pos ) );
    PROTECT(sig_pos0 = R_tryEval(e1, R_GlobalEnv, &errorOccurred) );

    SEXP sigpos, e2;
    PROTECT(e2 = lang2( install("unique"), sig_pos0 ) );
    PROTECT(sigpos = R_tryEval(e2, R_GlobalEnv, &errorOccurred) );
    int* sigpos_r = INTEGER(sigpos);

    //R:-----------------------------------
    //R: all_SNPs <-c();
    //   for (i in 1:par$simu_p)
    //   {
    //		fi <- which(sig_pos==i);
    //		if (length(fi)>0)
    //			all_SNPs <- cbind(all_SNPs, x[,fi[[1]]]);
    //		else
    //			all_SNPs <- cbind(all_SNPs, rnorm(par$simu_n, mean=0, sd=1) );
    //	}
    SEXP all_snps;
    PROTECT( all_snps = allocMatrix( REALSXP, simu_n, simu_p ) );
    double* all_snps_r = REAL( all_snps );

    GetRNGstate();
    for (int j=0; j<simu_p; j++)
        for (int i=0; i<simu_n; i++)
        {
            bool bFound = false;
            for(int k=0; k<5; k++)
            {
                if (j==sigpos_r[k])
                {
                    bFound=true;
                    all_snps_r [ MI(simu_n, simu_p, i, j) ] = x_r[  MI(simu_n, ncol, i, k)];
                }
            }
            if (!bFound)
                all_snps_r [MI(simu_n, simu_p, i, j) ] = rnorm( 0, 1 );
        }
    PutRNGstate();

    //R:-----------------------------------
    //R: for (i in 1:par$simu_n)
    //   {
    //		for(j in 1:par$simu_p)
    //		{
    //			tmp <- all_SNPs[i,j];
    //			all_SNPs[i,j] <- -1*(tmp<c) + (tmp>-c);
    //		}
    //	}
    SEXP ans;
    PROTECT( ans = allocMatrix( INTSXP, simu_n, simu_p ) );
    int* ans_r = INTEGER( ans );
    for (int i=0; i<simu_n; i++)
        for (int j=0; j<simu_p; j++)
            ans_r [ MI(simu_n, simu_p, i, j) ] = -1*(all_snps_r[MI(simu_n, simu_p, i, j)] < c ) + ( all_snps_r[MI(simu_n, simu_p, i, j)]>-c);

    //R:-----------------------------------
    //R:colnames(ans) <- paste("snp1", c(1:par$simu_p), sep="");
    //R:rownames(ans) <- paste("id", c(1:par$simu_n), sep="");
    SEXP cols, rows, dimnames;
    PROTECT(cols = allocVector(STRSXP, simu_p ));
    for (int i=0; i<simu_p; i++)
    {
        char buf[256];
        sprintf(buf, "%s-%d", "snp", i+1);
        SET_STRING_ELT(cols, i, mkChar( buf ));
    }

    PROTECT(rows = allocVector(STRSXP, simu_n ));
    for (int i=0; i<simu_n; i++)
    {
        char buf[256];
        sprintf(buf, "%s-%d", "id", i+1);
        SET_STRING_ELT(rows, i, mkChar( buf ));
    }

    PROTECT(dimnames = allocVector(VECSXP, 2));
    // setrownames
    SET_VECTOR_ELT(dimnames, 0, rows);
    // setcolnames
    SET_VECTOR_ELT(dimnames, 1, cols);
    setAttrib(ans, R_DimNamesSymbol, dimnames);

    //R:-----------------------------------
    //R: ret <- t(ans);
    SEXP ret;
    PROTECT(e1 = lang2(install("t.default"), ans ) );
    PROTECT( ret = R_tryEval(e1, R_GlobalEnv, &errorOccurred) );
    int* ret_r = INTEGER(ret);

	//Add Missing SNP
	if(simu_snpmiss!=0)
    for (int i=0; i<simu_p; i++)
    {
        GetRNGstate();
		int n_miss = round( runif( 0, simu_n * simu_snpmiss) );
		for (int k=0; k<n_miss ; k++)
		{
			int j = round( runif(1, simu_n) ) - 1;
        	ret_r[ MI(simu_p, simu_n, i, j) ] = -9;
		}
        PutRNGstate();
	}

    //R:-----------------------------------
    //gen_r <- ret;
    for (int i=0; i<simu_p; i++)
        for (int j=0; j<simu_n ; j++)
            gen_r[ MI(simu_p, simu_n, i, j) ] = ret_r[ MI(simu_p, simu_n, i, j) ];

    UNPROTECT(18);

    return(0);
}


int CFmSimulate::Simu_pheno_longdt( double* snp, double* vct_X_r, double* mat_Y_r, double* mat_Z_r, double* vct_Z_sample )
{
    int simu_p = m_par_longdt->simu_p;
    int simu_n = m_par_longdt->simu_n;
    int simu_a_len = m_par_longdt->simu_a_len;
    int simu_d_len = m_par_longdt->simu_d_len;
    int* simu_a_pos = m_par_longdt->simu_a_pos;
    int* simu_d_pos = m_par_longdt->simu_d_pos;
    double simu_rho = m_par_longdt->simu_rho;
    double simu_sigma2 = m_par_longdt->simu_sigma2;
    double* simu_mu = m_par_longdt->simu_mu;
    int* simu_z_range = m_par_longdt->simu_z_range;
    int* simu_z_count = m_par_longdt->simu_z_count;
    int simu_covar_len = m_par_longdt->simu_covar_len;
    double* simu_covar_range = m_par_longdt->simu_covar_range;

    _log_debug(_HI_, "P1, simu_p:%d, simu_n:%d.", simu_p, simu_n);

    //R:-----------------------------------
    //R: gen_a <- t(snp);
    SEXP gen_a;
    PROTECT( gen_a = allocMatrix( INTSXP, simu_n, simu_p ) );
    int* gen_a_r = INTEGER(gen_a);
    for (int i=0; i<simu_n; i++)
        for (int j=0; j<simu_p; j++)
            gen_a_r[ MI(simu_n, simu_p, i, j) ] = (int)snp[ MI(simu_p, simu_n,  j, i)];

    //R:-----------------------------------
    //R: gen_d <- 1-abs(gen_a);
    SEXP gen_d;
    PROTECT( gen_d = allocMatrix( INTSXP, simu_n, simu_p ) );
    int* gen_d_r = INTEGER( gen_d);
    for (int i=0; i<simu_n; i++)
        for (int j=0; j<simu_p; j++)
            gen_d_r [MI(simu_n, simu_p, i, j) ] = 1-abs(gen_a_r [MI(simu_n, simu_p, i, j) ]);

    //R:-----------------------------------
    //R: add_mat<-array(0, dim=c(par$simu_p, length(par$simu_a_effect[1,]));
    //   for (i in 1:length(par$simu_a_pos) )
    //		add_mat[ par$simu_a_pos[i], ] <- par$simu_a_effect[i];
    SEXP add_mat;
    PROTECT( add_mat = allocMatrix( REALSXP, simu_p, 4 ) );
    double* add_mat_r = REAL( add_mat );
    for (int k=0; k<simu_p*4 ; k++) add_mat_r[k]=0;
    for (int n=0; n<simu_a_len; n++)
    for (int j=0; j<4; j++)
    {
        int i = simu_a_pos[n]-1;
        add_mat_r [ MI(simu_p, 4, i, j) ] = m_par_longdt->simu_a_effect[n][j];
    }

    //R:-----------------------------------
    //R: dom_mat<-array(0, dim=c(par$simu_p, length(par$simu_d_effect[1,]));
    //   for (i in 1:length(par$simu_d_pos) )
    //		dom_mat[ par$simu_d_pos[i], ] <- par$simu_d_effect[i];
    SEXP dom_mat;
    PROTECT( dom_mat = allocMatrix( REALSXP, simu_p, 4 ) );
    double* dom_mat_r = REAL( dom_mat );
    for (int k=0; k<simu_p*4 ; k++) dom_mat_r[k]=0;
    for (int n=0; n<simu_d_len; n++)
    for (int j=0; j<4; j++)
    {
        int i = simu_d_pos[n]-1;
        dom_mat_r [MI(simu_p, 4, i, j)] = m_par_longdt->simu_d_effect[n][j];
    }

    _log_debug(_HI_, "P2, simu_z_count:%d,%d.", simu_z_count[0], simu_z_count[1]);

    //R:-----------------------------------
    //R: ind.mesu <- round(runif(par$simu_n, min=par$simu_z_count[1], max=par$simu_z_count[2]) );
    SEXP ind_mesu;
    PROTECT( ind_mesu = allocVector( REALSXP, simu_n ) );
    double* ind_mesu_r = REAL(ind_mesu);
    GetRNGstate();
    for (int i=0; i<simu_n; i++)
        ind_mesu_r[i] = round(runif( simu_z_count[0], simu_z_count[1] ) );
    PutRNGstate();

    _log_debug(_HI_, "P3, ind.Z[%d,%d]", simu_n, simu_z_count[1]);

    //R:-----------------------------------
    //R: ind.Z <- array(0, dim=c(par$simu_n, simu_z_count[2]) );
    SEXP ind_Z;
    PROTECT( ind_Z = allocMatrix( REALSXP, simu_n, simu_z_count[1] ) );
    double* ind_Z_r = REAL(ind_Z);

    //R:-----------------------------------
    //R: for(i in 1:par$simu_n)
    //		ind.Z[i,1:ind.mesu[i]] <- sort( sample(par$simu_age[1]:par$simu_age[2], ind.mesu[i]) );
    bool bAlikeReal=false;
    SEXP z_range;
    if (bAlikeReal)
    {
        PROTECT( z_range = allocVector( INTSXP, simu_z_count[1] ) );
        int* z_range_r = INTEGER(z_range);

        //R:-----------------------------------
        //R: c(runif(0,1)*(par$simu_age[2] - par$simu_age[1]) + par$simu_age[1];
        SEXP ages;
        PROTECT( ages = allocVector( INTSXP, simu_z_range[1]- simu_z_range[0]-1 ) );
        int* ages_r = INTEGER(ages);
        for (int j=0; j<simu_z_range[1] - simu_z_range[0]-1; j++)
            ages_r[j] = simu_z_range[0] + j+1;

        SEXP size;
        PROTECT(size = allocVector(INTSXP, 1));
        INTEGER(size)[0] = simu_z_count[1]-2;

        SEXP replace;
        PROTECT(replace = allocVector(LGLSXP, 1));
        INTEGER(replace)[0] = FALSE;

        SEXP prob = R_NilValue;

        int errorOccurred;
        SEXP sample_ret, e1;
        PROTECT(e1 = lang5( install("sample"), ages, size, replace, prob ) );
        PROTECT(sample_ret = R_tryEval(e1, R_GlobalEnv, &errorOccurred) );

        SEXP sample_mm;
        PROTECT(sample_mm = allocVector(INTSXP, simu_z_count[1]));
        for (int j=0;j<simu_z_count[1]-2;j++)
            INTEGER(sample_mm)[j] = INTEGER(sample_ret)[j];
        INTEGER(sample_mm)[ simu_z_count[1]-2 ] = simu_z_range[1];
        INTEGER(sample_mm)[ simu_z_count[1]-1 ] = simu_z_range[0];

        SEXP decreasing;
        PROTECT(decreasing = allocVector(LGLSXP, 1));
        INTEGER(decreasing)[0] = FALSE;

        SEXP e2, sort_ret;
        PROTECT(e2 = lang3( install("sort"), sample_mm,  decreasing) );
        PROTECT(sort_ret = R_tryEval(e2, R_GlobalEnv, &errorOccurred) );

        for (int j=0; j<simu_z_count[1]; j++)
            z_range_r[j] = INTEGER(sort_ret)[j];

        UNPROTECT(10);
    }
    else
    {
        //R:-----------------------------------
        //R: c(par$simu_age[1]:par$simu_age[2]);
        PROTECT( z_range = allocVector( INTSXP, simu_z_range[1]- simu_z_range[0]+1 ) );
        int* z_range_r = INTEGER(z_range);
        for (int j=0; j<simu_z_range[1] - simu_z_range[0]+1; j++)
            z_range_r[j] = simu_z_range[0] + j;
        UNPROTECT(1);
    }

    for (int i=0; i<LENGTH(z_range); i++)
        vct_Z_sample[i] = (double)(INTEGER(z_range)[i]);

    for (int i=0; i<simu_n; i++)
    {
        SEXP size;
        PROTECT(size = allocVector(INTSXP, 1));
        INTEGER(size)[0] = (int)ind_mesu_r[i];

        SEXP replace;
        PROTECT(replace = allocVector(LGLSXP, 1));
        INTEGER(replace)[0] = FALSE;

        SEXP prob = R_NilValue;

        int errorOccurred;
        SEXP age_sample, e1;
        PROTECT(e1 = lang5( install("sample"), z_range, size, replace, prob ) );
        PROTECT(age_sample = R_tryEval(e1, R_GlobalEnv, &errorOccurred) );

        SEXP decreasing;
        PROTECT(decreasing = allocVector(LGLSXP, 1));
        INTEGER(decreasing)[0] = FALSE;

        SEXP age_sort;
        PROTECT(e1 = lang3( install("sort"), age_sample,  decreasing) );
        PROTECT(age_sort = R_tryEval(e1, R_GlobalEnv, &errorOccurred) );
        int* age_sort_r = INTEGER(age_sort);

        for (int j=0; j<simu_z_count[1]; j++)
            if (j<ind_mesu_r[i])
                ind_Z_r[MI(simu_n, simu_z_count[1], i, j) ] = age_sort_r[j];
            else
                ind_Z_r[MI(simu_n, simu_z_count[1], i, j) ] = R_NaN;

        UNPROTECT(7);
    }

    _log_debug(_HI_, "P4, ind.X[%d]", simu_n);

    //R:-----------------------------------
    //R: ind.X <- round(runif(par$simu_n, min=1, max=2));
    SEXP ind_X;
    PROTECT( ind_X = allocMatrix( REALSXP, simu_n, simu_covar_len ) );
    double* ind_X_r = REAL(ind_X);
    GetRNGstate();
    for (int nCov=0; nCov<simu_covar_len; nCov++)
    for (int i=0; i<simu_n; i++)
    {
		if ( nCov==0)
			ind_X_r[MI(simu_n, simu_covar_len, i, nCov)] = round(runif( 0, 1) )*1.0;
        else
            ind_X_r[MI(simu_n, simu_covar_len, i, nCov)] = runif( simu_covar_range[0], simu_covar_range[1] ) ;
    }
    PutRNGstate();

    //R:-----------------------------------
    //R: ind.y<-c()
    SEXP ind_Y;
    PROTECT( ind_Y = allocMatrix( REALSXP, simu_n, simu_z_count[1]) );
    double* ind_Y_r = REAL(ind_Y);

    for (int i=0; i<simu_n; i++)
    {
        //R:-----------------------------------
        //R: tmp_age <-ind.Z[i,];
        int nozero_len =0;
        for (int j=0; j<simu_z_count[1]; j++)
            if (ind_Z_r[MI( simu_n, simu_z_count[1], i, j) ]>0 )
                 nozero_len++;

        SEXP tmp_age;
        PROTECT( tmp_age = allocVector( INTSXP, nozero_len ) );
        int* tmp_age_r = INTEGER(tmp_age);
        int k=0;
        for (int j=0; j<simu_z_count[1]; j++)
            if (ind_Z_r[MI( simu_n, simu_z_count[1], i, j) ]>0 )
            {
                tmp_age_r[k] = (int)ind_Z_r[MI( simu_n, simu_z_count[1], i, j) ];
                k++;
            }

        SEXP tp;
        PROTECT( tp = allocVector( REALSXP, nozero_len ) );
        double* tp_r = REAL(tp );
        for (int j=0; j<nozero_len; j++)
            tp_r[j] = - 1.0 + (2.0* (tmp_age_r[j]-simu_z_range[0]) )/(simu_z_range[1]-simu_z_range[0]);

        SEXP ui;
        PROTECT( ui = allocMatrix( REALSXP, nozero_len, 4  ) );
        double* ui_r = REAL(ui );
        for (int nRow=0; nRow<nozero_len; nRow++)
        {
            ui_r[MI( nozero_len, 4, nRow, 0)]	 = 1;
            ui_r[MI( nozero_len, 4, nRow, 1)]	 = tp_r[nRow];
            ui_r[MI( nozero_len, 4, nRow, 2)]	 = (3*pow(tp_r[nRow],2)-1)/2;
            ui_r[MI( nozero_len, 4, nRow, 3)]	 = (5*pow(tp_r[nRow],3)-3*tp_r[nRow])/2;
        }

        _log_debug(_HI_, "P5, gen_effect,%d,[%d]", i,nozero_len);

        SEXP gen_effect;
        PROTECT( gen_effect = allocVector( REALSXP, nozero_len) );
        double* gen_effect_r = REAL(gen_effect);

        SEXP tmp_prod;
        PROTECT( tmp_prod = allocVector( REALSXP, nozero_len) );
        double* tmp_prod_r = REAL( tmp_prod );

        matprod(ui_r, nozero_len, 4, 'N', simu_mu, 4, 1, 'N', gen_effect_r);

		for( int nCov=0; nCov<simu_covar_len; nCov++)
		{
			matprod(ui_r, nozero_len, 4, 'N', m_par_longdt->simu_covar_effect[nCov], 4, 1, 'N', tmp_prod_r);
			for (int k=0; k< nozero_len; k++)
				 gen_effect_r[k] = gen_effect_r[k]+	tmp_prod_r[k]*ind_X_r[ MI(simu_n, simu_covar_len, i, nCov) ];
		}

        for (int j=0; j< simu_p; j++)
        {
            if (gen_a_r[ MI(simu_n, simu_p, i, j)]!=0)
            {
                double add_mat_rj[4]={0};
                add_mat_rj[0]  = add_mat_r[ MI(simu_p, 4, j, 0) ];
                add_mat_rj[1]  = add_mat_r[ MI(simu_p, 4, j, 1) ];
                add_mat_rj[2]  = add_mat_r[ MI(simu_p, 4, j, 2) ];
                add_mat_rj[3]  = add_mat_r[ MI(simu_p, 4, j, 3) ];

                matprod(ui_r, nozero_len, 4, 'N', add_mat_rj, 4, 1, 'N', tmp_prod_r);

                for (int k=0; k< nozero_len; k++)
                    gen_effect_r[k] = gen_effect_r[k]+	tmp_prod_r[k]*gen_a_r[MI(simu_n, simu_p, i, j)];
            }


            if (gen_d_r[MI(simu_n, simu_p, i, j)]!=0)
            {
                double dom_mat_rj[4]={0};
                dom_mat_rj[0]  = dom_mat_r[MI(simu_p, 4, j, 0)];
                dom_mat_rj[1]  = dom_mat_r[MI(simu_p, 4, j, 1)];
                dom_mat_rj[2]  = dom_mat_r[MI(simu_p, 4, j, 2)];
                dom_mat_rj[3]  = dom_mat_r[MI(simu_p, 4, j, 3)];

                matprod(ui_r, nozero_len, 4, 'N', dom_mat_rj, 4, 1, 'N', tmp_prod_r);
                for (int k=0; k< nozero_len; k++)
                    gen_effect_r[k] = gen_effect_r[k]+	tmp_prod_r[k]*gen_d_r[MI(simu_n, simu_p, i, j)];
            }

        }

        _log_debug(_HI_, "P6, sigma,%d,[%d,%d]", i,nozero_len,nozero_len);

        SEXP sigma;
        PROTECT( sigma = allocMatrix( REALSXP, nozero_len, nozero_len) );
        double* sigma_r = REAL(sigma);
        for (int ii=0; ii<nozero_len; ii++)
            for (int jj=0; jj<nozero_len; jj++)
                sigma_r[ MI(nozero_len, nozero_len, ii, jj) ] = simu_sigma2 * pow( simu_rho, abs(tmp_age_r[ii]-tmp_age_r[jj] ) );

        SEXP mu;
        PROTECT( mu = allocVector( REALSXP, nozero_len) );
        double* mu_r = REAL(mu);
        for (int j=0; j<nozero_len; j++) mu_r[j] = 0;

        SEXP rdata;
        PROTECT( rdata = allocVector( REALSXP, nozero_len) );
        double* rdata_r = REAL(rdata);
        for (int k=0;k<nozero_len;k++) rdata_r[k]=0;

        _log_debug(_HI_, "P7, ind.Y,%d,[%d]", i,nozero_len );

        GetRNGstate();
        int ret = rmultnorm( 1, mu_r, sigma_r, nozero_len, rdata_r);
        PutRNGstate();
        if (ret!=0)
            return(ret);

        //R:-----------------------------------
        //R: yi<-t(gen_effect) + rmultnorm..
        for(int k=0; k<nozero_len; k++)
            gen_effect_r[k] += rdata_r[k];

        for (int j=0; j<simu_z_count[1];j++)
            if (j< nozero_len )
                ind_Y_r[ MI(simu_n, simu_z_count[1], i, j)] = gen_effect_r[j];
            else
                ind_Y_r[ MI(simu_n, simu_z_count[1], i, j)] = R_NaN;

        UNPROTECT(8);


        _log_debug(_HI_, "P8, %d", i );
    }

    for(int i= 0; i<simu_n * simu_z_count[1]; i++)	mat_Y_r[i] = ind_Y_r[i];
    for(int i= 0; i<simu_n * simu_z_count[1]; i++)	mat_Z_r[i] = ind_Z_r[i];
    for(int i= 0; i<simu_n * simu_covar_len; i++)	vct_X_r[i] = ind_X_r[i];

    UNPROTECT(8);

    return(0);
}

int CFmSimulate::Simu_pheno( double* snp, CFmMatrix* pMat_cov, CFmVector* pVct_Y )
{
    int simu_p = m_par->simu_p;
    int simu_n = m_par->simu_n;
    int simu_a_len = m_par->simu_a_len;
    int simu_d_len = m_par->simu_d_len;
    int simu_cov_len = m_par->simu_covar_len;
    int* simu_a_pos = m_par->simu_a_ppos;
    int* simu_d_pos = m_par->simu_d_ppos;
    double simu_sigma2 = m_par->simu_sigma2;
    double simu_mu = m_par->simu_mu;
    double* simu_covar_pcoefs = m_par->simu_covar_pcoefs;
    double* simu_t_range = m_par->simu_t_prange;
    double* simu_covar_range = m_par->simu_covar_prange;

    _log_debug(_HI_, "P1, simu_p:%d, simu_n:%d, sim_a_len:%d, sim_d_len:%d.", simu_p, simu_n, simu_a_len, simu_d_len );

    //R:-----------------------------------
    //R: gen_a <- t(snp);
    SEXP gen_a;
    PROTECT( gen_a = allocMatrix( INTSXP, simu_n, simu_p ) );
    int* gen_a_r = INTEGER(gen_a);
    for (int i=0; i<simu_n; i++)
        for (int j=0; j<simu_p; j++)
            gen_a_r[ MI(simu_n, simu_p, i, j) ] = (int)snp[ MI(simu_p, simu_n,  j, i)];

    //R:-----------------------------------
    //R: gen_d <- 1-abs(gen_a);
    SEXP gen_d;
    PROTECT( gen_d = allocMatrix( INTSXP, simu_n, simu_p ) );
    int* gen_d_r = INTEGER( gen_d);
    for (int i=0; i<simu_n; i++)
        for (int j=0; j<simu_p; j++)
            gen_d_r [MI(simu_n, simu_p, i, j) ] = 1-abs(gen_a_r [MI(simu_n, simu_p, i, j) ]);

    //R:-----------------------------------
    //R: add_mat<-array(0, dim=c(par$simu_p, length(par$simu_a_effect[1,]));
    //   for (i in 1:length(par$simu_a_pos) )
    //		add_mat[ par$simu_a_pos[i], ] <- par$simu_a_effect[i];
    SEXP add_vct;
    PROTECT( add_vct = allocVector( REALSXP, simu_p  ) );
    double* add_vct_r = REAL( add_vct );
    for (int k=0; k<simu_p ; k++) add_vct_r[k]=0;
    for (int n=0; n<simu_a_len; n++)
    {
        int i = simu_a_pos[n]-1;
        add_vct_r [ i ] = m_par->simu_a_peffect[n];
    }

    //R:-----------------------------------
    //R: dom_mat<-array(0, dim=c(par$simu_p, length(par$simu_d_effect[1,]));
    //   for (i in 1:length(par$simu_d_pos) )
    //		dom_mat[ par$simu_d_pos[i], ] <- par$simu_d_effect[i];
    SEXP dom_vct;
    PROTECT( dom_vct = allocVector( REALSXP, simu_p ) );
    double* dom_vct_r = REAL( dom_vct );
    for (int k=0; k<simu_p ; k++) dom_vct_r[k]=0;
    for (int n=0; n<simu_d_len; n++)
    {
        int i = simu_d_pos[n]-1;
        dom_vct_r [i] = m_par->simu_d_peffect[n];
    }

    GetRNGstate();
    for (int k=0; k<simu_cov_len; k++)
    {
	    for (int i=0; i<simu_n; i++)
        {
			if (k==0)
				pMat_cov->Set( i, k, round( runif( 1, 2 ) ) - 1 );
			else
				pMat_cov->Set( i, k, runif( -1, 1 ) );
		}
	}

    //R:-----------------------------------
    //R:simu_Y <- rnorm(par$simu_n, mean=0, sd=sqrt(par$simu_sigma2) )
    //  + par$simu_X*par$simu_alpha
    //  + par$simu_Z*par$simu_beta
    //  + Xi%*%additive
    //  + Zeta%*%dominant
    //  + par$simu_mu;

    CFmVector vct_Y( simu_n, simu_mu );
    for (int k=0; k<simu_cov_len; k++)
	    vct_Y = vct_Y + pMat_cov->GetCol(k) * simu_covar_pcoefs[k];

    for (int i=0; i< simu_n; i++)
    {
        double add_eff = 0;
        for (int j=0; j<simu_p; j++)
            add_eff += gen_a_r[ MI(simu_n, simu_p, i, j)]* add_vct_r[j];

        double dom_eff = 0;
        for (int j=0; j<simu_p; j++)
            dom_eff += gen_d_r[ MI(simu_n, simu_p, i, j)]* dom_vct_r[j];

        double e = rnorm( 0, sqrt( simu_sigma2 ) );
        vct_Y.Set(i, vct_Y.Get(i) + add_eff + dom_eff + e );
    }

    PutRNGstate();

 	*pVct_Y = vct_Y;

    UNPROTECT(4);

    return(0);
}

int CFmSimulate::SaveSnpFile( char* szSnpoutFile )
{
    int nGrps = 1;
    if (m_par)
        nGrps = m_par->simu_grps;
    else
        nGrps = m_par_longdt->simu_grps;

    FILE* fp = fopen(szSnpoutFile, "wt");
    if (!fp)
    {
        _log_error(_HI_, "SaveGenoFile: Failed to create genotype file(%s)", szSnpoutFile);
        return( ERR_CREATE_FILE );
    }

    fprintf(fp, "CHR,POS");
    for(int i=0; i<m_pSimuSnps->GetNumCols(); i++)
    {
        fprintf(fp, ",Sub%d", i+1);
    }
    fprintf(fp, "\n");

    for(int k=0; k<nGrps; k++)
    {
        _log_debug(_HI_, "SaveGenoFile: Start to output genotype data into the file%s.(Group:%d)", szSnpoutFile, k+1);

        int nSnpWidth = m_pSimuSnps->GetNumRows()/nGrps;
        int nStartSNP = k*nSnpWidth;
        int nStopSNP = (k+1)*nSnpWidth;
        if (k==nGrps-1)
            nStopSNP = m_pSimuSnps->GetNumRows();

        for(int i=nStartSNP; i < nStopSNP; i++)
        {
            fprintf(fp, "%s", m_pSimuSnps->GetRowName(i) );
            // Chrmosome is from 1 to ...
            fprintf(fp, ",%d,%d", k+1, i+1 );

            for(int j=0; j < m_pSimuSnps->GetNumCols(); j++)
            {
            	if ( m_pSimuSnps->Get(i,j) ==-9 )
            		fprintf(fp, ",NA");
				else
               		fprintf(fp, ",%d", (int)(m_pSimuSnps->Get(i,j)+1) );
			}
            fprintf(fp, "\n");
        }

        _log_debug(_HI_, "SaveGenoFile: The genotype data is saved into the file(%s)", szSnpoutFile);
    }

    fclose(fp);

    return(0);
}


int CFmSimulate::SavePhenoFile( char* szPheoutFile )
{
    FILE* fp = fopen(szPheoutFile, "wt");
    if (!fp)
    {
        _log_error(_HI_, "SavePhenoFile: Failed to create phenotype file(%s)", szPheoutFile);
        return( ERR_CREATE_FILE );
    }

    _log_debug(_HI_, "SavePhenoFile: Start to write phenotype data into the file(%s)", szPheoutFile);

    fprintf(fp, "ID");

    if ( m_pCovarX->GetNumCols()==1)
	    fprintf(fp, ",X");

    if ( m_pCovarX->GetNumCols()>1)
    {
		fprintf(fp, ",X_1");
        for(int i=1; i<m_pCovarX->GetNumCols(); i++)
        {
            fprintf(fp, ",X_%d", i+1);
        }
	}

	if(m_pCovarZ)
	{
        if ( m_pCovarZ->GetNumCols()==1)
            fprintf(fp, ",Z");
        else if ( m_pCovarZ->GetNumCols()>1)
            for(int i=0; i<m_pCovarZ->GetNumCols(); i++)
            {
                fprintf(fp, ",Z_%d", i+1);
            }
	}

    if ( m_pPhenoY->GetNumCols()==1)
        fprintf(fp, ",Y");
    else if ( m_pPhenoY->GetNumCols()>1)
        for(int i=0; i<m_pPhenoY->GetNumCols(); i++)
        {
            fprintf(fp, ",Y_%d", i+1);
        }
    fprintf(fp, "\n");

    for(int i=0; i<m_pPhenoY->GetNumRows(); i++)
    {
        fprintf(fp, "Sub%d", i+1);

        for(int j=0; j<m_pCovarX->GetNumCols(); j++)
        {
			fprintf(fp, ",%.3f", m_pCovarX->Get(i,j));
        }

	    if(m_pCovarZ)
	    {
            for(int j=0; j<m_pCovarZ->GetNumCols(); j++)
            {
                if ( !isnan(m_pCovarZ->Get(i,j)))
                {
                    fprintf(fp, ",%.3f", m_pCovarZ->Get(i,j));
                }
                else
                {
                    fprintf(fp, ", ");
                }
            }
	    }

        for(int j=0; j<m_pPhenoY->GetNumCols(); j++)
        {
            if (!isnan(m_pPhenoY->Get(i,j)))
            {
                fprintf(fp, ",%.3f", m_pPhenoY->Get(i,j));
            }
            else
            {
                fprintf(fp, ", ");
            }
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\n");
    fclose(fp);

    _log_debug(_HI_, "SavePhenoFile: Successful to save phenotype data");

    return(0);
}

void destroy(CFmSimulate* p)
{
	CFmNewTemp  fmRef;
	p->~CFmSimulate();
	operator delete(p, fmRef);
}
