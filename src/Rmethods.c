/* Rmethods.c  -	the interface for some R functions
 *
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#include "Rmethods.h"

SEXP getListElement(SEXP list, const char *str)
{
	SEXP elmt = R_NilValue;
	SEXP names = getAttrib( list, R_NamesSymbol);

	for (int i=0;i<length(list); i++)
		if (strcmp(CHAR(STRING_ELT(names, i)), str)==0)
		{
			elmt = VECTOR_ELT(list, i);
			break;
		}

	return elmt;
}

/*extern void F77_NAME(dchdc)(double*, int*, int*, double*, int*, int*, int*);
extern void F77_NAME(chol)(double *a, int *lda, int *n, double *v, int *info);*/

/*
//function rmvnorm(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
//    method = c("eigen", "svd", "chol"))
//{
//    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
//        stop("sigma must be a symmetric matrix")
//    }
//    if (length(mean) != nrow(sigma)) {
//        stop("mean and sigma have non-conforming size")
//    }
//    sigma1 <- sigma
//    dimnames(sigma1) <- NULL
//    if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
//        warning("sigma is numerically not symmetric")
//    }
//    method <- match.arg(method)
//    if (method == "eigen") {
//        ev <- eigen(sigma, symmetric = TRUE)
//        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
//            warning("sigma is numerically not positive definite")
//        }
//        retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*%
//            t(ev$vectors)
//    }
//    else if (method == "svd") {
//        sigsvd <- svd(sigma)
//        if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
//            warning("sigma is numerically not positive definite")
//        }
//        retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
//    }
//    else if (method == "chol") {
//        retval <- chol(sigma, pivot = TRUE)
//        o <- order(attr(retval, "pivot"))
//        retval <- retval[, o]
//    }
//    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
//    retval <- sweep(retval, 2, mean, "+")
//    colnames(retval) <- names(mean)
//    retval
//}
*/

/*chol <- function(x, ...) UseMethod("chol")

chol.default <- function(x, pivot = FALSE, LINPACK = pivot, ...)
{
    if (is.complex(x))
        stop("complex matrices not permitted at present")
    else if(!is.numeric(x))
	stop("non-numeric argument to 'chol'")

    if(is.matrix(x)) {
	if(nrow(x) != ncol(x))
	    stop("non-square matrix in 'chol'")
	n <- nrow(x)
    } else {
	if(length(x) != 1L)
	    stop("non-matrix argument to 'chol'")
	n <- 1L
    }
    if(!pivot && !LINPACK)
        return(.Call("La_chol", as.matrix(x), PACKAGE = "base"))

    if(!is.double(x)) storage.mode(x) <- "double"

    if(pivot) { ## code could be used in the other case too
        xx <- x
        xx[lower.tri(xx)] <- 0
        z <- .Fortran("dchdc",
                      x = xx,
                      n,
                      n,
                      double(n),
                      piv = integer(n),
                      as.integer(pivot),
                      rank = integer(1L),
                      DUP = FALSE, PACKAGE = "base")
        if (z$rank < n)
            if(!pivot) stop("matrix not positive definite")
            else warning("matrix not positive definite")
        robj <- z$x
        if (pivot) {
            attr(robj, "pivot") <- z$piv
            attr(robj, "rank") <- z$rank
            if(!is.null(cn <- colnames(x)))
                colnames(robj) <- cn[z$piv]
        }
        robj
    } else {
        z <- .Fortran("chol",
                      x = x,
                      n,
                      n,
                      v = matrix(0, nrow=n, ncol=n),
                      info = integer(1L),
                      DUP = FALSE, PACKAGE = "base")
        if(z$info)
            stop("non-positive definite matrix in 'chol'")
        z$v
    }
}*/


SEXP f_chol(SEXP x)
{
	SEXP dim = getAttrib( x, R_DimSymbol );
	int ncol = INTEGER(dim)[0];
	double* x_r = REAL( x );

		//R: xx <- x;
		SEXP xx;
		PROTECT(xx = allocMatrix(REALSXP, ncol, ncol) );
		double* xx_r = REAL( xx );

        // R: xx[lower.tri(xx)] <- 0
		for (int i=0; i<ncol; i++) // row
			for (int j=0; j<ncol; j++) // col
                if (j<=i)
				{
                    xx_r [i* ncol + j] = x_r [i* ncol + j];
				}
				else
                    xx_r [i* ncol + j] = 0;

        //R: z <- .Fortran("dchdc", x = xx, n, n, double(n), piv = integer(n), as.integer(pivot),
        //         rank = integer(1L), DUP = FALSE, PACKAGE = "base")
		SEXP double_n;
		PROTECT(double_n = allocVector(REALSXP, ncol) );

		SEXP piv;
		PROTECT(piv = allocVector(INTSXP, ncol) );
		int* piv_n = INTEGER(piv);

		SEXP rank;
		PROTECT(rank = allocVector(INTSXP, 1) );
		int* rank_n = INTEGER(rank);
        rank_n[0] = 0;

		int pivot = TRUE;
		F77_CALL(dchdc)( xx_r, &ncol, &ncol, REAL(double_n), piv_n, &pivot, rank_n);
        if ( *rank_n < ncol )
        {
			if(!pivot)
				error("matrix not positive definite");
            else
            	warning("matrix not positive definite");
		}

		// .Fortran will return value in xx_r, piv and rank field.
        if (pivot)
        {
            setAttrib( xx, install("pivot"), piv );
            setAttrib( xx, install("rank"),  rank );
        }

		SEXP dim2;
		PROTECT(dim2=allocVector(INTSXP,2));
		INTEGER(dim2)[0] = ncol;
		INTEGER(dim2)[1] = ncol;
		setAttrib( xx, R_DimSymbol, dim2);

		UNPROTECT(5);

        return(xx);
}

int matprod(double *x, int nrx, int ncx, char x_transx,
                   double *y, int nry, int ncy, char y_transx, double *z)
{
    //char *transa = "N" or "T", *transb = "N" or "T";
	if ( x_transx!='N' && x_transx!='T')
	{
		Rprintf("x_transx:%c\n",x_transx );
		return -1;
	}
	if ( y_transx!='N' && y_transx!='T')
	{
		Rprintf("y_transx:%c\n",y_transx );
		return -1;
	}

    int i,  j, k;
    double one = 1.0, zero = 0.0;
    double sum;
    int have_na = 0;

    if (x_transx=='N' && y_transx=='N')
    {
        char *transa = "N", *transb = "N";
        if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0)
        {
            /* Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
             * The test is only O(n) here
             */
            for (i = 0; i < nrx*ncx; i++)
                if (ISNAN(x[i]))
                {
                    have_na = 1;
                    break;
                }

            if (!have_na)
                for (i = 0; i < nry*ncy; i++)
                    if (ISNAN(y[i]))
                        {have_na = 1; break;}

            if (have_na)
            {
                for (i = 0; i < nrx; i++)
                for (k = 0; k < ncy; k++)
                {
                    sum = 0.0;
                    for (j = 0; j < ncx; j++)
                        sum += x[i + j * nrx] * y[j + k * nry];

                    z[i + k * nrx] = sum;
                }
            }
            else
                F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
                        x, &nrx, y, &nry, &zero, z, &nrx);
        }
        else /* zero-extent operations should return zeroes */
            for(i = 0; i < nrx*ncy; i++) z[i] = 0;
    }

    if (x_transx=='T' && y_transx=='T')
    {
        char *transa = "T", *transb = "N";
        double one = 1.0, zero = 0.0;
        if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
            F77_CALL(dgemm)(transa, transb, &ncx, &ncy, &nrx, &one,
                x, &nrx, y, &nry, &zero, z, &ncx);
        }
        else
        { /* zero-extent operations should return zeroes */
            int i;
            for(i = 0; i < ncx*ncy; i++) z[i] = 0;
        }
    }

    if (x_transx=='N' && y_transx=='T')
    {
        char *transa = "N", *transb = "T";
        double one = 1.0, zero = 0.0;
        if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0)
        {
            F77_CALL(dgemm)(transa, transb, &nrx, &nry, &ncx, &one,
                x, &nrx, y, &nry, &zero, z, &nrx);
        }
        else
        { /* zero-extent operations should return zeroes */
            int i;
            for(i = 0; i < nrx*nry; i++) z[i] = 0;
        }
    }

	return(0);
}


SEXP rmvnorm_chol(int n, SEXP mean, SEXP sigma)
{
	SEXP retval = f_chol( sigma );

	SEXP dim = getAttrib( retval, R_DimSymbol );
	int ncol = INTEGER(dim)[0];
	double* retval_r = REAL( retval );
    SEXP pivot = getAttrib(retval, install("pivot") );
	SEXP rank = getAttrib(retval, install("rank") );

	//R: o <- order(attr(retval, "pivot"))
    int errorOccurred=0;
	SEXP r_o, e1;
    PROTECT(e1 = lang2(install("order"), pivot ) );
    PROTECT( r_o = R_tryEval(e1, R_GlobalEnv, &errorOccurred) );

	//R: retval_new <- retval[, o]
	//means : subset(x, subset=TRUE, select=o,drop=FALSE)
	SEXP subset;
	PROTECT(subset = allocVector(LGLSXP, 1) );
	int* subset_b = LOGICAL(subset);
	subset_b[0] = 1;


	SEXP retval_new, e2;
    PROTECT(e2 = lang4(install("subset"), retval, subset, r_o ) );
    PROTECT( retval_new = R_tryEval(e2, R_GlobalEnv, &errorOccurred) );
	double* retval_new_r = REAL( retval_new );

	//R: rn_m <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval_new
	SEXP rn_m, rn_m2;
	PROTECT( rn_m = allocMatrix(REALSXP, n, ncol ) );
	double* rn_m_r=REAL( rn_m );
    //GetRNGstate();
	for (int i=0; i<n; i++)
		for (int j=0; j<ncol; j++)
			rn_m_r[ncol*i +j ] = rnorm( 0, 1);
    //PutRNGstate();

	PROTECT( rn_m2 = allocMatrix(REALSXP, n, ncol ) );
	double* rn_m2_r=REAL( rn_m2 );
	matprod( rn_m_r, n, ncol, 'N', retval_new_r, ncol, ncol, 'N', rn_m2_r );

	//R:rn_m <- sweep(rn_m, 2, mean, "+")
	double* mean_r = REAL( mean );
	for (int i=0; i<n; i++)
		for (int j=0; j<ncol; j++)
			rn_m2_r[ncol*i +j ] = rn_m2_r[ncol*i +j ] + mean_r[j];

	//R:colnames(rn_m) <- names(mean)
	UNPROTECT(7);

	return(rn_m2);
}


/* R-code
rmultnorm<-function (n, mu, vmat, tol = 1e-10)
{
    p <- ncol(vmat)
    if (length(mu) != p)
        stop(paste("mu vector is the wrong length:", length(mu)))
    vs <- La.svd(vmat)
    vsqrt <- t(t(vs$vt) %*% (t(vs$u) * sqrt(vs$d)))
    ans <- matrix(rnorm(n * p), nrow = n) %*% vsqrt
    ans <- sweep(ans, 2, mu, "+")
    dimnames(ans) <- list(NULL, dimnames(vmat)[[2]])
    ans
}

n = nrow(x);.
p = ncol(x);

svd<-function(x, nu = min(n, p), nv = min(n, p))
{
	z <- .Fortran("dsvdc", as.double(x), n, n, p, d = double(mm),
        double(p), u = u, n, v = v, p, double(n), as.integer(job),
        info = integer(1L), DUP = FALSE, PACKAGE = "base")[c("d",
        "u", "v", "info")]
    if (z$info)
        stop(gettextf("error %d in 'dsvdc'", z$info), domain = NA)
    z$d <- z$d[1L:mn]
    if (nv && nv < p)
        z$v <- z$v[, 1L:nv, drop = FALSE]
    z[c("d", if (nu) "u", if (nv) "v")]
}

function (x, nu = min(n, p), nv = min(n, p))
{
    if (!is.numeric(x) && !is.complex(x))
        stop("argument to 'La.svd' must be numeric or complex")
    if (any(!is.finite(x)))
        stop("infinite or missing values in 'x'")
    x <- as.matrix(x)
    if (is.numeric(x))
        storage.mode(x) <- "double"
    n <- nrow(x)
    p <- ncol(x)
    if (!n || !p)
        stop("0 extent dimensions")

    if (nu || nv)
    {
        np <- min(n, p)
        if (nu <= np && nv <= np)
        {
            jobu <- "S"
            u <- matrix(0, n, np)
            v <- matrix(0, np, p)
        }
        else
        {
            jobu <- "A"
            u <- matrix(0, n, n)
            v <- matrix(0, p, p)
        }
    }
    else
    {
        jobu <- "N"
        u <- matrix(0, 1L, 1L)
        v <- matrix(0, 1L, 1L)
    }

    jobv <- ""
    res <- .Call("La_svd", jobu, jobv, x, double(min(n, p)),
            u, v, "dgsedd", PACKAGE = "base")
    res <- res[c("d", if (nu) "u", if (nv) "vt")]
    if (nu)
        res$u <- res$u[, 1L:min(n, nu), drop = FALSE]
    if (nv)
        res$vt <- res$vt[1L:min(p, nv), , drop = FALSE]
    return(res)
}


//TEST code
x = diag(1,4);
rho=0.6
for (ii in 1:4)
	for (jj in 1:4)
		 x[ii,jj] <- 4*rho^(abs(ii-jj));
r<-svd(x);

r$u%*%diag(r$d)%*%t(r$v);
x;

*/

double* rmt_vs_d = NULL;
double* rmt_vs_e = NULL;
double* rmt_vs_u = NULL;
double* rmt_vs_v = NULL;
double* rmt_work = NULL;
double* rmt_vmat = NULL;
double* rmt_z =NULL;
double* rmt_norm = NULL;
int     rmt_dim=0;
int     rmt_norm_n=0;

int rmultnorm(int n, double* mu, double* vmat0, int m, double* ans)
{
	int job=11;
	int info=0;

    if (rmt_vs_d!=NULL && (rmt_dim<m || rmt_norm_n<n) )
    {
        Free(rmt_vs_d); rmt_vs_d = NULL;
        Free(rmt_vs_e); rmt_vs_e = NULL;
        Free(rmt_vs_u); rmt_vs_u = NULL;
        Free(rmt_vs_v); rmt_vs_v = NULL;
        Free(rmt_work); rmt_work = NULL;
        Free(rmt_vmat); rmt_vmat = NULL;
        Free(rmt_z);    rmt_z = NULL;
   }

    if (rmt_vs_d==NULL)
    {
//**MEMTEST: Rprintf("ALLOC MEMORY in Rmethods 1\n");
        rmt_vs_d = (double*)Calloc( m,   double );
        rmt_vs_e = (double*)Calloc( m,   double );
        rmt_vs_u = (double*)Calloc( m*m, double );
        rmt_vs_v = (double*)Calloc( m*m, double );
        rmt_work = (double*)Calloc( m,   double );
        rmt_vmat = (double*)Calloc( m*m, double );
        rmt_z    = (double*)Calloc( m*m, double );
        rmt_norm = (double*)Calloc( n*m, double );
        rmt_dim  = m;
        rmt_norm_n = n;
    }

    for(int i=0; i<m; i++) rmt_vs_d[i] = 0;
    for(int i=0; i<m; i++) rmt_vs_e[i] = 0;
    for(int i=0; i<m*m; i++) rmt_vs_u[i] = 0;
    for(int i=0; i<m*m; i++) rmt_vs_v[i] = 0;
    for(int i=0; i<m; i++) rmt_work[i] = 0;
    for(int i=0; i<m*m; i++) rmt_z[i] = 0;
    for(int i=0; i<m*n; i++) rmt_norm[i] = 0;
    for(int i=0; i<m*m; i++) rmt_vmat[i] = vmat0[i];

    F77_CALL(dsvdc)(rmt_vmat, &m, &m,
            &m, /*d*/rmt_vs_d,
            /*e*/rmt_vs_e,
            /*u*/rmt_vs_u, &m,
            /*v*/rmt_vs_v, &m,
            /*work*/rmt_work,
			/*job*/&job,
			/*info*/&info);
	if(info!=0)
		return -1;

    //R: vsqrt <- t(vs$v %*% (t(vs$u) * sqrt(vs$d)))
	for (int i=0; i<m; i++)
        for(int j=0;j<m;j++)
            rmt_vs_u[ MI(m, m, i, j) ] *= sqrt(rmt_vs_d[j]);


    int ret = matprod(rmt_vs_v, m, m, 'N', rmt_vs_u, m, m, 'T', rmt_z);
	if (ret!=0)
		return(ret);

	for(int i=0; i<m; i++)
		for(int j=i+1; j<m; j++)
		{
            double swp = rmt_z[MI(m,m,i,j)];
            rmt_z[MI(m,m,i,j)] = rmt_z[MI(m,m,j,i)];
            rmt_z[MI(m,m,j,i)] = swp;
		}

    //R: ans <- matrix(rnorm(n * p), nrow = n) %*% vsqrt
    //GetRNGstate();
    for(int i=0; i<n*m; i++) rmt_norm[i] = rnorm(0,1);
    //PutRNGstate();

    ret = matprod(rmt_norm, n, m, 'N', rmt_z, m, m, 'N', ans);
    if (ret)
		return(ret);

    //R: ans <- sweep(ans, 2, mu, "+")
	for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            ans[MI(n,m,i,j)] = ans[MI(n,m,i,j)] + mu[j];

    return(0);
}

int* sol_ipiv = NULL;
double* sol_avals = NULL;
int sol_an = 0;

// A * X = B;
int solve(double* a, int an, double* b, int bx, int by)
{
    /*
	//MKL_INT ipiv[N];
	double a[LDA*N] = {
		6.80, -2.11,  5.66,  5.97,  8.23,
		-6.05, -3.30,  5.36, -4.44,  1.08,
		-0.45,  2.58, -2.70,  0.27,  9.04,
		8.32,  2.71,  4.35, -7.17,  2.14,
           -9.67, -5.14, -7.26,  6.08, -6.87
        };
        double b[LDB*NRHS] = {
            4.02,  6.19, -8.22, -7.57, -3.03,
           -1.56,  4.00, -8.67,  1.75,  2.86,
            9.81, -4.09, -4.57, -8.61,  8.99
        };
        info = LAPACKE_dgesv( LAPACK_COL_MAJOR, n, nrhs, a, lda, ipiv, b, ldb );*/

    if(an == 0)
    	error("'a' is 0-diml");

    if(by/*p*/ == 0)
    	error("no right-hand side in 'b'");

    if( bx != an)
		error("'b' (%d x %d) must be compatible with 'a' (%d x %d)", bx, by, an, an);

    if(sol_ipiv!=NULL && an>sol_an)
    {
        Free(sol_avals);sol_avals=NULL;
        Free(sol_ipiv);sol_ipiv=NULL;
        sol_an = 0;
    }

    if (sol_ipiv==NULL)
    {
//**MEMTEST: Rprintf("ALLOC MEMORY in Rmethods 2 \n");
        sol_ipiv = (int *) Calloc( an , int );
        sol_avals = (double *) Calloc (an * an , double );
        sol_an = an;
    }

    memset(sol_ipiv, 0, an * sizeof(int));
    memcpy( sol_avals, a, an * an * sizeof(double));

	int info;
    F77_CALL(dgesv)(&an, &by, sol_avals, &an, sol_ipiv, b, &an, &info);

    if (info < 0)
    {
        //error( "argument %d of Lapack routine %s had invalid value",  -info, "dgesv");
        Rprintf( "argument %d of Lapack routine %s had invalid value",  -info, "dgesv");
        return(-1);
    }

    if (info > 0)
    {
        //error( "Lapack routine dgesv: system is exactly singular" );
        Rprintf( "Lapack routine dgesv: system is exactly singular" );
        return(-1);
    }

    /*The following is checking the correction of answer???*/
    /*double* anorm = F77_CALL(dlange)("1", &an, &an, a, &an, (double*) NULL);
    double* work = (double *) Calloc (4*an,  double );
    double rcond;
    F77_CALL(dgecon)("1", &an, a, &an, &anorm, &rcond, work, sol_ipiv, &info);
    double tol=0.0001;
    Free(work);
    if (rcond < tol)
        error("system is computationally singular: reciprocal condition number = %g",rcond);*/

    return(0);
}

void show_object( SEXP obj )
{
	SEXP show, e1;
	int errorOccurred;
    PROTECT(e1 = lang2(install("print.default"), obj ) );
    PROTECT(show= R_tryEval(e1, R_GlobalEnv, &errorOccurred) );
    UNPROTECT(2);
}

void show_sexp( SEXP obj, const char* szName )
{
    Rprintf("\nName: %s\n", szName);
    show_object(obj);
}

void show_matrix2( double* mat, int nRows, int nCols, const char* szName)
{
    Rprintf("\nMatrix: %s[%d,%d]\n", szName, nRows, nCols);

	SEXP matrix;

	PROTECT( matrix = allocMatrix(REALSXP, nRows, nCols) );
	double* matrix_r = REAL(matrix);

	for (int i=0;i<nRows*nCols;i++) matrix_r[i] = mat[i];

    show_object(matrix);

	UNPROTECT(1);
}

void show_vector2( double* vct, int nLen, const char* szName )
{
    Rprintf("\nVector: %s[%d]\n", szName, nLen );

    SEXP vector;

	PROTECT( vector = allocVector(REALSXP, nLen  ) );
	double* vector_r = REAL(vector);
    for (int i=0;i<nLen;i++) vector_r[i] = vct[i];

    show_object( vector );

	UNPROTECT(1);
}

int* det_jpvt=NULL;
double* det_pData=NULL;
int det_length=0;

double GetDeterminant( double* mat, int n )
{
    int sign = 1;
    double modulus = 0.0; /* -Wall */

    if (det_jpvt!=NULL && det_length<n )
    {
        Free(det_jpvt); det_jpvt = NULL;
        Free(det_pData); det_pData = NULL;
        det_length = 0;
   }

    if(det_jpvt==NULL)
    {
//MEMTEST: Rprintf("ALLOC MEMORY in Rmethods 3 \n");
        det_jpvt = (int *) Calloc( n, int );
        det_pData = (double*) Calloc( n*n, double );
        det_length = n;
    }

    memset(det_jpvt, 0, n*sizeof(int) );
    memcpy( det_pData, mat, n*n*sizeof(double) );

    int info;
    F77_CALL(dgetrf)(&n, &n, det_pData, &n, det_jpvt, &info);

    if (info < 0)
		error( "error code %d from Lapack routine '%s'", info, "dgetrf");
    else if (info > 0)
    {
		/* Singular matrix:  U[i,i] (i := info) is 0 */
		/*warning("Lapack dgetrf(): singular matrix: U[%d,%d]=0", info,info);*/
		modulus = 0;
		return(0);
    }
    else
    {
		for (int i = 0; i < n; i++)
            if (det_jpvt[i] != (i + 1))
	    		sign = -sign;

		modulus = 1.0;
		for (int i = 0; i < n; i++)
            modulus *= det_pData[i*(n + 1)];

		if (modulus < 0)
		{
			modulus = -modulus;
			sign = -sign;
		}

		return( sign * modulus );
    }
}

static int R_DefaultSaveFormatVersion = 2;

#define R_MAGIC_ASCII_V2   2001
#define R_MAGIC_BINARY_V2  2002
#define R_MAGIC_XDR_V2     2003

static int R_WriteMagic(FILE *fp, int number)
{
    unsigned char buf[5];
    size_t res;

    number = abs(number);
    switch (number)
    {
    case R_MAGIC_ASCII_V2:
        /* Version >=2 - R Data, ASCII Format */
        strcpy((char*)buf, "RDA2");
        break;
    case R_MAGIC_BINARY_V2:
        /* Version >=2 - R Data, Binary Format */
        strcpy((char*)buf, "RDB2");
        break;
    case R_MAGIC_XDR_V2:
        /* Version >=2 - R Data, XDR Binary Format */
        strcpy((char*)buf, "RDX2");
        break;
    default:
        buf[0] = (number/1000) % 10 + '0';
        buf[1] = (number/100) % 10 + '0';
        buf[2] = (number/10) % 10 + '0';
        buf[3] = number % 10 + '0';
    }

    buf[4] = '\n';
    res = fwrite((char*)buf, sizeof(char), 5, fp);
    if(res != 5)
        return(-1);
    else
        return(0);
}

void R_SaveToFileV(SEXP obj, FILE *fp, int ascii)
{
    int version = R_DefaultSaveFormatVersion;
    struct R_outpstream_st out;
    R_pstream_format_t type;

    int magic;
    if (ascii)
    {
        magic = R_MAGIC_ASCII_V2;
        type = R_pstream_ascii_format;
    }
    else
    {
        magic = R_MAGIC_XDR_V2;
        type = R_pstream_xdr_format;
    }
    R_WriteMagic(fp, magic);

    R_InitFileOutPStream(&out, fp, type, version, NULL, NULL);

    R_Serialize(obj, &out);
}

// R: set.seed( dat$taskId^2 + dat$taskId^3 - 1
int R_set_seed(int seed_num)
{
    SEXP seeds, e1;

    PROTECT( seeds = allocVector( INTSXP, 1));
    int* seeds_r = INTEGER(seeds);
    seeds_r[0] = seed_num;
    //seeds_r[0] = 12;

    int errorOccurred;
    PROTECT(e1 = lang4(install("set.seed"), seeds, R_NilValue, R_NilValue ) );
    PROTECT(R_tryEval(e1, R_GlobalEnv, &errorOccurred) );

    UNPROTECT(3);

    return(errorOccurred);
}

