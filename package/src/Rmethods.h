// Rmethods.h: some functions provided by R.
//
//////////////////////////////////////////////////////////////////////

#ifndef _RMETHODS_H_
#define _RMETHODS_H_

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MI(rows, cols, nrow, ncol) (nrow+ncol*rows)

void show_sexp( SEXP obj, const char* szName );
void show_matrix2( double* mat, int nRows, int nCols, const char* szName);
void show_vector2( double* vct, int nLen, const char* szName );

void R_SaveToFileV(SEXP obj, FILE *fp, int ascii);

SEXP getListElement(SEXP list, const char *str);

SEXP rmvnorm_chol(int n, SEXP mean, SEXP sigma);

int matprod(double *x, int nrx, int ncx, char x_transx,
                   double *y, int nry, int ncy, char y_transx, double *z);

int rmultnorm(int n, double* mu, double* vmat, int m, double* ans);

int solve(double* a, int an, double* b, int bx, int by );

double GetDeterminant( double* mat, int n );

int R_set_seed(int seed_num);

#ifdef __cplusplus
}
#endif

#endif
