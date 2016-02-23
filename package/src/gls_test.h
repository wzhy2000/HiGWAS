// test.h: interface for the Test Data
//
//////////////////////////////////////////////////////////////////////

#ifndef GLS_TEST_H
#define GLS_TEST_H

#include "fm_matrix.h"
#include "fm_vector.h"

#ifdef __cplusplus
extern "C" {
#endif

void proc_gls_simu_geno_test( CFmMatrix* pSnps);
void proc_gls_simu_pheno_test( CFmMatrix* pPhenoY,
        CFmMatrix* pPhenoZ, CFmMatrix* pCovarX,CFmVector* pZRange);

#ifdef __cplusplus
}
#endif

#endif // GLS_TEST_H
