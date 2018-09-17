// gpu_sharep.h: interface for the GPU memory
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_GPU_H_)
#define _FM_GPU_H_

#include <curand_kernel.h>

#include "fm_matrix.h"
#include "fm_vector.h"

#ifdef __cplusplus
extern "C" {
#endif

#define LG 4

#define THREADS_PER_BLOCK 1024

#define _M(a, i, j) a[ (i) + (j) * (int)(a[1]) + 3]

#define _V(a, i ) a[ (i) + 3]

#define PERR(call) \
  if (call) {\
   fprintf(stderr, "%s:%d Error [%d=%s] on "#call"\n", __FILE__, __LINE__,\
      call, cudaGetErrorString(cudaGetLastError()));\
   exit(1);\
  }

#define ERRCHECK \
  if (cudaPeekAtLastError()) { \
    fprintf(stderr, "%s:%d Error [%s]\n", __FILE__, __LINE__,\
       cudaGetErrorString(cudaGetLastError()));\
    exit(1);\
  }

clock_t startTimer();
double stopTimer(clock_t begin);

int _copy_fmMatrix_list( double** cudaList, int size, CFmMatrix** matList);
int _copy_fmVector_list( double** cudaList, int size, CFmVector** matList);
int _copy_fmVector_Device( double* gCudaDstVec, CFmVector& fmSrcVec );
int _copy_fmMatrix_Device( double* gCudaDstMat, CFmMatrix& fmSrcMat, bool bShow=false );
int _copyback_fmVector_Host( CFmVector& fmDstVec, double* gCudaSrcVec );
int _copyback_fmMatrix_Host( CFmMatrix& fmDstMat, double* gCudaSrcMat );

double** make_vector_list( double** pList, unsigned int N, unsigned int nLen, bool verbose=FALSE  );
double** make_matrix_list( double** pList, unsigned int N, unsigned int nRow, unsigned int nCol, bool verbose=FALSE );
double* make_matrix_ongpu( unsigned int nRow, unsigned int nCol );

#ifdef __cplusplus
}
#endif

#endif