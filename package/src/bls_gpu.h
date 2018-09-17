// bls_gpu.h: interface for the GPU memory
//
//////////////////////////////////////////////////////////////////////

#if !defined(_BLS_GPU_H_)
#define _BLS_GPU_H_

#include <curand_kernel.h>

#include "fm_matrix.h"
#include "fm_vector.h"

#ifdef __cplusplus
extern "C" {
#endif

__device__ struct blsGPUobj{
	int nSize;

    //# cudaMatrix; (dim: N, P)
    double* gen_pA;
    //# cudaMatrix; (dim: N, p)
    double* gen_pD;
    //# cudaVector; (dim: P, 1)
    double* ad;
    //# cudaVector; (dim: P, 1)
    double* Var;
	//# cudaVector  (dim: P, 1)
	double* tau2;
	//# cudaVector  (dim: P, 1)
	double* tmp;
	//# cudaVector  (dim: N, 1)
	double* spdY;
	//# cudaVector  (dim: N, 1)
	double* tmp3;

	void* pNext;
};

int _cuda_bpart1( struct blsGPUobj* gCuda, struct blsGPUobj* gCpuObj, int N, int P, CFmVector& spdY, CFmVector& aVar, CFmVector& tmp, CFmVector& tmp3, CFmVector& a);
int _cuda_bpart2( struct blsGPUobj* gCuda, struct blsGPUobj* gCpuObj, int P, double lambda2, double sigma2, CFmVector& a, CFmVector& tau2 );
int _cuda_bpart3( struct blsGPUobj* gCuda, struct blsGPUobj* gCpuObj, int N, int P, CFmVector& spdY, CFmVector& dVar, CFmVector& tmp, CFmVector& tmp3, CFmVector& d);
int _cuda_bpart4( struct blsGPUobj* gCuda, struct blsGPUobj* gCpuObj, int P, double lambda_st2, double sigma2, CFmVector& d, CFmVector& tau2_st );

int Init_blsGPU(struct blsGPUobj** pCpuObj, struct blsGPUobj** pGpuObj, int N, int P);
int Free_blsGPU(struct blsGPUobj* pGpuObj, struct blsGPUobj* pCpuObj, int N);

#ifdef __cplusplus
}
#endif

#endif