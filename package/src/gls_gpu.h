// gls_gpu.h: interface for the GPU memory
//
//////////////////////////////////////////////////////////////////////

#if !defined(_GLS_GPU_H_)
#define _GLS_GPU_H_

#include <curand_kernel.h>

#include "fm_matrix.h"
#include "fm_vector.h"

#ifdef __cplusplus
extern "C" {
#endif


__device__ struct GPUobj{
	int nSize;

    //# array of cudaMatrix; ( size:N, Q*Q)
	double** all_corMat;
    //# array of cudaMatrix; ( size:N, Q*Q)
	double** all_corMat_MH;
    //# array of cudaMatrix; ( size:N, Q*Q)
    double** all_corTimes;
    //# array of cudaMatrix; ( size:N, Q*Q)
    double** all_corMat_Inv;
    //# array of cudaMatrix; ( size:N, Q*Q)
    double** all_corMat_MH_Inv;
    //# array of cudaVector; ( size:N, Q)
    double** all_yi;
    //# array of cudaVector; ( size:N, Q)
    double** all_rd;
    //# array of cudaMatrix; ( size:N, Q*LG)
    double** all_ui;

    //# cudaMatrix; (size:N*Q)
    double* X;
	//# cudaMatrix; (size:N*Q)
	double* Z;
	//# cudaMatrix; (size:N*Q)
	double* Z0;
	//# cudaVector; (size:N)
    double* mInZ;

    //# cudaVector; (size: LG)
    double* mu;
    //# cudaMatrix; (size, NC+1*LG )
    double* alpha;
    //# cudaMatrix; (size: P*LG)
    double* a;
    //# cudaMatrix; (size: P*LG)
    double* a_old;
    //# cudaMatrix; (size: P*LG)
    double* d;
    //# cudaMatrix; (size: P*LG)
    double* d_old;
	//# cudaVector  (size: P)
	double* vctP;
    //# cudaMatrix; (size: P*N)
    double* gen_a;
    //# cudaMatrix; (size: P*N)
    double* gen_d;
	//# cudaVector  (size: P)
	double* tau2;
	//# cudaVector  (size: P)
	double* tau2_x;

	//# cudaVector for temporary use(size:N)!
	double* tVN0;
	//# cudaVector for temporary use(size:N)!
	double* tVN1;
	//# cudaVector for temporary use(size:N)!
	double* tVN2;

	//# array of cudaMatrix for temporary use! (size:N, Q*Q)
	double** tempMat;
	// size: N, Q*Q
	double** tmp2;
	// size: N, Q*Q
	double** tmp3;

	void* pNext;
};


__device__ struct GPUShare{
    //# array of cudaMatrix for temporary use! (size: Q*Q)
    double* tempMat0;
    //# array of cudaMatrix for temporary use! (size: Q*Q)
    double* tempMat1;
    //# array of cudaMatrix for temporary use! (size: Q*Q)
    double* tempMat2;
    //# array of cudaMatrix for temporary use! (size: Q*Q)
    double* tempMat3;
    //# array of cudaMatrix for temporary use! (size: Q*Q)
    double* tempMat4;
    //# array of cudaMatrix for temporary use! (size: Q*Q)
    double* tMA;
    //# array of cudaMatrix for temporary use! (size: Q*Q)
    double* tMB;
    //# array of cudaMatrix for temporary use! (size: Q*Q)
    double* tMC;
    //# array of cudaMatrix for temporary use! (size: Q*Q)
    double* tMD;
    //# array of cudaMatrix for temporary use! (size: Q*Q)
    double* tME;

    double* pNext;
};

int _cuda_gpart0( int N, double rho, double tmp5 );
int _cuda_gpart1( struct GPUobj* gCuda, struct GPUobj* gCpuObj, int N, double rho, double tmp5 );
int _cuda_gpart2( struct GPUobj* gCuda, struct GPUobj* gCpuObj, struct GPUobj* gGpuMap, int N, int Q, int nC, double sigma2, CFmMatrix& alpha, CFmMatrix& tmp2, CFmMatrix& tmp3 );
int _cuda_gpart3( struct GPUobj* gCuda, struct GPUobj* gCpuObj, int N, int Q, int P, CFmMatrix& a, CFmMatrix& d );
int _cuda_gpart4( struct GPUobj* gCuda, struct GPUobj* gCpuObj, struct GPUobj* gGpuMap, int N, int Q, int j, int nC, double sigma2, CFmVector& mu, CFmMatrix& alpha, CFmMatrix& a, CFmMatrix& tmp2, CFmMatrix& tmp3);
int _cuda_gpart5( struct GPUobj* gCuda, struct GPUobj* gCpuObj, int N, int Q, int j, CFmMatrix& a, CFmMatrix& a_old );
int _cuda_gpart6( struct GPUobj* gCuda, struct GPUobj* gCpuObj, int P, double sigma2, double lambda2, double lambda2_x, CFmVector& vctP, CFmMatrix& a, CFmVector& tau2, CFmVector& tau2_x );
int _cuda_gpart7( struct GPUobj* gCuda, struct GPUobj* gCpuObj, struct GPUobj* gGpuMap, int N, int Q, int j, int nC, double sigma2, CFmMatrix& alpha, CFmVector& mu, CFmMatrix& d, CFmMatrix& tmp2, CFmMatrix& tmp3);
int _cuda_gpart8( struct GPUobj* gCuda, struct GPUobj* gCpuObj, int N, int Q, int j, CFmMatrix& d, CFmMatrix& d_old);
int _cuda_gpart9( struct GPUobj* gCuda, struct GPUobj* gCpuObj, int P, double lambda_st2, double lambda_st2_x, double sigma2, CFmVector& vctP, CFmMatrix& d, CFmVector& tau_st2, CFmVector& tau_st2_x);
int _cuda_gpart10( struct GPUobj* gCuda, struct GPUobj* gCpuObj, struct GPUobj* gGpuMap, int N, int Q, int nC, int nX, double sigma2, CFmMatrix& alpha, CFmVector& mu, CFmMatrix& tmp2, CFmMatrix& tmp3 );
int _cuda_gpart11( struct GPUobj* gCuda, struct GPUobj* gCpuObj, struct GPUobj* gGpuMap, int N, int Q, int nC, CFmMatrix& alpha, CFmVector& mu, double* sigma2_scale);
int _cuda_gpart12( struct GPUobj* gCuda, struct GPUobj* gCpuObj, struct GPUobj* gGpuMap, int N, int Q, int nC, double sigma2, CFmMatrix& alpha, CFmVector& mu, double* exp_diff);

int Init_GPUobj(struct GPUobj** pCpuObj, struct GPUobj** pGpuObj, struct GPUobj** pGpuMap, int N, int P, int Q, int NC);
int Free_GPUobj(struct GPUobj* pGpuObj, struct GPUobj* pCpuObj, struct GPUobj* pGpuMap, int N);

#ifdef __cplusplus
}
#endif

#endif