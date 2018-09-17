#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

#include <stdio.h>

#include "bls_gpu.h"

#include "gpu_share.h" 

__device__ double _cuda_invGau_same(int seed, double theta, double chi)
{
    curandState s ;

    // seed a random number generator
    curand_init ( (unsigned int) seed, 0, 0, &s) ;
    double _rn = curand_normal_double (&s);   
    double _ru = curand_uniform_double (&s);

#ifdef USECUDA
#ifdef MONTIME
    _rn = 0.67;
    _ru = 0.45;
#endif
#endif

    //squared normal, i.e., chi-square with df=1
    double chisq1 = _rn * _rn;
    double y1    = theta + 0.5*theta/chi * ( theta*chisq1 - sqrt(4*theta*chi*chisq1 + theta*theta*chisq1*chisq1) );
    double y2    = theta*theta/y1;
    double out_1 = _ru < (theta/(theta+y1));
    double value = out_1*y1+(1-out_1)*y2;

    return(value);
}


void _initBlsGPU(struct blsGPUobj* gCpu, struct blsGPUobj* gGpuObj, unsigned int nsize, unsigned int N, unsigned int P )
{
    memset(gCpu, 0, nsize);
    gCpu->nSize = nsize;
    gCpu->pNext = &(gCpu->pNext);

    gCpu->gen_pA  = make_matrix_ongpu( N, P );
    gCpu->gen_pD  = make_matrix_ongpu( N, P );
    gCpu->ad      = make_matrix_ongpu( P, 1 );
    gCpu->Var     = make_matrix_ongpu( P, 1 );
    gCpu->tmp     = make_matrix_ongpu( P, 1 );
    gCpu->tau2    = make_matrix_ongpu( P, 1 );
    gCpu->tmp3    = make_matrix_ongpu( N, 1 );
    gCpu->spdY    = make_matrix_ongpu( N, 1 );
}

int Init_blsGPU(struct blsGPUobj** pCpuObj, struct blsGPUobj** pGpuObj, int N, int P )
{
    int ngCudaSize = sizeof(struct blsGPUobj) ;

    struct blsGPUobj* gGpuObj;
    PERR( cudaMalloc( (void **)&gGpuObj, ngCudaSize ) );
    *pGpuObj = gGpuObj;

    struct blsGPUobj* gCpuCopy = (struct blsGPUobj*)Calloc( ngCudaSize, char);
    *pCpuObj = gCpuCopy;

    _initBlsGPU( gCpuCopy, gGpuObj, ngCudaSize, N, P );

    PERR( cudaMemcpy( gGpuObj, gCpuCopy, ngCudaSize, cudaMemcpyHostToDevice ) );
    
printf("CPU copy=%p\n", *pCpuObj);
printf("GPU addr=%p\n", *pGpuObj);
    
    return(0);
}

void _freeBlsGPU(struct blsGPUobj* pCpuObj, int N )
{
    PERR( cudaFree( pCpuObj->gen_pA ) );
    PERR( cudaFree( pCpuObj->gen_pD ) );
    PERR( cudaFree( pCpuObj->ad ) );
    PERR( cudaFree( pCpuObj->Var ) );
    PERR( cudaFree( pCpuObj->tau2 ) );
    PERR( cudaFree( pCpuObj->tmp ) );
    PERR( cudaFree( pCpuObj->spdY ) );
    PERR( cudaFree( pCpuObj->tmp3 ) );
}
    
int Free_blsGPU(struct blsGPUobj* pGpuObj, struct blsGPUobj* pCpuObj, int N)
{
    _freeBlsGPU( pCpuObj, N );

    PERR( cudaFree(pGpuObj) );

    Free(pCpuObj); 

    return(0);
}

__global__ void b_part1(struct blsGPUobj* gCuda, int N, int P )
{    
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if(j < P)
    {
        double a_j = _V(gCuda->ad, j );
        double aMu_j = 0.0;
        int N_AA = 0;
        int N_aa = 0;

        for (int k=0;k<N;k++)
        {
             double pA_kj = _M(gCuda->gen_pA, k, j);
             double spdY_k = _V(gCuda->spdY, k );
             
             double x0 = _V( gCuda->tmp3, k) - spdY_k + pA_kj * a_j;
             aMu_j = aMu_j + x0 * pA_kj;

             if (pA_kj>0) N_AA++;
             if (pA_kj<0) N_aa++;
        }     
        
        aMu_j = aMu_j / _V( gCuda->tmp, j);
        double aVar_j = _V(gCuda->Var, j);
        
        //repalce rnorm with curand_normal_double
        //double new_a = rnorm( aMu_j, sqrt(aVar_j) );
        
        curandState s ;
        // seed a random number generator
        curand_init ( (unsigned int) j*100, 0, 0, &s) ;
        double _rn = curand_normal_double (&s);   
        double new_a = _rn * sqrt(aVar_j) + aMu_j;
        
        if (N_AA==0 || N_aa==0) 
            new_a = 0.0;

        _V(gCuda->ad, j) = new_a;

   }
}


int _cuda_bpart1( struct blsGPUobj* gCuda, struct blsGPUobj* gCpuObj, int N, int P, CFmVector& spdY, CFmVector& aVar, CFmVector& tmp, CFmVector& tmp3, CFmVector& a)
{
    _copy_fmVector_Device( gCpuObj->spdY, spdY);
    _copy_fmVector_Device( gCpuObj->tmp, tmp );
    _copy_fmVector_Device( gCpuObj->tmp3, tmp3 );
    _copy_fmVector_Device( gCpuObj->Var, aVar);
    _copy_fmVector_Device( gCpuObj->ad, a );
    
    b_part1<<< (N + (THREADS_PER_BLOCK-1)) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>(gCuda, N, P);
    cudaDeviceSynchronize();
    ERRCHECK;

    _copyback_fmVector_Host( a, gCpuObj->ad);
    
//printf("End of part1\n");
    return (0);
}

__global__ void b_part2(struct blsGPUobj* gCuda, int P, double lambda2, double sigma2 )
{    
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if(j < P)
    {
        double _tau2 = 0.0;
        double aj = _V(gCuda->ad, j);
        if ( aj != 0.0 )
        {
            double InvTau2_1 = sqrt( lambda2 * sigma2)/fabs(aj);
            _tau2 = 1/ _cuda_invGau_same(j*15,  InvTau2_1, lambda2);
        }   
        _V(gCuda->tau2, j) = _tau2;
   }
}

int _cuda_bpart2( struct blsGPUobj* gCuda, struct blsGPUobj* gCpuObj, int P, double lambda2, double sigma2, CFmVector& a, CFmVector& tau2 )
{
    _copy_fmVector_Device( gCpuObj->ad, a );
    
    b_part2<<< (P + (THREADS_PER_BLOCK-1)) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>( gCuda, P, lambda2, sigma2  );
    cudaDeviceSynchronize();
    ERRCHECK;

    _copyback_fmVector_Host( tau2, gCpuObj->tau2);
    
//printf("End of part2\n");
    return (0);

}


__global__ void b_part3(struct blsGPUobj* gCuda, int N, int P )
{    
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if(j < P)
    {
        double d_j = _V(gCuda->ad, j );
        double dMu_j = 0.0;
        int N_Aa = 0;

        for (int k=0;k<N;k++)
        {
             double pD_kj = _M(gCuda->gen_pD, k, j);
             double spdY_k = _V(gCuda->spdY, k );
             
             double x0 = _V( gCuda->tmp3, k) - spdY_k + pD_kj * d_j;
             dMu_j = dMu_j + x0 * pD_kj;

             if (pD_kj == 0) N_Aa++;
        }     
        
        dMu_j = dMu_j/ _V( gCuda->tmp, j);
             
        double dVar_j =  _V(gCuda->Var, j);
        
        //repalce rnorm with curand_normal_double
        //double new_d = rnorm( dMu_j, sqrt(dVar_j) );

        curandState s ;
        // seed a random number generator
        curand_init ( (unsigned int) j*10, 0, 0, &s) ;
        double _rn = curand_normal_double (&s);   
        double new_d = _rn * sqrt(dVar_j) + dMu_j;
        
        if (N_Aa == 0) 
            new_d = 0.0;

        _V(gCuda->ad, j) = new_d;
   }
}

int _cuda_bpart3( struct blsGPUobj* gCuda, struct blsGPUobj* gCpuObj, int N, int P, CFmVector& spdY, CFmVector& dVar, CFmVector& tmp, CFmVector& tmp3, CFmVector& d)
{
    _copy_fmVector_Device( gCpuObj->spdY, spdY);
    _copy_fmVector_Device( gCpuObj->tmp, tmp );
    _copy_fmVector_Device( gCpuObj->tmp3, tmp3 );
    _copy_fmVector_Device( gCpuObj->Var, dVar);
    _copy_fmVector_Device( gCpuObj->ad, d );
    
    b_part3<<< (N + (THREADS_PER_BLOCK-1)) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>(gCuda, N, P);
    cudaDeviceSynchronize();
    ERRCHECK;

    _copyback_fmVector_Host( d, gCpuObj->ad);
    
//printf("End of part3\n");
    return (0);

}

int _cuda_bpart4( struct blsGPUobj* gCuda, struct blsGPUobj* gCpuObj, int P, double lambda_st2, double sigma2, CFmVector& d, CFmVector& tau2_st )
{
    _copy_fmVector_Device( gCpuObj->ad, d );
    
    b_part2<<< (P + (THREADS_PER_BLOCK-1)) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>( gCuda, P, lambda_st2, sigma2  );
    cudaDeviceSynchronize();
    ERRCHECK;

    _copyback_fmVector_Host( tau2_st, gCpuObj->tau2);
    
//printf("End of part4\n");
    return (0);
}
