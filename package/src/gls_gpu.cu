#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

#include <stdio.h>

#include "gls_gpu.h"
#include "gpu_share.h" 
#include "gpu_matrix.cux"
#include "gpu_reduction.cux"

int _CheckCuda()
{
	int nCount = 0;
	if ( cudaGetDeviceCount ( &nCount ) != cudaSuccess )
    {
		Rprintf("Failed to call CUDA library, the GLS model may not run on GPU nodes.\n");
		return(0);
    }

    Rprintf("%d GPU card(s) is/are available for the GLS model.\n", nCount);
	return(nCount);
}

__device__ double _cuda_invGau(int seed, double theta, double chi)
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

 
__global__ void g_part0(int N, double rho, double tmp5 )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < N)
    {
        printf("hello %d\n", i);   
    } 
}

int _cuda_gpart0( int N, double rho, double tmp5)
{   
    g_part0<<< 1, 5>>>(N, rho, tmp5 );
    cudaDeviceSynchronize();
    ERRCHECK;
    
    return(0);
}
 
__device__ GPUShare* ConvetSharePoint( double* smb, int Q )
{
    GPUShare* p = (GPUShare*)smb; 
    p->tempMat0 = (double*)(&(p->pNext));
    p->tempMat1 = p->tempMat0 + (Q*Q+3);
    p->tempMat2 = p->tempMat1 + (Q*Q+3);
    p->tempMat3 = p->tempMat2 + (Q*Q+3);
    p->tempMat4 = p->tempMat3 + (Q*Q+3);
    p->tMA  = p->tempMat4 + (Q*Q+3);
    p->tMB  = p->tMA + (Q*Q+3);
    p->tMC  = p->tMB + (Q*Q+3);
    p->tMD  = p->tMC + (Q*Q+3);
    p->tME  = p->tMD + (Q*Q+3);

    p->tempMat0[0]= Q*Q;
    p->tempMat1[0]= Q*Q;
    p->tempMat2[0]= Q*Q;
    p->tempMat3[0]= Q*Q;
    p->tempMat4[0]= Q*Q;
    p->tMA[0]= Q*Q;
    p->tMB[0]= Q*Q;
    p->tMC[0]= Q*Q;
    p->tMD[0]= Q*Q;
    p->tME[0]= Q*Q;

    return(p);
}
        
void _initGPUobj(struct GPUobj* gCpu, struct GPUobj* gGpuObj, struct GPUobj* gGpuMap, unsigned int nsize, unsigned int N, unsigned int P, unsigned int Q, unsigned int NC )
{
printf("%p, %d, %d, %d, %d \n", nsize, N, P, Q, NC );

    memset(gCpu, 0, nsize);
    gCpu->nSize = nsize;
    gCpu->pNext = &(gCpu->pNext);

    gCpu->X     = make_matrix_ongpu( N, NC+1 );
    gCpu->Z     = make_matrix_ongpu( N, Q );
    gCpu->Z0    = make_matrix_ongpu( N, Q );
    gCpu->mInZ  = make_matrix_ongpu( N, 1 );
    gCpu->mu    = make_matrix_ongpu( LG, 1 );
    gCpu->alpha = make_matrix_ongpu( NC+1, LG );
    gCpu->a     = make_matrix_ongpu( P, LG );
    gCpu->a_old = make_matrix_ongpu( P, LG );
    gCpu->d     = make_matrix_ongpu( P, LG );
    gCpu->d_old = make_matrix_ongpu( P, LG );
    gCpu->vctP  = make_matrix_ongpu( P, 1 );
    gCpu->tau2  = make_matrix_ongpu( P, 1 );
    gCpu->tau2_x= make_matrix_ongpu( P, 1 );
    gCpu->gen_a = make_matrix_ongpu( P, N );
    gCpu->gen_d = make_matrix_ongpu( P, N );
    gCpu->tVN0  = make_matrix_ongpu( N, 1 );
    gCpu->tVN1  = make_matrix_ongpu( N, 1 );
    gCpu->tVN2  = make_matrix_ongpu( N, 1 );

printf("gCpu->X=%p\n", gCpu->X );
printf("gCpu->gen_a=%p \n", gCpu->gen_a );
printf("gCpu->gen_d=%p \n", gCpu->gen_d );
printf("Test %p +1 %p +2=%p\n", &(gCpu->pNext), &(gCpu->pNext) + 1, &(gCpu->pNext) + N  );

    gCpu->all_corTimes = make_matrix_list( (double**)(&(gCpu->pNext)), N, Q, Q );
    gCpu->all_corMat = make_matrix_list( (double**)( gCpu->all_corTimes + N), N, Q, Q );
    gCpu->all_corMat_MH = make_matrix_list( (double**)( gCpu->all_corMat + N), N, Q, Q );
    gCpu->all_corMat_Inv = make_matrix_list( (double**)( gCpu->all_corMat_MH + N), N, Q, Q );
    gCpu->all_corMat_MH_Inv = make_matrix_list( (double**)( gCpu->all_corMat_Inv + N), N, Q, Q );

printf("Test p1=%p p2=%p  p3==%p\n", gCpu->all_corTimes, gCpu->all_corMat, gCpu->all_corMat_MH  );

    gCpu->all_yi  = make_vector_list( (double**)( gCpu->all_corMat_MH_Inv + N), N, Q );
    gCpu->all_rd  = make_vector_list( (double**)( gCpu->all_yi + N), N, Q );
    gCpu->all_ui  = make_matrix_list( (double**)( gCpu->all_rd + N), N, Q, LG  );

    gCpu->tempMat  = make_matrix_list( (double**)(gCpu->all_ui + N), N, Q, Q  );
    gCpu->tmp2     = make_matrix_list( (double**)(gCpu->tempMat + N), N, LG, LG );
    gCpu->tmp3     = make_matrix_list( (double**)(gCpu->tmp2 + N), N, 1, LG );
    
    memcpy(gGpuMap, gCpu, nsize); 

printf("gCpu            =%p\n", gCpu);
printf("gGpuObj         =%p\n", gGpuObj);

    #define MAP_ADDR(x)  (double**)((char*)gGpuObj + (unsigned int)((char*)(x) - (char*)gCpu) )

    gGpuMap->all_corTimes = MAP_ADDR( gCpu->all_corTimes );
    gGpuMap->all_corMat = MAP_ADDR( gCpu->all_corMat );
    gGpuMap->all_corMat_MH = MAP_ADDR( gCpu->all_corMat_MH );
    gGpuMap->all_corMat_Inv = MAP_ADDR( gCpu->all_corMat_Inv );
    gGpuMap->all_corMat_MH_Inv = MAP_ADDR( gCpu->all_corMat_MH_Inv );
    gGpuMap->all_yi   = MAP_ADDR( gCpu->all_yi );
    gGpuMap->all_rd   = MAP_ADDR( gCpu->all_rd );
    gGpuMap->all_ui   = MAP_ADDR( gCpu->all_ui );

    gGpuMap->tempMat  = MAP_ADDR( gCpu->tempMat );
    gGpuMap->tmp2     = MAP_ADDR( gCpu->tmp2 );
    gGpuMap->tmp3     = MAP_ADDR( gCpu->tmp3 );

}

int Init_GPUobj(struct GPUobj** pCpuObj, struct GPUobj** pGpuObj, struct GPUobj** pGpuMap, int N, int P, int Q, int NC)
{
    int ngCudaSize = sizeof(struct GPUobj) + 50* N * sizeof(double*);

    //cublasHandle_t h;
    //cublasCreate(&h);
    //cublasSetPointerMode(h, CUBLAS_POINTER_MODE_DEVICE);
 
    struct GPUobj* gGpuObj;
    PERR( cudaMalloc( (void **)&gGpuObj, ngCudaSize ) );
    *pGpuObj = gGpuObj;

    struct GPUobj* gCpuCopy = (struct GPUobj*)Calloc( ngCudaSize, char);
    *pCpuObj = gCpuCopy;

    struct GPUobj* gGpuMap = (struct GPUobj*)Calloc( ngCudaSize, char);
    *pGpuMap = gGpuMap;

    _initGPUobj( gCpuCopy, gGpuObj, gGpuMap, ngCudaSize, N, P, Q, NC );

    PERR( cudaMemcpy( gGpuObj, gGpuMap, ngCudaSize, cudaMemcpyHostToDevice ) );
    
printf("CPU copy=%p\n", *pCpuObj);
printf("GPU addr=%p\n", *pGpuObj);
printf("GPU map =%p\n", *gGpuMap);
    
    return(0);
}

void _freeGPUobj(struct GPUobj* pCpuObj, int N )
{
    PERR( cudaFree( pCpuObj->X ) );
    PERR( cudaFree( pCpuObj->Z ) );
    PERR( cudaFree( pCpuObj->Z0 ) );
    PERR( cudaFree( pCpuObj->mInZ ) );
    PERR( cudaFree( pCpuObj->mu ) );
    PERR( cudaFree( pCpuObj->alpha ) );
    PERR( cudaFree( pCpuObj->a ) );
    PERR( cudaFree( pCpuObj->d ) );
    PERR( cudaFree( pCpuObj->a_old ) );
    PERR( cudaFree( pCpuObj->d_old ) );
    PERR( cudaFree( pCpuObj->vctP ) );
    PERR( cudaFree( pCpuObj->tau2 ) );
    PERR( cudaFree( pCpuObj->tau2_x ) );
    PERR( cudaFree( pCpuObj->gen_a ) );
    PERR( cudaFree( pCpuObj->gen_d ) );
    PERR( cudaFree( pCpuObj->tVN0 ) );
    PERR( cudaFree( pCpuObj->tVN1 ) );
    PERR( cudaFree( pCpuObj->tVN2 ) );

    for(int i=0;i<N;i++)
    {
        PERR( cudaFree( pCpuObj->all_corTimes[i] ) );
        PERR( cudaFree( pCpuObj->all_corMat[i] ) );
        PERR( cudaFree( pCpuObj->all_corMat_MH[i] ) );
        PERR( cudaFree( pCpuObj->all_corMat_Inv[i] ) );
        PERR( cudaFree( pCpuObj->all_corMat_MH_Inv[i] ) );
        PERR( cudaFree( pCpuObj->all_yi[i] ) );
        PERR( cudaFree( pCpuObj->all_rd[i] ) );
        PERR( cudaFree( pCpuObj->all_ui[i] ) );
        PERR( cudaFree( pCpuObj->tmp2[i] ) );
        PERR( cudaFree( pCpuObj->tmp3[i] ) );
        PERR( cudaFree( pCpuObj->tempMat[i] ) );
    }
    
}
    
int Free_GPUobj(struct GPUobj* pGpuObj, struct GPUobj* pCpuObj, struct GPUobj* gGpuMap, int N)
{
    _freeGPUobj( pCpuObj, N );

    PERR( cudaFree(pGpuObj) );

    Free(gGpuMap); 

    Free(pCpuObj); 

    return(0);
}

__global__ void g_part1(struct GPUobj* gCuda, int N, double rho, double tmp5 )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < N)
    {
        int m = (int)_V( gCuda->mInZ, i) ;

        double* Z0 = gCuda->all_corTimes[i];

        for (int ii=0; ii<m; ii++)
          for (int jj=00; jj<m; jj++)
          {
              _M( gCuda->all_corMat[i], ii, jj ) =  pow(rho, _M(Z0, ii, jj));  
              _M( gCuda->all_corMat_MH[i], ii, jj ) =  pow(tmp5, _M(Z0, ii, jj));
          }    
    } 
}

int _cuda_gpart1( struct GPUobj* gCuda, struct GPUobj* gCpuObj, int N, double rho, double tmp5)
{   
//printf("part1 %p, %d, %f, %f %d, %d \n", gCuda, N, rho, tmp5, (N + (THREADS_PER_BLOCK-1)) / THREADS_PER_BLOCK, THREADS_PER_BLOCK);

    g_part1<<< (N + (THREADS_PER_BLOCK-1)) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>(gCuda, N, rho, tmp5 );
    cudaDeviceSynchronize();
    ERRCHECK;

//printf("End of part1\n");
    // dont need copy all_corMat and all_corMat_MH back to CPU
    return(0);
}

__global__ void g_part2(struct GPUobj* gCuda, int N, int Q, int nC, double sigma2 )
{    
    extern __shared__ double smb[];

    //int i = blockIdx.x * blockDim.x + threadIdx.x;
    int col = threadIdx.x;
    int idx = blockIdx.x;

    if(idx < N)
    {
#define tmp          gShare->tMA
#define tmp4         gShare->tMB
#define tui          gShare->tMC
#define tempMat0     gShare->tempMat0
#define tempMat1     gShare->tempMat1
#define tempMat2     gShare->tempMat2

        GPUShare* gShare = (GPUShare*)smb;
        if(col==0) gShare = ConvetSharePoint( smb, Q );
        __syncthreads();

        Matrix_Transpose_Col( col, gCuda->all_ui[idx], tui);
        __syncthreads();

        Matrix_mult_Double_Col( col, gCuda->all_corMat_Inv[idx], 1/sigma2, tempMat0 );
        __syncthreads();

        Matrix_mult_Matrix_Col( col, tempMat0, gCuda->all_ui[idx], tmp );
         __syncthreads();

        Matrix_Resize_Col( col, tmp4, 1.0,  Vector_GetLength(gCuda->all_rd[idx]), TRUE );
        __syncthreads();

       if(nC > 0)
            for( int nX = 0; nX < 1; nX++ )
            {
               Matrix_GetRow_Col( col, gCuda->alpha, nX, tempMat0 ); 
               __syncthreads();

               Matrix_mult_Matrix_Col( col, tempMat0, tui, tempMat1);
               __syncthreads();

               Matrix_mult_Double_Col( col, tempMat1, _M(gCuda->X, idx, nX+1), tempMat1 );
               __syncthreads();

               Matrix_add_Matrix_Col( col, tempMat1, tmp4, tmp4 );
               __syncthreads();
            }
            
        __syncthreads();

        Matrix_Transpose_Col( col, gCuda->all_rd[idx], tempMat1 ); 
        __syncthreads();

        Matrix_sub_Matrix_Col( col, tempMat1, tmp4, tempMat1 );
        __syncthreads();

        //tmp3
        Matrix_mult_Matrix_Col( col, tempMat1, tmp, tempMat0);
        __syncthreads();
        Matrix_Copy_Col( col, gCuda->tmp3[idx], tempMat0);
        __syncthreads();

        //tmp2
        Matrix_mult_Matrix_Col( col, tui, tmp, tempMat0);
        __syncthreads();
        Matrix_Copy_Col( col, gCuda->tmp2[idx], tempMat0);
        __syncthreads();

#undef tempMat0 
#undef tempMat1 
#undef tempMat2 
#undef tmp
#undef tui
#undef tmp4
   }
}

int _cuda_gpart2( struct GPUobj* gCuda, struct GPUobj* gCpuObj, struct GPUobj* gGpuMap, int N, int Q, int nC, double sigma2, 
                 CFmMatrix& alpha, CFmMatrix& tmp2, CFmMatrix& tmp3 )
{
//printf("part2 %p, %d, %f, %f %d, %d \n", gCuda, N, nC, sigma2 );
    _copy_fmMatrix_Device( gCpuObj->alpha, alpha );
    
    int nShareSize = ((Q*Q+3)*10 + (Q+3)*10 ) * sizeof(double);
    
    g_part2<<< N,  16, nShareSize >>>(gCuda, N, Q, nC, sigma2 );
    cudaDeviceSynchronize();
    ERRCHECK;

    // tmp2 
    g_reduce_matrix( gGpuMap->tmp2, gGpuMap->tempMat, N );
    _copyback_fmMatrix_Host( tmp2, gCpuObj->tempMat[0]);

    // tmp3  
    g_reduce_matrix( gGpuMap->tmp3, gGpuMap->tempMat, N );
    _copyback_fmMatrix_Host( tmp3, gCpuObj->tempMat[0]);
    
//printf("End of part2\n");

    return(0);
}
 
__global__ void g_part3(struct GPUobj* gCuda, int N, int Q, int P )
{    
    extern __shared__ double smb[];

    //int i = blockIdx.x * blockDim.x + threadIdx.x;
    int col = threadIdx.x;
    int idx = blockIdx.x;

    if(idx < N)
    {
#define ui           gShare->tMA
#define mean_effect  gShare->tMB
#define tempMat0     gShare->tempMat0
#define tempMat1     gShare->tempMat1
#define tempMat2     gShare->tempMat2

        GPUShare* gShare = (GPUShare*)smb;
        if(col==0) gShare = ConvetSharePoint( smb, Q );
        __syncthreads();

         Matrix_Resize_Col( col, mean_effect, Vector_GetLength( gCuda->all_yi[idx]), 1, TRUE );
        __syncthreads();
         
        Matrix_Copy_Col( col, ui, gCuda->all_ui[idx]);
        __syncthreads();
         
         for(int jj=0; jj < P; jj++)
         {
             double a0 = _M(gCuda->gen_a, jj, idx);
             if(  a0 != 0)
             {
                 Matrix_GetRow_Col( col, gCuda->a, jj, tempMat0 );
                 __syncthreads();

                 Matrix_Transpose_Col( col, tempMat0, tempMat1 );
                 __syncthreads();

                 Matrix_mult_Matrix_Col( col, ui, tempMat1, tempMat0 );
                 __syncthreads();

                 Matrix_mult_Double_Col( col, tempMat0, a0, tempMat0 );
                 __syncthreads();

                 Matrix_add_Matrix_Col( col, mean_effect, tempMat0, mean_effect); 
                 __syncthreads();
             }
             
             double d0 = _M(gCuda->gen_d, jj, idx);
             if( d0 != 0)
             {
                 Matrix_GetRow_Col( col, gCuda->d, jj, tempMat0 );
                 __syncthreads();

                 Matrix_Transpose_Col( col, tempMat0, tempMat1 );
                 __syncthreads();

                 Matrix_mult_Matrix_Col( col, ui, tempMat1, tempMat0 );
                 __syncthreads();

                 Matrix_mult_Double_Col( col, tempMat0, d0, tempMat0 );
                 __syncthreads();

                 Matrix_add_Matrix_Col( col, mean_effect, tempMat0, mean_effect ); 
                 __syncthreads();
             }   
         }
         
         Matrix_sub_Matrix_Col( col, gCuda->all_yi[idx], mean_effect, gCuda->all_rd[idx]);
         __syncthreads(); 

#undef ui         
#undef mean_effect
#undef tempMat0
#undef tempMat1
#undef tempMat2
    }
}
         
int _cuda_gpart3( struct GPUobj* gCuda, struct GPUobj* gCpuObj, int N, int Q, int P, CFmMatrix& a, CFmMatrix& d)
{
//printf("part3 %p, %d, %d\n", gCuda, N, P );
    _copy_fmMatrix_Device( gCpuObj->a, a );
    _copy_fmMatrix_Device( gCpuObj->d, d );
    
    int nShareSize = ((Q*Q+3)*10 + (Q+3)*10 ) * sizeof(double);    
    g_part3<<< N, 16, nShareSize >>>( gCuda, N, Q, P );
    cudaDeviceSynchronize();
    ERRCHECK;
    
//printf("End of part3\n");
    return (0);
}


__global__ void g_part4(struct GPUobj* gCuda, int N, int Q, int j, int nC, double sigma2 )
{    
    extern __shared__ double smb[];

    //int i = blockIdx.x * blockDim.x + threadIdx.x;
    int col = threadIdx.x;
    int idx = blockIdx.x;
    
    if(idx < N)
    {
#define tmp          gShare->tMA
#define tui          gShare->tMC
#define tmp4         gShare->tMB
#define ui           gShare->tMD
#define tempMat0     gShare->tempMat0
#define tempMat1     gShare->tempMat1
#define tempMat2     gShare->tempMat2

        GPUShare* gShare = (GPUShare*)smb;
        if(col==0) gShare = ConvetSharePoint( smb, Q );
        __syncthreads();

        Matrix_Resize_Col( col, tmp4, 1, Vector_GetLength(gCuda->all_rd[idx]), TRUE );
        __syncthreads();

        Matrix_Transpose_Col( col, gCuda->all_ui[idx], tui );
        __syncthreads();
        Matrix_Transpose_Col( col, tui, ui );
        __syncthreads();

        double a0 = _M(gCuda->gen_a, j, idx);
        if( a0 != 0.0)
        {
            if(nC > 0)
            {
                for(int nX = 0; nX < nC; nX++)
                {
                    Matrix_GetRow_Col( col, gCuda->alpha, nX, tempMat0 ); 
                    __syncthreads();
                    Matrix_mult_Matrix_Col( col, tempMat0, tui, tempMat1 );
                    __syncthreads();
                    Matrix_mult_Double_Col( col, tempMat1, _M(gCuda->X, idx, nX+1), tempMat1 );
                    __syncthreads();
                    Matrix_add_Matrix_Col( col, tempMat1, tmp4, tmp4 );
                    __syncthreads();
                }
             }
             
             Matrix_mult_Double_Col( col, gCuda->all_corMat_Inv[idx], a0/sigma2 , tempMat0);
             __syncthreads();
             Matrix_mult_Matrix_Col( col, tempMat0, ui, tmp);
             __syncthreads();
                 
             Matrix_Transpose_Col( col, gCuda->mu, tempMat0);
             __syncthreads();
             Matrix_mult_Matrix_Col( col, tempMat0, tui, tempMat1 );
             __syncthreads();

             
             Matrix_Transpose_Col( col, gCuda->all_rd[idx], tempMat2); 
             __syncthreads();
             Matrix_sub_Matrix_Col( col, tempMat2, tempMat1, tempMat1);
             __syncthreads();
             Matrix_sub_Matrix_Col( col, tempMat1, tmp4, tempMat0);
             __syncthreads();

             Matrix_GetRow_Col( col, gCuda->a, j, tempMat1 );
             __syncthreads();
             Matrix_mult_Matrix_Col( col, tempMat1, tui, tempMat2);
             __syncthreads();
             Matrix_mult_Double_Col(col, tempMat2, a0, tempMat1);
             __syncthreads();

             Matrix_add_Matrix_Col( col, tempMat0, tempMat1, tempMat0 );
             __syncthreads();
             Matrix_mult_Matrix_Col( col, tempMat0, tmp , gCuda->tmp3[idx]);
             __syncthreads();
             
             Matrix_mult_Matrix_Col( col, tui, tmp, tempMat0);
             __syncthreads();
             Matrix_mult_Double_Col( col, tempMat0, a0, gCuda->tmp2[idx]); 
             __syncthreads();
        }
        else
        {
             Matrix_Resize_Col( col, gCuda->tmp2[idx], LG, LG, TRUE);
             __syncthreads();
             Matrix_Resize_Col( col, gCuda->tmp3[idx], 1, LG, TRUE);
             __syncthreads();
        }   

#undef tempMat0 
#undef tempMat1 
#undef tempMat2 
#undef tmp4
#undef tui         
#undef ui 
#undef tmp 
        
    }
}


int _cuda_gpart4( struct GPUobj* gCuda, struct GPUobj* gCpuObj, struct GPUobj* gGpuMap, int N, int Q, int j, int nC, double sigma2, 
                 CFmVector& mu, CFmMatrix& alpha, CFmMatrix& a, CFmMatrix& tmp2, CFmMatrix& tmp3)
{
//printf("part4 %p, %d, %d, %d \n", gCuda, N, j, nC );
    if (j==0)
    {
        _copy_fmMatrix_Device( gCpuObj->a, a);
        _copy_fmVector_Device( gCpuObj->mu, mu);
        _copy_fmMatrix_Device( gCpuObj->alpha, alpha);
    }
    
    int nShareSize = ((Q*Q+3)*10 + (Q+3)*10 ) * sizeof(double);    
    g_part4<<< N, 16, nShareSize >>>(gCuda, N, Q, j, nC,sigma2 );
    cudaDeviceSynchronize();
    ERRCHECK;

    // tmp2 and tmp 3 
    g_reduce_matrix( gGpuMap->tmp2, gGpuMap->tempMat, N );
    _copyback_fmMatrix_Host( tmp2, gCpuObj->tempMat[0] );

    g_reduce_matrix( gGpuMap->tmp3, gGpuMap->tempMat, N );
    _copyback_fmMatrix_Host( tmp3, gCpuObj->tempMat[0] );

//printf("End of part4\n");

    return(0);
}


__global__ void g_part5(struct GPUobj* gCuda, int N, int Q, int j )
{    
    extern __shared__ double smb[];

    //int i = blockIdx.x * blockDim.x + threadIdx.x;
    int col = threadIdx.x;
    int idx = blockIdx.x;
    
    if(idx < N)
    {
#define ui       gShare->tMA
#define rd       gShare->tMB
#define tempMat0     gShare->tempMat0
#define tempMat1     gShare->tempMat1
#define tempMat2     gShare->tempMat2

        GPUShare* gShare = (GPUShare*)smb;
        if(col==0) gShare = ConvetSharePoint( smb, Q );
        __syncthreads();

        Matrix_Copy_Col( col, ui, gCuda->all_ui[idx] );
        Matrix_Copy_Col( col, rd, gCuda->all_rd[idx] );
        __syncthreads();

        double gen_a = _M(gCuda->gen_a, j, idx);  
        if(  gen_a != 0.0 )
        {
             Matrix_GetRow_Col( col, gCuda->a_old, j, tempMat0);
             __syncthreads();
             Matrix_Transpose_Col( col, tempMat0, tempMat1);
             __syncthreads();
             Matrix_mult_Matrix_Col( col, ui, tempMat1, tempMat2);
             __syncthreads();
             
             //to save up transfer time, here only transfer one vector into VN0
             //Matrix_GetRow_Col( col, gCuda->a, j, tempMat0);
             //__syncthreads();
             //Matrix_Transpose_Col( col, tempMat0, tempMat1 );
             //__syncthreads();
             
             Matrix_Copy_Col( col, tempMat1, gCuda->tVN0);
             __syncthreads();
             Matrix_mult_Matrix_Col( col, ui, tempMat1, tempMat0);
             __syncthreads();
             Matrix_sub_Matrix_Col( col, tempMat2, tempMat0, tempMat0);
             __syncthreads();
             
             Matrix_mult_Double_Col( col, tempMat0, gen_a, tempMat0);
             __syncthreads();
         
             Matrix_add_Matrix_Col( col, rd, tempMat0, gCuda->all_rd[idx] );
             __syncthreads();
        }

#undef ui       
#undef rd       
#undef tempMat0     
#undef tempMat1     
#undef tempMat2     

    }
}

int _cuda_gpart5( struct GPUobj* gCuda, struct GPUobj* gCpuObj, int N, int Q, int j, CFmMatrix& a, CFmMatrix& a_old)
{
    //only fully update a and a_old at the first time to save up computational times.
    if(j==0)
        _copy_fmMatrix_Device( gCpuObj->a_old, a_old);

    CFmVector new_aj(LG, 0); 
    new_aj = a.GetRow(j);
    _copy_fmVector_Device(gCpuObj->tVN0, new_aj );
    
    int nShareSize = ((Q*Q+3)*10 + (Q+3)*10 ) * sizeof(double);

    g_part5<<< N, 16, nShareSize >>>(gCuda, N, Q, j );
    cudaDeviceSynchronize();
    ERRCHECK;

//printf("End of part5\n");

    return(0);
} 



__global__ void g_part6(struct GPUobj* gCuda, int P, double sigma2, double lambda2, double lambda2_x)
{    
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < P )
    {
#define tau2    gCuda->tau2
#define tau2_x  gCuda->tau2_x

        double vp = _V( gCuda->vctP, i ); 
        if( vp == 2.0)
        {
            double InvTau2_1 = sqrt( LG * lambda2 * sigma2/ Matrix_RowProd( gCuda->a, i, i));
            _V( tau2, i) = 1/_cuda_invGau( i*10, InvTau2_1, LG * lambda2);
            _V( tau2_x, i) = 0.0;
         }
         else if(vp == 1.0)
         {
             double InvTau2_1 = sqrt( LG * lambda2_x * sigma2/Matrix_RowProd( gCuda->a, i, i) );
             _V( tau2_x, i) = 1/_cuda_invGau( i*10+1, InvTau2_1, LG*lambda2_x);
             _V( tau2, i) = 0.0;
         }
         else
         {
             _V( tau2, i) = 0.0;
             _V( tau2_x, i) = 0.0;
         }

#undef tau2
#undef tau2_x
    }
}

int _cuda_gpart6( struct GPUobj* gCuda, struct GPUobj* gCpuObj, int P, double sigma2, double lambda2, double lambda2_x, 
                 CFmVector& vctP, CFmMatrix& a, CFmVector& tau2, CFmVector& tau2_x )
{
    _copy_fmVector_Device( gCpuObj->vctP, vctP);
    _copy_fmMatrix_Device( gCpuObj->a, a);
    
    g_part6<<< (P + (THREADS_PER_BLOCK-1)) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>( gCuda, P, sigma2, lambda2, lambda2_x );
    cudaDeviceSynchronize();
    ERRCHECK;
    
    // tau2 and tau2_x
    _copyback_fmVector_Host(tau2, gCpuObj->tau2 );
    _copyback_fmVector_Host(tau2_x, gCpuObj->tau2_x );
    
    return(0);
}

__global__ void g_part7(struct GPUobj* gCuda, int N, int Q, int j, int nC, double sigma2)
{    
    extern __shared__ double smb[];

    //int i = blockIdx.x * blockDim.x + threadIdx.x;
    int col = threadIdx.x;
    int idx = blockIdx.x;
    
    if(idx < N)
    {
#define tmp          gShare->tMA
#define tmp4         gShare->tMB
#define tui          gShare->tMC
#define gtmp2        gShare->tMD 
#define gtmp3        gShare->tME
#define tempMat0     gShare->tempMat0
#define tempMat1     gShare->tempMat1
#define tempMat2     gShare->tempMat2

        GPUShare* gShare = (GPUShare*)smb;
        if(col==0) gShare = ConvetSharePoint( smb, Q );
        __syncthreads();

        Matrix_Resize_Col( col, tmp4, 1, Vector_GetLength(gCuda->all_rd[idx]), TRUE );
        Matrix_Resize_Col( col, gtmp2, LG, LG, TRUE );
        Matrix_Resize_Col( col, gtmp3, 1, LG, TRUE );
        Matrix_Transpose_Col( col, gCuda->all_ui[idx], tui );
        __syncthreads();

        double d0 = _M(gCuda->gen_d, j, idx);
        if( d0 != 0.0 )
        {
             if(nC > 0)
             {
                 for(int nX = 0; nX < nC; nX++)
                 {
                     Matrix_GetRow_Col( col, gCuda->alpha, nX, tempMat0); 
                     __syncthreads();
                     Matrix_mult_Matrix_Col( col, tempMat0, tui, tempMat1);
                     __syncthreads();
                     Matrix_mult_Double_Col( col, tempMat1, _M( gCuda->X, idx, nX+1), tempMat1 );
                     __syncthreads();
                     Matrix_add_Matrix_Col( col, tempMat1, tmp4, tmp4);
                     __syncthreads();
                 }
             }
                 
             Matrix_mult_Double_Col( col, gCuda->all_corMat_Inv[idx], d0/sigma2, tempMat0 );
             __syncthreads();
             Matrix_Transpose_Col( col, tui, tempMat1 );
             __syncthreads();
             Matrix_mult_Matrix_Col( col, tempMat0, tempMat1, tmp);
             __syncthreads();
                 
             Matrix_Transpose_Col( col, gCuda->mu, tempMat0 );
             __syncthreads();
             Matrix_mult_Matrix_Col( col, tempMat0, tui, tempMat1 );
             __syncthreads();
             Matrix_Transpose_Col( col, tempMat1, tempMat0 );
             __syncthreads();

             Matrix_Transpose_Col( col, tmp4, tempMat1 );
             __syncthreads();

             Matrix_sub_Matrix_Col( col, gCuda->all_rd[idx], tempMat1, tempMat1 );
             __syncthreads();
             Matrix_sub_Matrix_Col( col, tempMat1, tempMat0, tempMat0 );
             __syncthreads();

             Matrix_GetRow_Col( col, gCuda->d, j, tempMat2 );
             __syncthreads();
             Matrix_mult_Matrix_Col( col, tempMat2, tui, tempMat1 );
             __syncthreads();
             Matrix_mult_Double_Col( col, tempMat1, d0, tempMat1 );
             __syncthreads();
             Matrix_Transpose_Col( col, tempMat1, tempMat2 );
             __syncthreads();
    
             Matrix_add_Matrix_Col( col, tempMat0, tempMat2, tempMat1 );
             __syncthreads();
             Matrix_Transpose_Col( col, tempMat1, tempMat0 );
             __syncthreads();
             Matrix_mult_Matrix_Col( col, tempMat0, tmp, gtmp3);
             __syncthreads();
             
             Matrix_mult_Matrix_Col( col, tui, tmp, tempMat0 );
             __syncthreads();
             Matrix_mult_Double_Col( col, tempMat0, d0,  gtmp2);
             __syncthreads();
        }

        Matrix_Copy_Col( col, gCuda->tmp2[idx], gtmp2);
        __syncthreads();
        Matrix_Copy_Col( col, gCuda->tmp3[idx], gtmp3);
        __syncthreads();


#undef tmp          
#undef tmp4       
#undef tui      
#undef tempMat0     
#undef tempMat1
#undef tempMat2
#undef gtmp2
#undef gtmp3

    }
}


int _cuda_gpart7( struct GPUobj* gCuda, struct GPUobj* gCpuObj, struct GPUobj* gGpuMap, int N, int Q, int j, int nC, double sigma2, 
                 CFmMatrix& alpha, CFmVector& mu, CFmMatrix& d, CFmMatrix& tmp2, CFmMatrix& tmp3 )
{
    if (j==0)
    {
        _copy_fmMatrix_Device( gCpuObj->alpha, alpha);
        _copy_fmMatrix_Device( gCpuObj->d, d);
        _copy_fmVector_Device( gCpuObj->mu, mu);
    }
    
    int nShareSize = ((Q*Q+3)*10 + (Q+3)*10 ) * sizeof(double);
    g_part7<<< N, 16, nShareSize>>>(gCuda, N, Q, j, nC, sigma2 );
    cudaDeviceSynchronize();
    ERRCHECK;

    // tmp2 and tmp3 
    g_reduce_matrix( gGpuMap->tmp2, gGpuMap->tempMat, N );
    _copyback_fmMatrix_Host( tmp2, gCpuObj->tempMat[0] );

    g_reduce_matrix( gGpuMap->tmp3, gGpuMap->tempMat, N );
    _copyback_fmMatrix_Host( tmp3, gCpuObj->tempMat[0]);
    
    return(0);
}


__global__ void g_part8( struct GPUobj* gCuda, int N, int Q, int j )
{    
    extern __shared__ double smb[];

    //int i = blockIdx.x * blockDim.x + threadIdx.x;
    int col = threadIdx.x;
    int idx = blockIdx.x;
    
    if(idx < N)
    {
#define ui           gShare->tMA
#define tempMat0     gShare->tempMat0
#define tempMat1     gShare->tempMat1
#define tempMat2     gShare->tempMat2

        GPUShare* gShare = (GPUShare*)smb;
        if(col==0) gShare = ConvetSharePoint( smb, Q );
        __syncthreads();
        
        Matrix_Copy_Col( col, ui, gCuda->all_ui[idx] );
        __syncthreads();

        double d0 = _M(gCuda->gen_d, j, idx );
        if( d0!= 0)
        {
             Matrix_GetRow_Col( col, gCuda->d_old, j, tempMat0 );
             __syncthreads();
             Matrix_Transpose_Col( col, tempMat0, tempMat1 );
            __syncthreads();
             Matrix_mult_Matrix_Col( col, ui, tempMat1, tempMat0);
            __syncthreads();

             // Save time of data transfer.
             //Matrix_GetRow_Col( col, gCuda->d, j, tempMat1 );
             //__syncthreads();
             //Matrix_Transpose_Col( col, tempMat1, tempMat2 );
             //__syncthreads();

             Matrix_Copy_Col( col, tempMat2, gCuda->tVN0);
             __syncthreads();
             Matrix_mult_Matrix_Col( col, ui, tempMat2, tempMat1 );
             __syncthreads();
             
             Matrix_sub_Matrix_Col( col, tempMat0, tempMat1, tempMat0 );
             __syncthreads();
             Matrix_Transpose_Col( col, tempMat0, tempMat1);
            __syncthreads();
             
             Matrix_mult_Double_Col( col, tempMat1, d0, tempMat1 );
             __syncthreads();

             Matrix_Transpose_Col( col, tempMat1, tempMat0 );
             __syncthreads();
         
             Matrix_add_Matrix_Col( col, gCuda->all_rd[idx], tempMat0, gCuda->all_rd[idx]);
             __syncthreads();
        }
#undef ui      
#undef tempMat0     
#undef tempMat1
#undef tempMat2
    }
}

int _cuda_gpart8( struct GPUobj* gCuda, struct GPUobj* gCpuObj, int N, int Q, int j, CFmMatrix& d, CFmMatrix& d_old)
{
    if(j==0)
        _copy_fmMatrix_Device( gCpuObj->d_old, d_old );
        

    //_copy_fmMatrix_Device( gCpuObj->d, d );
    CFmVector new_dj(LG, 0); 
    new_dj = d.GetRow(j);
    _copy_fmVector_Device(gCpuObj->tVN0, new_dj );
    
    
    int nShareSize = ((Q*Q+3)*10 + (Q+3)*10 ) * sizeof(double);
    g_part8<<< N, 16, nShareSize >>>( gCuda, N, Q, j );
    cudaDeviceSynchronize();
    ERRCHECK;

    return(0);
}


__global__ void g_part9(struct GPUobj* gCuda, int P, double lambda2_st, double lambda2_st_x, double sigma2 )
{    
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i < P )
    {
#define tau2_st    gCuda->tau2
#define tau2_st_x  gCuda->tau2_x

        double vp = _V( gCuda->vctP, i ); 
        if( vp == 2.0)
        {
            double InvTau2_1 = sqrt( LG * lambda2_st * sigma2/ Matrix_RowProd(gCuda->d, i, i));
            _V( tau2_st, i) = 1/_cuda_invGau( i*10, InvTau2_1, LG * lambda2_st);
            _V( tau2_st_x, i) = 0.0;
        }
        else if( vp == 1.0)
        {
            double InvTau2_1 = sqrt( LG * lambda2_st_x * sigma2/Matrix_RowProd(gCuda->d, i, i) );
            _V( tau2_st_x, i) = 1/_cuda_invGau( i*10+1, InvTau2_1, LG*lambda2_st_x);
            _V( tau2_st, i) = 0;
        }
        else
        {
            _V( tau2_st, i) = 0;
            _V( tau2_st_x, i) = 0;
        }
#undef tau2_st
#undef tau2_st_x 
    }
}


int _cuda_gpart9( struct GPUobj* gCuda, struct GPUobj* gCpuObj, int P, double lambda_st2, double lambda_st2_x, double sigma2, 
                 CFmVector& vctP, CFmMatrix& d, CFmVector& tau2_st, CFmVector& tau2_st_x  )
{
    _copy_fmVector_Device( gCpuObj->vctP, vctP );
    _copy_fmMatrix_Device( gCpuObj->d, d );

    g_part9<<< (P + (THREADS_PER_BLOCK-1)) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>( gCuda, P, lambda_st2, lambda_st2_x, sigma2 );
    cudaDeviceSynchronize();
    ERRCHECK;
    
    _copyback_fmVector_Host( tau2_st, gCpuObj->tau2);
    _copyback_fmVector_Host( tau2_st_x, gCpuObj->tau2_x);
   
    return(0);
}


__global__ void g_part10( struct GPUobj* gCuda, int N, int Q, int nC, int nX, double sigma2 )
{    
    extern __shared__ double smb[];

    //int i = blockIdx.x * blockDim.x + threadIdx.x;
    int col = threadIdx.x;
    int idx = blockIdx.x;
    
    if(idx < N)
    {
#define tui      gShare->tMA
#define ui       gShare->tMB
#define tmp      gShare->tMC
#define tmp4     gShare->tMD 
#define tempMat0     gShare->tempMat0
#define tempMat1     gShare->tempMat1
#define tempMat2     gShare->tempMat2

        GPUShare* gShare = (GPUShare*)smb;
        if(col==0) gShare = ConvetSharePoint( smb, Q );
        __syncthreads();

         Matrix_Resize_Col( col, tmp4, 1.0, Vector_GetLength( gCuda->all_rd[idx] ), TRUE );
         Matrix_Copy_Col( col, ui, gCuda->all_ui[idx]);
        __syncthreads();

         Matrix_Transpose_Col( col, ui, tui );
        __syncthreads();

         if(nC > 0)
            for( int nX2 = 0; nX2 < nC; nX2++ )
            {
                if(nX2 != nX)
                {
                    Matrix_GetRow_Col( col, gCuda->alpha, nX2, tempMat0 ); 
                   __syncthreads();
                    
                    Matrix_mult_Matrix_Col( col, tempMat0, tui, tempMat1 );
                    __syncthreads();
                    
                    Matrix_mult_Double_Col( col, tempMat1, _M( gCuda->X, idx, nX2 + 1), tempMat1 );
                    __syncthreads();
                    
                    Matrix_add_Matrix_Col( col, tmp4, tempMat1, tmp4 );
                    __syncthreads();
                }  
            }
            
         double x0 = _M( gCuda->X, idx, nX+1);   
         
         Matrix_mult_Double_Col( col, gCuda->all_corMat_Inv[idx], 1/sigma2, tempMat0 );
         __syncthreads();
         Matrix_mult_Double_Col( col, tempMat0, x0, tempMat1 );
         __syncthreads();
         Matrix_mult_Matrix_Col( col, tempMat1, ui, tmp );
         __syncthreads();

         Matrix_Transpose_Col( col, gCuda->mu, tempMat0 );
         __syncthreads();
         Matrix_mult_Matrix_Col( col, tempMat0, tui, tempMat1 );
         __syncthreads();

         Matrix_Transpose_Col( col, gCuda->all_rd[idx], tempMat0 );
         __syncthreads();

         Matrix_sub_Matrix_Col( col, tempMat0, tempMat1, tempMat0 );
         __syncthreads();
         Matrix_sub_Matrix_Col( col, tempMat0, tmp4, tempMat0 );
         __syncthreads();
         Matrix_mult_Matrix_Col( col, tempMat0, tmp, gCuda->tmp3[idx] );
         __syncthreads();
  
         Matrix_mult_Matrix_Col( col, tui, tmp, tempMat0 );
         __syncthreads();
         Matrix_mult_Double_Col( col, tempMat0, x0, gCuda->tmp2[idx]);
         __syncthreads();

#undef tui         
#undef ui         
#undef rd       
#undef tmp  
#undef tmp4 
#undef tempMat0        
#undef tempMat1        
#undef tempMat2
     }
}

int _cuda_gpart10( struct GPUobj* gCuda, struct GPUobj* gCpuObj, GPUobj* gGpuMap, int N, int Q, int nC, int nX, double sigma2, 
                  CFmMatrix& alpha, CFmVector& mu, CFmMatrix& tmp2, CFmMatrix& tmp3 )
{
    _copy_fmMatrix_Device( gCpuObj->alpha, alpha);
    _copy_fmVector_Device( gCpuObj->mu, mu);
    
    int nShareSize = ((Q*Q+3)*10 + (Q+3)*10 ) * sizeof(double);
    g_part10<<< N, 16, nShareSize >>>( gCuda, N, Q, nC, nX, sigma2 );
    cudaDeviceSynchronize();
    ERRCHECK;

    // tmp2 
    g_reduce_matrix( gGpuMap->tmp2, gGpuMap->tempMat, N );
    _copyback_fmMatrix_Host( tmp2, gCpuObj->tempMat[0]);
    
    // tmp3
    g_reduce_matrix( gGpuMap->tmp3, gGpuMap->tempMat, N );
    _copyback_fmMatrix_Host( tmp3, gCpuObj->tempMat[0]);

    return(0);
}


__global__ void g_part11( struct GPUobj* gCuda, int N, int Q, int nC )
{    
    extern __shared__ double smb[];

    //int i = blockIdx.x * blockDim.x + threadIdx.x;
    int col = threadIdx.x;
    int idx = blockIdx.x;
    
    if(idx < N)
    {
#define sigma2_scale gCuda->tVN0
#define tmp4         gShare->tMA
#define tmpv_sigma2  gShare->tMB
#define tmpv_get2    gShare->tMC
#define tui          gShare->tMD 
#define tempMat0     gShare->tempMat0
#define tempMat1     gShare->tempMat1
#define tempMat2     gShare->tempMat2


        GPUShare* gShare = (GPUShare*)smb;
        if(col==0) gShare = ConvetSharePoint( smb, Q );
        __syncthreads();

        Matrix_Resize_Col( col, tmp4, 1, Vector_GetLength(gCuda->all_rd[idx]), TRUE);
        Matrix_Transpose_Col( col, gCuda->all_ui[idx], tui );
        __syncthreads();
        
        if(nC > 0)
        {
            for(int nX = 0; nX < nC; nX++)
            {
                Matrix_GetRow_Col( col, gCuda->alpha, nX, tempMat0 );
                __syncthreads();
                Matrix_mult_Matrix_Col( col, tempMat0, tui, tempMat1 );
                __syncthreads();
                Matrix_mult_Double_Col( col, tempMat1, _M(gCuda->X, idx, nX+1), tempMat1 );
                __syncthreads();
                Matrix_add_Matrix_Col( col, tempMat1, tmp4, tmp4);
                __syncthreads();
            }
        }
        
        Matrix_Transpose_Col( col, gCuda->mu, tempMat0 );
        __syncthreads();
        Matrix_mult_Matrix_Col( col, tempMat0, tui, tempMat1 );
        __syncthreads();
        
        Matrix_Transpose_Col( col,  gCuda->all_rd[idx], tempMat0 );
        __syncthreads();
        Matrix_sub_Matrix_Col( col,  tempMat0, tempMat1, tempMat0 );
        __syncthreads();
        Matrix_sub_Matrix_Col( col,  tempMat0, tmp4, tmpv_sigma2 );
        __syncthreads();
        
        Matrix_mult_Matrix_Col( col, tmpv_sigma2, gCuda->all_corMat_Inv[idx], tempMat0 );
        __syncthreads();
        Matrix_Transpose_Col( col,  tmpv_sigma2, tempMat1 );
        __syncthreads();
        Matrix_mult_Matrix_Col( col, tempMat0, tempMat1, tmpv_get2 );
        __syncthreads();
        
        _V(sigma2_scale, idx ) = _M( tmpv_get2, 0, 0);

#undef sigma2_scale 
#undef tmp4         
#undef tmpv_sigma2  
#undef tmpv_get2    
#undef tui          
#undef tempMat0     
#undef tempMat1     
#undef tempMat2     

    }
}

int _cuda_gpart11( struct GPUobj* gCuda, struct GPUobj* gCpuObj, struct GPUobj* gGpuMap, int N, int Q, int nC, 
                  CFmMatrix& alpha, CFmVector& mu, double* sigma2_scale)
{
    _copy_fmMatrix_Device( gCpuObj->alpha, alpha);
    _copy_fmVector_Device( gCpuObj->mu, mu);
    
    int nShareSize = ((Q*Q+3)*10 + (Q+3)*10 ) * sizeof(double);
    g_part11<<< N, 16, nShareSize >>>( gCuda, N, Q, nC );
    cudaDeviceSynchronize();
    ERRCHECK;
    
    //gCuda->tVA <--> sigma2_scale
    g_reduce_double( gGpuMap->tVN0, gGpuMap->tVN1, N, sigma2_scale);

    if(1) 
    { 
        double dRet = g_reduce_double_test( gGpuMap->tVN0, N );
        if( (dRet - *sigma2_scale)> 1e-6 )
           printf("Part11 Failed to reduce(sigma2_scale C:%.5f, G:%.5f)\n", dRet, *sigma2_scale );
    }

    return(0);
}

__global__ void g_part12(struct GPUobj* gCuda, int N, int Q, int nC, double sigma2 )
{
    extern __shared__ double smb[];

    //int i = blockIdx.x * blockDim.x + threadIdx.x;
    int col = threadIdx.x;
    int idx = blockIdx.x;
    
    if(idx < N)
    {
#define exp_diff  gCuda->tVN0
#define tmp4      gShare->tMA
#define tui       gShare->tMB
#define exp_mat   gShare->tMC
#define tempMat0  gShare->tempMat0
#define tempMat1  gShare->tempMat1
#define tempMat2  gShare->tempMat2
#define allRd_i   gShare->tempMat4


        GPUShare* gShare = (GPUShare*)smb;
        if(col==0) gShare = ConvetSharePoint( smb, Q );
        __syncthreads();

        Matrix_Resize_Col( col, tmp4, 1, Vector_GetLength(gCuda->all_rd[idx]), TRUE);
        Matrix_Transpose_Col( col, gCuda->all_ui[idx], tui );
        __syncthreads();

        if(nC > 0)
        {
            for(int nX = 0; nX < nC; nX++)
            {
                Matrix_GetRow_Col( col, gCuda->alpha, nX, tempMat0);
                __syncthreads();
                Matrix_mult_Matrix_Col( col, tempMat0, tui, tempMat1);
                __syncthreads();
                Matrix_mult_Double_Col( col, tempMat1, _M( gCuda->X, idx, nX+1), tempMat1 );
                __syncthreads();
                Matrix_add_Matrix_Col( col, tmp4, tempMat1, tmp4 );
                __syncthreads();
            }
        }
        
        Matrix_Transpose_Col( col, gCuda->mu, tempMat0 );
        Matrix_mult_Matrix_Col( col, tempMat0, tui, tempMat2);

        Matrix_Transpose_Col( col, gCuda->all_rd[idx], tempMat0);
        Matrix_sub_Matrix_Col( col, tempMat0, tempMat2, tempMat0);
        Matrix_sub_Matrix_Col( col, tempMat0, tmp4, allRd_i );

        Matrix_mult_Double_Col( col, gCuda->all_corMat_MH_Inv[idx], 1/sigma2, tempMat0);
        Matrix_mult_Matrix_Col( col, allRd_i, tempMat0, tempMat1 );
        Matrix_Transpose_Col(col, allRd_i, tempMat2); 
        Matrix_mult_Matrix_Col( col, tempMat1, tempMat2, exp_mat );

        Matrix_mult_Double_Col( col, gCuda->all_corMat_Inv[idx], 1/sigma2, tempMat0);
        Matrix_mult_Matrix_Col( col, allRd_i, tempMat0, tempMat1 );
        Matrix_Transpose_Col( col, allRd_i, tempMat2); 
        Matrix_mult_Matrix_Col( col, tempMat1, tempMat2, tempMat0 );

        Matrix_sub_Matrix_Col( col, exp_mat, tempMat0, exp_mat );

        _V( exp_diff, idx ) = _M( exp_mat, 0, 0 );


#undef exp_diff  
#undef tmp4      
#undef tui       
#undef exp_mat   
#undef tempMat0  
#undef tempMat1  
#undef tempMat2  
#undef allRd_i   

    }
}

int _cuda_gpart12( struct GPUobj* gCuda, struct GPUobj* gCpuObj, struct GPUobj* gGpuMap, int N, int Q, int nC, double sigma2, 
                  CFmMatrix& alpha, CFmVector& mu, double* exp_diff)
{
    _copy_fmMatrix_Device( gCpuObj->alpha, alpha);
    _copy_fmVector_Device( gCpuObj->mu, mu);

    int nShareSize = ((Q*Q+3)*10 + (Q+3)*10 ) * sizeof(double);
    g_part12<<< N, 16, nShareSize >>>( gCuda, N, Q, nC, sigma2 );
    cudaDeviceSynchronize();
    ERRCHECK;

    g_reduce_double( gGpuMap->tVN0, gGpuMap->tVN1, N, exp_diff);
    
    return(0);
}

__global__ void g_showTmp3(struct GPUobj* gCuda, int N )
{    
    for(int i=0; i<N; i++)
    {
        Print_Matrix( gCuda->tmp3[i] );
    }
}
//g_showTmp3<<< 1,  1>>> (gCuda, N);

