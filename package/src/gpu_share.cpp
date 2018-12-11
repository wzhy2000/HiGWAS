#ifdef USECUDA

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

#include "gpu_share.h"

/*
__global__ void g_cuda_show(double* gMat )
{
    double* tx = gMat;
    printf("GPUADDR=%p  (%.1f %.1f %.1f)\n", tx,  tx[0], tx[1],  tx[2]);
    if(round(tx[2])!=1)
       for( int i=0; i<5; i++ )
            printf("[%f %f %f %f %f ]\n", _M(tx, i, 0), _M(tx, i, 1), _M(tx, i, 2 ), _M(tx, i, 3 ), _M(tx, i, 4 ) );
    else
        printf("[%f %f %f %f %f ]\n", _V(tx, 0), _V(tx, 1), _V(tx, 2 ), _V(tx, 3 ), _V(tx, 4 ) );
}

int _cuda_show( double* gMat)
{
    g_cuda_show<<< 1, 1>>>(gMat);
    cudaDeviceSynchronize();
    ERRCHECK;

    return(0);
}*/

//*run on CPU;
int _copy_fmMatrix_list( double** cudaList, int size, CFmMatrix** matList)
{
//printf("_copy_fmMatrix_list %p-->%p\n", &matList, cudaList);
    for(int i=0;i<size;i++)
    {
        double* pMat = cudaList[i];
        CFmMatrix* pFm = matList[i];
        double dRow=(double)(pFm->GetNumRows());
        double dCol=(double)(pFm->GetNumCols());

        PERR( cudaMemcpy( pMat+3, pFm->GetData(), (int)(dRow*dCol*sizeof(double)), cudaMemcpyHostToDevice ) );
        PERR( cudaMemcpy( pMat+2, &dCol, sizeof(double), cudaMemcpyHostToDevice ) );
        PERR( cudaMemcpy( pMat+1, &dRow, sizeof(double), cudaMemcpyHostToDevice ) );
    }

    return(0);
}

//*run on CPU;
int _copy_fmVector_list( double** cudaList, int size, CFmVector** vecList)
{
//printf("_copy_fmVector_list %p-->%p\n", &vecList, cudaList);

    for(int i=0;i<size;i++)
    {
        double* pMat = cudaList[i];
        CFmVector* pFm = vecList[i];
        double dRow=(double)(pFm->GetLength());
        double dCol=1.0;

        PERR( cudaMemcpy( pMat+3, pFm->GetData(), (int)(dRow*dCol*sizeof(double)), cudaMemcpyHostToDevice ) );
        PERR( cudaMemcpy( pMat+2, &dCol, sizeof(double), cudaMemcpyHostToDevice ) );
        PERR( cudaMemcpy( pMat+1, &dRow, sizeof(double), cudaMemcpyHostToDevice ) );
    }

    return(0);
}

//*run on CPU;
int _copy_fmVector_Device(double* gcudaDstVec, CFmVector& fmSrcVec )
{
//printf("_copy_fmVector_Device %p-->%p\n", &fmSrcVec, gcudaDstVec);

    double dSize=0.0;
    PERR( cudaMemcpy( &dSize, gcudaDstVec, sizeof(double), cudaMemcpyDeviceToHost ) );

    if(fmSrcVec.GetLength() <= round(dSize) )
    {
        double fLen = fmSrcVec.GetLength()*1.0;
        double fLen2 = 1.0;
        double *pData = fmSrcVec.GetData();

        PERR( cudaMemcpy( gcudaDstVec+3, pData, fmSrcVec.GetLength()*sizeof(double), cudaMemcpyHostToDevice ) );
        PERR( cudaMemcpy( gcudaDstVec+2, &fLen2, sizeof(double), cudaMemcpyHostToDevice ) );
        PERR( cudaMemcpy( gcudaDstVec+1, &fLen, sizeof(double), cudaMemcpyHostToDevice ) );

//_cuda_show( gcudaDstVec );

        return (0);
    }
    else
    {
    	Rprintf("cudaSize(%d) < vector size(%d)\n", round(dSize), fmSrcVec.GetLength());
	    return (-1);
    }
}

//*run on CPU;
int _copy_fmMatrix_Device(double* gcudaDstMat, CFmMatrix& fmSrcMat, bool bShow)
{
//printf("_copy_fmMatrix_Device %p-->%p\n", &fmSrcMat, gcudaDstMat);
    double dSize = 0.0;
    PERR(  cudaMemcpy( &dSize, gcudaDstMat, sizeof(double), cudaMemcpyDeviceToHost ) );

    if( fmSrcMat.GetNumRows() * fmSrcMat.GetNumCols() <= (int)round(dSize) )
    {
        double dRow = fmSrcMat.GetNumRows();
        double dCol = fmSrcMat.GetNumCols();
        double *pData = fmSrcMat.GetData();

        PERR( cudaMemcpy( gcudaDstMat+1, &dRow, sizeof(double), cudaMemcpyHostToDevice ) );
        PERR( cudaMemcpy( gcudaDstMat+2, &dCol, sizeof(double), cudaMemcpyHostToDevice ) );
        PERR( cudaMemcpy( gcudaDstMat+3, pData, round(dRow*dCol)* sizeof(double), cudaMemcpyHostToDevice ) );

//if (bShow) _cuda_show( gcudaDstMat );

        return (0);
    }
    else
    {
    	Rprintf("cudaSize(%d) < matrix size(%d*%d)\n", round(dSize), fmSrcMat.GetNumRows(), fmSrcMat.GetNumCols());
	    return (-1);
    }
}

//*run on CPU;
int _copyback_fmVector_Host(CFmVector& fmDstVec, double* gcudaSrcVec )
{
//printf("_copyback_fmVector_Host %p-->%p\n", gcudaSrcVec, &fmDstVec);

    double dSize=0;
    PERR( cudaMemcpy( &dSize, gcudaSrcVec+1, sizeof(double), cudaMemcpyDeviceToHost ) );

//printf("size=%.1f\n", dSize);

    fmDstVec.Resize((int)dSize, TRUE);
    PERR( cudaMemcpy( fmDstVec.GetData(), gcudaSrcVec + 3, dSize * sizeof(double), cudaMemcpyDeviceToHost ) );
    return (0);
}

//*run on CPU;
int _copyback_fmMatrix_Host(CFmMatrix& fmDstMat, double* gcudaSrcMat )
{
//printf("_copyback_fmMatrix_Host %p-->%p\n", gcudaSrcMat, &fmDstMat);

    double dRow, dCol;
    PERR( cudaMemcpy( &dRow, gcudaSrcMat+1, sizeof(double), cudaMemcpyDeviceToHost ) );
    PERR( cudaMemcpy( &dCol, gcudaSrcMat+2, sizeof(double), cudaMemcpyDeviceToHost ) );

    fmDstMat.Resize( round(dRow), round(dCol), TRUE);
    PERR( cudaMemcpy( fmDstMat.GetData(), gcudaSrcMat + 3, round(dRow * dCol * sizeof(double)), cudaMemcpyDeviceToHost ) );
    return (0);
}

double** make_vector_list( double** pList, unsigned int N, unsigned int nLen, bool verbose  )
{
    printf("make_vector_list, N=%d, %p\n", N, (void*)pList);

    for(unsigned int i=0;i<N;i++)
    {
        double* p;
        PERR( cudaMalloc( &p, (nLen + 3 )*sizeof(double) ) );
        pList[i] = p;

        double nMatHead[3];
        nMatHead[0] = ( nLen * 1.0 );
        nMatHead[1] = ( nLen * 1.0 );
        nMatHead[2] = ( 1.0 );
        PERR( cudaMemcpy( p, &nMatHead, 3*sizeof(double), cudaMemcpyHostToDevice ) );
if(verbose) printf("I=%d, %p\n", i, p);
    }

    return(pList);
}

double** make_matrix_list( double** pList, unsigned int N, unsigned int nRow, unsigned int nCol, bool verbose )
{
    printf("make_matrix_list, N=%d, %p\n", N, (void*)pList);

    for(unsigned int i=0;i<N;i++)
    {
        double* p;
        PERR( cudaMalloc( &p, (nRow*nCol+3)*sizeof(double) ) );
        pList[i] = p;

        double nMatHead[3];
        nMatHead[0] = ( nRow * nCol * 1.0 );
        nMatHead[1] = ( nRow * 1.0 );
        nMatHead[2] = ( nCol * 1.0 );
        PERR( cudaMemcpy( p, &nMatHead, 3*sizeof(double), cudaMemcpyHostToDevice ) );
if(verbose) printf("I=%d, %p\n", i, p);
    }

    return(pList);
}

double* make_matrix_ongpu( unsigned int nRow, unsigned int nCol )
{
    double* p;
    PERR( cudaMalloc( &p, (nRow*nCol+3)*sizeof(double) ) );

    double nMatHead[3];
    nMatHead[0] = ( nRow * nCol * 1.0 );
    nMatHead[1] = ( nRow * 1.0 );
    nMatHead[2] = ( nCol * 1.0 );
    PERR( cudaMemcpy( p, &nMatHead, 3*sizeof(double), cudaMemcpyHostToDevice ) );

    return(p);
}

#endif

