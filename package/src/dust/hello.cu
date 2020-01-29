#include <curand_kernel.h>

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>


__global__ void tau2_calc(int P, double *a, double *tau2, double lambda2, double sigma2)
{
    double _tau2=0;
    /* insert code to calculate the index properly using blockIdx.x, blockDim.x, threadIdx.x */
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if( j<P)
    {
         curandState s ;
         // seed a random number generator
         curand_init ( (unsigned int)j , 0, 0, &s) ;
	            
         if (a[j] != 0.0 )
         {
              double InvTau2_1 = sqrt( lambda2 * sigma2)/fabs(a[j]);
             // double _tau2 = 1/func_invGau(InvTau2_1, lambda2);
	    
	          double y1 = curand_normal_double (&s);   
              double y = y1*y1;
              double x = InvTau2_1 + 0.5*InvTau2_1/lambda2 * ( InvTau2_1*y - sqrt(4*InvTau2_1*lambda2*y + InvTau2_1*InvTau2_1*y*y) );
	          if (  curand_uniform_double (&s) <= (InvTau2_1/(InvTau2_1+x)))
                  _tau2 = 1/x;
	          else
                  _tau2 = 1/((InvTau2_1*InvTau2_1)/x);
            
              tau2[j] = _tau2;
        }
        else
        {
            tau2[j] = 0.0;
        }        
    }  
}

#define THREADS_PER_BLOCK 64

void calculate_tau(double *a, double *tau2, double lambda2, double sigma2, int P)
{
    double *d_a, *d_tau2;

    cudaMalloc( (void **) &d_a, P*sizeof(double) );
    cudaMalloc( (void **) &d_tau2, P*sizeof(double) );
    cudaMemcpy( d_a, a, P*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy( d_tau2, tau2, P*sizeof(double), cudaMemcpyHostToDevice );

    tau2_calc<<< (P + (THREADS_PER_BLOCK-1)) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>( P, d_a, d_tau2, lambda2, sigma2);
    
/*{
#cudaError_t cudaerr = cudaDeviceSynchronize();
#if ( cudaerr!=CUDA_SUCCESS)
#Rprintf("kernel launch failed with err \"%s\".\n",cudaGetErrorString(cudaerr));
#else
#Rprintf("Kernel launch good\n");
#}
*/

    cudaMemcpy( tau2, d_tau2, P*sizeof(double), cudaMemcpyDeviceToHost );

    cudaFree( d_a );
    cudaFree( d_tau2 );
}
