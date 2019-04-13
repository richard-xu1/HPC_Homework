#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <omp.h>
//#include "utils.h"

#define BLOCK_DIM 34


double norm(double *err, long N) {
  double sum = 0;
  #pragma omp parallel for reduction (+:sum)
  for (long i = 0; i < N*N ; i+=1) {
     sum += err[i]*err[i];
     //printf("sum = %f\n",sum);
  }
  return sqrt(sum);
} 

/*  
__global__ void jacobi_kernel(double* u, double* f, long N){
  __shared__ double smem[BLOCK_DIM][BLOCK_DIM];
 long grid_dim = (N-2)/(BLOCK_DIM - 2) ;
 int idx = (blockIdx.x)*(BLOCK_DIM-2) + threadIdx.x;
 int idy = (blockIdx.y)*(BLOCK_DIM-2) + threadIdx.y;
 smem[threadIdx.x+1][threadIdx.y+1] = u[idx][idy];
 
 if (blockIdx.x == 0) smem[0][threadIdx.y] = 0;
 if (blockIdx.y == 0) smem[threadIdx.x][0] = 0;
 if (blockIdx.x == grid_dim) smem[grid_dim][threadIdx.y] = 0;
 if (blockIdx.y == grid_dim) smem[threadIdx.x][grid_dim]= 0;

}
*/

__global__ void jacobi_kernel_nsmem(double* u, double* f, double * err, double * temp, long N, double h){
 int idx = (blockIdx.x)*(BLOCK_DIM-2) + threadIdx.x + 1;
 int idy = (blockIdx.y)*(BLOCK_DIM-2) + threadIdx.y + 1;
 temp[idx*N+idy] = (h*h*f[idx*N+idy] + u[(idx-1)*N+idy] + u[idx*N+idy-1] + u[(idx+1)*N+idy] + u[idx*N+idy+1])/4;
 //if (threadIdx.x == 0) {
   //printf("temp = %f", temp[idx*N+idy]); 
 //}
 __syncthreads();
 u[idx*N + idy] = temp[idx*N+idy];
 //printf("f, idx, idy = %f %d %d\n ",f[idx*N+idy], idx, idy);
 err[idx*N+idy] = (-u[(idx-1)*N+idy] - u[idx*N+idy-1] + 4*u[idx*N+idy] - u[(idx+1)*N+idy] - u[idx*N+idy+1] )/(h*h) - f[idx*N+idy];
}


int main(){
  long N = 6402; // dimension of 2D space
  double h = 1.0/(N+1); // size of update step
  double *u, *f, *err; // u is the solution, f is the forcing function
  cudaMallocHost((void**)&u, (N)*(N)*sizeof(double));
  cudaMallocHost((void**)&f, (N)*(N)*sizeof(double));
  cudaMallocHost((void**)&err, N*N*sizeof(double));  

  for(long i = 0; i < N; i++){
    for(long j =0; j < N; j++) {
       u[i*N+j] = 0;
       f[i*N+j] = 1;
       //temp[i*N +j] = 0;
       err[i*N+j] = 0;
   }
  }

  dim3 blockDim(BLOCK_DIM-2,BLOCK_DIM-2);
  dim3 gridDim( (N-2)/(BLOCK_DIM-2),(N-2)/(BLOCK_DIM-2));
  double *u_d, *f_d, *temp_d, *err_d;
  cudaMalloc(&u_d, N*N*sizeof(double));
  cudaMalloc(&f_d, N*N*sizeof(double));
  cudaMalloc(&temp_d, N*N*sizeof(double));
  cudaMalloc(&err_d, N*N*sizeof(double));

  cudaMemcpyAsync(u_d, u, N*N*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(f_d, f, N*N*sizeof(double), cudaMemcpyHostToDevice);
  cudaDeviceSynchronize();
  long t = 1;
  double tt = omp_get_wtime();
  while (t < 10000) {
   jacobi_kernel_nsmem<<<gridDim, blockDim, 0>>>(u_d,f_d,err_d, temp_d, N, h );
  if ((t % 1000) == 0) { 
    printf("time per 1000 iterations = %f s \n ", (omp_get_wtime()-tt) );
    //printf("Bandwidth = %f\n", 1000*7*(N-2)*(N-2)/(omp_get_wtime()-tt)/1e9);
    cudaMemcpyAsync(err,err_d,N*N*sizeof(double), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    printf("Err = %f\n", norm(err, N));
    tt = omp_get_wtime();
  }
  t+=1;
}
 cudaFree(u_d);
 cudaFree(f_d);
 cudaFree(temp_d);
 cudaFree(err_d);

 cudaFreeHost(u);
 cudaFreeHost(f);
 cudaFreeHost(err);

}

