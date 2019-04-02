/******************************************************************************
* FILE: omp_bug4.c
* DESCRIPTION:
*   This very simple program causes a segmentation fault.
* AUTHOR: Blaise Barney  01/09/04
* LAST REVISED: 04/06/05
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define N 1048

int main (int argc, char *argv[]) 
{
int nthreads, tid, i, j;
double *a;

/* Fork a team of threads with explicit variable scoping */
#pragma omp parallel shared(nthreads) private(i,j,tid,a) 
  {
  a = (double*) malloc(N * N * sizeof(double));
  // the matrix was too big to fit on the stack so we used malloc for 
  // dynamic memory allocation

  /* Obtain/print thread info */
  tid = omp_get_thread_num();
  if (tid == 0) 
    {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }
  printf("Thread %d starting...\n", tid);
  //#pragma omp for schedule(dynamic)
  /* Each thread works on its own private copy of the array */
  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      //printf("tid, i, j %d %d %d \n",tid,i,j );
      a[i+j*N] = tid + i + j;
      //printf("a[i+j*N] is %f\n", a[i+j*N] );
    }
  }

  /* For confirmation */
  printf("Thread %d done. Last element= %f\n",tid,a[N*N-1]);
  free(a);

  }  /* All threads join master thread and disband */
  
}

