#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <iostream>

// Scan A array and write result into prefix_sum array;
// use long data type to avoid overflow
void scan_seq(long* prefix_sum, const long* A, long n) {
  if (n == 0) return;
  prefix_sum[0] = A[0];
  for (long i = 1; i < n; i++) {
    prefix_sum[i] = prefix_sum[i-1] + A[i];
  }
}

long min(long A, long B){
  if (A < B){
    return A;
  }
  return B;
}

void scan_omp(long* prefix_sum, const long* A, long n) {
  // TODO: implement multi-threaded OpenMP scan
  const int num_threads = 4;
  const long n_chunk = n/num_threads;
  omp_set_num_threads(num_threads);
  #pragma omp parallel  
  { int thread_num = omp_get_thread_num();
    scan_seq(&prefix_sum[thread_num*n_chunk],&A[thread_num*n_chunk],n_chunk);
    if (thread_num == 0){
      scan_seq(&prefix_sum[num_threads*n_chunk],&A[num_threads*n_chunk],n-num_threads*n_chunk);
    }
  }
  for(long i = 1; i < num_threads; i++){
    for(long j = 0; j < n_chunk; j++){
      prefix_sum[j + i*n_chunk] = prefix_sum[j + i*n_chunk] + prefix_sum[i*n_chunk-1];
    }
  }
  for(long i = 0; i < n-num_threads*n_chunk; i++){
    prefix_sum[i+num_threads*n_chunk] = prefix_sum[i+num_threads*n_chunk] + prefix_sum[num_threads*n_chunk-1];
  }
}

int main() {
  long N = 100000000;
  long* A = (long*) malloc(N * sizeof(long));
  long* B0 = (long*) malloc(N * sizeof(long));
  long* B1 = (long*) malloc(N * sizeof(long));
  for (long i = 0; i < N; i++) {
    A[i] = rand() % 10;
    B0[i] = 0;
    B1[i] = 0;
  }

  double tt = omp_get_wtime();
  scan_seq(B0, A, N);
  // for(short i = 0; i < N; i++){
  //   std::cout<< B0[i] <<  " " << i << "\n";
  // }
  printf("sequential-scan = %fs\n", omp_get_wtime() - tt);

  tt = omp_get_wtime();
  scan_omp(B1, A, N);
  // for(short i = 0; i < N; i++){
  //   std::cout<< B1[i] << " " << i << "\n";
  // }
  printf("parallel-scan   = %fs\n", omp_get_wtime() - tt);

  long err = 0;
  for (long i = 0; i < N; i++) err = std::max(err, std::abs(B0[i] - B1[i]));
  printf("error = %ld\n", err);

  free(A);
  free(B0);
  free(B1);
  return 0;
}
