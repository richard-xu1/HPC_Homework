#include <stdio.h>
#include <array>
#include <iostream>
#include <cmath>
#include "utils.h"
#include <omp.h>
#define NUM_THREADS 4

const long N = 5000;
const double h = (double) 1/(N+1);


double norm(double f[], double Au[]){
	double sum = 0;
	for (long i = 0; i < (N+2)*(N+2); i += 1){
		sum = sum + (Au[i] - f[i])*(Au[i] - f[i]);
	}
	return sqrt(sum);
}

double compute_error(double u[], double f[], double Au[]){
	for (long i = 1; i < N+1; i++){
		for (long j = 1; j < N+1; j++){
			Au[i + j*(N+2)] = (-u[i-1 + j*(N+2)] - u[i + (j-1)*(N+2)] + 4*u[i + j*(N+2)] - u[i+1 + j*(N+2)] - u[i + (j+1)*(N+2)] )/(h*h);
			//std::cout<< "Auij = " << Au[i + j*(N+2)] << "\n";
		}
	} 
	double error = norm(Au, f);
	return error;
}
double gauss_seidel(double u[], double f[], double temp[], double Au[]){
	//double init_error = compute_error(A,u,f);
	//double err_tol = 1e-6;
	//double error = init_error;
	double iter = 0;
	std::cout << "Working on iter: " << iter << std::endl;
	omp_set_num_threads(NUM_THREADS);
	while( (iter < 3)){ //(error > init_error*err_tol) &&
		#pragma omp parallel for 
		for (long i = 1; i < N+1; i += 1){
			for(long j = 1; j < N+1; j += 1 ){
				if((i + j) % 2 == 0){
					temp[i + j*(N+2)] = (f[i + j*(N+2)]*h*h + u[i-1 + j*(N+2)] + u[i + (j-1)*(N+2)] + u[i+1 + j*(N+2)] + u[i + (j+1)*(N+2)])/4;
					//std::cout << "uk+1: " << temp[i+j*(N+2)] << "\n";
				}
			}	
		}
		//#pragma omp barrier
		
		#pragma omp parallel for
		for (long i = 1; i < N+1; i += 1){
			for(long j = 1; j < N+1; j += 1 ){
				if((i+j) % 2 == 1){
					temp[i + j*(N+2)] = (f[i + j*(N+2)]*h*h + temp[i-1 + j*(N+2)] + temp[i + (j-1)*(N+2)] + temp[i+1 + j*(N+2)] + temp[i + (j+1)*(N+2)])/4;
				}

			}	
		}
		//u = temp;
		for (long i = 1; i < N+1; i += 1){
			for(long j = 1; j < N+1; j += 1 ){
				u[i + j*(N+2)] = temp[i + j*(N+2)];
			}	
		}
		double error = compute_error(u,f, Au);
		std::cout << "Iteration: "<< iter << " Error: " << error <<  "\n";
		iter += 1;		
	}
	//std::cout << "Initial Error: "<< init_error << std::endl;
	//std::cout << "Iter: " << iter << std::endl;
	//std::cout << "Error: "<< error << std::endl;
	return *u;	
}

int main(int argc, char** argv) {
  //double * A = (double*) malloc( (N+2) * (N+2) * (N+2)* sizeof(double)); 
  double * u = (double*) malloc((N+2) * (N+2) * sizeof(double));
  double * f = (double*) malloc((N+2) * (N+2) * sizeof(double));
  double * temp = (double*) malloc((N+2) * (N+2) * sizeof(double));
  double * Au = (double*) malloc((N+2) * (N+2) * sizeof(double)); 


  for (long k = 0; k < N+2; k += 1) {
  	for (long l = 0; l < N+2; l += 1){
  	u[k + l*(N+2)] = 0;
	f[k + l*(N+2)] = 1;
	temp[k + l*(N+2)] = 0;
	Au[k + l*(N+2)] = 1;
		}
	}

  // for (long i = 0; i < N; i += 1) {
  // 	for (long j = 0; j < N; j += 1) {
  // 		if (i == j){
  // 			*(A + i*N + j) = 2/(double) N;
  // 		} else if ((i == j+1)||( i == j-1)) {
  // 			*(A + i*N + j)  = -1/(double) N;
  // 		} else {
  // 			*(A + i*N + j)  = 0;
  // 		}
  // 	}
  // }

  Timer t;
  t.tic();
  gauss_seidel(u,f,temp, Au);
  free(temp);
  double time = t.toc();
  std::cout << "Total time: " << time << std::endl;
  //free(A);
  free(u);
  free(f);
  free(Au);
  return 0;
}
