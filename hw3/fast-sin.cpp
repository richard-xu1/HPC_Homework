#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include "utils.h"
#include "intrin-wrapper.h"

// Headers for intrinsics
#ifdef __SSE__
#include <xmmintrin.h>
#endif
#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __AVX__
#include <immintrin.h>
#endif


// coefficients in the Taylor series expansion of sin(x)
static constexpr double c3  = -1/(((double)2)*3);
static constexpr double c5  =  1/(((double)2)*3*4*5);
static constexpr double c7  = -1/(((double)2)*3*4*5*6*7);
static constexpr double c9  =  1/(((double)2)*3*4*5*6*7*8*9);
static constexpr double c11 = -1/(((double)2)*3*4*5*6*7*8*9*10*11);
// sin(x) = x + c3*x^3 + c5*x^5 + c7*x^7 + c9*x^9 + c11*x^11

static constexpr double c2  = -1/(((double)2));
static constexpr double c4  =  1/(((double)2)*3*4);
static constexpr double c6  = -1/(((double)2)*3*4*5*6);
static constexpr double c8  =  1/(((double)2)*3*4*5*6*7*8);
static constexpr double c10 = -1/(((double)2)*3*4*5*6*7*8*9*10);
static constexpr double c12 = 1/(((double)2)*3*4*5*6*7*8*9*10*11*12);

void sin4_reference(double* sinx, const double* x) {
  for (long i = 0; i < 4; i++) sinx[i] = sin(x[i]);
}

double sin_taylor(double x) {
  double x1  = x;
  double x2  = x1 * x1;
  double x3  = x1 * x2;
  double x5  = x3 * x2;
  double x7  = x5 * x2;
  double x9  = x7 * x2;
  double x11 = x9 * x2;

  double s = x1;
  s += x3  * c3;
  s += x5  * c5;
  s += x7  * c7;
  s += x9  * c9;
  s += x11 * c11;
  return s;
}

double cos_taylor(double x){
  double x1  = x;
  double x2  = x1 * x1;
  double x4  = x2 * x2;
  double x6  = x4 * x2;
  double x8  = x6 * x2;
  double x10  = x8 * x2;
  double x12 = x10 * x2;

  double s = 1.0;
  s += x2  * c2;
  s += x4  * c4;
  s += x6  * c6;
  s += x8  * c8;
  s += x10 * c10;
  s += x12 * c12;
  return s; 
}

void sin4_mod_taylor(double* sinx, const double* x){
  for (long i = 0; i < 4; i++){
  	double xi = x[i];
  	xi = fmod(xi, M_PI*2);
  	double s;
    if (abs(xi) <= M_PI/4){
      s = sin_taylor(xi);
    }
    else if ((xi > M_PI/4)&& (xi <= M_PI/2) ){
      s = cos_taylor(M_PI/2 - xi);
    }
    else if ((xi < -M_PI/4) && (xi >= -M_PI/2) ){
      s = -cos_taylor(M_PI/2 + xi);
    }
    else if (xi > M_PI/2) {
      s = -sin_taylor(xi - M_PI);
    }
    else if (xi < - M_PI/2) {
      s = -sin_taylor(xi + M_PI);
    }
    sinx[i] = s;
    //std::cout << sinx[i] << "\n";
  }
}

void sin4_taylor(double* sinx, const double* x) {
  for (int i = 0; i < 4; i++) {
    double x1  = x[i];
    double x2  = x1 * x1;
    double x3  = x1 * x2;
    double x5  = x3 * x2;
    double x7  = x5 * x2;
    double x9  = x7 * x2;
    double x11 = x9 * x2;

    double s = x1;
    s += x3  * c3;
    s += x5  * c5;
    s += x7  * c7;
    s += x9  * c9;
    s += x11 * c11;
    sinx[i] = s;
  }
}

void sin4_vector(double* sinx, const double* x) {
  // The Vec class is defined in the file intrin-wrapper.h
  typedef Vec<double,4> Vec4;
  Vec4 x1, x2, x3, x5, x7, x9, x11;
  x1  = Vec4::LoadAligned(x);
  x2  = x1 * x1;
  x3  = x1 * x2;
  x5  = x3 * x2;
  x7  = x5 * x2;
  x9  = x7 * x2;
  x11 = x9 * x2;

  Vec4 s = x1;
  s += x3  * c3 ;
  s += x5  * c5;
  s += x7  * c7;
  s += x9  * c9;
  s += x11 * c11;
  s.StoreAligned(sinx);
}

double err(double* x, double* y, long N) {
  double error = 0;
  for (long i = 0; i < N; i++) error = std::max(error, fabs(x[i]-y[i]));
  return error;
}

int main() {
  Timer tt;
  long N = 1000000;
  double* x = (double*) aligned_malloc(N*sizeof(double));
  double* xm = (double*) aligned_malloc(N*sizeof(double));
  double* sinx_ref = (double*) aligned_malloc(N*sizeof(double));
  double* sinx_taylor = (double*) aligned_malloc(N*sizeof(double));
  double* sinx_mod = (double*) aligned_malloc(N*sizeof(double));
  double* sinx_vector = (double*) aligned_malloc(N*sizeof(double));
  for (long i = 0; i < N; i++) {
    x[i] =  (drand48() - 0.5) * M_PI/2; // [-pi/4,pi/4]
    xm[i] = (drand48() - 0.5) * M_PI*10;   // random 
    sinx_ref[i] = 0;
    sinx_taylor[i] = 0;
    sinx_mod[i] = 0;
    sinx_vector[i] = 0;
  }

  tt.tic();
  for (long rep = 0; rep < 1000; rep++) {
    for (long i = 0; i < N; i+=4) {
      sin4_reference(sinx_ref+i, x+i);
    }
  }
  printf("Reference time: %6.4f\n", tt.toc());

  tt.tic();
  for (long rep = 0; rep < 1000; rep++) {
    for (long i = 0; i < N; i+=4) {
      sin4_taylor(sinx_taylor+i, x+i);
    }
  }

  printf("Original Taylor time:    %6.4f      Error: %e\n", tt.toc(), err(sinx_ref, sinx_taylor, N));
  tt.tic();
  for (long rep = 0; rep < 1000; rep++) {
    for (long i = 0; i < N; i+=4) {
      sin4_mod_taylor(sinx_mod+i, xm+i);
    }
  }
  printf("Modified Taylor time:    %6.4f      Error: %e\n", tt.toc(), err(sinx_ref, sinx_taylor, N));

  tt.tic();
  for (long rep = 0; rep < 1000; rep++) {
    for (long i = 0; i < N; i+=4) {
      sin4_vector(sinx_vector+i, x+i);
    }
  }
  printf("Vector time:    %6.4f      Error: %e\n", tt.toc(), err(sinx_ref, sinx_vector, N));
  printf("The Modified Taylor time calculates the sin of \n a variable x in an arbitrary range by replacing sin(x) with -cos(pi/2 + x) \n if x < -pi/4 or cos(pi/2-x) if x > pi/4");
  aligned_free(x);
  aligned_free(sinx_ref);
  aligned_free(sinx_taylor);
  aligned_free(sinx_mod);
  aligned_free(sinx_vector);
}

