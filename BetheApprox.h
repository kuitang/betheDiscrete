#pragma once
#include <cmath>
#include "MinSum.h"


inline double sigmoid(double x) {
  return 1 / (1 + exp(-x));
}

// Returns xi and stores marginals in out.
inline double marginalize(double alpha, double q_i, double q_j, double *out) {
  double xi, beta, R;

  if (alpha == 0) {
    xi = q_i * q_j;
  } else {
    beta = 1 / alpha;
    R = beta + q_i + q_j;
    xi = 0.5*(R - copysign(sqrt(R*R - 4*(1 + beta)*q_i*q_j), beta));
  }

  out[0] = 1 + xi - q_i - q_j;
  out[1] = q_i - xi;
  out[2] = q_j - xi;
  out[3] = xi;

  return xi;
}

inline double binaryEntropy(double p) {
  if (p == 0) {
    return 0;
  } else {
    return -p * log(p) - (1 - p) * log(1 - p);
  }
}

template<int N>
inline double entropy(const double *p) {
  double ret = 0;
  // expect this to unroll
  for (int i = 0; i < N; i++) {
    if (p[i] > 0) {
      ret = ret - p[i]*log(p[i]);
    }
  }
  return ret;
}

inline bool oneNormConverged(size_t N, double *A, double *B, double thresh) {
  for (size_t i = 0; i < N; i++) {
    if (fabs(A[i] - B[i]) >= thresh) { return false; }
  }
  return true;
}

// Compressed sparse column sparse matrix (MATLAB's format)
struct cscMatrix {
  size_t N, M;
  size_t nzMax;
  double *pr;
  size_t *ir, *jc;
};

inline int degree(const cscMatrix &W, int j) {
  return W.jc[j+1] - W.jc[j];
}


// Compressed sparse column layout. (MATLAB's format)
//
// jc[j] is the linear index of the first nonzero element in column j and
// jc[j+1]-1 is the linear index of the last nonzero element in column j.
//
// ir[i] is the ROW of the node at linear index i.
inline cscMatrix extractCSC(const mxArray *M) {
  if (!mxIsDouble(M) || !mxIsSparse(M)) {
    mexErrMsgIdAndTxt("extractCSC:M", "matrix M must be double and sparse");
  }

  return { .N = mxGetN(M),
           .M = mxGetM(M),
           .nzMax = mxGetNzmax(M),
           .pr = mxGetPr(M),
           .ir = mxGetIr(M),
           .jc = mxGetJc(M) };
}

inline void deleteCscMatrix(cscMatrix &M) {
  delete[] M.pr;
  delete[] M.ir;
  delete[] M.jc;
}

bool propogateBetheBound(size_t nNodes,
                         const double *theta,
                         const cscMatrix &W,
                         double thresh,
                         int maxIter,
                         double *A, // outputs
                         double *B,
                         double *alpha);

bool mooijBound(size_t nNodes,
                const double *theta,
                const cscMatrix &W,
                double thresh,
                int maxIter,
                double *A, // outputs
                double *B,
                double *alpha);

double getIntervalSz(size_t nNodes,
                     const double *A,
                     const double *B,
                     const cscMatrix &W,
                     const double *alpha,
                     double epsilon);

void makeBetheMinSum(size_t nNodes,
                     const double *theta,
                     const cscMatrix &W,
                     const double *A,
                     const double *B,
                     const double *alpha,
                     double intervalSz,
                     MinSum &m);
