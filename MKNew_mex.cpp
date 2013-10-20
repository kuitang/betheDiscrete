#include "mexcpp.h"
#include "tictoc.h"
#include <algorithm>
#include <cmath>

using namespace mexcpp;

enum {
  oA,
  oB
};

void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  clock_t mexBegin = tic();

  const char *usage = "Usage: [A, B] = MKNew_mex(theta, W, thresh). W must be sparse and symmetric.";
  if (nIn < 2 || nOut != 2) {
    mexErrMsgIdAndTxt("MKNew_mex:args", usage);
  }

  /////////////////////////////////////////////////
  // Extract arguments and outputs
  /////////////////////////////////////////////////
  Mat<double> theta                      (pIn[0]);
  cscMatrix   W = extractCSC             (pIn[1]);


  double thresh;
  if (nIn == 3) {
    thresh = scalar<double>(pIn[2]);
  } else {
    // Default
    thresh = 1e6;
  }

  int N = theta.length;
  int E = W.nzMax;

  Mat<double> A(N, 1);
  Mat<double> B(N, 1);
  pOut[oA] = A;
  pOut[oB] = B;

  double J[E];
  double etaLo[E];
  double etaHi[E];

  std::fill(etaLo, etaLo + E, -HUGE_VAL);
  std::fill(etaHi, etaHi + E, HUGE_VAL);

  double etaLoNew[E];
  double etaHiNew[E];
  double eta[N];

  mwIndex i;

  for (mwIndex idx = 0; idx < E; idx++) {
    J[idx] = 0.25 * W.pr[idx];
  }

  // double check
  /*
  for (mwIndex j = 0; j < N; j++) {
    for (mwIndex ijIdx = W.jc[j]; ijIdx < W.jc[j+1]; ijIdx++) {
      i = W.ir[ijIdx];
      mexPrintf("W[%d,%d] = %g ; J[%d,%d] = %g\n",
                i, j, W.pr[ijIdx], i, j, J[ijIdx]);
    }
  }
  */

  // This loop implements
  // eta = 0.5 * theta + sum(J, 2)
//  for (int j = 0; j < N; j++) {
//  }

  for (mwIndex j = 0; j < N; j++) {
    eta[j] = 0.5 * theta[j];
    for (mwIndex ijIdx = W.jc[j]; ijIdx < W.jc[j+1]; ijIdx++) {
      i = W.ir[ijIdx];
      eta[j] += J[ijIdx];
    }
    //mexPrintf("eta[%d] = %g\n", j, eta[j]);
  }

  int nIter = 0;
  double a, b, dLo, dHi;
  /*
  double b;
  double dLo;
  double dHi;
  */
  mwIndex k;

  while (true) {
    for (mwIndex j = 0; j < N; j++) {
      for (mwIndex ijIdx = W.jc[j]; ijIdx < W.jc[j+1]; ijIdx++) {
        i = W.ir[ijIdx];

        etaLoNew[ijIdx] = eta[i];
        etaHiNew[ijIdx] = eta[i];

        // find(J(:,i) ~= 0)
        for (mwIndex kiIdx = W.jc[i]; kiIdx < W.jc[i+1]; kiIdx++) {
          k = W.ir[kiIdx];

          if (k != j) {
            //mexPrintf("%s:%d (i,j) = (%d,%d), k = %d\n", __FILE__, __LINE__, i + 1, j + 1, k + 1);
            a = atanh(tanh(J[kiIdx]) * tanh(etaLo[kiIdx]));
            b = atanh(tanh(J[kiIdx]) * tanh(etaHi[kiIdx]));
            etaLoNew[ijIdx] += std::min(a, b);
            etaHiNew[ijIdx] += std::max(a, b);
          }
        }
      }
    }

    // Column-uniform convergence check
    for (mwIndex j = 0; j < N; j++) {
      dLo = 0;
      dHi = 0;

      for (mwIndex ijIdx = W.jc[j]; ijIdx < W.jc[j+1]; ijIdx++) {
        dLo += fabs(etaLo[ijIdx] - etaLoNew[ijIdx]);
        dHi += fabs(etaHi[ijIdx] - etaHiNew[ijIdx]);
      }

      //mexPrintf("col %d conv %g\n", j, dLo + dHi);

      if (dLo + dHi >= thresh) {
        // Column j already failed the uniform convergence test
        nIter++;

        std::copy(etaLoNew, etaLoNew + E, etaLo);
        std::copy(etaHiNew, etaHiNew + E, etaHi);

        goto nextIter;
      }
    }
    // If we got here, we converged
    break;

nextIter:
    ; // no-op
  }

  // intermediate check
  //mexPrintf("etaLo, etaHi\n");
  for (mwIndex idx = 0; idx < E; idx++) {
    //mexPrintf("%g %g\n", etaLo[idx], etaHi[idx]);
  }

  double betaLo[N];
  double betaHi[N];

  // Beliefs in the tanh parameterization.
  std::copy(eta, eta + N, betaLo);
  std::copy(eta, eta + N, betaHi);

  for (mwIndex j = 0; j < N; j++) {
    for (mwIndex ijIdx = W.jc[j]; ijIdx < W.jc[j+1]; ijIdx++) {
      i = W.ir[ijIdx];

      a = atanh(tanh(J[ijIdx]) * tanh(etaLo[ijIdx]));
      b = atanh(tanh(J[ijIdx]) * tanh(etaHi[ijIdx]));

      //mexPrintf("ijIdx = %d; a = %g; b = %g\n", ijIdx, a, b);

      betaLo[j] += std::min(a, b);
      betaHi[j] += std::max(a, b);
    }
  }


  //mexPrintf("betaLo betaHi\n");
  for (mwIndex j = 0; j < N; j++) {
    //mexPrintf("%g %g\n", betaLo, betaHi);
  }

  for (mwIndex j = 0; j < N; j++) {
    A[j] = (tanh(betaLo[j]) + 1) / 2.0;
    B[j] = 1 - ((tanh(betaHi[j]) + 1) / 2.0);
  }

}
