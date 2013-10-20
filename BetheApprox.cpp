#include "BetheApprox.h"
#include <algorithm>
#include <limits>
#include "mex.h"
#include "gsl_roots.h"

double TEN_EPS = 2.2204e-15;
void fixBounds(size_t nNodes, double *A, double *B) {
  // Fix small numerical issues where the bounds cross
  for (int n = 0; n < nNodes; n++) {
    double diff = 1 - B[n] - A[n];
    if (diff < -TEN_EPS) {
      mexErrMsgIdAndTxt("BetheApprox.cpp:fixBounds", "A, B bound crossed; invalid\n");
    } else if (diff < 0) {
      // Tie them to the average
      double m = 0.5 * (A[n] + 1 - B[n]);
      A[n] = m;
      B[n] = 1 - m;
    }
  }
}

// Return true if convergence threshold met
// Output are A, B, and alpha, which must be preallocated!
bool propogateBetheBound(size_t nNodes,
                         const double *theta,
                         const cscMatrix &W,
                         double thresh,
                         int maxIter,
                         double *A, // outputs
                         double *B,
                         double *alpha) {
  double posW[nNodes], negW[nNodes];
  for (size_t j = 0; j < nNodes; j++) {
    posW[j] = 0;
    negW[j] = 0;

    for (size_t idx = W.jc[j]; idx < W.jc[j+1]; idx++) {
      double w = W.pr[idx];
      if (w > 0) {
        posW[j] += w;
      } else {
        negW[j] -= w;
      }

      alpha[idx] = exp(fabs(w)) - 1;
    }

    A[j] = sigmoid(theta[j] - negW[j]);
    B[j] = 1 - sigmoid(theta[j] + posW[j]);
  }

  bool converged = false;
  double oldA[nNodes], oldB[nNodes];
  double L, U;

  int t; // scope it outside for debug
  for (t = 0; !converged && t < maxIter; t++) {
    std::copy(A, A + nNodes, oldA);
    std::copy(B, B + nNodes, oldB);

    for (size_t j = 0; j < nNodes; j++) {
      L = 1.0;
      U = 1.0;
      for (int idx = W.jc[j]; idx < W.jc[j+1]; idx++) {
        int i = W.ir[idx];
        double a = alpha[idx];

        if (W.pr[idx] > 0) {
          L = L * (1 + a*A[i] / (1 + a*(1 - B[j])*(1 - A[i])));
          U = U * (1 + a*B[i] / (1 + a*(1 - A[j])*(1 - B[i])));
        } else {
          L = L * (1 + a*B[i] / (1 + a*(1 - B[j])*(1 - B[i])));
          U = U * (1 + a*A[i] / (1 + a*(1 - A[j])*(1 - A[i])));
        }
      }

      A[j] = 1 / (1 + exp(-theta[j] + negW[j]) / L);
      B[j] = 1 / (1 + exp(theta[j]  + posW[j]) / U);
    }

    converged = oneNormConverged(nNodes, A, oldA, thresh) &&
                oneNormConverged(nNodes, B, oldB, thresh);
  }

#ifndef NDEBUG
  mexPrintf("DEBUG: Ran for %d iterations; maxIter = %d\n", t, maxIter);
#endif

  fixBounds(nNodes, A, B);
  return converged;
}

bool mooijBound(size_t nNodes,
                const double *theta,
                const cscMatrix &W,
                double thresh,
                int maxIter,
                double *A, // outputs
                double *B,
                double *alpha) {

  // copy-pasted alpha calculation (rest of code relies on it)
  double posW[nNodes], negW[nNodes];
  for (size_t j = 0; j < nNodes; j++) {
    posW[j] = 0;
    negW[j] = 0;

    for (size_t idx = W.jc[j]; idx < W.jc[j+1]; idx++) {
      double w = W.pr[idx];
      if (w > 0) {
        posW[j] += w;
      } else {
        negW[j] -= w;
      }

      alpha[idx] = exp(fabs(w)) - 1;
    }
  }

  ////////////////////////////////////////////////////////
  // Convert to {-1, +1} format.
  ////////////////////////////////////////////////////////
  // Jpr is the entries of the same sparsity pattern as W.
  double Jpr[W.nzMax];
  std::transform(W.pr, W.pr + W.nzMax, Jpr, [](double x) { return 0.25 * x;});

  double eta[nNodes];
  for (int n = 0; n < nNodes; n++) {
    eta[n] = 0.5 * theta[n];
    for (size_t idx = W.jc[n]; idx < W.jc[n+1]; idx++) {
      eta[n] += Jpr[idx];
      //mexPrintf("%s:%d eta[%d] = %g; Jpr[%d] = %g\n", __FILE__, __LINE__, n, eta[n], idx, Jpr[idx]);
    }
  }

  ////////////////////////////////////////////////////////
  // Propagate bounds on the cavity fields
  ////////////////////////////////////////////////////////
  double etaLo[W.nzMax];
  double etaHi[W.nzMax];
  std::fill_n(etaLo, W.nzMax, std::numeric_limits<double>::lowest());
  std::fill_n(etaHi, W.nzMax, std::numeric_limits<double>::max());

  double etaLoNew[W.nzMax];
  double etaHiNew[W.nzMax];

  double a, b;
  bool converged = false;
  int t;
  for (t = 0; t < maxIter; t++) {
    for (size_t j = 0; j < nNodes; j++) {
      for (size_t outIdx = W.jc[j]; outIdx < W.jc[j+1]; outIdx++) {
        size_t i = W.ir[outIdx];

        etaLoNew[outIdx] = eta[i];
        etaHiNew[outIdx] = eta[i];

        for (size_t inIdx = W.jc[i]; inIdx < W.jc[i+1]; inIdx++) {
          size_t k = W.ir[inIdx];
          if (k == j) { continue; }


          a = atanh(tanh(Jpr[inIdx]) * tanh(etaLo[inIdx]));
          b = atanh(tanh(Jpr[inIdx]) * tanh(etaHi[inIdx]));

          etaLoNew[outIdx] += std::min(a, b);
          etaHiNew[outIdx] += std::max(a, b);

          /*
          if (t < 3) {
            mexPrintf("%s:%d i = %d, j = %d, k = %d\n", __FILE__, __LINE__, i, j, k);
            mexPrintf("%s:%d a = %g, b = %g, etaLo[%d] = %g, etaHi[%d] = %g, etaLoNew[%d] = %g, etaHiNew[%d] = %g\n", __FILE__, __LINE__, a, b, inIdx, etaLo[inIdx], inIdx, etaHi[inIdx], outIdx, etaLoNew[outIdx], outIdx, etaHiNew[outIdx]);
          }
          */
        }
      }
    }

    // Check for convergence before the obligatory copy
    // (This might be a bit slow)
    if (oneNormConverged(W.nzMax, etaLoNew, etaLo, thresh) &&
        oneNormConverged(W.nzMax, etaHiNew, etaHi, thresh)) {
      converged = true;
      break;
    }

    std::copy(etaLoNew, etaLoNew + W.nzMax, etaLo);
    std::copy(etaHiNew, etaHiNew + W.nzMax, etaHi);

    /*
    if (t < 3) {
      for (int n = 0; n < W.nzMax; n++) {
        mexPrintf("%s:%d %g %g %g %g\n", __FILE__, __LINE__, etaLo[n], etaHi[n], etaLoNew[n], etaHiNew[n]);
      }
    }
    */
  }


  // debug
  /*
  for (size_t idx = 0; idx < W.nzMax; idx++) {
    mexPrintf("%s:%d etaLo[%d] = %g; etaHi[%d] = %g\n", __FILE__, __LINE__, idx, etaLo[idx], idx, etaHi[idx]);
  }
  */

  ////////////////////////////////////////////////////////
  // Calculate beliefs in tanh parameterization
  ////////////////////////////////////////////////////////
  double betaLo[nNodes];
  double betaHi[nNodes];
  std::copy(eta, eta + nNodes, betaLo);
  std::copy(eta, eta + nNodes, betaHi);

  for (size_t j = 0; j < nNodes; j++) {
    for (size_t idx = W.jc[j]; idx < W.jc[j+1]; idx++) {
      size_t i = W.ir[idx];
      //mexPrintf("%s:%d i = %d, j = %d\n", __FILE__, __LINE__, i, j);

      a = atanh(tanh(Jpr[idx]) * tanh(etaLo[idx]));
      b = atanh(tanh(Jpr[idx]) * tanh(etaHi[idx]));
      betaLo[j] += std::min(a, b);
      betaHi[j] += std::max(a, b);
    }
  }

  /*
  for (size_t j = 0; j < nNodes; j++) {
    mexPrintf("%s:%d betaLo[%d] = %g; betaHi[%d] = %g\n", __FILE__, __LINE__, j, betaLo[j], j, betaHi[j]);
  }
  */


  ////////////////////////////////////////////////////////
  // Calculate A,B bounds on pseudomargina
  ////////////////////////////////////////////////////////
  std::transform(betaLo, betaLo + nNodes, A, [](double x) { return 0.5 * (tanh(x) + 1); });
  std::transform(betaHi, betaHi + nNodes, B, [](double x) { return 1 - (0.5 * (tanh(x) + 1)); });

  fixBounds(nNodes, A, B);
  return converged;
}



// For node n, find points in the interval [A[n], 1 - B[n]] such that the
// distance between two consecutive points is at most intervalSz.
//
// Returns a cscMatrix. Remember to free[] the pointers yourself.
cscMatrix calcIntervals(size_t nNodes, const double *A, const double *B, double intervalSz) {
  // Precalculate my number of intervals
  size_t totPoints = 0;
  size_t maxPoints = 0;

  mxAssert(intervalSz > 0, "intervalSz must be positive!");

  // Number of endpoints is ceil(intervalLength / intervalSz) + 1 (for endpoint)
  //
  // i.e.
  // -----------
  // |  |  |  ||
  //
  // With floating point numbers, we will assume that the endpoint is never
  // reached by stepping with intervalSz.
  //
  // We might have intervalLength < intervalSz. In that case, we still need
  // the left endpoint, hence the max.
  for (size_t n = 0; n < nNodes; n++) {
    double intervalLength = 1 - B[n] - A[n];
    mxAssert(intervalLength >= 0, "Bounds failed.");

    double dPoints = intervalLength / intervalSz;
    if (dPoints > std::numeric_limits<int>::max()) {
      mexPrintf("%s:%d -- intervalLength = %g ; intervalSz = %g. But dPoints = intervalLength  / intervalSz = %g overflows int.\n",
                __FILE__, __LINE__, intervalLength, intervalSz, dPoints);
      mexErrMsgIdAndTxt("calcIntervals:tooManyPoints", "Interval size too small; discretization points too many. Aborting.");
    }

    size_t points = std::max((int) std::ceil(intervalLength / intervalSz), 1) + 1;
    mxAssert(points >= 2, "Did not create at least two states.");

    if (points > maxPoints) { maxPoints = points; }
    totPoints += points;
    //mexPrintf("%s:%d intervalLength = %g, intervalSz = %g, / = %g, points = %d, totPoints = %d\n",
    //          __FILE__, __LINE__, intervalLength, intervalSz, intervalLength / intervalSz, points, totPoints);
  }

  // maxPoints x nNodes sparse matrix
  cscMatrix m;
  m.M = maxPoints;
  m.N = nNodes;
  m.nzMax = totPoints;
  m.pr = new double[totPoints];
  m.ir = new size_t[totPoints];
  m.jc = new size_t[m.N+1];

  size_t ind = 0;
  for (size_t n = 0; n < nNodes; n++) {
    m.jc[n] = ind;
    size_t k = 0;

    for (double q = A[n]; q < 1 - B[n]; q += intervalSz) {
      if (ind >= totPoints) {
        mexPrintf("%s:%d totPoints fuckup message: q = %g, A[n] = %g, 1 - B[n] = %g, intervalSz = %g, ind = %d, totPoints = %d\n",
                  __FILE__, __LINE__, q, A[n], 1 - B[n], intervalSz, ind, totPoints);
      }
      mxAssert(ind < totPoints, "You fucked up calculating totPoints!");
      m.ir[ind] = k;
      m.pr[ind] = q;
      ind++;
      k++;
    }

    // Add the endpoint
    m.ir[ind] = k;
    m.pr[ind] = 1 - B[n];
    ind++;
  }
  mxAssert(ind == totPoints, "Point calculation failed.");
  // Last column; one-past-end
  // We don't need to add 1, because the ind++ above means that at the last
  // iterations, we were already one past the end.
  m.jc[nNodes] = ind;

  return m;
}

double getIntervalSz(size_t nNodes,
                     const double *A,
                     const double *B,
                     const cscMatrix &W,
                     const double *alpha,
                     double epsilon) {
  // (Eq 16)
  double eta[nNodes];
  for (size_t n = 0; n < nNodes; n++) {
    eta[n] = std::min(A[n], B[n]);
  }

  double aMax = -std::numeric_limits<double>::min();
  double aa, aij;
  for (size_t j = 0; j < nNodes; j++) {
    for (int idx = W.jc[j]; idx < W.jc[j+1]; idx++) {
      int i = W.ir[idx];

      // Upper triangular
      if (j > i) {
        aij  = alpha[idx];
        aa = aij*(aij + 1) / (4*(2*aij + 1)*eta[i]*eta[j]*(1 - eta[i])*(1 - eta[j]));
        if (aa > aMax) {
          aMax = aa;
        }
      }
    }
  }

  double bb, bInner, alphaSum;
  double bMax = -std::numeric_limits<double>::min();
  for (size_t j = 0; j < nNodes; j++) {
    alphaSum = 0.0;
    // loop over the neighbors
    for (size_t idx = W.jc[j]; idx < W.jc[j+1]; idx++) {
      double a = alpha[idx];
      alphaSum += (a + 1)*(a + 1) / (2*a + 1);
    }

    bInner = 1 - degree(W, j) + alphaSum;
    bb = 1 / (eta[j] * (1 - eta[j])) * bInner;
    if (bb > bMax) {
      bMax = bb;
    }
  }

  double Omega = std::max(aMax, bMax);
  // cast inside to force double division
  double density = (W.nzMax + nNodes) / double(W.N * W.M);

  return sqrt(2*epsilon / (nNodes * Omega * sqrt(density)));
}



cscMatrix adaptiveMinSumIntervals(size_t nNodes,
                                  const double *theta,
                                  const double *A,
                                  const double *B,
                                  double epsilon) {

  // Compute a single (half) BBP iteration
  double L[nNodes];
  double U[nNodes];
  for (size_t j = 0; j < nNodes; j++) {
    L[j] = 1.0;
    U[j] = 1.0;
    for (int idx = W.jc[j]; idx < W.jc[j+1]; idx++) {
      int i = W.ir[idx];
      double a = alpha[idx];

      if (W.pr[idx] > 0) {
        L[j] = L[j] * (1 + a*A[i] / (1 + a*(1 - B[j])*(1 - A[i])));
        U[j] = U[j] * (1 + a*B[i] / (1 + a*(1 - A[j])*(1 - B[i])));
      } else {
        L[j] = L[j] * (1 + a*B[i] / (1 + a*(1 - B[j])*(1 - B[i])));
        U[j] = U[j] * (1 + a*A[i] / (1 + a*(1 - A[j])*(1 - A[i])));
      }
    }
  }

  // Compute negative/positive weights.
  double posW[nNodes], negW[nNodes]
  for (size_t j = 0; j < nNodes; j++) {
    posW[j] = 0;
    negW[j] = 0;

    for (size_t idx = W.jc[j]; idx < W.jc[j+1]; idx++) {
      double w = W.pr[idx];
      if (w > 0) {
        posW[j] += w;
      } else {
        negW[j] -= w;
      }
    }
  }

  double D[nNodes];
  double lb[nNodes], ub[nNodes];
  double ll, ur;

  double sqrtSD[nNodes];
  double sumSqrtSD = 0;
  for (size_t n = 0; n < nNodes; n++) {
    lb[n] = -theta[n] - Wpos[n] + log(U[n])
    D[n]  = std::max(fabs(ll[n]), fabs(ur[n]));
    S[n]  = 1 - B[n] - A[n];
    if (S[n] < 0) { S[n] = 0; }

    sqrtSD[n] = sqrt(S[n] * D[n]);
    sumSqrtSD += sqrtSD[n];
  }

  double keps[nNodes];
  for (size_t n = 0; n < nNodes; n++) {
    keps[n] = epsilon / sumSqrtSD * sqrtSD[n];
  }

  // Now, actually do the adaptation!
  double Uconst, Lconst, prevr;
  // but, you _don't know_ totPoints in advance this time.
  //
  // refer to the MATLAB code for reference.

  cscMatrix m;
  //m.M = maxPoints;
  m.N = nNodes;
  //m.nzMax = totPoints;
  //m.pr = new double[totPoints];
  //m.ir = new size_t[totPoints];

  m.jc = new size_t[m.N+1];
  std::vector<double> pr;
  std::vector<size_t> ir;

  size_t currIdx = 0;

  double m, prevr, nextr;
  m.jc[0] = 0;
  for (size_t n = 0; n < nNodes; n++) {
    int nStates = 0;

    if (A[n] < 1 - B[n]) {
      while (prevr < 1 - B[n]) {
        adaptNextPoint(prevr, ub[n], lb[n], keps[n], A[n], B[n], &m, &nextr);
        prevr = nextr;

        nStates++;
        wVec.append(m);
        ir.append(nStates);
      }
      mxAssert(wVec.size() == ir.size(), "wVec and ir had different sizes.");
    } else {
      // just put a point on A[n]
      wVec.append(A[n]);
      nStates = 1;
      ir.append(0);
    }

    m.jc[n+1] = m.jc[n] + nStates; // cf. line after comment "Unaries"
  }

  m.nzMax = wVec.size();
  m.pr = new double[m.nzMax];
  m.ir = new size_t[m.nzMax];

  std::copy(wVec.data(), wVec.data() + m.nzMax, m.pr);
  std::copy(ir.data(), ir.data() + m.nzMax, m.ir);

  return m;
}

struct integralParams {

};

inline double xlogx(double x) {
  if (x == 0) { return 0; }
  return log(x);
}

void boundIntegral(double upperLimit, void *params, double *f, double *df) {
  *df = log(upperLimit) - log1p(-upperLimit);

  integralParams *p = (integralParams *) params;
  *f = p->Uconst * (upperLimit - p->lowerLimit) + xlogx(upperLimit) + xlogx(1 - upperLimit)
       -xlogx(p->lowerLimit) - xlogx(1 - p->lowerLimit);
}

void adaptNextPoint(double prevr,
                    const double *ub,
                    const double *lb,
                    const double *keps,
                    const double *A,
                    const double *B,
                    double *m,
                    double *nextr) {
  double p = prevr;

  //
}

void makeBetheMinSum(size_t nNodes,
                     const double *theta,
                     const cscMatrix &W,
                     const double *A,
                     const double *B,
                     const double *alpha,
                     double intervalSz,
                     MinSum &m) {
  mxAssert(m.nNodes == nNodes, "nNodes was a lie.");

  //cscMatrix intervals = calcIntervals(nNodes, A, B, intervalSz);
  cscMatrix intervals = adaptiveMinSumIntervals(nNodes, A, B, epsilon);

  // Unaries: Term two of (Eq 4)
  for (size_t j = 0; j < nNodes; j++) {
    int nStates = intervals.jc[j+1] - intervals.jc[j];
    mxAssert(nStates >= 2, "states did not exceed two!");
    Node &nj = m.addNode(j, nStates);

    int degMinusOne = degree(W, j) - 1;

    for (size_t k = 0; k < nStates; k++) {
      double q = intervals.pr[intervals.jc[j] + k];
      nj(k) = -theta[j] * q + degMinusOne * binaryEntropy(q);
      mexPrintf("%s:%d -- Node %d[%d] is %g\n", __FILE__, __LINE__, j, k, nj(k));
    }
  }

  // We know ahead of time there will be nnz distinct potentials.
  m.potentials.reserve(W.nzMax / 2);

  // Pairwise: (Eq 5)
  for (size_t hi = 0; hi < nNodes; hi++) {
    int nHiStates = m.nodes[hi].nStates;
    mxAssert(nHiStates >= 0, "nHiStates cannot be negative.");

    for (size_t nodeIdx = W.jc[hi]; nodeIdx < W.jc[hi+1]; nodeIdx++) {
      int lo = W.ir[nodeIdx];

      if (hi > lo) {
        int nLoStates = m.nodes[lo].nStates;
        mxAssert(nLoStates >= 0, "nLoStates cannot be negative.");

        // Only look at the upper triangular (minus diagonal)
        double w = W.pr[nodeIdx];
        double aij = alpha[nodeIdx];

        mxAssert(aij != 0, "alpha cannot be zero.");

        Potential &potential = m.addPotential(nLoStates, nHiStates);
        m.addEdge(lo, hi, 1.0, &potential);

        for (size_t iqLo = 0; iqLo < nLoStates; iqLo++) {
          double qLo = intervals.pr[intervals.jc[lo] + iqLo];

          for (size_t iqHi = 0; iqHi < nHiStates; iqHi++) {
            double qHi = intervals.pr[intervals.jc[hi] + iqHi];
            double marginals[4];
            double xi = marginalize(aij, qLo, qHi, marginals);
            potential(iqLo, iqHi) = -w*xi - entropy<4>(marginals);

            //mexPrintf("%s:%d -- Potential at %lx entry (%d, %d) is %g\n",
            //          __FILE__, __LINE__, &potential, iqLo, iqHi, potential(iqLo, iqHi));
          }
        }

      }
    }
  }

  deleteCscMatrix(intervals);
}

