#include "mexcpp.h"
#include "matrix.h"
#include "tictoc.h"
#include "BetheApprox.h"
#include "mex.h"

using namespace mexcpp;

void makeBetheGamsMinSum(size_t nNodes,
                     const double *theta,
                     const cscMatrix &W,
                     CellMat<Mat<double> > &gams,
                     MinSum &m,
                     double *alpha) {
  mxAssert(m.nNodes == nNodes, "nNodes was a lie.");
  mxAssert(m.nNodes == gams.length, "gams.length was a lie.");

  // Unaries: Term two of (Eq 4)
  for (size_t j = 0; j < nNodes; j++) {
    int nStates = gams[j].length;
    mxAssert(nStates >= 2, "states did not exceed two!");
    Node &nj = m.addNode(j, nStates);

    int degMinusOne = degree(W, j) - 1;
    double q;

    for (size_t k = 0; k < nStates; k++) {
      q = gams[j][k];
      nj(k) = -theta[j] * q + degMinusOne * binaryEntropy(q);
      //mexPrintf("%s:%d -- Node %d[%d] is %g\n", __FILE__, __LINE__, j, k, nj(k));
    }
  }

  for (size_t j = 0; j < nNodes; j++) {
    for (size_t idx = W.jc[j]; idx < W.jc[j+1]; idx++) {
      alpha[idx] = exp(fabs(W.pr[idx])) - 1;
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
          double qLo = gams[lo][iqLo];

          for (size_t iqHi = 0; iqHi < nHiStates; iqHi++) {
            double qHi = gams[hi][iqHi];
            double marginals[4];
            double xi = marginalize(aij, qLo, qHi, marginals);
            potential(iqLo, iqHi) = -w*xi - entropy<4>(marginals);

//            mexPrintf("%s:%d -- Potential at %lx entry (%d, %d) is %g\n",
//                      __FILE__, __LINE__, &potential, iqLo, iqHi, potential(iqLo, iqHi));
          }
        }

      }
    }
  }
}

enum {
  oLogZ,           /* 0 */
  oOneMarginals,   /* 1 */
  oTwoMarginals,   /* 2 */
  oMisc            /* 3 */
};


void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  clock_t mexBegin = tic();

  const char *usage = "Usage: [logZ, oneMarginals, twoMarginals, misc] = BetheGams_mex(theta, W, gams)";
  if (nIn != 3 || nOut != 4) {
    mexErrMsgIdAndTxt("BetheGams_mex:args", usage);
  }

  /////////////////////////////////////////////////
  // Extract arguments and outputs
  /////////////////////////////////////////////////
  Mat<double> theta                      (pIn[0]);
  cscMatrix   W = extractCSC             (pIn[1]);
  CellMat<Mat<double> > gams             (pIn[2]);

  int nNodes = theta.length;
  mwSize nEdges = W.nzMax / 2; // symmetric mat

  auto misc = StructMat(1, 1, {"Vm", "elMat", "e",
      "maxFlow", "x", "nBKNodes", "nBKEdgeBound", "nBKEdges", "nZInfEdges",
      "nSTEdges", "nPairEdges", "BKConstructTime", "BKMaxFlowTime", "computeEnergyTime",
      "makeMinSumTime", "mexTotTime"});

  pOut[oMisc] = misc;

  /////////////////////////////////////////////////
  // Prepare the MinSum problem
  /////////////////////////////////////////////////

  // TODO: Principled estimation of memory; not just 10MM.
  // 2 gigs
  //const size_t MAX_DOUBLES = 200000000;
  // 16 gigs. Diminishing marginals returns to using more.
  const size_t MAX_DOUBLES = 2147483648;
  MinSum m(nNodes, mexErrMsgTxt, mexPrintf, 10000000, MAX_DOUBLES);

  double alpha[W.nzMax];
  makeBetheGamsMinSum(nNodes, theta, W, gams, m, alpha);

#ifndef NDEBUG
  auto fc = m.validate();
  if (fc != MinSum::FailCode::SUCCESS) {
    mexErrMsgIdAndTxt("MinSum:validate", "Construction of MinSUm problem failed validation: fail code %d; defined in MinSum.h\n", fc);
  }
#endif

  std::vector<int> x;
  double energy[3];
  double maxFlow;

  /////////////////////////////////////////////////
  // Compute the MinSum solution
  /////////////////////////////////////////////////
  MinStats stats = m.minimize(x, energy, maxFlow);

  misc.set("e", Mat<double>(1, 3, energy));
  misc.set("maxFlow", scalar<double>(maxFlow));
  misc.set("x", Mat<int>(x, /*colMat = */ false));

  /////////////////////////////////////////////////
  // Output the time and problem size statistics
  /////////////////////////////////////////////////
#define setStat(field) misc.setS(#field, stats.field)
  setStat(nBKNodes);
  setStat(nBKEdgeBound);
  setStat(nBKEdges);
  setStat(nZInfEdges);
  setStat(nSTEdges);
  setStat(nPairEdges);
  setStat(BKConstructTime);
  setStat(BKMaxFlowTime);
  setStat(computeEnergyTime);

  /////////////////////////////////////////////////
  // Compute marginals from the MinSum solution
  /////////////////////////////////////////////////

  pOut[oLogZ] = scalar<double>(-energy[0]);

  Mat<double> oneMarginals(nNodes, 1);
  pOut[oOneMarginals] = oneMarginals;
  for (mwIndex n = 0; n < nNodes; n++) {
    oneMarginals[n] = gams[n][ x[n] ];
  }

  //mexPrintf("%s:%d -- oneMarginals success\n", __FILE__, __LINE__);

  // Recover marginals
  mwSize dims[] = {2, 2, nEdges};
  // TODO: Write wrapper into mexcpp
  pOut[oTwoMarginals] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  double *pTwoMarg = mxGetPr(pOut[oTwoMarginals]);

  mwIndex i, off;
  mwIndex mIdx = 0;
  // is there a cache-friendlier way to write this? does it even matter?
  for (mwIndex j = 0; j < nNodes; j++) {
    for (mwIndex wIdx = W.jc[j]; wIdx < W.jc[j+1]; wIdx++) {
      i = W.ir[wIdx];

      // Upper triangular only
      if (j > i) {
        off = 4 * mIdx; // we vary along the 3rd dimension.
        marginalize(alpha[wIdx], oneMarginals[i], oneMarginals[j], pTwoMarg + off);
        mIdx++;
      }
    }
  }

  /////////////////////////////////////////////////
  // Output more intermediate values if we're debugging
  /////////////////////////////////////////////////
#ifndef NDEBUG
  Mat<double> elMat(m.debugEdges.size(), 4);
  for (int i = 0; i < m.debugEdges.size(); i++) {
    const auto &de = m.debugEdges[i];
    elMat(i,0) = de.src;
    elMat(i,1) = de.dst;
    elMat(i,2) = de.fw;
    elMat(i,3) = de.rw;
  }
  misc.set("elMat", elMat);

  // Output the potentials
  CellMat<Mat<double>> Vm(m.potentials.size(), 1);
  for (int i = 0; i < m.potentials.size(); i++) {
    Potential &p = m.potentials[i];
    Mat<double> potMat(p.nLo, p.nHi);
    std::copy(&p(0,0), &p(0,0) + p.nLo * p.nHi, potMat.re);
    Vm.set(i, static_cast<mxArray *>(potMat));
  }
  misc.set("Vm", Vm);


#endif

  misc.set("mexTotTime", scalar<double>(toc(mexBegin)));
  //mexPrintf("%s:%d -- twoMarginals success\n", __FILE__, __LINE__);
}

