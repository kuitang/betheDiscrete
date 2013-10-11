#include "mexcpp.h"
#include "tictoc.h"
#include "graph.h"
#include "BetheApprox.h"

using namespace mexcpp;

enum {
  oEnergy,           /* 0 */
  oX,   /* 1 */
  oMisc            /* 2 */
};

#if !defined(MAX)
#define    MAX(A, B)    ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define    MIN(A, B)    ((A) < (B) ? (A) : (B))
#endif

void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  clock_t mexBegin = tic();

  const char *usage = "Usage: [energy, xmap, misc] = submodularMAP_mex(theta, W)";
  if (nIn != 2 || nOut != 3) {
    mexErrMsgIdAndTxt("submodularMAP_mex:args", usage);
  }

  /////////////////////////////////////////////////
  // Extract arguments and outputs
  /////////////////////////////////////////////////
  Mat<double> theta                      (pIn[0]);
  cscMatrix   W = extractCSC             (pIn[1]);

  int nNodes = theta.length;
  mwSize nEdges = W.nzMax / 2; // symmetric mat

  auto misc = StructMat(1, 1, {"maxFlow", "mexTotTime"});
  pOut[oMisc] = misc;

  std::vector<int> x(nNodes);

  // Construct flow graph. We follow notation of Greig 1989, which
  // requires symmetric potentials (00 and 11). This is given by
  // a reparameterization
  // [W/2 0 ; 0 W/2] = [W/2 0]' + [0 0 ; 0 W] + [0 -W/2]
  // so that
  // [0 0 ; 0 W] = [-W/2 0]' + [W/2 0 ; 0 W/2] + [0 W/2]

	typedef Graph<double,double,double> GraphType;
	GraphType *g = new GraphType(nNodes, nEdges + 2 * nNodes);
  g->add_node(nNodes);

  // Pairwise are normal (bidirectional) edges.
  mwIndex i, off;
  mwIndex mIdx = 0;
  double w;

  // Use Kolmogorov's submodular construction, Kolmogorov04 pp. 151 with
  // A = B = C = 0, D = -W, E_i(0) = 0, E_i(1) = -\theta_i

  // is there a cache-friendlier way to write this? does it even matter?

  double tAug;
  double s;
  double t;
  for (mwIndex j = 0; j < nNodes; j++) {
    tAug = 0.0;

    for (mwIndex wIdx = W.jc[j]; wIdx < W.jc[j+1]; wIdx++) {
      // Upper triangular only
      i = W.ir[wIdx];
      w = W.pr[wIdx];

      if (j > i) {
        mxAssert(w > 0, "W had a negative entry.");
        tAug += w;
        g->add_edge(i, j, w, 0);
        //mexPrintf("W[%d,%d] = %g\n", i, j, w);
      }
    }

    s = MAX(0, -theta[j]);
    t = MAX(0, theta[j]) + tAug;
    g->add_tweights(j, s, t);

    //mexPrintf("s->%d = %g ; %d->t = %g\n", j, s, j, t);
  }

  double maxFlow = g->maxflow();

  for (int i = 0; i < nNodes; i++) {
    x[i] = g->what_segment(i) == GraphType::SINK;
    //mexPrintf("x[%d] = %d, ws = %d\n", i, x[i], g->what_segment(i));
  }
  misc.set("maxFlow", scalar<double>(maxFlow));

  // Compute final energies
  double energy = 0.0;
  for (mwIndex j = 0; j < nNodes; j++) {
    energy -= theta[j] * x[j];

    for (mwIndex wIdx = W.jc[j]; wIdx < W.jc[j+1]; wIdx++) {
      // Upper triangular only
      i = W.ir[wIdx];
      if (j > i) {

        energy -= W.pr[wIdx] * x[i]*x[j];
      }
    }
  }
  pOut[oX]      = Mat<int>(x);
  pOut[oEnergy] = scalar<double>(energy);

  delete g;

  misc.set("mexTotTime", scalar<double>(toc(mexBegin)));
}

