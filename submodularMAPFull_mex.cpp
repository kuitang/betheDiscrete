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

  const char *usage = "[energy, xmap, misc] = submodularMAPFull_mex(theta, pW, pots).\n"
                      " Solve MAP problem with overcomplete parameterization.\n"
                      " theta is nNodes x 2, first column is 0 potential, second is 1.\n"
                      " W is nNodes x nNodes sparse *indexing* a potential in pots, 1-based.\n"
                      " pots is nPots cell array, each 2x2.\n";

  if (nIn != 3) {
    mexErrMsgIdAndTxt("submodularMAP_mex:args", usage);
  }

  /////////////////////////////////////////////////
  // Extract arguments and outputs
  /////////////////////////////////////////////////
  Mat<double> theta         (pIn[0]);
  cscMatrix   W = extractCSC(pIn[1]);
  CellMat<Mat<double> > pots(pIn[2]);

  int nNodes = theta.M;
  mwSize nEdges = W.nzMax / 2; // symmetric mat
  int nPots  = pots.length;

  auto misc = StructMat(1, 1, {"maxFlow", "mexTotTime"});
  pOut[oMisc] = misc;

  std::vector<int> x(nNodes);

	typedef Graph<double,double,double> GraphType;
  GraphType *g = new GraphType(nNodes, nEdges + 2 * nNodes);
  g->add_node(nNodes);

  // Pairwise are normal (bidirectional) edges.
  mwIndex i, off;
  // Potential idx
  int pidx;
  double w;

  // Use Kolmogorov's submodular construction, Kolmogorov04 pp. 151
  // is there a cache-friendlier way to write this? does it even matter?
  std::vector<double> sAug(nNodes);
  std::vector<double> tAug(nNodes);

  for (mwIndex j = 0; j < nNodes; j++) {
    for (mwIndex wIdx = W.jc[j]; wIdx < W.jc[j+1]; wIdx++) {
      // Upper triangular only
      i = W.ir[wIdx];

      // some dirty trickery since we represent an int as a double.
      if (j > i) {
        if (pidx >= 0 && pidx < nPots) {
          pidx = W.pr[wIdx] - 1;
          // We may have weird numerical issues.
          if (pidx >= 0 && pidx < nPots) {
            sAug[i] += pots[pidx](1,0) - pots[pidx](0,0); // C - A
            tAug[j] += pots[pidx](1,0) - pots[pidx](1,1); // C - D
            w = pots[pidx](0,1) + pots[pidx](1,0) - pots[pidx](0,0) - pots[pidx](1,1); // B + C - A - D
            g->add_edge(i, j, w, 0);

            mexPrintf("%d->%d = %g\n", i, j, w);
          }
        }
      }
    }
  }

  double dTheta;
  double s;
  double t;
  for (int j = 0; j < nNodes; j++) {
    dTheta = theta(j,1) - theta(j,0);
    s      = MAX(0, dTheta)  + sAug[j];
    t      = MAX(0, -dTheta) + tAug[j];

    g->add_tweights(j, s, t);

    mexPrintf("s->%d = %g\t; %d->t = %g\t; sAug = %g\t; tAug = %g\n", j, s, j, t, sAug[j], tAug[j]);
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
    energy += theta(j,x[j]);

    for (mwIndex wIdx = W.jc[j]; wIdx < W.jc[j+1]; wIdx++) {
      // Upper triangular only
      i = W.ir[wIdx];
      if (j > i) {
        if (pidx >= 0 && pidx < nPots) {
          pidx    = W.pr[wIdx] - 1;
          energy += pots[pidx](x[i],x[j]);
        }
      }
    }
  }
  pOut[oX]      = Mat<int>(x);
  pOut[oEnergy] = scalar<double>(energy);

  delete g;

  misc.set("mexTotTime", scalar<double>(toc(mexBegin)));
}

