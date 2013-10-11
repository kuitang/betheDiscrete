#include <cmath>
#include <cstdio>
#include <numeric>
#include "MinSum.h"
#include "tictoc.h"

// TODO: errMsg function

const double TEN_EPS = 2.2204e-15;
const double INF = 1e100;

// Checks that potentials are submodular and that the neighbors matrix
// is symmetric and potentials are transposed when needed.

bool MinSum::isSubModular(Potential &q) {
  char errMsg[1000];
  for (int lo = 0; lo < q.nLo - 1; lo++) {
    for (int hi = 0; hi < q.nHi - 1; hi++) {
      double diff = q(lo,hi) + q(lo+1, hi+1) - q(lo, hi+1) - q(lo+1, hi);
      if (diff - 5 * TEN_EPS > 0) {
        snprintf(errMsg, 1000, "isSubModular: GEQ Potential %p entry (%d,%d) failed (%g + %g > %g + %g) (diff = %g)\n", &q, lo, hi, q(lo,hi), q(lo+1,hi+1), q(lo,hi+1), q(lo+1,hi), diff);

        errFunc(errMsg);
        return false;
      }
    }
  }

  return true;
}

Node &MinSum::addNode(size_t n, int nStates_) {
  // Recall that all nodes were preallocated
  mxAssert(n < nNodes, "No space for this node!");

  Node &node = nodes[n];
  node.nStates = nStates_;
  node.vals = alloc(nStates_);
  return node;
}

// Potentials are always indexed lo on the rows, hi on the columns.
void MinSum::addEdge(size_t src, size_t dst, double w, Potential *ep) {
  size_t lo = std::min(src, dst);
  size_t hi = std::max(src, dst);

  mxAssert(lo >= 0 && lo < nNodes && hi >= 0 && hi < nNodes, "addEdge index out of bounds");

  neighbors[lo].emplace_back(hi, /*srcHigher*/ false, w, ep);
  neighbors[hi].emplace_back(lo, /*srcHigher*/ true,  w, ep);

  nEdges++;
}

// should be fast because this happens a lot
Potential &MinSum::addPotential(int nLo, int nHi) {
  mxAssert(nLo > 0 && nHi > 0, "Potential was allocated without entries!");

  double *vals = alloc(nLo * nHi);
  potentials.emplace_back(nLo, nHi, vals);
  return potentials.back();
}

double *MinSum::alloc(size_t nDoubles) {
  size_t nextMemUsed = memUsed + nDoubles;
  if (nextMemUsed > MAX_MEM) { errFunc("MinSum::alloc: Exceeded memory cap!"); }
  if (nextMemUsed >= mem.size()) {
    double *oldBase = mem.data();
    // single resize
//    mem.resize(MAX_MEM);
    mem.resize(2 * nextMemUsed);
    double *newBase = mem.data();

    printFunc("%s:%d -- Realloc; old base was %lp; new base is %lp; new size (in doubles) is %d\n",
               __FILE__, __LINE__, oldBase, newBase, mem.size());

    if (newBase != oldBase) { // translate all old pointers
      for (Node &node : nodes) {
        if (node.vals != nullptr) {
          ptrdiff_t off = node.vals - oldBase;
          node.vals = newBase + off;
        }
      }
      for (Potential &pot : potentials) {
        if (pot.vals != nullptr) {
          ptrdiff_t off = pot.vals - oldBase;
          pot.vals = newBase + off;
        }
      }
    }
  }

  double *ret = &mem[memUsed];
  memUsed = nextMemUsed;
  return ret;
}


MinSum::FailCode MinSum::validate() {
  for (Potential &p : potentials) {
    if (!isSubModular(p)) { return FailCode::NOT_SUBMODULAR;  }
  }

  // TODO: Check for symmetry

  for (int src = 0; src < nNodes; src++) {
    for (const Edge &e : neighbors[src]) {
      int dst = e.dst;
      if ( (src > dst && !e.srcHigher) ||
           (src < dst && e.srcHigher) ) { return FailCode::NOT_TRANSPOSED_POTENTIAL; }
      if (src == dst) { return FailCode::SELF_LOOP; }
    }
  }

  for (const Node &n : nodes) {
    if (n.nStates < 2) {
      return FailCode::INSUFFICIENT_STATES;
    }
  }

  if (nNodes != nodes.size()) {
    return FailCode::INSUFFICIENT_NODES;
  }

  return FailCode::SUCCESS;
}

MinStats MinSum::minimize(std::vector<int> &x, double *energy, double &maxFlow) {
  // Compute offsets to convert (r, k) indices to linear

  clock_t constructBegin = tic();

  /////////////////////////////////////////////////
  // Compute the offset array. This translates between
  // the (node, state) index into a linear (bkNode) index.
  /////////////////////////////////////////////////
  std::vector<int> offsets;
  offsets.push_back(0);
  for (int r = 1; r < nNodes; r++) {
    offsets.push_back(offsets[r - 1] + nodes[r - 1].nStates - 1);
  }
  mxAssert(offsets.size() == nNodes, "Offsets and nodes epic fail");

  MinStats stats;
  stats.nBKNodes = offsets.back() + nodes.back().nStates - 1;

//  for (int r = 0; r < nNodes; r++) {
//    mexPrintf("%s:%d -- Node %d@%p had %d represented states. Offset is %d\n",
//              __FILE__, __LINE__, r, nodes[r], nodes[r]->nStates - 1, offsets[r]);
//  }

  /////////////////////////////////////////////////
  // Upper-bound the edge count so BK graph can preallocate.
  /////////////////////////////////////////////////
  stats.nBKEdgeBound = offsets.back() + nodes.back().nStates;

  // Count the pairwise edges
  for (int r = 0; r < nNodes; r++) {
    for (const Edge &e : neighbors[r]) {
      int rr = e.dst;
      int rStates  = nodes[r].nStates - 1;
      int rrStates = nodes[rr].nStates - 1;
      stats.nBKEdgeBound += rStates * rrStates;
    }
  }

  // Make our graph
  dGraph *gp = new dGraph(stats.nBKNodes, stats.nBKEdgeBound, errFunc);

  // Add nodes
  gp->add_node(stats.nBKNodes);

  /////////////////////////////////////////////////
  // Add unary source/sink edges.
  /////////////////////////////////////////////////
  double qrk;
  stats.nSTEdges = 0;
  for (int r = 0; r < nNodes; r++) {
    for (int k = 0; k < nodes[r].nStates - 1; k++) {
      qrk = 0;

      for (const Edge &e : neighbors[r]) {
        int end = nodes[e.dst].nStates - 1;
        qrk += e.w * (e.p(k,0) + e.p(k,end) - e.p(k+1,0) - e.p(k+1,end));
        //mexPrintf("%s:%d qrk for r = %d, k = %d, rr = %d is now %g\n",
        //          __FILE__, __LINE__, r, k, e.dst, qrk);
        //mexPrintf("%s:%d Components: %g %g %g %g\n", __FILE__, __LINE__, e.p(k,1), e.p(k,end), e.p(k+1,1), e.p(k+1,end));
      }

      qrk = qrk / 2.0;
      qrk = qrk + nodes[r](k) - nodes[r](k+1);
      //mexPrintf("%s:%d qrk for r = %d, k = %d, FINAL is now %g\n", __FILE__, __LINE__, r, k, qrk);
      //mexPrintf("%s:%d Unary term was %g\n", __FILE__, __LINE__, nr(k) - nr(k+1));
      int node = offsets[r] + k;

      stats.nSTEdges++;
      if (qrk > 0) {
        gp->add_tweights(node, qrk, 0);
#ifndef NDEBUG
        debugEdges.emplace_back(-1, node + 1, qrk, 0);
#endif
        //mexPrintf("%s:%d -- S EDGE s -> %d W = %g\n",
        //          __FILE__, __LINE__, node, qrk);
      } else {
        gp->add_tweights(node, 0, -qrk);
#ifndef NDEBUG
        debugEdges.emplace_back(node + 1, -2, -qrk, 0);
#endif
        //mexPrintf("%s:%d -- T EDGE %d -> t W = %g\n",
        //          __FILE__, __LINE__, node, -qrk);
      }
    }
  }

  /////////////////////////////////////////////////
  // Add zero/infinity edges between nodes representing
  // successive states.
  /////////////////////////////////////////////////
  int loNode, hiNode;
  stats.nZInfEdges = 0;
  for (int r = 0; r < nodes.size(); r++) {
    Node &currNode = nodes[r];
    //mexPrintf("%s:%d -- r = %d, currNode->nStates = %d\n", __FILE__, __LINE__, r, currNode->nStates);

    for (int k = 0; k < currNode.nStates - 2; k++) {
      //mexPrintf("%s:%d -- r = %d, k = %d, currNode->nStates = %d\n",
      //          __FILE__, __LINE__, r, k, currNode->nStates);
      loNode = offsets[r] + k;
      hiNode = offsets[r] + k + 1;

      mxAssert(loNode >= 0 && hiNode >= 0, "No negative nodes!");
      mxAssert(loNode < stats.nBKNodes && hiNode < stats.nBKNodes, "Node idx out of bounds!");
      gp->add_edge(loNode, hiNode, 0.0, INF);
#ifndef NDEBUG
      debugEdges.emplace_back(loNode + 1, hiNode + 1, 0.0, INF);
#endif
      stats.nZInfEdges++;

      //mexPrintf("%s:%d -- EDGE %d -> %d FW = %g RW = %g\n",
      //          __FILE__, __LINE__, loNode, hiNode, 0.0, INF);
    }
  }


  /////////////////////////////////////////////////
  // Add pairwise edges between all states for nodes
  // that were neighbors in the original graph.
  /////////////////////////////////////////////////
  stats.nPairEdges = 0;
  int rr, rNode, rrNode;
  double inside, arr;
  for (int r = 0; r < nNodes; r++) {
    for (const Edge &e : neighbors[r]) {
      rr = e.dst;

      for (int k = 0; k < nodes[r].nStates - 1; k++) {
        rNode = offsets[r] + k;

        for (int kk = 0; kk < nodes[rr].nStates - 1; kk++) {
          rrNode = offsets[rr] + kk;
          //mexPrintf("%s:%d -- r = %d, rr = %d, k = %d, kk = %d, nodes[r]->nStates = %d, nodes[rr]->nStates = %d", __FILE__, __LINE__, r, rr, k, kk, nodes[r]->nStates, nodes[rr]->nStates);
          inside = e.p(k,kk) + e.p(k+1,kk+1) - e.p(k+1,kk) - e.p(k,kk+1);

          // We may have small perturbations
          if (fabs(inside) < 10 * TEN_EPS) {
            inside = 0.0;
          }

          arr = -(e.w * inside) / 2.0;
          // Will be true if problem was submodular.
          if (arr < 0) {
            mexPrintf("%s:%d -- arr = %g\n", __FILE__, __LINE__, arr);
          }
          mxAssert(arr >= 0, "not submodular!");

          gp->add_edge(rNode, rrNode, arr, 0.0);
#ifndef NDEBUG
          debugEdges.emplace_back(rNode + 1, rrNode + 1, arr, 0.0);
#endif
          stats.nPairEdges++;
          //mexPrintf("%s:%d -- EDGE %d -> %d FW = %g RW = %g\n",
          //          __FILE__, __LINE__, rNode, rrNode, arr, 0.0);
        }
      }
    }
  }
  stats.nBKEdges = stats.nZInfEdges + stats.nPairEdges;
  stats.BKConstructTime = toc(constructBegin);

  /////////////////////////////////////////////////
  // Run maxflow
  /////////////////////////////////////////////////
  clock_t maxFlowBegin = tic();
  maxFlow = gp->maxflow();
  stats.BKMaxFlowTime = toc(maxFlowBegin);

  /////////////////////////////////////////////////
  // Recover energy
  /////////////////////////////////////////////////
  clock_t computeEnergyBegin = tic();
  std::vector<bool> cut(stats.nBKNodes);
  for (int n = 0; n < stats.nBKNodes; n++) {
    cut[n] = gp->what_segment(n) == dGraph::SINK;
  }

  //mexPrintf("%s:%d -- Before assigning cut\n", __FILE__, __LINE__);
  x.resize(nNodes);

  // The state of node r is the lowest state asssigned to the sink.
  for (int r = 0; r < nNodes; r++) {
    int nst = nodes[r].nStates;
    int k = 0;
    while (k < nst - 1 && !cut[offsets[r] + k]) {
      k++;
    }
    x[r] = k;
  }

  energy[1] = 0.0;
  energy[2] = 0.0;
  for (int src = 0; src < nNodes; src++) {
    for (Edge &e : neighbors[src]) {
      int dst = e.dst;
      if (src < dst) { // count just once
        energy[2] += e.w * e.p(x[src], x[dst]);
      }
    }

    energy[1] += nodes[src](x[src]);
  }
  energy[0] = energy[1] + energy[2];

  delete gp;
  stats.computeEnergyTime = toc(computeEnergyBegin);

  return stats;
}

