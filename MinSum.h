#pragma once
#include <algorithm>
#include <cstdio>
#include <vector>
#include <utility>
#include "graph.h"
#include "mex.h"
#include "matrix.h"

struct MinSum;

// Make really dumb to avoid memory trouble
struct Node {
  // Used a signed type to detection < subtraction fail.
  int nStates;
  double *vals;

  Node() : nStates(0), vals(nullptr) { }
  Node(int64_t nStates_, double *vals_) : nStates(nStates_), vals(vals_) { }

  Node(const Node &other) : nStates(other.nStates), vals(other.vals) { } 
  Node &operator=(const Node &other) {
    nStates = other.nStates;
    vals = other.vals;

    return *this;
  }

  double &operator()(int i) {
    mxAssert(i >= 0 && i < nStates, "Node operator() out of bounds");
    mxAssert(vals != nullptr, "Grrgh");
    return vals[i];
  }
};

//static_assert(sizeof(Node) % sizeof(double) == 0, "Node cannot be stored in double array");

struct Potential {
  int nLo;
  int nHi;
  double *vals;

  Potential() : nLo(0), nHi(0), vals(nullptr) { }
  Potential(int nLo_, int nHi_, double *vals_) : nLo(nLo_), nHi(nHi_), vals(vals_) { }

  Potential(const Potential &other) : nLo(other.nLo), nHi(other.nHi), vals(other.vals) { }
  Potential &operator=(const Potential &other) {
    nLo = other.nLo;
    nHi = other.nHi;
    vals = other.vals;

    return *this;
  }

  // Column major
  double &operator()(int lo, int hi) {
    mxAssert(lo >= 0 && lo < nLo && hi >= 0 && hi < nHi, "Potential operator() out of bounds");
    return vals[nLo*hi + lo];
  }
};

static_assert(sizeof(Potential) % sizeof(double) == 0, "Potential cannot be stored in double array");

struct Edge {
  // If E \in edges[i], then E.dst = i.
  int dst;
  bool srcHigher;
  double w;
  Potential *potential;

  Edge(int dst_, bool srcHigher_, double w_, Potential *potential_) :
    dst(dst_),
    srcHigher(srcHigher_),
    w(w_),
    potential(potential_) { }

  Edge(const Edge &other) : dst(other.dst), srcHigher(other.srcHigher), w(other.w), potential(other.potential) { }
  Edge &operator=(const Edge &other) {
    dst = other.dst;
    srcHigher = other.srcHigher;
    w = other.w;
    potential = other.potential;

    return *this;
  }

  // Evaluate the potential, transposing if necessary.
  double &p(size_t srcState, size_t dstState) const {
    // I'm hoping the branch predictor catches this
    if (srcHigher) {
      return (*potential)(dstState, srcState);
    } else {
      return (*potential)(srcState, dstState);
    }
  }
};

#ifndef NDEBUG
struct DebugEdge {
  int src;
  int dst;
  double fw;
  double rw;

  DebugEdge(int s, int d, double f, double r) : src(s), dst(d), fw(f), rw(r) { }
};
#endif

struct MinStats {
  int nBKNodes;
  int nBKEdgeBound;

  int nBKEdges;
  int nZInfEdges;
  int nSTEdges;
  int nPairEdges;

  double BKConstructTime;
  double BKMaxFlowTime;
  double computeEnergyTime;
};

struct MinSum {
  int nNodes;
  int nEdges;

  // Undirected edge
  // TODO: Replace with a cscMatrix (but that would require pointers)
  std::vector<std::vector<Edge>> neighbors;
  std::vector<Potential> potentials;
  std::vector<Node> nodes;
#ifndef NDEBUG
  std::vector<DebugEdge> debugEdges;
#endif

  typedef void (*TerrFunc)(const char *);
  TerrFunc errFunc;

  typedef int (*TprintFunc) (const char * format, ...);
  TprintFunc printFunc;

  // from BK
  typedef Graph<double, double, double> dGraph;

  // Fixed length allocation
  std::vector<double> mem;
  size_t memUsed;
  size_t MAX_MEM;
  double *alloc(size_t nDoubles);

  // Construct MinSum with an initial number of nodes, a initial memory size
  // to store nodes and potentials of initMem DOUBLES (not bytes) and an
  // error message. Note that nNodes must be fixed ahead of time.
  //
  // The fixed nNodes restriction is not difficult to remove, but incurs a
  // small runtime cost. Perhaps remove in future.
  MinSum(size_t nNodes_,
         TerrFunc errFunc_ = (void (*)(const char *)) puts,
         TprintFunc printFunc_ = printf,
         size_t initMem=10000,
         size_t maxMem=200000000) :
    nNodes(nNodes_),
    nodes(nNodes_),
    neighbors(nNodes_),
    errFunc(errFunc_),
    printFunc(printFunc_),
    mem(initMem),
    MAX_MEM(maxMem),
    memUsed(0)
  {
    potentials.reserve(0.5 * nNodes * nNodes);
  }

  Node &addNode(size_t n, int nStates_);
  void addEdge(size_t src, size_t dst, double w, Potential *ep);
  Potential &addPotential(int nLo, int nHi);

  enum class FailCode {
    SUCCESS = 0,
    NOT_SUBMODULAR,
    NOT_SYMMETRIC_W,
    NOT_TRANSPOSED_POTENTIAL,
    SELF_LOOP,
    INSUFFICIENT_STATES,
    INSUFFICIENT_NODES,
  };

  bool isSubModular(Potential &p);
  FailCode validate();
  // energy is 3-vector; [0] is total, [1] is unary, [2] is pairwise.
  MinStats minimize(std::vector<int> &x, double *energy, double &maxFlow);
};

