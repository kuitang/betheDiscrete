#include <cmath>
#include <sstream>
#include "mexcpp.h"
#include "BetheApprox.h"

using namespace mexcpp;

const char *usage = "Usage: [psi, [strPsi]] = makePotential_mex(w, qi, qj) where w is\n"
                    "weight of edge (i,j) and qi, qj are vectors of discretized\n"
                    "marginal values. psi is a length(qi) x length(qj) matrix,\n"
                    "where the rows index i and the columns index j.\n\n";

                    //"If strPsi is an output, prints a string format of the NEGATIVE potential (UAI form).";

void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  if (nIn != 3 || nOut > 2) {
    mexErrMsgIdAndTxt("makePotential_mex:args", usage);
  }

  if (nIn != 3 || nOut > 1) {
    mexErrMsgIdAndTxt("makePotential_mex:args", usage);
  }

  double w = scalar<double>(pIn[0]);
  Mat<double> qi(pIn[1]);
  Mat<double> qj(pIn[2]);

  Mat<double> psi(qi.length, qj.length);

  double alpha = exp(w) - 1;
  double marginals[4];
  double xi;

  for (size_t i = 0; i < qi.length; i++) {
    for (size_t j = 0; j < qj.length; j++) {
      xi = marginalize(alpha, qi[i], qj[j], marginals);
      psi(i,j) = -w * xi - entropy<4>(marginals);
    }
  }

  pOut[0] = psi;

  /*
  if (nOut == 2) {
    std::stringstream ss;
    for (size_t i = 0; i < qi.length; i++) {
      for (size_t j = 0; j < qi.length; j++) {
        ss << -j << " ";
      }
      ss << "\n";
    }
    pOut[1] = scalar(ss.str());
  }
  */
}
