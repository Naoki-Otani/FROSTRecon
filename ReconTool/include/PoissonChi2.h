#pragma once

#include <cmath>

#include "Mapping.h"
#include "FROSTConstants.h"

namespace FROST {

// Poisson (Baker-Cousins) chi-square used in the original code.
//
// chi2 = 2 * sum_i [ mu_i - n_i + n_i*log(n_i) - n_i*log(mu_i) ]
// with the convention that for n_i ~ 0 the log term is skipped.
inline double PoissonChi2(const SingleHitMap& map,
                          const double* datax,  // size kNfibX
                          const double* datay,  // size kNfibY
                          double x_mm,
                          double y_mm) {
  double chi2 = 0.0;

  for (int i = 0; i < kNfibX; ++i) {
    const double n = datax[i];
    double mu = map.ExpectedX(i, x_mm);
    if (mu <= 0.0) mu = 1e-9;

    if (n > 0.1) {
      chi2 += 2.0 * (mu - n + n * std::log(n) - n * std::log(mu));
    } else {
      chi2 += 2.0 * mu;
    }
  }

  for (int i = 0; i < kNfibY; ++i) {
    const double n = datay[i];
    double mu = map.ExpectedY(i, y_mm);
    if (mu <= 0.0) mu = 1e-9;

    if (n > 0.1) {
      chi2 += 2.0 * (mu - n + n * std::log(n) - n * std::log(mu));
    } else {
      chi2 += 2.0 * mu;
    }
  }

  return chi2;
}

}  // namespace FROST
