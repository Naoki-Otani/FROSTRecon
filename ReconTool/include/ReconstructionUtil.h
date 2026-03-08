#pragma once

#include <cmath>

#include "FROSTConstants.h"

namespace FROST {

inline double ComputeXg(const double* photonx, double alpha) {
  double sum = 0.0;
  double wsum = 0.0;
  for (int i = 0; i < kNfibX; ++i) {
    const double ly = (photonx[i] > 0.0) ? photonx[i] : 0.0; // ensure non-negative
    const double w = std::pow(ly, alpha);
    sum += w * FiberXmm(i);
    wsum += w;
  }
  if (wsum <= 0.0) return 0.0;
  return sum / wsum;
}

inline double ComputeYg(const double* photony, double alpha) {
  double sum = 0.0;
  double wsum = 0.0;
  for (int i = 0; i < kNfibY; ++i) {
    const double ly = (photony[i] > 0.0) ? photony[i] : 0.0; // ensure non-negative
    const double w = std::pow(ly, alpha);
    sum += w * FiberYmm(i);
    wsum += w;
  }
  if (wsum <= 0.0) return 0.0;
  return sum / wsum;
}

}  // namespace FROST
