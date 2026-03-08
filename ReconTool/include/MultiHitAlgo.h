#pragma once

#include <vector>

class TH1D;

namespace FROST {

struct MultiHitAlgoConfig {
  double alpha = 1.0;       // Weight exponent for the light yield centroid.
  double threshold1 = 10.0; // Threshold (p.e.) used in the first grouping stage.
  double threshold2 = 20.0; // Threshold (p.e.) used to derive the adaptive grouping threshold.
};

// Group bins whose content is above threshold1 into contiguous clusters.
//
// This is a direct port of the logic used in the original FROST_reconstruction.cc.
void DivideGroup1(const TH1D* h, std::vector<std::vector<int>>& group, const MultiHitAlgoConfig& cfg);

// Refine the grouping by re-evaluating each group with an adaptive threshold based on
// the local maximum. This may split a broad cluster into multiple clusters.
//
// This is a direct port of the logic used in the original FROST_reconstruction.cc.
void DivideGroup2(const TH1D* h, std::vector<std::vector<int>>& group, const MultiHitAlgoConfig& cfg);

// Compute candidate light-yield centroids (x_g) for multi-hit events.
//
// rectype encodes how the candidate should be mapped to a true position:
//   0 : "minus" candidate (use g1xm/g1ym in the two-hit map)
//   1 : "plus"  candidate (use g1xp/g1yp in the two-hit map)
//   2 : generic single-hit-style candidate (use gx/gy in the single-hit map)
//
// This is a direct port of the logic used in the original FROST_reconstruction.cc.
void MultiReconstructionX(const TH1D* h,
                          std::vector<double>& xgmulti,
                          std::vector<int>& rectype,
                          const std::vector<std::vector<int>>& groupx,
                          const MultiHitAlgoConfig& cfg);

void MultiReconstructionY(const TH1D* h,
                          std::vector<double>& ygmulti,
                          std::vector<int>& rectype,
                          const std::vector<std::vector<int>>& groupy,
                          const MultiHitAlgoConfig& cfg);

}  // namespace FROST
