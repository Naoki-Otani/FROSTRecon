#pragma once

#include <array>

namespace FROST {

// Geometry / binning constants used across the FROST analysis.
//
// Notes:
//  - kNfibX/kNfibY are the number of readout fibers (bins) in each direction.
//  - The coordinate convention follows the original analysis code:
//      x(mm) = -5*kNfibX + 5 + 10*i  (i = 0..kNfibX-1)
//      y(mm) = -5*kNfibY + 5 + 10*i  (i = 0..kNfibY-1)
//
// This repository focuses on normal incidence (tanx=tany=0), so angle-dependent
// constants are intentionally omitted.

constexpr int kNx = 93;
constexpr int kNy = 97;
constexpr int kNfibX = 132;
constexpr int kNfibY = 140;

// Number of MC points used in the original mapping tables.
// They correspond to a 1 mm grid spanning approximately [-659, +659] mm for X
// and [-699, +699] mm for Y (including endpoints).
constexpr int kNmcPX = 1319;
constexpr int kNmcPY = 1399;

// Labels (true positions in mm) used to generate mapping points.
extern const std::array<int, kNx> kXLabel;
extern const std::array<int, kNy> kYLabel;

// Helper: convert fiber index (0-based) to coordinate in mm.
inline double FiberXmm(int i) { return -5.0 * kNfibX + 5.0 + 10.0 * static_cast<double>(i); }
inline double FiberYmm(int i) { return -5.0 * kNfibY + 5.0 + 10.0 * static_cast<double>(i); }

}  // namespace FROST
