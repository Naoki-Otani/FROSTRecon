#pragma once

#include <array>
#include <cstring>
#include <memory>
#include <optional>
#include <string>

#include "TFile.h"
#include "TTree.h"

#include "FROSTConstants.h"

namespace FROST {

// Reader for MC/Data light-yield TTrees.
//
// Supported formats:
//
// MC format:
//   - tree name typically: "wls"
//   - branches:
//       lightyieldx[132]/D
//       lightyieldy[140]/D
//   - optional truth branches:
//       x_muon1, y_muon1, x_muon2, y_muon2
//
// Data format:
//   - tree name typically: "tree"
//   - branch:
//       lightyield[272][8]/D
//
// Internal convention for Data:
//   raw lightyield[0..139][b]   -> photony_bunch_[139..0][b]
//   raw lightyield[140..271][b] -> photonx_bunch_[131..0][b]
//
// PhotonX()/PhotonY() return the currently selected bunch view for Data,
// and the plain MC arrays for MC.
class TreeReader {
 public:
  explicit TreeReader(const std::string& path,
                     const std::string& tree_name = "wls",
                     bool read_x = true,
                     bool read_y = true);

  Long64_t Entries() const;
  bool GetEntry(Long64_t i);

  // Current selected bunch view (for Data) or direct MC view.
  const double* PhotonX() const { return photonx_.data(); }
  const double* PhotonY() const { return photony_.data(); }

  // Access a specific bunch explicitly (Data only; for MC returns the same arrays).
  const double* PhotonXAtBunch(int bunch) const;
  const double* PhotonYAtBunch(int bunch) const;

  void SetBunch(int bunch);
  int Bunch() const { return bunch_; }

  bool IsData() const { return is_data_; }
  bool IsMC() const { return !is_data_; }

  // Optional truth positions used in calibration samples (MC only).
  std::optional<double> TrueX() const { return truex_; }
  std::optional<double> TrueY() const { return truey_; }
  std::optional<double> TrueX2() const { return truex2_; }
  std::optional<double> TrueY2() const { return truey2_; }

 private:
  bool HasBranch(const char* name) const;
  void UpdateSelectedBunchView();

  std::unique_ptr<TFile> file_;
  TTree* tree_ = nullptr;

  bool read_x_ = true;
  bool read_y_ = true;
  bool is_data_ = false;
  int bunch_ = 0;  // used only for Data

  // Public-facing arrays for the currently selected bunch (or direct MC arrays).
  std::array<double, kNfibX> photonx_{};
  std::array<double, kNfibY> photony_{};

  // Raw MC storage.
  std::array<double, kNfibX> mc_lightyieldx_{};
  std::array<double, kNfibY> mc_lightyieldy_{};

  // Raw Data storage: ROOT branch lightyield[272][8]/D
  double data_lightyield_raw_[kNfibX + kNfibY][8] = {};

  // Converted Data storage.
  std::array<std::array<double, 8>, kNfibX> photonx_bunch_{};
  std::array<std::array<double, 8>, kNfibY> photony_bunch_{};

  // Optional truth branches.
  std::optional<double> truex_;
  std::optional<double> truey_;
  std::optional<double> truex2_;
  std::optional<double> truey2_;

  std::unique_ptr<double> truex_store_;
  std::unique_ptr<double> truey_store_;
  std::unique_ptr<double> truex2_store_;
  std::unique_ptr<double> truey2_store_;
};

}  // namespace FROST
