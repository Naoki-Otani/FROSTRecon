#include "TreeReader.h"

#include <stdexcept>
#include <algorithm>

namespace FROST {

TreeReader::TreeReader(const std::string& path,
                     const std::string& tree_name,
                     bool read_x,
                     bool read_y)
    : file_(TFile::Open(path.c_str(), "READ")),
      read_x_(read_x),
      read_y_(read_y) {
  if (!file_ || file_->IsZombie()) {
    throw std::runtime_error("Failed to open input ROOT file: " + path);
  }

  tree_ = dynamic_cast<TTree*>(file_->Get(tree_name.c_str()));
  if (!tree_) {
    throw std::runtime_error("Failed to find TTree '" + tree_name + "' in: " + path);
  }

  // Decide input format.
  //
  // MC:
  //   lightyieldx[132], lightyieldy[140]
  //
  // Data:
  //   lightyield[272][8]
  //
  if (HasBranch("lightyield")) {
    is_data_ = true;

    // Data format
    tree_->SetBranchAddress("lightyield", data_lightyield_raw_);
  } else {
    is_data_ = false;

    // MC format
    if (read_x_) {
      if (!HasBranch("lightyieldx")) {
        throw std::runtime_error("Missing branch 'lightyieldx' in MC input.");
      }
      tree_->SetBranchAddress("lightyieldx", mc_lightyieldx_.data());
    }

    if (read_y_) {
      if (!HasBranch("lightyieldy")) {
        throw std::runtime_error("Missing branch 'lightyieldy' in MC input.");
      }
      tree_->SetBranchAddress("lightyieldy", mc_lightyieldy_.data());
    }

    // Optional truth position branches used in calibration samples.
    if (HasBranch("x_muon1")) {
      truex_store_ = std::make_unique<double>(0.0);
      tree_->SetBranchAddress("x_muon1", truex_store_.get());
    }
    if (HasBranch("y_muon1")) {
      truey_store_ = std::make_unique<double>(0.0);
      tree_->SetBranchAddress("y_muon1", truey_store_.get());
    }
    if (HasBranch("x_muon2")) {
      truex2_store_ = std::make_unique<double>(0.0);
      tree_->SetBranchAddress("x_muon2", truex2_store_.get());
    }
    if (HasBranch("y_muon2")) {
      truey2_store_ = std::make_unique<double>(0.0);
      tree_->SetBranchAddress("y_muon2", truey2_store_.get());
    }
  }
}

bool TreeReader::HasBranch(const char* name) const {
  return tree_->GetBranch(name) != nullptr;
}

Long64_t TreeReader::Entries() const {
  return tree_->GetEntries();
}

void TreeReader::SetBunch(int bunch) {
  if (bunch < 0 || bunch >= 8) {
    throw std::runtime_error("Bunch index out of range: " + std::to_string(bunch));
  }
  bunch_ = bunch;
  if (is_data_) UpdateSelectedBunchView();
}

const double* TreeReader::PhotonXAtBunch(int bunch) const {
  if (!is_data_) return photonx_.data();
  if (bunch < 0 || bunch >= 8) {
    throw std::runtime_error("Bunch index out of range: " + std::to_string(bunch));
  }
  return &(photonx_bunch_[0][bunch]);
}

const double* TreeReader::PhotonYAtBunch(int bunch) const {
  if (!is_data_) return photony_.data();
  if (bunch < 0 || bunch >= 8) {
    throw std::runtime_error("Bunch index out of range: " + std::to_string(bunch));
  }
  return &(photony_bunch_[0][bunch]);
}

void TreeReader::UpdateSelectedBunchView() {
  if (!is_data_) return;

  if (read_x_) {
    for (int i = 0; i < kNfibX; ++i) {
      photonx_[i] = photonx_bunch_[i][bunch_];
    }
  }

  if (read_y_) {
    for (int i = 0; i < kNfibY; ++i) {
      photony_[i] = photony_bunch_[i][bunch_];
    }
  }
}

bool TreeReader::GetEntry(Long64_t i) {
  if (i < 0 || i >= Entries()) return false;
  tree_->GetEntry(i);

  if (is_data_) {
    // Convert raw Data lightyield[272][8] to MC-like views:
    //
    // raw [0..139]   corresponds to MC lightyieldy[139..0]
    // raw [140..271] corresponds to MC lightyieldx[131..0]
    //
    // i.e. both are reversed.
    //
    if (read_y_) {
      for (int raw = 0; raw < kNfibY; ++raw) {
        const int iy = kNfibY - 1 - raw;  // reverse order
        for (int b = 0; b < 8; ++b) {
          photony_bunch_[iy][b] = std::max(0.0, data_lightyield_raw_[raw][b]);
        }
      }
    }

    if (read_x_) {
      for (int raw = 0; raw < kNfibX; ++raw) {
        const int ix = kNfibX - 1 - raw;  // reverse order
        for (int b = 0; b < 8; ++b) {
          photonx_bunch_[ix][b] = std::max(0.0, data_lightyield_raw_[kNfibY + raw][b]);
        }
      }
    }

    UpdateSelectedBunchView();

    // No truth positions for Data.
    truex_.reset();
    truey_.reset();
    truex2_.reset();
    truey2_.reset();

  } else {
    // MC: just expose the arrays directly through photonx_/photony_.
    if (read_x_) {
      for (int ixf = 0; ixf < kNfibX; ++ixf) {
        photonx_[ixf] = std::max(0.0, mc_lightyieldx_[ixf]);
      }
    }

    if (read_y_) {
      for (int iyf = 0; iyf < kNfibY; ++iyf) {
        photony_[iyf] = std::max(0.0, mc_lightyieldy_[iyf]);
      }
    }

    if (truex_store_) truex_ = *truex_store_;
    else truex_.reset();

    if (truey_store_) truey_ = *truey_store_;
    else truey_.reset();

    if (truex2_store_) truex2_ = *truex2_store_;
    else truex2_.reset();

    if (truey2_store_) truey2_ = *truey2_store_;
    else truey2_.reset();
  }

  return true;
}

}  // namespace FROST
