#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "FROSTConstants.h"

namespace FROST {

// Helper for "is this a completed 1-mm grid graph?"
inline void RequireMinN(const TObject* obj, const std::string& name, int minN, const std::string& path) {
  const auto* g = dynamic_cast<const TGraph*>(obj);
  if (!g || g->GetN() < minN) throw std::runtime_error(name + " is not completed (1-mm grid). Regenerate map with mapfunction_tool. File: " + path);
}

class SingleHitMap {
 public:
  explicit SingleHitMap(const std::string& path) : file_(TFile::Open(path.c_str())) {
    if (!file_ || file_->IsZombie()) {
      throw std::runtime_error("Failed to open single-hit map file: " + path);
    }

    gx_ = dynamic_cast<TGraphErrors*>(file_->Get("gx"));
    gy_ = dynamic_cast<TGraphErrors*>(file_->Get("gy"));
    if (!gx_ || !gy_) {
      throw std::runtime_error("Missing gx/gy graphs in single-hit map file: " + path);
    }
    // Expect completed (1-mm grid) graphs for normal incidence.
    if (gx_->GetN() < kNfibX * 10 + 1 || gy_->GetN() < kNfibY * 10 + 1) {
      throw std::runtime_error("gx/gy are not completed (1-mm grid). Regenerate map with mapfunction_tool.");
    }

    gmcx_.resize(kNfibX, nullptr);
    for (int i = 0; i < kNfibX; ++i) {
      auto name = std::string("gmcx_") + std::to_string(i);
      gmcx_[i] = dynamic_cast<TGraph*>(file_->Get(name.c_str()));
      if (!gmcx_[i]) {
        throw std::runtime_error("Missing response graph: " + name + " in " + path);
      }
      if (gmcx_[i]->GetN() < kNmcPX) {
        throw std::runtime_error("gmcx graphs are not completed (1-mm grid). Regenerate map with mapfunction_tool.");
      }
    }

    gmcy_.resize(kNfibY, nullptr);
    for (int i = 0; i < kNfibY; ++i) {
      auto name = std::string("gmcy_") + std::to_string(i);
      gmcy_[i] = dynamic_cast<TGraph*>(file_->Get(name.c_str()));
      if (!gmcy_[i]) {
        throw std::runtime_error("Missing response graph: " + name + " in " + path);
      }
      if (gmcy_[i]->GetN() < kNmcPY) {
        throw std::runtime_error("gmcy graphs are not completed (1-mm grid). Regenerate map with mapfunction_tool.");
      }
    }
  }

  double XFromXg(double xg) const { return gx_->Eval(xg); }
  double YFromYg(double yg) const { return gy_->Eval(yg); }

  double ExpectedX(int fiber, double x_mm) const {
    if (fiber < 0 || fiber >= kNfibX) return 0.0;
    return gmcx_[fiber]->Eval(x_mm);
  }
  double ExpectedY(int fiber, double y_mm) const {
    if (fiber < 0 || fiber >= kNfibY) return 0.0;
    return gmcy_[fiber]->Eval(y_mm);
  }

 private:
  std::unique_ptr<TFile> file_;
  TGraphErrors* gx_ = nullptr;
  TGraphErrors* gy_ = nullptr;
  std::vector<TGraph*> gmcx_;
  std::vector<TGraph*> gmcy_;
};

class TwoHitMap {
 public:
  explicit TwoHitMap(const std::string& path) : file_(TFile::Open(path.c_str())) {
    if (!file_ || file_->IsZombie()) {
      throw std::runtime_error("Failed to open two-hit map file: " + path);
    }

    // Expect completed (1-mm grid) graphs built after interpolation/symmetry/clamp.
    // X: nmcpx+2 points ([-660..+660] inclusive) => kNmcPX + 2
    // Y: nmcpy+2 points ([-700..+700] inclusive) => kNmcPY + 2
    g1xm_ = dynamic_cast<TGraphErrors*>(file_->Get("g1xm"));
    g1xp_ = dynamic_cast<TGraphErrors*>(file_->Get("g1xp"));
    g1ym_ = dynamic_cast<TGraphErrors*>(file_->Get("g1ym"));
    g1yp_ = dynamic_cast<TGraphErrors*>(file_->Get("g1yp"));
    if (!g1xm_ || !g1xp_ || !g1ym_ || !g1yp_) {
      throw std::runtime_error("Missing g1xm/g1xp/g1ym/g1yp in two-hit map file: " + path);
    }
    RequireMinN(g1xm_, "g1xm", kNmcPX + 2, path);
    RequireMinN(g1xp_, "g1xp", kNmcPX + 2, path);
    RequireMinN(g1ym_, "g1ym", kNmcPY + 2, path);
    RequireMinN(g1yp_, "g1yp", kNmcPY + 2, path);
  }

  double XFromXgMinus(double xg) const { return g1xm_->Eval(xg); }
  double XFromXgPlus(double xg) const { return g1xp_->Eval(xg); }
  double YFromYgMinus(double yg) const { return g1ym_->Eval(yg); }
  double YFromYgPlus(double yg) const { return g1yp_->Eval(yg); }

 private:
  std::unique_ptr<TFile> file_;
  TGraphErrors* g1xm_ = nullptr;
  TGraphErrors* g1xp_ = nullptr;
  TGraphErrors* g1ym_ = nullptr;
  TGraphErrors* g1yp_ = nullptr;
};

}  // namespace FROST
