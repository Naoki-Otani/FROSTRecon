// DrawFrostExternalFitDiagnostics.C
//
// Usage in ROOT:
//   MC:
//     root -l -b -q 'DrawFrostExternalFitDiagnostics.C("input_dir", "output.pdf", "mc")'
//   data:
//     root -l -b -q 'DrawFrostExternalFitDiagnostics.C("input_dir", "output.pdf", "data")'
//
// The macro reads all .root files directly under input_dir, loops over the
// match_info tree, and saves a multi-page PDF for FROST external-track
// diagnostics.
//
// The median, mean, and sigma_68 values written on the plots and in the
// summary log are calculated from all selected residual values, including
// underflow and overflow values outside the displayed histogram range.
//
// In MC mode, the output contains:
//   - x_{ext}-x_{true} and y_{ext}-y_{true}
//   - x_{rec}-x_{ext} and y_{rec}-y_{ext}
//   - x_{PM-only}-x_{DWG-only} and y_{PM-only}-y_{DWG-only}
//   - x_{rec}-x_{true} and y_{rec}-y_{true}
//   - 2D correlation plots of x_{ext}-x_{true} versus x_{rec}-x_{true}
//     and y_{ext}-y_{true} versus y_{rec}-y_{true}
//
// In data mode, truth information is not available, so only
// x_{rec}-x_{ext}, y_{rec}-y_{ext}, and PM-only minus DWG-only split-fit
// residual pages are drawn.
//
// Selection:
//   common:
//      trackmatch_has_match == 1
//      trackmatch_ninja_track_type == 1
//      abs(frost_nearest_x) < ACCEPTANCE_X mm
//      abs(frost_nearest_y) < ACCEPTANCE_Y mm
//      trackmatch_external_num_planes_downstream_wagasci_x >= NHIT_DWG_X
//      trackmatch_external_num_planes_proton_module_x >= NHIT_PM_X
//      trackmatch_external_num_planes_downstream_wagasci_y >= NHIT_DWG_Y
//      trackmatch_external_num_planes_proton_module_y >= NHIT_PM_Y
//
//   data only:
//      bsd_good_spill_flag != 0
//      detector_flags[0] == 1  // Proton Module detector flag
//      detector_flags[2] == 1  // downstream WAGASCI detector flag
//
// Angle-binned x plots use theta_x = atan(|tan_x|).
// Angle-binned y plots use theta_y = atan(|tan_y|).

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TString.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TPad.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int NHIT_DWG_X = 4;
constexpr int NHIT_PM_X = 10;
constexpr int NHIT_DWG_Y = 4;
constexpr int NHIT_PM_Y = 10;

constexpr double ACCEPTANCE_X = 560.0; //mm
constexpr double ACCEPTANCE_Y = 600.0; //mm

constexpr double kNonInitializedThreshold = -9.0e7;

// Angle bins [deg]:
//   [0,5), [5,10), [10,15), [15,20), [20,25),
//   [25,30), [30,35), [35,40)
constexpr int kNAngleBins = 8;
const double kAngleBins[kNAngleBins + 1] = {
  0.0, 5.0, 10.0, 15.0, 20.0,
  25.0, 30.0, 35.0, 40.0
};

enum class SampleType {
  kMonteCarlo,
  kData
};

SampleType ParseSampleType(const char *sample_type_arg) {
  std::string sample_type = sample_type_arg ? sample_type_arg : "mc";
  std::transform(sample_type.begin(), sample_type.end(), sample_type.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

  if (sample_type == "mc" ||
      sample_type == "montecarlo" ||
      sample_type == "monte_carlo") {
    return SampleType::kMonteCarlo;
  }
  if (sample_type == "data" ||
      sample_type == "realdata" ||
      sample_type == "real_data") {
    return SampleType::kData;
  }

  throw std::invalid_argument(
    "Unknown sample type: " + sample_type + " (use mc or data)");
}

bool IsMonteCarlo(SampleType sample_type) {
  return sample_type == SampleType::kMonteCarlo;
}

bool IsData(SampleType sample_type) {
  return sample_type == SampleType::kData;
}

struct CorrelationStats {
  Long64_t n = 0;
  double sum_e = 0.;
  double sum_r = 0.;
  double sum_ee = 0.;
  double sum_rr = 0.;
  double sum_er = 0.;

  void Fill(double external_minus_true, double reco_minus_true) {
    ++n;
    sum_e += external_minus_true;
    sum_r += reco_minus_true;
    sum_ee += external_minus_true * external_minus_true;
    sum_rr += reco_minus_true * reco_minus_true;
    sum_er += external_minus_true * reco_minus_true;
  }

  bool Calculate(double &mean_e,
                 double &mean_r,
                 double &rms_e,
                 double &rms_r,
                 double &covariance,
                 double &correlation) const {
    if (n < 2) {
      return false;
    }

    mean_e = sum_e / static_cast<double>(n);
    mean_r = sum_r / static_cast<double>(n);

    const double variance_e =
      sum_ee / static_cast<double>(n) - mean_e * mean_e;
    const double variance_r =
      sum_rr / static_cast<double>(n) - mean_r * mean_r;

    if (variance_e <= 0. || variance_r <= 0.) {
      return false;
    }

    rms_e = std::sqrt(variance_e);
    rms_r = std::sqrt(variance_r);
    covariance = sum_er / static_cast<double>(n) - mean_e * mean_r;
    correlation = covariance / (rms_e * rms_r);

    return std::isfinite(correlation);
  }
};

struct ResidualSummary {
  bool valid = false;
  Long64_t n = 0;
  double mean = 0.;
  double median = 0.;
  double sigma68 = 0.;
  double q16 = 0.;
  double q84 = 0.;
};

using ResidualValueMap = std::map<const TH1D*, std::vector<double>>;

double CalculateQuantileFromSortedValues(const std::vector<double> &sorted_values,
                                         double probability) {
  if (sorted_values.empty()) {
    return 0.;
  }

  if (probability <= 0.) {
    return sorted_values.front();
  }
  if (probability >= 1.) {
    return sorted_values.back();
  }

  const double position = probability *
    static_cast<double>(sorted_values.size() - 1);
  const std::size_t lower_index = static_cast<std::size_t>(std::floor(position));
  const std::size_t upper_index = static_cast<std::size_t>(std::ceil(position));
  const double fraction = position - static_cast<double>(lower_index);

  if (lower_index == upper_index) {
    return sorted_values.at(lower_index);
  }

  return sorted_values.at(lower_index) * (1. - fraction) +
         sorted_values.at(upper_index) * fraction;
}

ResidualSummary CalculateResidualSummary(const std::vector<double> &values) {
  ResidualSummary summary;
  if (values.empty()) {
    return summary;
  }

  std::vector<double> sorted_values = values;
  std::sort(sorted_values.begin(), sorted_values.end());

  summary.n = static_cast<Long64_t>(sorted_values.size());
  summary.mean = std::accumulate(sorted_values.begin(),
                                 sorted_values.end(),
                                 0.0) / static_cast<double>(summary.n);
  summary.q16 = CalculateQuantileFromSortedValues(sorted_values, 0.16);
  summary.median = CalculateQuantileFromSortedValues(sorted_values, 0.50);
  summary.q84 = CalculateQuantileFromSortedValues(sorted_values, 0.84);
  summary.sigma68 = 0.5 * (summary.q84 - summary.q16);
  summary.valid = std::isfinite(summary.mean) &&
                  std::isfinite(summary.median) &&
                  std::isfinite(summary.sigma68);
  return summary;
}

const std::vector<double> &GetResidualValues(
    const ResidualValueMap &residual_values,
    const TH1D *hist) {
  static const std::vector<double> empty_values;
  const auto it = residual_values.find(hist);
  if (it == residual_values.end()) {
    return empty_values;
  }
  return it->second;
}

void FillResidual(TH1D *hist,
                  ResidualValueMap &residual_values,
                  double value) {
  if (!hist || !std::isfinite(value)) {
    return;
  }

  hist->Fill(value);
  residual_values[hist].push_back(value);
}

bool IsValidValue(double value) {
  return std::isfinite(value) && value > kNonInitializedThreshold;
}

std::vector<std::string> FindRootFiles(const std::string &input_dir) {
  std::vector<std::string> root_files;

  TSystemDirectory directory("input_directory", input_dir.c_str());
  std::unique_ptr<TList> files(directory.GetListOfFiles());

  if (!files) {
    std::cerr << "Error: cannot list directory: " << input_dir << std::endl;
    return root_files;
  }

  TIter next(files.get());
  while (auto *object = next()) {
    auto *file = dynamic_cast<TSystemFile *>(object);
    if (!file) {
      continue;
    }

    TString name = file->GetName();
    if (file->IsDirectory()) {
      continue;
    }
    if (!name.EndsWith(".root")) {
      continue;
    }

    TString full_path = input_dir.c_str();
    if (!full_path.EndsWith("/")) {
      full_path += "/";
    }
    full_path += name;
    root_files.emplace_back(full_path.Data());
  }

  std::sort(root_files.begin(), root_files.end());
  return root_files;
}

template <typename T>
bool SetVectorBranchAddress(TTree *tree,
                            const char *branch_name,
                            std::vector<T> **branch_ptr) {
  if (!tree->GetBranch(branch_name)) {
    std::cerr << "Error: required branch is missing: "
              << branch_name << std::endl;
    return false;
  }

  tree->SetBranchAddress(branch_name, branch_ptr);
  return true;
}

template <typename T>
bool HasIndex(const std::vector<T> *values, std::size_t index) {
  return values && index < values->size();
}

int FindAngleBin(double angle_deg) {
  for (int i = 0; i < kNAngleBins; ++i) {
    if (angle_deg >= kAngleBins[i] && angle_deg < kAngleBins[i + 1]) {
      return i;
    }
  }
  return -1;
}

std::vector<TH1D*> MakeAngleHistograms(const char *name_prefix,
                                       const char *theta_label,
                                       const char *x_axis_title,
                                       int nbins,
                                       double xmin_mm,
                                       double xmax_mm) {
  std::vector<TH1D*> histograms;
  histograms.reserve(kNAngleBins);

  for (int i = 0; i < kNAngleBins; ++i) {
    histograms.push_back(new TH1D(
      Form("%s_angle_%02d", name_prefix, i),
      Form("%.0f #leq #theta_{%s} < %.0f deg;%s;Number of tracks",
           kAngleBins[i], theta_label, kAngleBins[i + 1], x_axis_title),
      nbins,
      xmin_mm,
      xmax_mm));
  }

  return histograms;
}

std::vector<TH2D*> MakeAngleCorrelationHistograms(const char *name_prefix,
                                                  const char *theta_label,
                                                  const char *x_axis_title,
                                                  const char *y_axis_title,
                                                  int nbins,
                                                  double xmin_mm,
                                                  double xmax_mm) {
  std::vector<TH2D*> histograms;
  histograms.reserve(kNAngleBins);

  for (int i = 0; i < kNAngleBins; ++i) {
    histograms.push_back(new TH2D(
      Form("%s_angle_%02d", name_prefix, i),
      Form("%.0f #leq #theta_{%s} < %.0f deg;%s;%s",
           kAngleBins[i], theta_label, kAngleBins[i + 1],
           x_axis_title, y_axis_title),
      nbins,
      xmin_mm,
      xmax_mm,
      nbins,
      xmin_mm,
      xmax_mm));
  }

  return histograms;
}

void ProcessOneFile(const std::string &file_path,
                    SampleType sample_type,
                    TH1D *hist_x,
                    TH1D *hist_y,
                    TH1D *hist_reco_external_x,
                    TH1D *hist_reco_external_y,
                    TH1D *hist_pm_dwg_split_x,
                    TH1D *hist_pm_dwg_split_y,
                    TH1D *hist_reco_truth_x,
                    TH1D *hist_reco_truth_y,
                    TH2D *hist_correlation_x,
                    TH2D *hist_correlation_y,
                    std::vector<TH1D*> &hist_external_truth_x_by_angle,
                    std::vector<TH1D*> &hist_external_truth_y_by_angle,
                    std::vector<TH1D*> &hist_reco_external_x_by_angle,
                    std::vector<TH1D*> &hist_reco_external_y_by_angle,
                    std::vector<TH1D*> &hist_pm_dwg_split_x_by_angle,
                    std::vector<TH1D*> &hist_pm_dwg_split_y_by_angle,
                    std::vector<TH1D*> &hist_reco_truth_x_by_angle,
                    std::vector<TH1D*> &hist_reco_truth_y_by_angle,
                    std::vector<TH2D*> &hist_correlation_x_by_angle,
                    std::vector<TH2D*> &hist_correlation_y_by_angle,
                    ResidualValueMap &residual_values,
                    CorrelationStats &correlation_stats_x,
                    CorrelationStats &correlation_stats_y,
                    std::vector<CorrelationStats> &correlation_stats_x_by_angle,
                    std::vector<CorrelationStats> &correlation_stats_y_by_angle,
                    Long64_t &n_tracks_x,
                    Long64_t &n_tracks_y,
                    Long64_t &n_tracks_pm_dwg_split_x,
                    Long64_t &n_tracks_pm_dwg_split_y) {
  std::unique_ptr<TFile> file(TFile::Open(file_path.c_str(), "READ"));
  if (!file || file->IsZombie()) {
    std::cerr << "Warning: cannot open file: " << file_path << std::endl;
    return;
  }

  auto *tree = dynamic_cast<TTree *>(file->Get("match_info"));
  if (!tree) {
    std::cerr << "Warning: match_info tree is missing in "
              << file_path << std::endl;
    return;
  }

  std::vector<int> *has_match = nullptr;
  std::vector<int> *n_dwg_x = nullptr;
  std::vector<int> *n_dwg_y = nullptr;
  std::vector<int> *n_pm_x = nullptr;
  std::vector<int> *n_pm_y = nullptr;
  std::vector<double> *external_expected_x = nullptr;
  std::vector<double> *external_expected_y = nullptr;
  std::vector<double> *pm_only_expected_x = nullptr;
  std::vector<double> *pm_only_expected_y = nullptr;
  std::vector<double> *dwg_only_expected_x = nullptr;
  std::vector<double> *dwg_only_expected_y = nullptr;
  std::vector<double> *frost_nearest_x = nullptr;
  std::vector<double> *frost_nearest_y = nullptr;
  std::vector<double> *true_frost_x = nullptr;
  std::vector<double> *true_frost_y = nullptr;
  std::vector<double> *baby_mind_tangent_x = nullptr;
  std::vector<double> *baby_mind_tangent_y = nullptr;
  std::vector<int> *ninja_track_type = nullptr;

  Int_t bsd_good_spill_flag = 0;
  Int_t detector_flags[8] = {};

  bool ok = true;
  ok &= SetVectorBranchAddress(tree, "trackmatch_has_match", &has_match);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_external_num_planes_downstream_wagasci_x", &n_dwg_x);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_external_num_planes_downstream_wagasci_y", &n_dwg_y);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_external_num_planes_proton_module_x", &n_pm_x);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_external_num_planes_proton_module_y", &n_pm_y);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_external_expected_x", &external_expected_x);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_external_expected_y", &external_expected_y);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_external_pm_only_expected_x", &pm_only_expected_x);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_external_pm_only_expected_y", &pm_only_expected_y);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_external_dwg_only_expected_x", &dwg_only_expected_x);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_external_dwg_only_expected_y", &dwg_only_expected_y);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_frost_nearest_x", &frost_nearest_x);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_frost_nearest_y", &frost_nearest_y);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_baby_mind_tangent_x", &baby_mind_tangent_x);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_baby_mind_tangent_y", &baby_mind_tangent_y);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_ninja_track_type", &ninja_track_type);

  if (IsMonteCarlo(sample_type)) {
    ok &= SetVectorBranchAddress(
      tree, "trackmatch_true_frost_nearest_position_x", &true_frost_x);
    ok &= SetVectorBranchAddress(
      tree, "trackmatch_true_frost_nearest_position_y", &true_frost_y);
  }

  if (IsData(sample_type)) {
    if (!tree->GetBranch("bsd_good_spill_flag")) {
      std::cerr << "Error: required branch is missing: bsd_good_spill_flag"
                << std::endl;
      ok = false;
    } else {
      tree->SetBranchAddress("bsd_good_spill_flag", &bsd_good_spill_flag);
    }

    if (!tree->GetBranch("detector_flags")) {
      std::cerr << "Error: required branch is missing: detector_flags"
                << std::endl;
      ok = false;
    } else {
      tree->SetBranchAddress("detector_flags", detector_flags);
    }
  }

  if (!ok) {
    std::cerr << "Warning: skip file because required branches are missing: "
              << file_path << std::endl;
    return;
  }

  const Long64_t nentries = tree->GetEntries();
  for (Long64_t entry = 0; entry < nentries; ++entry) {
    tree->GetEntry(entry);

    if (IsData(sample_type)) {
      if (bsd_good_spill_flag == 0) {
        continue;
      }
      if (detector_flags[0] != 1 || detector_flags[2] != 1) {
        continue;
      }
    }

    if (!has_match) {
      continue;
    }

    const std::size_t ntracks = has_match->size();
    for (std::size_t itrack = 0; itrack < ntracks; ++itrack) {
      if (has_match->at(itrack) != 1) {
        continue;
      }

      if (!HasIndex(ninja_track_type, itrack) ||
          !HasIndex(frost_nearest_x, itrack) ||
          !HasIndex(frost_nearest_y, itrack) ||
          !IsValidValue(frost_nearest_x->at(itrack)) ||
          !IsValidValue(frost_nearest_y->at(itrack))) {
        continue;
      }

      const bool pass_track_type =
        (ninja_track_type->at(itrack) == 1);

      const bool pass_position =
        (std::abs(frost_nearest_x->at(itrack)) < ACCEPTANCE_X &&
         std::abs(frost_nearest_y->at(itrack)) < ACCEPTANCE_Y);

      if (!pass_track_type) {
        continue;
      }

      if (!pass_position) {
        continue;
      }
      if (!HasIndex(baby_mind_tangent_x, itrack) ||
          !HasIndex(baby_mind_tangent_y, itrack) ||
          !IsValidValue(baby_mind_tangent_x->at(itrack)) ||
          !IsValidValue(baby_mind_tangent_y->at(itrack))) {
        continue;
      }

      const double tangent_x = baby_mind_tangent_x->at(itrack);
      const double tangent_y = baby_mind_tangent_y->at(itrack);
      const double theta_x_deg =
        std::atan(std::abs(tangent_x)) * 180.0 / TMath::Pi();
      const double theta_y_deg =
        std::atan(std::abs(tangent_y)) * 180.0 / TMath::Pi();

      const int angle_bin_x = FindAngleBin(theta_x_deg);
      const int angle_bin_y = FindAngleBin(theta_y_deg);
      const bool has_angle_bin_x = (angle_bin_x >= 0);
      const bool has_angle_bin_y = (angle_bin_y >= 0);

      const bool common_external_plane_selected =
        HasIndex(n_dwg_x, itrack) &&
        HasIndex(n_pm_x, itrack) &&
        HasIndex(n_dwg_y, itrack) &&
        HasIndex(n_pm_y, itrack) &&
        n_dwg_x->at(itrack) >= NHIT_DWG_X &&
        n_pm_x->at(itrack) >= NHIT_PM_X &&
        n_dwg_y->at(itrack) >= NHIT_DWG_Y &&
        n_pm_y->at(itrack) >= NHIT_PM_Y;

      if (!common_external_plane_selected) {
        continue;
      }

      const bool split_x_selected =
        HasIndex(pm_only_expected_x, itrack) &&
        HasIndex(dwg_only_expected_x, itrack) &&
        IsValidValue(pm_only_expected_x->at(itrack)) &&
        IsValidValue(dwg_only_expected_x->at(itrack));

      if (split_x_selected) {
        const double pm_minus_dwg_x =
          pm_only_expected_x->at(itrack) - dwg_only_expected_x->at(itrack);
        FillResidual(hist_pm_dwg_split_x,
                     residual_values,
                     pm_minus_dwg_x);
        ++n_tracks_pm_dwg_split_x;

        if (has_angle_bin_x) {
          FillResidual(hist_pm_dwg_split_x_by_angle.at(angle_bin_x),
                       residual_values,
                       pm_minus_dwg_x);
        }
      }

      const bool x_base_selected =
        HasIndex(external_expected_x, itrack) &&
        IsValidValue(external_expected_x->at(itrack));

      if (x_base_selected &&
          HasIndex(frost_nearest_x, itrack) &&
          IsValidValue(frost_nearest_x->at(itrack))) {
        FillResidual(hist_reco_external_x,
                     residual_values,
                     frost_nearest_x->at(itrack) -
                     external_expected_x->at(itrack));
        ++n_tracks_x;

        if (has_angle_bin_x) {
          FillResidual(hist_reco_external_x_by_angle.at(angle_bin_x),
                       residual_values,
                       frost_nearest_x->at(itrack) -
                       external_expected_x->at(itrack));
        }
      }

      if (IsMonteCarlo(sample_type) &&
          x_base_selected &&
          HasIndex(true_frost_x, itrack) &&
          IsValidValue(true_frost_x->at(itrack))) {
        FillResidual(hist_x,
                     residual_values,
                     external_expected_x->at(itrack) -
                     true_frost_x->at(itrack));
        if (has_angle_bin_x) {
          FillResidual(hist_external_truth_x_by_angle.at(angle_bin_x),
                       residual_values,
                       external_expected_x->at(itrack) -
                       true_frost_x->at(itrack));
        }
      }

      if (IsMonteCarlo(sample_type) &&
          x_base_selected &&
          HasIndex(frost_nearest_x, itrack) &&
          HasIndex(true_frost_x, itrack) &&
          IsValidValue(frost_nearest_x->at(itrack)) &&
          IsValidValue(true_frost_x->at(itrack))) {
        const double external_minus_true_x =
          external_expected_x->at(itrack) - true_frost_x->at(itrack);
        const double reco_minus_true_x =
          frost_nearest_x->at(itrack) - true_frost_x->at(itrack);

        FillResidual(hist_reco_truth_x,
                     residual_values,
                     reco_minus_true_x);
        hist_correlation_x->Fill(external_minus_true_x, reco_minus_true_x);
        correlation_stats_x.Fill(external_minus_true_x, reco_minus_true_x);

        if (has_angle_bin_x) {
          FillResidual(hist_reco_truth_x_by_angle.at(angle_bin_x),
                       residual_values,
                       reco_minus_true_x);
          hist_correlation_x_by_angle.at(angle_bin_x)->Fill(
            external_minus_true_x, reco_minus_true_x);
          correlation_stats_x_by_angle.at(angle_bin_x).Fill(
            external_minus_true_x, reco_minus_true_x);
        }
      }

      const bool split_y_selected =
        HasIndex(pm_only_expected_y, itrack) &&
        HasIndex(dwg_only_expected_y, itrack) &&
        IsValidValue(pm_only_expected_y->at(itrack)) &&
        IsValidValue(dwg_only_expected_y->at(itrack));

      if (split_y_selected) {
        const double pm_minus_dwg_y =
          pm_only_expected_y->at(itrack) - dwg_only_expected_y->at(itrack);
        FillResidual(hist_pm_dwg_split_y,
                     residual_values,
                     pm_minus_dwg_y);
        ++n_tracks_pm_dwg_split_y;

        if (has_angle_bin_y) {
          FillResidual(hist_pm_dwg_split_y_by_angle.at(angle_bin_y),
                       residual_values,
                       pm_minus_dwg_y);
        }
      }

      const bool y_base_selected =
        HasIndex(external_expected_y, itrack) &&
        IsValidValue(external_expected_y->at(itrack));

      if (y_base_selected &&
          HasIndex(frost_nearest_y, itrack) &&
          IsValidValue(frost_nearest_y->at(itrack))) {
        FillResidual(hist_reco_external_y,
                     residual_values,
                     frost_nearest_y->at(itrack) -
                     external_expected_y->at(itrack));
        ++n_tracks_y;

        if (has_angle_bin_y) {
          FillResidual(hist_reco_external_y_by_angle.at(angle_bin_y),
                       residual_values,
                       frost_nearest_y->at(itrack) -
                       external_expected_y->at(itrack));
        }
      }

      if (IsMonteCarlo(sample_type) &&
          y_base_selected &&
          HasIndex(true_frost_y, itrack) &&
          IsValidValue(true_frost_y->at(itrack))) {
        FillResidual(hist_y,
                     residual_values,
                     external_expected_y->at(itrack) -
                     true_frost_y->at(itrack));

        if (has_angle_bin_y) {
          FillResidual(hist_external_truth_y_by_angle.at(angle_bin_y),
                       residual_values,
                       external_expected_y->at(itrack) -
                       true_frost_y->at(itrack));
        }
      }

      if (IsMonteCarlo(sample_type) &&
          y_base_selected &&
          HasIndex(frost_nearest_y, itrack) &&
          HasIndex(true_frost_y, itrack) &&
          IsValidValue(frost_nearest_y->at(itrack)) &&
          IsValidValue(true_frost_y->at(itrack))) {
        const double external_minus_true_y =
          external_expected_y->at(itrack) - true_frost_y->at(itrack);
        const double reco_minus_true_y =
          frost_nearest_y->at(itrack) - true_frost_y->at(itrack);

        FillResidual(hist_reco_truth_y,
                     residual_values,
                     reco_minus_true_y);
        hist_correlation_y->Fill(external_minus_true_y, reco_minus_true_y);
        correlation_stats_y.Fill(external_minus_true_y, reco_minus_true_y);

        if (has_angle_bin_y) {
          FillResidual(hist_reco_truth_y_by_angle.at(angle_bin_y),
                       residual_values,
                       reco_minus_true_y);
          hist_correlation_y_by_angle.at(angle_bin_y)->Fill(
            external_minus_true_y, reco_minus_true_y);
          correlation_stats_y_by_angle.at(angle_bin_y).Fill(
            external_minus_true_y, reco_minus_true_y);
        }
      }
    }
  }
}

void DrawCentral68Text(const std::vector<double> &values) {
  const ResidualSummary summary = CalculateResidualSummary(values);
  if (!summary.valid) {
    return;
  }

  auto *text = new TPaveText(0.6, 0.40, 0.85, 0.50, "NDC");
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->SetTextSize(0.035);
  text->AddText(Form("median = %.3f mm", summary.median));
  text->AddText(Form("#sigma_{68} = %.3f mm", summary.sigma68));
  text->Draw("same");
}

void DrawHistogram(TH1D *hist,
                   const ResidualValueMap &residual_values) {
  gPad->SetLeftMargin(0.17);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.13);

  hist->SetLineWidth(2);
  hist->GetXaxis()->SetTitleOffset(1.15);
  hist->GetYaxis()->SetTitleOffset(1.65);
  hist->Draw("hist");

  if (hist->GetEntries() > 0) {
    TString fit_name = "fit_";
    fit_name += hist->GetName();

    auto *fit = new TF1(
      fit_name,
      "gaus",
      hist->GetXaxis()->GetXmin(),
      hist->GetXaxis()->GetXmax());

    hist->Fit(fit, "Q");
    fit->SetLineWidth(2);
    fit->Draw("same");
  }

  DrawCentral68Text(GetResidualValues(residual_values, hist));
}

void DrawTwoPanelPage(TCanvas &canvas,
                      TH1D *hist_left,
                      TH1D *hist_right,
                      const ResidualValueMap &residual_values) {
  canvas.Clear();
  canvas.Divide(2, 1);

  canvas.cd(1);
  DrawHistogram(hist_left, residual_values);

  canvas.cd(2);
  DrawHistogram(hist_right, residual_values);
}

void DrawCorrelationText(const CorrelationStats &stats) {
  double mean_e = 0.;
  double mean_r = 0.;
  double rms_e = 0.;
  double rms_r = 0.;
  double covariance = 0.;
  double correlation = 0.;

  if (!stats.Calculate(mean_e, mean_r, rms_e, rms_r,
                       covariance, correlation)) {
    return;
  }

  auto *text = new TPaveText(0.15, 0.68, 0.55, 0.88, "NDC");
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->SetTextSize(0.032);
  text->AddText(Form("N = %lld", stats.n));
  text->AddText(Form("Cov(E,R) = %.3f mm^{2}", covariance));
  text->AddText(Form("#rho(E,R) = %.3f", correlation));
  text->Draw("same");
}

void DrawCorrelationHistogram(TH2D *hist, const CorrelationStats &stats) {
  gPad->SetLeftMargin(0.17);
  gPad->SetRightMargin(0.16);
  gPad->SetBottomMargin(0.13);

  hist->GetXaxis()->SetTitleOffset(1.15);
  hist->GetYaxis()->SetTitleOffset(1.65);
  hist->Draw("colz");
  DrawCorrelationText(stats);
}

void DrawCorrelationTwoPanelPage(TCanvas &canvas,
                                 TH2D *hist_x,
                                 TH2D *hist_y,
                                 const CorrelationStats &stats_x,
                                 const CorrelationStats &stats_y) {
  canvas.Clear();
  canvas.Divide(2, 1);

  canvas.cd(1);
  DrawCorrelationHistogram(hist_x, stats_x);

  canvas.cd(2);
  DrawCorrelationHistogram(hist_y, stats_y);
}

void PrintAnglePages(TCanvas &canvas,
                     const char *output_pdf,
                     const std::vector<TH1D*> &hist_x_by_angle,
                     const std::vector<TH1D*> &hist_y_by_angle,
                     const ResidualValueMap &residual_values,
                     bool close_pdf = false) {
  TString pdf_close = output_pdf;
  pdf_close += ")";

  for (int i = 0; i < kNAngleBins; ++i) {
    DrawTwoPanelPage(canvas,
                     hist_x_by_angle.at(i),
                     hist_y_by_angle.at(i),
                     residual_values);

    if (close_pdf && i == kNAngleBins - 1) {
      canvas.Print(pdf_close);
    } else {
      canvas.Print(output_pdf);
    }
  }
}

void PrintAngleCorrelationPages(
    TCanvas &canvas,
    const char *output_pdf,
    const std::vector<TH2D*> &hist_x_by_angle,
    const std::vector<TH2D*> &hist_y_by_angle,
    const std::vector<CorrelationStats> &stats_x_by_angle,
    const std::vector<CorrelationStats> &stats_y_by_angle,
    bool close_pdf = false) {
  TString pdf_close = output_pdf;
  pdf_close += ")";

  for (int i = 0; i < kNAngleBins; ++i) {
    DrawCorrelationTwoPanelPage(canvas,
                                hist_x_by_angle.at(i),
                                hist_y_by_angle.at(i),
                                stats_x_by_angle.at(i),
                                stats_y_by_angle.at(i));

    if (close_pdf && i == kNAngleBins - 1) {
      canvas.Print(pdf_close);
    } else {
      canvas.Print(output_pdf);
    }
  }
}

const char *SampleTypeLabel(SampleType sample_type) {
  return IsMonteCarlo(sample_type) ? "mc" : "data";
}

TString MakeSummaryLogFileName(const char *output_pdf,
                               const char *output_log) {
  if (output_log && std::string(output_log).size() > 0) {
    return TString(output_log);
  }

  TString log_file = output_pdf;
  if (log_file.EndsWith(".pdf")) {
    log_file.Remove(log_file.Length() - 4);
  }
  log_file += ".log";
  return log_file;
}

void WriteResidualSummaryRow(std::ostream &out,
                             SampleType sample_type,
                             const char *quantity_name,
                             const char *angle_bin_label,
                             const std::vector<double> &values) {
  const ResidualSummary summary = CalculateResidualSummary(values);

  out << SampleTypeLabel(sample_type) << ','
      << quantity_name << ','
      << angle_bin_label << ',';

  if (!summary.valid) {
    out << 0 << ",,,,," << '\n';
    return;
  }

  out << summary.n << ','
      << std::setprecision(10) << summary.mean << ','
      << summary.median << ','
      << summary.sigma68 << ','
      << summary.q16 << ','
      << summary.q84 << '\n';
}

void WriteResidualSummaryRowsForAngles(
    std::ostream &out,
    SampleType sample_type,
    const char *quantity_name,
    const char *theta_label,
    const std::vector<TH1D*> &hist_by_angle,
    const ResidualValueMap &residual_values) {
  for (int i = 0; i < kNAngleBins; ++i) {
    const TString angle_bin_label =
      Form("theta_%s_%.0f_%.0f_deg",
           theta_label,
           kAngleBins[i],
           kAngleBins[i + 1]);
    WriteResidualSummaryRow(out,
                            sample_type,
                            quantity_name,
                            angle_bin_label.Data(),
                            GetResidualValues(residual_values,
                                              hist_by_angle.at(i)));
  }
}

void WriteResidualSummaryLog(
    const char *output_log,
    SampleType sample_type,
    const ResidualValueMap &residual_values,
    TH1D *hist_external_truth_x,
    TH1D *hist_external_truth_y,
    TH1D *hist_reco_external_x,
    TH1D *hist_reco_external_y,
    TH1D *hist_pm_dwg_split_x,
    TH1D *hist_pm_dwg_split_y,
    TH1D *hist_reco_truth_x,
    TH1D *hist_reco_truth_y,
    const std::vector<TH1D*> &hist_external_truth_x_by_angle,
    const std::vector<TH1D*> &hist_external_truth_y_by_angle,
    const std::vector<TH1D*> &hist_reco_external_x_by_angle,
    const std::vector<TH1D*> &hist_reco_external_y_by_angle,
    const std::vector<TH1D*> &hist_pm_dwg_split_x_by_angle,
    const std::vector<TH1D*> &hist_pm_dwg_split_y_by_angle,
    const std::vector<TH1D*> &hist_reco_truth_x_by_angle,
    const std::vector<TH1D*> &hist_reco_truth_y_by_angle) {
  std::ofstream log(output_log);
  if (!log) {
    std::cerr << "Warning: cannot open summary log file: "
              << output_log << std::endl;
    return;
  }

  log << "sample_type,quantity,angle_bin,n,mean_mm,median_mm,"
      << "sigma68_mm,q16_mm,q84_mm" << '\n';

  if (IsMonteCarlo(sample_type)) {
    WriteResidualSummaryRow(log, sample_type, "x_ext_minus_x_true", "all",
                            GetResidualValues(residual_values,
                                              hist_external_truth_x));
    WriteResidualSummaryRow(log, sample_type, "y_ext_minus_y_true", "all",
                            GetResidualValues(residual_values,
                                              hist_external_truth_y));
  }

  WriteResidualSummaryRow(log, sample_type, "x_rec_minus_x_ext", "all",
                          GetResidualValues(residual_values,
                                            hist_reco_external_x));
  WriteResidualSummaryRow(log, sample_type, "y_rec_minus_y_ext", "all",
                          GetResidualValues(residual_values,
                                            hist_reco_external_y));

  WriteResidualSummaryRow(log,
                          sample_type,
                          "x_pm_only_minus_x_dwg_only",
                          "all",
                          GetResidualValues(residual_values,
                                            hist_pm_dwg_split_x));
  WriteResidualSummaryRow(log,
                          sample_type,
                          "y_pm_only_minus_y_dwg_only",
                          "all",
                          GetResidualValues(residual_values,
                                            hist_pm_dwg_split_y));

  if (IsMonteCarlo(sample_type)) {
    WriteResidualSummaryRow(log, sample_type, "x_rec_minus_x_true", "all",
                            GetResidualValues(residual_values,
                                              hist_reco_truth_x));
    WriteResidualSummaryRow(log, sample_type, "y_rec_minus_y_true", "all",
                            GetResidualValues(residual_values,
                                              hist_reco_truth_y));

    WriteResidualSummaryRowsForAngles(log,
                                      sample_type,
                                      "x_ext_minus_x_true",
                                      "x",
                                      hist_external_truth_x_by_angle,
                                      residual_values);
    WriteResidualSummaryRowsForAngles(log,
                                      sample_type,
                                      "y_ext_minus_y_true",
                                      "y",
                                      hist_external_truth_y_by_angle,
                                      residual_values);
  }

  WriteResidualSummaryRowsForAngles(log,
                                    sample_type,
                                    "x_rec_minus_x_ext",
                                    "x",
                                    hist_reco_external_x_by_angle,
                                    residual_values);
  WriteResidualSummaryRowsForAngles(log,
                                    sample_type,
                                    "y_rec_minus_y_ext",
                                    "y",
                                    hist_reco_external_y_by_angle,
                                    residual_values);

  WriteResidualSummaryRowsForAngles(log,
                                    sample_type,
                                    "x_pm_only_minus_x_dwg_only",
                                    "x",
                                    hist_pm_dwg_split_x_by_angle,
                                    residual_values);
  WriteResidualSummaryRowsForAngles(log,
                                    sample_type,
                                    "y_pm_only_minus_y_dwg_only",
                                    "y",
                                    hist_pm_dwg_split_y_by_angle,
                                    residual_values);

  if (IsMonteCarlo(sample_type)) {
    WriteResidualSummaryRowsForAngles(log,
                                      sample_type,
                                      "x_rec_minus_x_true",
                                      "x",
                                      hist_reco_truth_x_by_angle,
                                      residual_values);
    WriteResidualSummaryRowsForAngles(log,
                                      sample_type,
                                      "y_rec_minus_y_true",
                                      "y",
                                      hist_reco_truth_y_by_angle,
                                      residual_values);
  }
}

}  // namespace

// void DrawExternalFrostTruthResidual(const char *input_dir="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/6-TrackMatch_externalfit_PMandDWG",
//                                     const char *output_pdf="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/6-TrackMatch_externalfit_PMandDWG/external_fit_diagnostics.pdf",
//                                     const char *sample_type_arg="mc",
//                                     int nbins = 200,
//                                     double xmin_mm = -50.,
//                                     double xmax_mm = 50.,
//                                     const char *output_log="") {
void DrawExternalFrostTruthResidual(const char *input_dir="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM/2-rootfile_after_TrackMatch_externalfit_PMandDWG",
                                    const char *output_pdf="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM/2-rootfile_after_TrackMatch_externalfit_PMandDWG/external_fit_diagnostics.pdf",
                                    const char *sample_type_arg="data",
                                    int nbins = 200,
                                    double xmin_mm = -50.,
                                    double xmax_mm = 50.,
                                    const char *output_log="") {
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1111);

  const SampleType sample_type = ParseSampleType(sample_type_arg);

  const std::vector<std::string> root_files = FindRootFiles(input_dir);
  if (root_files.empty()) {
    std::cerr << "Error: no .root files found in directory: "
              << input_dir << std::endl;
    return;
  }

  auto *hist_x = new TH1D(
    "hist_external_truth_residual_x",
    "All angles;x_{ext}-x_{true} [mm];Number of tracks",
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_y = new TH1D(
    "hist_external_truth_residual_y",
    "All angles;y_{ext}-y_{true} [mm];Number of tracks",
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_reco_external_x = new TH1D(
    "hist_reco_external_residual_x",
    "All angles;x_{rec}-x_{ext} [mm];Number of tracks",
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_reco_external_y = new TH1D(
    "hist_reco_external_residual_y",
    "All angles;y_{rec}-y_{ext} [mm];Number of tracks",
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_pm_dwg_split_x = new TH1D(
    "hist_pm_dwg_split_residual_x",
    "All angles;x_{PM-only}-x_{DWG-only} [mm];Number of tracks",
    nbins,
    2*xmin_mm,
    2*xmax_mm);

  auto *hist_pm_dwg_split_y = new TH1D(
    "hist_pm_dwg_split_residual_y",
    "All angles;y_{PM-only}-y_{DWG-only} [mm];Number of tracks",
    nbins,
    2*xmin_mm,
    2*xmax_mm);

  auto *hist_reco_truth_x = new TH1D(
    "hist_reco_truth_residual_x",
    "All angles;x_{rec}-x_{true} [mm];Number of tracks",
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_reco_truth_y = new TH1D(
    "hist_reco_truth_residual_y",
    "All angles;y_{rec}-y_{true} [mm];Number of tracks",
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_correlation_x = new TH2D(
    "hist_correlation_x",
    "All angles;x_{ext}-x_{true} [mm];x_{rec}-x_{true} [mm]",
    nbins,
    xmin_mm,
    xmax_mm,
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_correlation_y = new TH2D(
    "hist_correlation_y",
    "All angles;y_{ext}-y_{true} [mm];y_{rec}-y_{true} [mm]",
    nbins,
    xmin_mm,
    xmax_mm,
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_correlation_x_by_angle = MakeAngleCorrelationHistograms(
    "hist_correlation_x",
    "x",
    "x_{ext}-x_{true} [mm]",
    "x_{rec}-x_{true} [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_correlation_y_by_angle = MakeAngleCorrelationHistograms(
    "hist_correlation_y",
    "y",
    "y_{ext}-y_{true} [mm]",
    "y_{rec}-y_{true} [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_external_truth_x_by_angle = MakeAngleHistograms(
    "hist_external_truth_residual_x",
    "x",
    "x_{ext}-x_{true} [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_external_truth_y_by_angle = MakeAngleHistograms(
    "hist_external_truth_residual_y",
    "y",
    "y_{ext}-y_{true} [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_reco_external_x_by_angle = MakeAngleHistograms(
    "hist_reco_external_residual_x",
    "x",
    "x_{rec}-x_{ext} [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_reco_external_y_by_angle = MakeAngleHistograms(
    "hist_reco_external_residual_y",
    "y",
    "y_{rec}-y_{ext} [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_pm_dwg_split_x_by_angle = MakeAngleHistograms(
    "hist_pm_dwg_split_residual_x",
    "x",
    "x_{PM-only}-x_{DWG-only} [mm]",
    nbins,
    2*xmin_mm,
    2*xmax_mm);

  auto hist_pm_dwg_split_y_by_angle = MakeAngleHistograms(
    "hist_pm_dwg_split_residual_y",
    "y",
    "y_{PM-only}-y_{DWG-only} [mm]",
    nbins,
    2*xmin_mm,
    2*xmax_mm);

  auto hist_reco_truth_x_by_angle = MakeAngleHistograms(
    "hist_reco_truth_residual_x",
    "x",
    "x_{rec}-x_{true} [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_reco_truth_y_by_angle = MakeAngleHistograms(
    "hist_reco_truth_residual_y",
    "y",
    "y_{rec}-y_{true} [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  ResidualValueMap residual_values;

  CorrelationStats correlation_stats_x;
  CorrelationStats correlation_stats_y;
  std::vector<CorrelationStats> correlation_stats_x_by_angle(kNAngleBins);
  std::vector<CorrelationStats> correlation_stats_y_by_angle(kNAngleBins);

  Long64_t n_tracks_x = 0;
  Long64_t n_tracks_y = 0;
  Long64_t n_tracks_pm_dwg_split_x = 0;
  Long64_t n_tracks_pm_dwg_split_y = 0;

  for (const auto &file_path : root_files) {
    std::cout << "Processing: " << file_path << std::endl;
    ProcessOneFile(file_path,
                   sample_type,
                   hist_x,
                   hist_y,
                   hist_reco_external_x,
                   hist_reco_external_y,
                   hist_pm_dwg_split_x,
                   hist_pm_dwg_split_y,
                   hist_reco_truth_x,
                   hist_reco_truth_y,
                   hist_correlation_x,
                   hist_correlation_y,
                   hist_external_truth_x_by_angle,
                   hist_external_truth_y_by_angle,
                   hist_reco_external_x_by_angle,
                   hist_reco_external_y_by_angle,
                   hist_pm_dwg_split_x_by_angle,
                   hist_pm_dwg_split_y_by_angle,
                   hist_reco_truth_x_by_angle,
                   hist_reco_truth_y_by_angle,
                   hist_correlation_x_by_angle,
                   hist_correlation_y_by_angle,
                   residual_values,
                   correlation_stats_x,
                   correlation_stats_y,
                   correlation_stats_x_by_angle,
                   correlation_stats_y_by_angle,
                   n_tracks_x,
                   n_tracks_y,
                   n_tracks_pm_dwg_split_x,
                   n_tracks_pm_dwg_split_y);
  }

  std::cout << "Selected tracks for x: " << n_tracks_x << std::endl;
  std::cout << "Selected tracks for y: " << n_tracks_y << std::endl;
  std::cout << "Selected tracks for PM-only minus DWG-only x: "
            << n_tracks_pm_dwg_split_x << std::endl;
  std::cout << "Selected tracks for PM-only minus DWG-only y: "
            << n_tracks_pm_dwg_split_y << std::endl;

  if (IsMonteCarlo(sample_type)) {
    double mean_e = 0.;
    double mean_r = 0.;
    double rms_e = 0.;
    double rms_r = 0.;
    double covariance = 0.;
    double correlation = 0.;

    if (correlation_stats_x.Calculate(mean_e, mean_r, rms_e, rms_r,
                                      covariance, correlation)) {
      std::cout << "All-angle x correlation:"
                << " N=" << correlation_stats_x.n
                << ", Cov(E,R)=" << covariance << " mm^2"
                << ", rho=" << correlation << std::endl;
    }
    if (correlation_stats_y.Calculate(mean_e, mean_r, rms_e, rms_r,
                                      covariance, correlation)) {
      std::cout << "All-angle y correlation:"
                << " N=" << correlation_stats_y.n
                << ", Cov(E,R)=" << covariance << " mm^2"
                << ", rho=" << correlation << std::endl;
    }

    for (int i = 0; i < kNAngleBins; ++i) {
      if (correlation_stats_x_by_angle.at(i).Calculate(
            mean_e, mean_r, rms_e, rms_r, covariance, correlation)) {
        std::cout << Form("Angle %.0f-%.0f deg theta_x correlation: ",
                          kAngleBins[i], kAngleBins[i + 1])
                  << "N=" << correlation_stats_x_by_angle.at(i).n
                  << ", Cov(E,R)=" << covariance << " mm^2"
                  << ", rho=" << correlation << std::endl;
      }
      if (correlation_stats_y_by_angle.at(i).Calculate(
            mean_e, mean_r, rms_e, rms_r, covariance, correlation)) {
        std::cout << Form("Angle %.0f-%.0f deg theta_y correlation: ",
                          kAngleBins[i], kAngleBins[i + 1])
                  << "N=" << correlation_stats_y_by_angle.at(i).n
                  << ", Cov(E,R)=" << covariance << " mm^2"
                  << ", rho=" << correlation << std::endl;
      }
    }
  }

  const TString output_log_file = MakeSummaryLogFileName(output_pdf, output_log);
  WriteResidualSummaryLog(output_log_file.Data(),
                          sample_type,
                          residual_values,
                          hist_x,
                          hist_y,
                          hist_reco_external_x,
                          hist_reco_external_y,
                          hist_pm_dwg_split_x,
                          hist_pm_dwg_split_y,
                          hist_reco_truth_x,
                          hist_reco_truth_y,
                          hist_external_truth_x_by_angle,
                          hist_external_truth_y_by_angle,
                          hist_reco_external_x_by_angle,
                          hist_reco_external_y_by_angle,
                          hist_pm_dwg_split_x_by_angle,
                          hist_pm_dwg_split_y_by_angle,
                          hist_reco_truth_x_by_angle,
                          hist_reco_truth_y_by_angle);
  std::cout << "Wrote residual summary log: "
            << output_log_file.Data() << std::endl;

  TCanvas canvas("canvas", "FROST external-fit diagnostics", 1200, 600);
  TString pdf_open = output_pdf;
  pdf_open += "(";
  TString pdf_close = output_pdf;
  pdf_close += ")";

  if (IsMonteCarlo(sample_type)) {
    DrawTwoPanelPage(canvas, hist_x, hist_y, residual_values);
    canvas.Print(pdf_open);

    DrawTwoPanelPage(canvas, hist_reco_external_x, hist_reco_external_y, residual_values);
    canvas.Print(output_pdf);

    DrawTwoPanelPage(canvas, hist_pm_dwg_split_x, hist_pm_dwg_split_y, residual_values);
    canvas.Print(output_pdf);

    DrawTwoPanelPage(canvas, hist_reco_truth_x, hist_reco_truth_y, residual_values);
    canvas.Print(output_pdf);

    PrintAnglePages(canvas,
                    output_pdf,
                    hist_external_truth_x_by_angle,
                    hist_external_truth_y_by_angle,
                    residual_values);

    PrintAnglePages(canvas,
                    output_pdf,
                    hist_reco_external_x_by_angle,
                    hist_reco_external_y_by_angle,
                    residual_values);

    PrintAnglePages(canvas,
                    output_pdf,
                    hist_pm_dwg_split_x_by_angle,
                    hist_pm_dwg_split_y_by_angle,
                    residual_values);

    PrintAnglePages(canvas,
                    output_pdf,
                    hist_reco_truth_x_by_angle,
                    hist_reco_truth_y_by_angle,
                    residual_values);

    DrawCorrelationTwoPanelPage(canvas,
                                hist_correlation_x,
                                hist_correlation_y,
                                correlation_stats_x,
                                correlation_stats_y);
    canvas.Print(output_pdf);

    PrintAngleCorrelationPages(canvas,
                               output_pdf,
                               hist_correlation_x_by_angle,
                               hist_correlation_y_by_angle,
                               correlation_stats_x_by_angle,
                               correlation_stats_y_by_angle,
                               true);
  } else {
    DrawTwoPanelPage(canvas, hist_reco_external_x, hist_reco_external_y, residual_values);
    canvas.Print(pdf_open);

    DrawTwoPanelPage(canvas, hist_pm_dwg_split_x, hist_pm_dwg_split_y, residual_values);
    canvas.Print(output_pdf);

    PrintAnglePages(canvas,
                    output_pdf,
                    hist_reco_external_x_by_angle,
                    hist_reco_external_y_by_angle,
                    residual_values);

    PrintAnglePages(canvas,
                    output_pdf,
                    hist_pm_dwg_split_x_by_angle,
                    hist_pm_dwg_split_y_by_angle,
                    residual_values,
                    true);
  }
}
