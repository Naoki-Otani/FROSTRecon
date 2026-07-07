// DrawFrostExternalFitDiagnostics.C
//
// Usage in ROOT:
//   MC:
//     root -l -b -q 'DrawFrostExternalFitDiagnostics.C("input_dir", "output.pdf", "mc")'
//   data:
//     root -l -b -q 'DrawFrostExternalFitDiagnostics.C("input_dir", "output.pdf", "data")'
//
//   Multiple input directories can be given as a comma-separated list:
//     root -l -b -q 'DrawFrostExternalFitDiagnostics.C("input_dir1,input_dir2", "output.pdf", "mc")'
//
//   To disable the position reweighting:
//     root -l -b -q 'DrawFrostExternalFitDiagnostics.C("input_dir", "output.pdf", "mc", 200, -50, 50, "", false)'
//
// The macro reads all .root files directly under the input directory or
// directories, loops over the match_info tree, and saves a multi-page PDF for
// FROST external-track diagnostics.
//
// The median, mean, and sigma_68 values written on the plots and in the
// summary log are calculated directly from selected residual values within the
// displayed histogram range. Underflow and overflow values are not used for
// these residual statistics, but they are still filled into the histograms as
// ROOT underflow/overflow entries.
//
// By default, the selected tracks are reweighted to reduce biases from the
// FROST hit-position distribution. The accepted FROST area,
//   abs(frost_nearest_x) < ACCEPTANCE_X
//   abs(frost_nearest_y) < ACCEPTANCE_Y
// is divided into 10 x 10 bins in the x-y plane. The event weight is assigned
// so that each occupied position bin has the same total weight. The histogram
// contents and the residual statistics in the summary log are calculated with
// these position weights. Empty position bins are ignored.
// In the summary log, n is the number of unweighted residual values inside the
// displayed histogram range, while sum_weight is the total weight used for the
// weighted statistics.
//
// In MC mode, the output contains:
//   - x_{ext}-x_{true} and y_{ext}-y_{true}
//   - x_{rec}-x_{ext} and y_{rec}-y_{ext}
//   - x_{PM-only}-x_{DWG-only} and y_{PM-only}-y_{DWG-only}
//   - x_{rec}-x_{true} and y_{rec}-y_{true}
//   - 2D correlation plots of x_{ext}-x_{true} versus x_{rec}-x_{true}
//     and y_{ext}-y_{true} versus y_{rec}-y_{true}
//   - final pages with plane-count histograms for downstream WAGASCI and
//     Proton Module in x and y, for all angles and for each angle bin
//
// In data mode, truth information is not available, so only
// x_{rec}-x_{ext}, y_{rec}-y_{ext}, and PM-only minus DWG-only split-fit
// residual pages are drawn. The final pages also include plane-count
// histograms for downstream WAGASCI and Proton Module in x and y, for
// all angles and for each angle bin.
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
#include <TH1.h>
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
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr int NHIT_DWG_X = 4;
constexpr int NHIT_PM_X = 10;
constexpr int NHIT_DWG_Y = 4;
constexpr int NHIT_PM_Y = 10;

constexpr double ACCEPTANCE_X = 560.0; //mm
constexpr double ACCEPTANCE_Y = 600.0; //mm

constexpr int kNPositionBinsX = 10;
constexpr int kNPositionBinsY = 10;
constexpr int kNPositionBins = kNPositionBinsX * kNPositionBinsY;

constexpr double kNonInitializedThreshold = -9.0e7;

// Angle bins [deg]:
//   [0,5), [5,10), [10,15), [15,20), [20,30),
//   [30,50)
constexpr int kNAngleBins = 6;
const double kAngleBins[kNAngleBins + 1] = {
  0.0, 5.0, 10.0, 15.0, 20.0,
  30.0, 50.0
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
  double sum_weight = 0.;
  double sum_e = 0.;
  double sum_r = 0.;
  double sum_ee = 0.;
  double sum_rr = 0.;
  double sum_er = 0.;

  void Fill(double external_minus_true,
            double reco_minus_true,
            double weight) {
    if (!std::isfinite(weight) || weight <= 0.) {
      return;
    }

    ++n;
    sum_weight += weight;
    sum_e += weight * external_minus_true;
    sum_r += weight * reco_minus_true;
    sum_ee += weight * external_minus_true * external_minus_true;
    sum_rr += weight * reco_minus_true * reco_minus_true;
    sum_er += weight * external_minus_true * reco_minus_true;
  }

  bool Calculate(double &mean_e,
                 double &mean_r,
                 double &rms_e,
                 double &rms_r,
                 double &covariance,
                 double &correlation) const {
    if (n < 2 || sum_weight <= 0.) {
      return false;
    }

    mean_e = sum_e / sum_weight;
    mean_r = sum_r / sum_weight;

    const double variance_e =
      sum_ee / sum_weight - mean_e * mean_e;
    const double variance_r =
      sum_rr / sum_weight - mean_r * mean_r;

    if (variance_e <= 0. || variance_r <= 0.) {
      return false;
    }

    rms_e = std::sqrt(variance_e);
    rms_r = std::sqrt(variance_r);
    covariance = sum_er / sum_weight - mean_e * mean_r;
    correlation = covariance / (rms_e * rms_r);

    return std::isfinite(correlation);
  }
};
struct WeightedResidualValue {
  double value = 0.;
  double weight = 1.;
};

struct ResidualSummary {
  bool valid = false;
  Long64_t n = 0;
  double sum_weight = 0.;
  double mean = 0.;
  double median = 0.;
  double sigma68 = 0.;
  double q16 = 0.;
  double q84 = 0.;
};

using ResidualValueMap = std::map<const TH1D*, std::vector<WeightedResidualValue>>;

double CalculateWeightedQuantileFromSortedValues(
    const std::vector<WeightedResidualValue> &sorted_values,
    double probability,
    double sum_weight) {
  if (sorted_values.empty() || sum_weight <= 0.) {
    return 0.;
  }

  if (probability <= 0.) {
    return sorted_values.front().value;
  }
  if (probability >= 1.) {
    return sorted_values.back().value;
  }

  const double target_weight = probability * sum_weight;
  double cumulative_weight = 0.;

  for (const auto &entry : sorted_values) {
    cumulative_weight += entry.weight;
    if (cumulative_weight >= target_weight) {
      return entry.value;
    }
  }

  return sorted_values.back().value;
}

ResidualSummary CalculateResidualSummary(
    const std::vector<WeightedResidualValue> &values) {
  ResidualSummary summary;
  if (values.empty()) {
    return summary;
  }

  std::vector<WeightedResidualValue> sorted_values;
  sorted_values.reserve(values.size());
  for (const auto &entry : values) {
    if (std::isfinite(entry.value) &&
        std::isfinite(entry.weight) &&
        entry.weight > 0.) {
      sorted_values.push_back(entry);
    }
  }

  if (sorted_values.empty()) {
    return summary;
  }

  std::sort(sorted_values.begin(),
            sorted_values.end(),
            [](const WeightedResidualValue &lhs,
               const WeightedResidualValue &rhs) {
              return lhs.value < rhs.value;
            });

  summary.n = static_cast<Long64_t>(sorted_values.size());
  double weighted_sum = 0.;
  for (const auto &entry : sorted_values) {
    summary.sum_weight += entry.weight;
    weighted_sum += entry.weight * entry.value;
  }

  if (summary.sum_weight <= 0.) {
    return summary;
  }

  summary.mean = weighted_sum / summary.sum_weight;
  summary.q16 =
    CalculateWeightedQuantileFromSortedValues(sorted_values,
                                              0.16,
                                              summary.sum_weight);
  summary.median =
    CalculateWeightedQuantileFromSortedValues(sorted_values,
                                              0.50,
                                              summary.sum_weight);
  summary.q84 =
    CalculateWeightedQuantileFromSortedValues(sorted_values,
                                              0.84,
                                              summary.sum_weight);
  summary.sigma68 = 0.5 * (summary.q84 - summary.q16);
  summary.valid = std::isfinite(summary.mean) &&
                  std::isfinite(summary.median) &&
                  std::isfinite(summary.sigma68);
  return summary;
}

const std::vector<WeightedResidualValue> &GetResidualValues(
    const ResidualValueMap &residual_values,
    const TH1D *hist) {
  static const std::vector<WeightedResidualValue> empty_values;
  const auto it = residual_values.find(hist);
  if (it == residual_values.end()) {
    return empty_values;
  }
  return it->second;
}

bool IsInsideHistogramRange(const TH1D *hist, double value) {
  if (!hist || !std::isfinite(value)) {
    return false;
  }

  const TAxis *axis = hist->GetXaxis();
  if (!axis) {
    return false;
  }

  return value >= axis->GetXmin() && value < axis->GetXmax();
}

void FillResidual(TH1D *hist,
                  ResidualValueMap &residual_values,
                  double value,
                  double weight) {
  if (!hist ||
      !std::isfinite(value) ||
      !std::isfinite(weight) ||
      weight <= 0.) {
    return;
  }

  hist->Fill(value, weight);
  // Keep the raw residual value for mean/median/sigma_68 calculations only
  // when it is inside the displayed histogram range. This avoids TH1 binning
  // effects while excluding underflow and overflow outliers from the summary
  // statistics.
  if (IsInsideHistogramRange(hist, value)) {
    residual_values[hist].push_back({value, weight});
  }
}

bool IsValidValue(double value) {
  return std::isfinite(value) && value > kNonInitializedThreshold;
}

std::string TrimString(const std::string &input) {
  const std::string whitespace = " \t\n\r";
  const std::size_t first = input.find_first_not_of(whitespace);
  if (first == std::string::npos) {
    return "";
  }

  const std::size_t last = input.find_last_not_of(whitespace);
  return input.substr(first, last - first + 1);
}

std::vector<std::string> ParseInputDirectories(const char *input_dirs_arg) {
  std::vector<std::string> input_dirs;
  if (!input_dirs_arg) {
    return input_dirs;
  }

  std::string arg = input_dirs_arg;
  for (char &character : arg) {
    if (character == ';' || character == '\n') {
      character = ',';
    }
  }

  std::stringstream stream(arg);
  std::string item;
  while (std::getline(stream, item, ',')) {
    const std::string directory = TrimString(item);
    if (!directory.empty()) {
      input_dirs.push_back(directory);
    }
  }

  return input_dirs;
}

std::vector<std::string> FindRootFilesInDirectory(
    const std::string &input_dir) {
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

std::vector<std::string> FindRootFiles(
    const std::vector<std::string> &input_dirs) {
  std::vector<std::string> root_files;

  for (const auto &input_dir : input_dirs) {
    const std::vector<std::string> files =
      FindRootFilesInDirectory(input_dir);
    root_files.insert(root_files.end(), files.begin(), files.end());
  }

  std::sort(root_files.begin(), root_files.end());
  root_files.erase(std::unique(root_files.begin(), root_files.end()),
                   root_files.end());
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

int FindPositionBin(double value, double acceptance, int n_bins) {
  if (!std::isfinite(value) || std::abs(value) >= acceptance) {
    return -1;
  }

  const double fraction = (value + acceptance) / (2. * acceptance);
  int bin = static_cast<int>(std::floor(fraction * n_bins));
  if (bin < 0) {
    bin = 0;
  }
  if (bin >= n_bins) {
    bin = n_bins - 1;
  }
  return bin;
}

int FindPositionBinIndex(double frost_x, double frost_y) {
  const int bin_x = FindPositionBin(frost_x,
                                    ACCEPTANCE_X,
                                    kNPositionBinsX);
  const int bin_y = FindPositionBin(frost_y,
                                    ACCEPTANCE_Y,
                                    kNPositionBinsY);
  if (bin_x < 0 || bin_y < 0) {
    return -1;
  }
  return bin_y * kNPositionBinsX + bin_x;
}

struct PositionBinWeights {
  std::array<Long64_t, kNPositionBins> counts{};
  std::array<double, kNPositionBins> weights{};
  Long64_t total_tracks = 0;
  int occupied_bins = 0;
  double target_count = 0.;
  bool enabled = true;

  void Count(double frost_x, double frost_y) {
    const int bin_index = FindPositionBinIndex(frost_x, frost_y);
    if (bin_index < 0) {
      return;
    }

    ++counts.at(bin_index);
    ++total_tracks;
  }

  void Disable() {
    enabled = false;
    weights.fill(1.);
  }

  void CalculateWeights() {
    enabled = true;
    occupied_bins = 0;
    for (const auto count : counts) {
      if (count > 0) {
        ++occupied_bins;
      }
    }

    weights.fill(0.);
    if (occupied_bins == 0) {
      return;
    }

    target_count =
      static_cast<double>(total_tracks) / static_cast<double>(occupied_bins);

    for (std::size_t i = 0; i < counts.size(); ++i) {
      if (counts.at(i) > 0) {
        weights.at(i) = target_count / static_cast<double>(counts.at(i));
      }
    }
  }

  double GetWeight(double frost_x, double frost_y) const {
    const int bin_index = FindPositionBinIndex(frost_x, frost_y);
    if (bin_index < 0) {
      return 0.;
    }
    if (!enabled) {
      return 1.;
    }
    return weights.at(bin_index);
  }
};

void CountPositionBinsOneFile(const std::string &file_path,
                              SampleType sample_type,
                              PositionBinWeights &position_weights) {
  std::unique_ptr<TFile> file(TFile::Open(file_path.c_str(), "READ"));
  if (!file || file->IsZombie()) {
    std::cerr << "Warning: cannot open file for position-bin counting: "
              << file_path << std::endl;
    return;
  }

  auto *tree = dynamic_cast<TTree *>(file->Get("match_info"));
  if (!tree) {
    std::cerr << "Warning: match_info tree is missing in "
              << file_path << std::endl;
    return;
  }

  std::vector<int> *has_match = nullptr;
  std::vector<int> *ninja_track_type = nullptr;
  std::vector<double> *frost_nearest_x = nullptr;
  std::vector<double> *frost_nearest_y = nullptr;
  std::vector<double> *external_tangent_x = nullptr;
  std::vector<double> *external_tangent_y = nullptr;

  Int_t bsd_good_spill_flag = 0;
  Int_t detector_flags[8] = {};

  bool ok = true;
  ok &= SetVectorBranchAddress(tree, "trackmatch_has_match", &has_match);
  ok &= SetVectorBranchAddress(tree,
                               "trackmatch_ninja_track_type",
                               &ninja_track_type);
  ok &= SetVectorBranchAddress(tree,
                               "trackmatch_frost_nearest_x",
                               &frost_nearest_x);
  ok &= SetVectorBranchAddress(tree,
                               "trackmatch_frost_nearest_y",
                               &frost_nearest_y);
  ok &= SetVectorBranchAddress(tree,
                               "trackmatch_external_tangent_x",
                               &external_tangent_x);
  ok &= SetVectorBranchAddress(tree,
                               "trackmatch_external_tangent_y",
                               &external_tangent_y);

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
    std::cerr << "Warning: skip file for position-bin counting because "
              << "required branches are missing: "
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
          !HasIndex(external_tangent_x, itrack) ||
          !HasIndex(external_tangent_y, itrack) ||
          !IsValidValue(frost_nearest_x->at(itrack)) ||
          !IsValidValue(frost_nearest_y->at(itrack)) ||
          !IsValidValue(external_tangent_x->at(itrack)) ||
          !IsValidValue(external_tangent_y->at(itrack))) {
        continue;
      }

      if (ninja_track_type->at(itrack) != 1) {
        continue;
      }

      const double frost_x = frost_nearest_x->at(itrack);
      const double frost_y = frost_nearest_y->at(itrack);
      if (std::abs(frost_x) >= ACCEPTANCE_X ||
          std::abs(frost_y) >= ACCEPTANCE_Y) {
        continue;
      }

      position_weights.Count(frost_x, frost_y);
    }
  }
}

void PrintPositionWeightSummary(const PositionBinWeights &position_weights) {
  std::cout << "Position reweighting: "
            << kNPositionBinsX << " x " << kNPositionBinsY
            << " bins in |x| < " << ACCEPTANCE_X
            << " mm, |y| < " << ACCEPTANCE_Y << " mm" << std::endl;
  std::cout << "  total selected tracks for weighting: "
            << position_weights.total_tracks << std::endl;
  std::cout << "  occupied position bins: "
            << position_weights.occupied_bins << " / "
            << kNPositionBins << std::endl;
  std::cout << "  target weighted entries per occupied bin: "
            << position_weights.target_count << std::endl;
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

std::vector<TH1D*> MakeAngleCountHistograms(const char *name_prefix,
                                            const char *theta_label,
                                            const char *x_axis_title,
                                            int min_value,
                                            int max_value) {
  std::vector<TH1D*> histograms;
  histograms.reserve(kNAngleBins);

  const int nbins = max_value - min_value + 1;
  const double xmin = min_value - 0.5;
  const double xmax = max_value + 0.5;

  for (int i = 0; i < kNAngleBins; ++i) {
    histograms.push_back(new TH1D(
      Form("%s_angle_%02d", name_prefix, i),
      Form("%.0f #leq #theta_{%s} < %.0f deg;%s;Number of tracks",
           kAngleBins[i], theta_label, kAngleBins[i + 1], x_axis_title),
      nbins,
      xmin,
      xmax));
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
                    TH1D *hist_num_planes_dwg_x_all,
                    TH1D *hist_num_planes_pm_x_all,
                    TH1D *hist_num_planes_dwg_y_all,
                    TH1D *hist_num_planes_pm_y_all,
                    std::vector<TH1D*> &hist_num_planes_dwg_x_by_angle,
                    std::vector<TH1D*> &hist_num_planes_pm_x_by_angle,
                    std::vector<TH1D*> &hist_num_planes_dwg_y_by_angle,
                    std::vector<TH1D*> &hist_num_planes_pm_y_by_angle,
                    const PositionBinWeights &position_weights,
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
  std::vector<double> *external_tangent_x = nullptr;
  std::vector<double> *external_tangent_y = nullptr;
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
    tree, "trackmatch_external_tangent_x", &external_tangent_x);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_external_tangent_y", &external_tangent_y);
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
      if (!HasIndex(external_tangent_x, itrack) ||
          !HasIndex(external_tangent_y, itrack) ||
          !IsValidValue(external_tangent_x->at(itrack)) ||
          !IsValidValue(external_tangent_y->at(itrack))) {
        continue;
      }

      const double position_weight =
        position_weights.GetWeight(frost_nearest_x->at(itrack),
                                   frost_nearest_y->at(itrack));
      if (position_weight <= 0.) {
        continue;
      }

      const double tangent_x = external_tangent_x->at(itrack);
      const double tangent_y = external_tangent_y->at(itrack);
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

      hist_num_planes_dwg_x_all->Fill(n_dwg_x->at(itrack), position_weight);
      hist_num_planes_pm_x_all->Fill(n_pm_x->at(itrack), position_weight);
      hist_num_planes_dwg_y_all->Fill(n_dwg_y->at(itrack), position_weight);
      hist_num_planes_pm_y_all->Fill(n_pm_y->at(itrack), position_weight);

      if (has_angle_bin_x) {
        hist_num_planes_dwg_x_by_angle.at(angle_bin_x)->Fill(
          n_dwg_x->at(itrack), position_weight);
        hist_num_planes_pm_x_by_angle.at(angle_bin_x)->Fill(
          n_pm_x->at(itrack), position_weight);
      }

      if (has_angle_bin_y) {
        hist_num_planes_dwg_y_by_angle.at(angle_bin_y)->Fill(
          n_dwg_y->at(itrack), position_weight);
        hist_num_planes_pm_y_by_angle.at(angle_bin_y)->Fill(
          n_pm_y->at(itrack), position_weight);
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
                     pm_minus_dwg_x,
                     position_weight);
        ++n_tracks_pm_dwg_split_x;

        if (has_angle_bin_x) {
          FillResidual(hist_pm_dwg_split_x_by_angle.at(angle_bin_x),
                       residual_values,
                       pm_minus_dwg_x,
                     position_weight);
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
                     external_expected_x->at(itrack),
                     position_weight);
        ++n_tracks_x;

        if (has_angle_bin_x) {
          FillResidual(hist_reco_external_x_by_angle.at(angle_bin_x),
                       residual_values,
                       frost_nearest_x->at(itrack) -
                       external_expected_x->at(itrack),
                     position_weight);
        }
      }

      if (IsMonteCarlo(sample_type) &&
          x_base_selected &&
          HasIndex(true_frost_x, itrack) &&
          IsValidValue(true_frost_x->at(itrack))) {
        FillResidual(hist_x,
                     residual_values,
                     external_expected_x->at(itrack) -
                     true_frost_x->at(itrack),
                     position_weight);
        if (has_angle_bin_x) {
          FillResidual(hist_external_truth_x_by_angle.at(angle_bin_x),
                       residual_values,
                       external_expected_x->at(itrack) -
                       true_frost_x->at(itrack),
                     position_weight);
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
                     reco_minus_true_x,
                     position_weight);
        hist_correlation_x->Fill(external_minus_true_x,
                                 reco_minus_true_x,
                                 position_weight);
        correlation_stats_x.Fill(external_minus_true_x,
                                 reco_minus_true_x,
                                 position_weight);

        if (has_angle_bin_x) {
          FillResidual(hist_reco_truth_x_by_angle.at(angle_bin_x),
                       residual_values,
                       reco_minus_true_x,
                     position_weight);
          hist_correlation_x_by_angle.at(angle_bin_x)->Fill(
            external_minus_true_x,
            reco_minus_true_x,
            position_weight);
          correlation_stats_x_by_angle.at(angle_bin_x).Fill(
            external_minus_true_x,
            reco_minus_true_x,
            position_weight);
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
                     pm_minus_dwg_y,
                     position_weight);
        ++n_tracks_pm_dwg_split_y;

        if (has_angle_bin_y) {
          FillResidual(hist_pm_dwg_split_y_by_angle.at(angle_bin_y),
                       residual_values,
                       pm_minus_dwg_y,
                     position_weight);
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
                     external_expected_y->at(itrack),
                     position_weight);
        ++n_tracks_y;

        if (has_angle_bin_y) {
          FillResidual(hist_reco_external_y_by_angle.at(angle_bin_y),
                       residual_values,
                       frost_nearest_y->at(itrack) -
                       external_expected_y->at(itrack),
                     position_weight);
        }
      }

      if (IsMonteCarlo(sample_type) &&
          y_base_selected &&
          HasIndex(true_frost_y, itrack) &&
          IsValidValue(true_frost_y->at(itrack))) {
        FillResidual(hist_y,
                     residual_values,
                     external_expected_y->at(itrack) -
                     true_frost_y->at(itrack),
                     position_weight);

        if (has_angle_bin_y) {
          FillResidual(hist_external_truth_y_by_angle.at(angle_bin_y),
                       residual_values,
                       external_expected_y->at(itrack) -
                       true_frost_y->at(itrack),
                     position_weight);
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
                     reco_minus_true_y,
                     position_weight);
        hist_correlation_y->Fill(external_minus_true_y,
                                 reco_minus_true_y,
                                 position_weight);
        correlation_stats_y.Fill(external_minus_true_y,
                                 reco_minus_true_y,
                                 position_weight);

        if (has_angle_bin_y) {
          FillResidual(hist_reco_truth_y_by_angle.at(angle_bin_y),
                       residual_values,
                       reco_minus_true_y,
                     position_weight);
          hist_correlation_y_by_angle.at(angle_bin_y)->Fill(
            external_minus_true_y,
            reco_minus_true_y,
            position_weight);
          correlation_stats_y_by_angle.at(angle_bin_y).Fill(
            external_minus_true_y,
            reco_minus_true_y,
            position_weight);
        }
      }
    }
  }
}

void DrawCentral68Text(const std::vector<WeightedResidualValue> &values) {
  const ResidualSummary summary = CalculateResidualSummary(values);
  if (!summary.valid) {
    return;
  }

  auto *text = new TPaveText(0.7, 0.40, 0.95, 0.50, "NDC");
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->SetTextSize(0.03);
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

void DrawPlaneCountHistogram(TH1D *hist) {
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.13);

  hist->SetLineWidth(2);
  hist->GetXaxis()->SetTitleOffset(1.15);
  hist->GetYaxis()->SetTitleOffset(1.55);
  hist->Draw("hist");
}

void DrawPlaneCountPage(TCanvas &canvas,
                        TH1D *hist_dwg_x,
                        TH1D *hist_pm_x,
                        TH1D *hist_dwg_y,
                        TH1D *hist_pm_y) {
  canvas.Clear();
  canvas.Divide(2, 2);

  canvas.cd(1);
  DrawPlaneCountHistogram(hist_dwg_x);

  canvas.cd(2);
  DrawPlaneCountHistogram(hist_dwg_y);

  canvas.cd(3);
  DrawPlaneCountHistogram(hist_pm_x);

  canvas.cd(4);
  DrawPlaneCountHistogram(hist_pm_y);
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
  text->AddText(Form("#Sigma w = %.1f", stats.sum_weight));
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

void PrintPlaneCountPages(
    TCanvas &canvas,
    const char *output_pdf,
    TH1D *hist_num_planes_dwg_x_all,
    TH1D *hist_num_planes_pm_x_all,
    TH1D *hist_num_planes_dwg_y_all,
    TH1D *hist_num_planes_pm_y_all,
    const std::vector<TH1D*> &hist_num_planes_dwg_x_by_angle,
    const std::vector<TH1D*> &hist_num_planes_pm_x_by_angle,
    const std::vector<TH1D*> &hist_num_planes_dwg_y_by_angle,
    const std::vector<TH1D*> &hist_num_planes_pm_y_by_angle,
    bool close_pdf = false) {
  TString pdf_close = output_pdf;
  pdf_close += ")";

  DrawPlaneCountPage(canvas,
                     hist_num_planes_dwg_x_all,
                     hist_num_planes_pm_x_all,
                     hist_num_planes_dwg_y_all,
                     hist_num_planes_pm_y_all);

  if (close_pdf && kNAngleBins == 0) {
    canvas.Print(pdf_close);
  } else {
    canvas.Print(output_pdf);
  }

  for (int i = 0; i < kNAngleBins; ++i) {
    DrawPlaneCountPage(canvas,
                       hist_num_planes_dwg_x_by_angle.at(i),
                       hist_num_planes_pm_x_by_angle.at(i),
                       hist_num_planes_dwg_y_by_angle.at(i),
                       hist_num_planes_pm_y_by_angle.at(i));

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
                             const std::vector<WeightedResidualValue> &values) {
  const ResidualSummary summary = CalculateResidualSummary(values);

  out << SampleTypeLabel(sample_type) << ','
      << quantity_name << ','
      << angle_bin_label << ',';

  if (!summary.valid) {
    out << 0 << ",,,,,," << '\n';
    return;
  }

  out << summary.n << ','
      << std::setprecision(10) << summary.sum_weight << ','
      << summary.mean << ','
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

  log << "sample_type,quantity,angle_bin,n,sum_weight,mean_mm,median_mm,"
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

// void DrawFrostExternalFitDiagnostics(const char *input_dirs="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC_retry/6-TrackMatch_externalfit_PMandDWG,/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/6-TrackMatch_externalfit_PMandDWG",
//                                     const char *output_pdf="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC_retry/6-TrackMatch_externalfit_PMandDWG/external_fit_diagnostics.pdf",
//                                     const char *sample_type_arg="mc",
//                                     int nbins = 200,
//                                     double xmin_mm = -50.,
//                                     double xmax_mm = 50.,
//                                     const char *output_log="",
//                                     bool use_position_reweighting = true) {
// void DrawFrostExternalFitDiagnostics(const char *input_dirs="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/6-TrackMatch_externalfit_PMandDWG",
//                                     const char *output_pdf="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/6-TrackMatch_externalfit_PMandDWG/external_fit_diagnostics.pdf",
//                                     const char *sample_type_arg="mc",
//                                     int nbins = 200,
//                                     double xmin_mm = -50.,
//                                     double xmax_mm = 50.,
//                                     const char *output_log="",
//                                     bool use_position_reweighting = true) {
void DrawFrostExternalFitDiagnostics(const char *input_dirs="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM/2-rootfile_after_TrackMatch_externalfit_PMandDWG_shift2",
                                    const char *output_pdf="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM/2-rootfile_after_TrackMatch_externalfit_PMandDWG_shift2/external_fit_diagnostics.pdf",
                                    const char *sample_type_arg="data",
                                    int nbins = 200,
                                    double xmin_mm = -50.,
                                    double xmax_mm = 50.,
                                    const char *output_log="",
                                    bool use_position_reweighting = true) {
  TH1::SetDefaultSumw2(kTRUE);

  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1111);

  const SampleType sample_type = ParseSampleType(sample_type_arg);

  const std::vector<std::string> input_dir_list =
    ParseInputDirectories(input_dirs);
  if (input_dir_list.empty()) {
    std::cerr << "Error: no input directory is specified" << std::endl;
    return;
  }

  std::cout << "Input directories:" << std::endl;
  for (const auto &input_dir : input_dir_list) {
    std::cout << "  " << input_dir << std::endl;
  }

  const std::vector<std::string> root_files = FindRootFiles(input_dir_list);
  if (root_files.empty()) {
    std::cerr << "Error: no .root files found in the input directories"
              << std::endl;
    return;
  }

  PositionBinWeights position_weights;
  if (use_position_reweighting) {
    for (const auto &file_path : root_files) {
      CountPositionBinsOneFile(file_path, sample_type, position_weights);
    }
    position_weights.CalculateWeights();
    PrintPositionWeightSummary(position_weights);
  } else {
    position_weights.Disable();
    std::cout << "Position reweighting: disabled" << std::endl;
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

  auto *hist_num_planes_dwg_x_all = new TH1D(
    "hist_num_planes_dwg_x_all",
    "All angles;N_{plane}^{DWG,x};Number of tracks",
    10 - NHIT_DWG_X + 1,
    NHIT_DWG_X - 0.5,
    10.5);

  auto *hist_num_planes_pm_x_all = new TH1D(
    "hist_num_planes_pm_x_all",
    "All angles;N_{plane}^{PM,x};Number of tracks",
    25 - NHIT_PM_X + 1,
    NHIT_PM_X - 0.5,
    25.5);

  auto *hist_num_planes_dwg_y_all = new TH1D(
    "hist_num_planes_dwg_y_all",
    "All angles;N_{plane}^{DWG,y};Number of tracks",
    10 - NHIT_DWG_Y + 1,
    NHIT_DWG_Y - 0.5,
    10.5);

  auto *hist_num_planes_pm_y_all = new TH1D(
    "hist_num_planes_pm_y_all",
    "All angles;N_{plane}^{PM,y};Number of tracks",
    25 - NHIT_PM_Y + 1,
    NHIT_PM_Y - 0.5,
    25.5);

  auto hist_num_planes_dwg_x_by_angle = MakeAngleCountHistograms(
    "hist_num_planes_dwg_x",
    "x",
    "N_{plane}^{DWG,x}",
    NHIT_DWG_X,
    10);

  auto hist_num_planes_pm_x_by_angle = MakeAngleCountHistograms(
    "hist_num_planes_pm_x",
    "x",
    "N_{plane}^{PM,x}",
    NHIT_PM_X,
    25);

  auto hist_num_planes_dwg_y_by_angle = MakeAngleCountHistograms(
    "hist_num_planes_dwg_y",
    "y",
    "N_{plane}^{DWG,y}",
    NHIT_DWG_Y,
    10);

  auto hist_num_planes_pm_y_by_angle = MakeAngleCountHistograms(
    "hist_num_planes_pm_y",
    "y",
    "N_{plane}^{PM,y}",
    NHIT_PM_Y,
    25);

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
                   hist_num_planes_dwg_x_all,
                   hist_num_planes_pm_x_all,
                   hist_num_planes_dwg_y_all,
                   hist_num_planes_pm_y_all,
                   hist_num_planes_dwg_x_by_angle,
                   hist_num_planes_pm_x_by_angle,
                   hist_num_planes_dwg_y_by_angle,
                   hist_num_planes_pm_y_by_angle,
                   position_weights,
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
                               false);

    PrintPlaneCountPages(canvas,
                         output_pdf,
                         hist_num_planes_dwg_x_all,
                         hist_num_planes_pm_x_all,
                         hist_num_planes_dwg_y_all,
                         hist_num_planes_pm_y_all,
                         hist_num_planes_dwg_x_by_angle,
                         hist_num_planes_pm_x_by_angle,
                         hist_num_planes_dwg_y_by_angle,
                         hist_num_planes_pm_y_by_angle,
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
                    false);

    PrintPlaneCountPages(canvas,
                         output_pdf,
                         hist_num_planes_dwg_x_all,
                         hist_num_planes_pm_x_all,
                         hist_num_planes_dwg_y_all,
                         hist_num_planes_pm_y_all,
                         hist_num_planes_dwg_x_by_angle,
                         hist_num_planes_pm_x_by_angle,
                         hist_num_planes_dwg_y_by_angle,
                         hist_num_planes_pm_y_by_angle,
                         true);
  }
}
