// CalculatePositionResolution.C
//
// Simplified diagnostic macro for the FROST position-resolution study.
//
// It produces a multi-page PDF, a CSV log and a ROOT log file.  The ROOT
// file also stores the histograms that are drawn in the PDF.
//
// Plots are produced for each angle bin, with x on the left pad and y on
// the right pad, in the following order:
//   1) PM-only - DWG-only, Data/MC comparison
//   2) x_rec - x_ext, Data/MC comparison
//   3) x_rec - x_true in MC
//   4) x_ext - x_true in MC
//
// The event selection, position acceptance and position reweighting follow
// the previous version of CalculatePositionResolution.C.  All median and
// sigma68 values are computed from weighted residual values within the
// displayed histogram range.
//
// Usage in ROOT:
//   root -l -b -q 'CalculatePositionResolution.C("mc_dir", "data_dir", "out.pdf")'
//
// Multiple input directories can be given as a comma-separated list.

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TList.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TString.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace {

// constexpr int NHIT_DWG_X = 4;
// constexpr int NHIT_PM_X = 10;
// constexpr int NHIT_DWG_Y = 4;
// constexpr int NHIT_PM_Y = 10;

constexpr int NHIT_DWG_X = 5;
constexpr int NHIT_PM_X = 10;
constexpr int NHIT_DWG_Y = 5;
constexpr int NHIT_PM_Y = 10;

constexpr double ACCEPTANCE_X = 560.0; // mm
constexpr double ACCEPTANCE_Y = 600.0; // mm

constexpr int kNPositionBinsX = 10;
constexpr int kNPositionBinsY = 10;
constexpr int kNPositionBins = kNPositionBinsX * kNPositionBinsY;

constexpr double kNonInitializedThreshold = -9.0e7;

// Angle bins [deg]: [0,10), [10,20), [20,50)
constexpr int kNAngleBins = 3;
const double kAngleBins[kNAngleBins + 1] = {0.0, 10.0, 20.0, 50.0};
constexpr int kNAnalysisBins = kNAngleBins + 1; // all + angle bins

constexpr double kRecTrueXMin = -10.0;
constexpr double kRecTrueXMax = 10.0;

constexpr double kTableTextSize = 0.017;

enum class SampleType {
  kMonteCarlo,
  kData
};

enum class Axis {
  kX,
  kY
};

enum class PlotKind {
  kSplit,
  kRecExt,
  kRecTrue,
  kExtTrue
};

bool IsMonteCarlo(SampleType sample_type) {
  return sample_type == SampleType::kMonteCarlo;
}

bool IsData(SampleType sample_type) {
  return sample_type == SampleType::kData;
}

const char *AxisLabel(Axis axis) {
  return axis == Axis::kX ? "x" : "y";
}

std::string AnalysisBinLabel(int analysis_bin) {
  if (analysis_bin == 0) {
    return "all";
  }

  const int angle_bin = analysis_bin - 1;
  return Form("%.0f_%.0f_deg", kAngleBins[angle_bin], kAngleBins[angle_bin + 1]);
}

std::string AnalysisBinTitle(int analysis_bin, Axis axis) {
  if (analysis_bin == 0) {
    return "All angles";
  }

  const int angle_bin = analysis_bin - 1;
  return Form("%.0f #leq #theta_{%s} < %.0f deg",
              kAngleBins[angle_bin],
              AxisLabel(axis),
              kAngleBins[angle_bin + 1]);
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

std::vector<std::string> FindRootFilesInDirectory(const std::string &input_dir) {
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
    if (file->IsDirectory() || !name.EndsWith(".root")) {
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

std::vector<std::string> FindRootFiles(const std::vector<std::string> &input_dirs) {
  std::vector<std::string> root_files;
  for (const auto &input_dir : input_dirs) {
    const auto files = FindRootFilesInDirectory(input_dir);
    root_files.insert(root_files.end(), files.begin(), files.end());
  }
  std::sort(root_files.begin(), root_files.end());
  root_files.erase(std::unique(root_files.begin(), root_files.end()), root_files.end());
  return root_files;
}

template <typename T>
bool SetVectorBranchAddress(TTree *tree,
                            const char *branch_name,
                            std::vector<T> **branch_ptr) {
  if (!tree->GetBranch(branch_name)) {
    std::cerr << "Error: required branch is missing: " << branch_name << std::endl;
    return false;
  }
  tree->SetBranchAddress(branch_name, branch_ptr);
  return true;
}

template <typename T>
bool HasIndex(const std::vector<T> *values, std::size_t index) {
  return values && index < values->size();
}

bool IsValidValue(double value) {
  return std::isfinite(value) && value > kNonInitializedThreshold;
}

int FindAngleBin(double angle_deg) {
  for (int i = 0; i < kNAngleBins; ++i) {
    if (angle_deg >= kAngleBins[i] && angle_deg < kAngleBins[i + 1]) {
      return i;
    }
  }
  return -1;
}

int FindPositionBin(double value, double acceptance, int n_bins) {
  if (!std::isfinite(value) || std::abs(value) >= acceptance) {
    return -1;
  }

  const double fraction = (value + acceptance) / (2. * acceptance);
  int bin = static_cast<int>(std::floor(fraction * n_bins));
  if (bin < 0) bin = 0;
  if (bin >= n_bins) bin = n_bins - 1;
  return bin;
}

int FindPositionBinIndex(double frost_x, double frost_y) {
  const int bin_x = FindPositionBin(frost_x, ACCEPTANCE_X, kNPositionBinsX);
  const int bin_y = FindPositionBin(frost_y, ACCEPTANCE_Y, kNPositionBinsY);
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

    target_count = static_cast<double>(total_tracks) /
                   static_cast<double>(occupied_bins);
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

struct WeightedValue {
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

struct Chi2Result {
  bool valid = false;
  int ndf = 0;
  int n_bins_used = 0;
  double chi2 = 0.;
  double chi2_ndf = 0.;
  double p_value = 0.;
  double q01 = 0.;
  double q99 = 0.;
};

double WeightedQuantileFromSortedValues(
    const std::vector<WeightedValue> &sorted_values,
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

ResidualSummary CalculateSummary(const std::vector<WeightedValue> &values,
                                 double xmin,
                                 double xmax,
                                 bool apply_range = true) {
  ResidualSummary summary;
  std::vector<WeightedValue> selected_values;
  selected_values.reserve(values.size());

  for (const auto &entry : values) {
    if (!std::isfinite(entry.value) ||
        !std::isfinite(entry.weight) ||
        entry.weight <= 0.) {
      continue;
    }
    if (apply_range && (entry.value < xmin || entry.value >= xmax)) {
      continue;
    }
    selected_values.push_back(entry);
  }

  if (selected_values.empty()) {
    return summary;
  }

  std::sort(selected_values.begin(), selected_values.end(),
            [](const WeightedValue &lhs, const WeightedValue &rhs) {
              return lhs.value < rhs.value;
            });

  double weighted_sum = 0.;
  summary.n = static_cast<Long64_t>(selected_values.size());
  for (const auto &entry : selected_values) {
    summary.sum_weight += entry.weight;
    weighted_sum += entry.value * entry.weight;
  }
  if (summary.sum_weight <= 0.) {
    return summary;
  }

  summary.mean = weighted_sum / summary.sum_weight;
  summary.q16 = WeightedQuantileFromSortedValues(selected_values, 0.16, summary.sum_weight);
  summary.median = WeightedQuantileFromSortedValues(selected_values, 0.50, summary.sum_weight);
  summary.q84 = WeightedQuantileFromSortedValues(selected_values, 0.84, summary.sum_weight);
  summary.sigma68 = 0.5 * (summary.q84 - summary.q16);
  summary.valid = std::isfinite(summary.mean) &&
                  std::isfinite(summary.median) &&
                  std::isfinite(summary.sigma68);
  return summary;
}

std::vector<WeightedValue> MakeAreaNormalizedValues(
    const std::vector<WeightedValue> &values) {
  std::vector<WeightedValue> normalized_values;
  normalized_values.reserve(values.size());

  double sum_weight = 0.;
  for (const auto &entry : values) {
    if (std::isfinite(entry.value) &&
        std::isfinite(entry.weight) &&
        entry.weight > 0.) {
      sum_weight += entry.weight;
    }
  }
  if (sum_weight <= 0.) {
    return normalized_values;
  }

  for (const auto &entry : values) {
    if (std::isfinite(entry.value) &&
        std::isfinite(entry.weight) &&
        entry.weight > 0.) {
      normalized_values.push_back({entry.value, entry.weight / sum_weight});
    }
  }
  return normalized_values;
}

std::vector<WeightedValue> MakeCombinedNormalizedValues(
    const std::vector<WeightedValue> &data_values,
    const std::vector<WeightedValue> &mc_values) {
  const auto normalized_data = MakeAreaNormalizedValues(data_values);
  const auto normalized_mc = MakeAreaNormalizedValues(mc_values);

  std::vector<WeightedValue> combined_values;
  combined_values.reserve(normalized_data.size() + normalized_mc.size());
  combined_values.insert(combined_values.end(), normalized_data.begin(), normalized_data.end());
  combined_values.insert(combined_values.end(), normalized_mc.begin(), normalized_mc.end());
  return combined_values;
}

double CalculateWeightedQuantile(std::vector<WeightedValue> values,
                                 double probability) {
  values.erase(std::remove_if(values.begin(), values.end(),
                              [](const WeightedValue &entry) {
                                return !std::isfinite(entry.value) ||
                                       !std::isfinite(entry.weight) ||
                                       entry.weight <= 0.;
                              }),
               values.end());
  if (values.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::sort(values.begin(), values.end(),
            [](const WeightedValue &lhs, const WeightedValue &rhs) {
              return lhs.value < rhs.value;
            });

  double sum_weight = 0.;
  for (const auto &entry : values) {
    sum_weight += entry.weight;
  }
  if (sum_weight <= 0.) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return WeightedQuantileFromSortedValues(values, probability, sum_weight);
}

Chi2Result CalculateChi2Q01Q99(const std::vector<WeightedValue> &data_values,
                               const std::vector<WeightedValue> &mc_values,
                               int nbins) {
  static int counter = 0;
  Chi2Result result;
  if (data_values.empty() || mc_values.empty() || nbins <= 0) {
    return result;
  }

  const auto combined_values = MakeCombinedNormalizedValues(data_values, mc_values);
  if (combined_values.empty()) {
    return result;
  }

  result.q01 = CalculateWeightedQuantile(combined_values, 0.01);
  result.q99 = CalculateWeightedQuantile(combined_values, 0.99);
  if (!std::isfinite(result.q01) ||
      !std::isfinite(result.q99) ||
      result.q01 >= result.q99) {
    return result;
  }

  TString data_name = Form("hist_chi2_data_%d", counter);
  TString mc_name = Form("hist_chi2_mc_%d", counter);
  ++counter;

  TH1D data_hist(data_name, "", nbins, result.q01, result.q99);
  TH1D mc_hist(mc_name, "", nbins, result.q01, result.q99);
  data_hist.Sumw2();
  mc_hist.Sumw2();

  for (const auto &entry : data_values) {
    if (!std::isfinite(entry.value) ||
        !std::isfinite(entry.weight) ||
        entry.weight <= 0. ||
        entry.value < result.q01 ||
        entry.value >= result.q99) {
      continue;
    }
    data_hist.Fill(entry.value, entry.weight);
  }

  for (const auto &entry : mc_values) {
    if (!std::isfinite(entry.value) ||
        !std::isfinite(entry.weight) ||
        entry.weight <= 0. ||
        entry.value < result.q01 ||
        entry.value >= result.q99) {
      continue;
    }
    mc_hist.Fill(entry.value, entry.weight);
  }

  const double data_integral = data_hist.Integral();
  const double mc_integral = mc_hist.Integral();
  if (data_integral <= 0. || mc_integral <= 0.) {
    return result;
  }

  data_hist.Scale(1. / data_integral);
  mc_hist.Scale(1. / mc_integral);

  for (int ibin = 1; ibin <= nbins; ++ibin) {
    const double data_content = data_hist.GetBinContent(ibin);
    const double mc_content = mc_hist.GetBinContent(ibin);
    const double data_error = data_hist.GetBinError(ibin);
    const double mc_error = mc_hist.GetBinError(ibin);
    const double error2 = data_error * data_error + mc_error * mc_error;
    if (error2 <= 0.) {
      continue;
    }
    const double diff = data_content - mc_content;
    result.chi2 += diff * diff / error2;
    ++result.n_bins_used;
  }

  result.ndf = result.n_bins_used - 1;
  if (result.ndf <= 0) {
    return result;
  }
  result.chi2_ndf = result.chi2 / static_cast<double>(result.ndf);
  result.p_value = TMath::Prob(result.chi2, result.ndf);
  result.valid = std::isfinite(result.chi2_ndf) &&
                 std::isfinite(result.p_value);
  return result;
}

struct McResidualEvent {
  double p = 0.; // PM-only - true
  double d = 0.; // DWG-only - true
  double e = 0.; // external combined - true
  double r = 0.; // reconstructed FROST - true
  double weight = 1.;
};

struct DataResidualEvent {
  double split = 0.; // PM-only - DWG-only
  double o = 0.;     // reconstructed FROST - external combined
  double weight = 1.;
  bool has_split = false;
  bool has_rec_ext = false;
};

struct ResidualSamples {
  std::array<std::vector<McResidualEvent>, kNAnalysisBins> mc_x;
  std::array<std::vector<McResidualEvent>, kNAnalysisBins> mc_y;
  std::array<std::vector<DataResidualEvent>, kNAnalysisBins> data_x;
  std::array<std::vector<DataResidualEvent>, kNAnalysisBins> data_y;
};

void AddMcEvent(std::array<std::vector<McResidualEvent>, kNAnalysisBins> &bins,
                int angle_bin,
                const McResidualEvent &event) {
  bins.at(0).push_back(event);
  if (angle_bin >= 0) {
    bins.at(angle_bin + 1).push_back(event);
  }
}

void AddDataEvent(std::array<std::vector<DataResidualEvent>, kNAnalysisBins> &bins,
                  int angle_bin,
                  const DataResidualEvent &event) {
  bins.at(0).push_back(event);
  if (angle_bin >= 0) {
    bins.at(angle_bin + 1).push_back(event);
  }
}

std::vector<WeightedValue> MakeMcSplitValues(
    const std::vector<McResidualEvent> &events) {
  std::vector<WeightedValue> values;
  values.reserve(events.size());
  for (const auto &event : events) {
    values.push_back({event.p - event.d, event.weight});
  }
  return values;
}

std::vector<WeightedValue> MakeMcRecExtValues(
    const std::vector<McResidualEvent> &events) {
  std::vector<WeightedValue> values;
  values.reserve(events.size());
  for (const auto &event : events) {
    values.push_back({event.r - event.e, event.weight});
  }
  return values;
}

std::vector<WeightedValue> MakeMcExtTrueValues(
    const std::vector<McResidualEvent> &events) {
  std::vector<WeightedValue> values;
  values.reserve(events.size());
  for (const auto &event : events) {
    values.push_back({event.e, event.weight});
  }
  return values;
}

std::vector<WeightedValue> MakeMcRecTrueValues(
    const std::vector<McResidualEvent> &events) {
  std::vector<WeightedValue> values;
  values.reserve(events.size());
  for (const auto &event : events) {
    values.push_back({event.r, event.weight});
  }
  return values;
}

std::vector<WeightedValue> MakeDataSplitValues(
    const std::vector<DataResidualEvent> &events) {
  std::vector<WeightedValue> values;
  values.reserve(events.size());
  for (const auto &event : events) {
    if (event.has_split) {
      values.push_back({event.split, event.weight});
    }
  }
  return values;
}

std::vector<WeightedValue> MakeDataRecExtValues(
    const std::vector<DataResidualEvent> &events) {
  std::vector<WeightedValue> values;
  values.reserve(events.size());
  for (const auto &event : events) {
    if (event.has_rec_ext) {
      values.push_back({event.o, event.weight});
    }
  }
  return values;
}

bool TrackPassesCommonSelection(
    std::size_t itrack,
    const std::vector<int> *n_dwg_x,
    const std::vector<int> *n_dwg_y,
    const std::vector<int> *n_pm_x,
    const std::vector<int> *n_pm_y) {
  return HasIndex(n_dwg_x, itrack) &&
         HasIndex(n_pm_x, itrack) &&
         HasIndex(n_dwg_y, itrack) &&
         HasIndex(n_pm_y, itrack) &&
         n_dwg_x->at(itrack) >= NHIT_DWG_X &&
         n_pm_x->at(itrack) >= NHIT_PM_X &&
         n_dwg_y->at(itrack) >= NHIT_DWG_Y &&
         n_pm_y->at(itrack) >= NHIT_PM_Y;
}

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
    std::cerr << "Warning: match_info tree is missing in " << file_path << std::endl;
    return;
  }

  auto *ntbm_tree = dynamic_cast<TTree *>(file->Get("ntbm"));
  if (!ntbm_tree) {
    std::cerr << "Warning: ntbm tree is missing in " << file_path << std::endl;
    return;
  }
  if (tree->GetEntries() != ntbm_tree->GetEntries()) {
    std::cerr << "Warning: match_info and ntbm have different numbers of entries in "
              << file_path << ": match_info=" << tree->GetEntries()
              << ", ntbm=" << ntbm_tree->GetEntries() << std::endl;
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
  Int_t number_of_tracks = 0;

  bool ok = true;
  if (!ntbm_tree->GetBranch("number_of_tracks_")) {
    std::cerr << "Error: required branch is missing in ntbm: number_of_tracks_"
              << std::endl;
    ok = false;
  } else {
    ntbm_tree->SetBranchAddress("number_of_tracks_", &number_of_tracks);
  }
  ok &= SetVectorBranchAddress(tree, "trackmatch_has_match", &has_match);
  ok &= SetVectorBranchAddress(tree, "trackmatch_ninja_track_type", &ninja_track_type);
  ok &= SetVectorBranchAddress(tree, "trackmatch_frost_nearest_x", &frost_nearest_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_frost_nearest_y", &frost_nearest_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_tangent_x", &external_tangent_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_tangent_y", &external_tangent_y);

  if (IsData(sample_type)) {
    if (!tree->GetBranch("bsd_good_spill_flag")) {
      std::cerr << "Error: required branch is missing: bsd_good_spill_flag" << std::endl;
      ok = false;
    } else {
      tree->SetBranchAddress("bsd_good_spill_flag", &bsd_good_spill_flag);
    }
    if (!tree->GetBranch("detector_flags")) {
      std::cerr << "Error: required branch is missing: detector_flags" << std::endl;
      ok = false;
    } else {
      tree->SetBranchAddress("detector_flags", detector_flags);
    }
  }

  if (!ok) {
    std::cerr << "Warning: skip file for position-bin counting because required branches are missing: "
              << file_path << std::endl;
    return;
  }

  const Long64_t nentries = tree->GetEntries();
  for (Long64_t entry = 0; entry < nentries; ++entry) {
    tree->GetEntry(entry);
    ntbm_tree->GetEntry(entry);

    if (number_of_tracks != 1) {
      continue;
    }

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

void PrintPositionWeightSummary(const char *label,
                                const PositionBinWeights &position_weights) {
  std::cout << label << " position reweighting: "
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

void LoadResidualSamplesOneFile(const std::string &file_path,
                                SampleType sample_type,
                                const PositionBinWeights &position_weights,
                                ResidualSamples &samples) {
  std::unique_ptr<TFile> file(TFile::Open(file_path.c_str(), "READ"));
  if (!file || file->IsZombie()) {
    std::cerr << "Warning: cannot open file: " << file_path << std::endl;
    return;
  }

  auto *tree = dynamic_cast<TTree *>(file->Get("match_info"));
  if (!tree) {
    std::cerr << "Warning: match_info tree is missing in " << file_path << std::endl;
    return;
  }

  auto *ntbm_tree = dynamic_cast<TTree *>(file->Get("ntbm"));
  if (!ntbm_tree) {
    std::cerr << "Warning: ntbm tree is missing in " << file_path << std::endl;
    return;
  }
  if (tree->GetEntries() != ntbm_tree->GetEntries()) {
    std::cerr << "Warning: match_info and ntbm have different numbers of entries in "
              << file_path << ": match_info=" << tree->GetEntries()
              << ", ntbm=" << ntbm_tree->GetEntries() << std::endl;
    return;
  }

  std::vector<int> *has_match = nullptr;
  std::vector<int> *n_dwg_x = nullptr;
  std::vector<int> *n_dwg_y = nullptr;
  std::vector<int> *n_pm_x = nullptr;
  std::vector<int> *n_pm_y = nullptr;
  std::vector<int> *ninja_track_type = nullptr;

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

  Int_t bsd_good_spill_flag = 0;
  Int_t detector_flags[8] = {};
  Int_t number_of_tracks = 0;

  bool ok = true;
  if (!ntbm_tree->GetBranch("number_of_tracks_")) {
    std::cerr << "Error: required branch is missing in ntbm: number_of_tracks_"
              << std::endl;
    ok = false;
  } else {
    ntbm_tree->SetBranchAddress("number_of_tracks_", &number_of_tracks);
  }
  ok &= SetVectorBranchAddress(tree, "trackmatch_has_match", &has_match);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_num_planes_downstream_wagasci_x", &n_dwg_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_num_planes_downstream_wagasci_y", &n_dwg_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_num_planes_proton_module_x", &n_pm_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_num_planes_proton_module_y", &n_pm_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_ninja_track_type", &ninja_track_type);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_expected_x", &external_expected_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_expected_y", &external_expected_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_pm_only_expected_x", &pm_only_expected_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_pm_only_expected_y", &pm_only_expected_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_dwg_only_expected_x", &dwg_only_expected_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_dwg_only_expected_y", &dwg_only_expected_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_frost_nearest_x", &frost_nearest_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_frost_nearest_y", &frost_nearest_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_tangent_x", &external_tangent_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_tangent_y", &external_tangent_y);

  if (IsMonteCarlo(sample_type)) {
    ok &= SetVectorBranchAddress(tree, "trackmatch_true_frost_nearest_position_x", &true_frost_x);
    ok &= SetVectorBranchAddress(tree, "trackmatch_true_frost_nearest_position_y", &true_frost_y);
  }

  if (IsData(sample_type)) {
    if (!tree->GetBranch("bsd_good_spill_flag")) {
      std::cerr << "Error: required branch is missing: bsd_good_spill_flag" << std::endl;
      ok = false;
    } else {
      tree->SetBranchAddress("bsd_good_spill_flag", &bsd_good_spill_flag);
    }
    if (!tree->GetBranch("detector_flags")) {
      std::cerr << "Error: required branch is missing: detector_flags" << std::endl;
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
    ntbm_tree->GetEntry(entry);

    if (number_of_tracks != 1) {
      continue;
    }

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
      if (!TrackPassesCommonSelection(itrack, n_dwg_x, n_dwg_y, n_pm_x, n_pm_y)) {
        continue;
      }

      const double weight = position_weights.GetWeight(frost_x, frost_y);
      if (!std::isfinite(weight) || weight <= 0.) {
        continue;
      }

      const double theta_x_deg = std::atan(std::abs(external_tangent_x->at(itrack))) *
                                 180.0 / TMath::Pi();
      const double theta_y_deg = std::atan(std::abs(external_tangent_y->at(itrack))) *
                                 180.0 / TMath::Pi();
      const int angle_bin_x = FindAngleBin(theta_x_deg);
      const int angle_bin_y = FindAngleBin(theta_y_deg);

      const bool common_x = HasIndex(pm_only_expected_x, itrack) &&
                            HasIndex(dwg_only_expected_x, itrack) &&
                            HasIndex(external_expected_x, itrack) &&
                            HasIndex(frost_nearest_x, itrack) &&
                            IsValidValue(pm_only_expected_x->at(itrack)) &&
                            IsValidValue(dwg_only_expected_x->at(itrack)) &&
                            IsValidValue(external_expected_x->at(itrack)) &&
                            IsValidValue(frost_nearest_x->at(itrack));

      const bool common_y = HasIndex(pm_only_expected_y, itrack) &&
                            HasIndex(dwg_only_expected_y, itrack) &&
                            HasIndex(external_expected_y, itrack) &&
                            HasIndex(frost_nearest_y, itrack) &&
                            IsValidValue(pm_only_expected_y->at(itrack)) &&
                            IsValidValue(dwg_only_expected_y->at(itrack)) &&
                            IsValidValue(external_expected_y->at(itrack)) &&
                            IsValidValue(frost_nearest_y->at(itrack));

      if (IsData(sample_type)) {
        if (common_x) {
          DataResidualEvent event;
          event.split = pm_only_expected_x->at(itrack) - dwg_only_expected_x->at(itrack);
          event.o = frost_nearest_x->at(itrack) - external_expected_x->at(itrack);
          event.weight = weight;
          event.has_split = true;
          event.has_rec_ext = true;
          AddDataEvent(samples.data_x, angle_bin_x, event);
        }
        if (common_y) {
          DataResidualEvent event;
          event.split = pm_only_expected_y->at(itrack) - dwg_only_expected_y->at(itrack);
          event.o = frost_nearest_y->at(itrack) - external_expected_y->at(itrack);
          event.weight = weight;
          event.has_split = true;
          event.has_rec_ext = true;
          AddDataEvent(samples.data_y, angle_bin_y, event);
        }
      } else {
        if (common_x &&
            HasIndex(true_frost_x, itrack) &&
            IsValidValue(true_frost_x->at(itrack))) {
          const double truth = true_frost_x->at(itrack);
          McResidualEvent event;
          event.p = pm_only_expected_x->at(itrack) - truth;
          event.d = dwg_only_expected_x->at(itrack) - truth;
          event.e = external_expected_x->at(itrack) - truth;
          event.r = frost_nearest_x->at(itrack) - truth;
          event.weight = weight;
          AddMcEvent(samples.mc_x, angle_bin_x, event);
        }
        if (common_y &&
            HasIndex(true_frost_y, itrack) &&
            IsValidValue(true_frost_y->at(itrack))) {
          const double truth = true_frost_y->at(itrack);
          McResidualEvent event;
          event.p = pm_only_expected_y->at(itrack) - truth;
          event.d = dwg_only_expected_y->at(itrack) - truth;
          event.e = external_expected_y->at(itrack) - truth;
          event.r = frost_nearest_y->at(itrack) - truth;
          event.weight = weight;
          AddMcEvent(samples.mc_y, angle_bin_y, event);
        }
      }
    }
  }
}

double Difference(double data_value, double mc_value) {
  if (!std::isfinite(data_value) || !std::isfinite(mc_value)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return data_value - mc_value;
}

double RelativeDifference(double data_value, double mc_value) {
  if (!std::isfinite(data_value) ||
      !std::isfinite(mc_value) ||
      mc_value == 0.) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return (data_value - mc_value) / mc_value;
}

std::string FormatDouble(double value, int precision = 6) {
  if (!std::isfinite(value)) {
    return "nan";
  }
  std::ostringstream stream;
  stream << std::fixed << std::setprecision(precision) << value;
  return stream.str();
}

std::string FormatPercent(double fraction, int precision = 1) {
  if (!std::isfinite(fraction)) {
    return "nan";
  }
  std::ostringstream stream;
  stream << std::fixed << std::setprecision(precision) << 100. * fraction << "%";
  return stream.str();
}

std::string SanitizeName(std::string name) {
  for (char &c : name) {
    const bool ok = (c >= 'a' && c <= 'z') ||
                    (c >= 'A' && c <= 'Z') ||
                    (c >= '0' && c <= '9');
    if (!ok) {
      c = '_';
    }
  }
  return name;
}

struct AnalysisResult {
  Axis axis = Axis::kX;
  int analysis_bin = 0;
  bool valid = false;
  Long64_t n_mc = 0;
  Long64_t n_data = 0;
  Long64_t n_data_split = 0;
  Long64_t n_data_rec_ext = 0;
  ResidualSummary split_data;
  ResidualSummary split_mc;
  ResidualSummary rec_ext_data;
  ResidualSummary rec_ext_mc;
  ResidualSummary rec_true_mc;
  ResidualSummary ext_true_mc;
  Chi2Result split_chi2;
  Chi2Result rec_ext_chi2;
  std::string warning;
};

AnalysisResult EvaluateAnalysisResult(Axis axis,
                                      int analysis_bin,
                                      const std::vector<McResidualEvent> &mc_events,
                                      const std::vector<DataResidualEvent> &data_events,
                                      int nbins,
                                      double residual_xmin,
                                      double residual_xmax,
                                      double split_xmin,
                                      double split_xmax) {
  AnalysisResult result;
  result.axis = axis;
  result.analysis_bin = analysis_bin;
  result.n_mc = static_cast<Long64_t>(mc_events.size());
  result.n_data = static_cast<Long64_t>(data_events.size());

  const auto mc_split = MakeMcSplitValues(mc_events);
  const auto mc_rec_ext = MakeMcRecExtValues(mc_events);
  const auto mc_rec_true = MakeMcRecTrueValues(mc_events);
  const auto mc_ext_true = MakeMcExtTrueValues(mc_events);
  const auto data_split = MakeDataSplitValues(data_events);
  const auto data_rec_ext = MakeDataRecExtValues(data_events);

  result.n_data_split = static_cast<Long64_t>(data_split.size());
  result.n_data_rec_ext = static_cast<Long64_t>(data_rec_ext.size());

  result.split_mc = CalculateSummary(mc_split, split_xmin, split_xmax, true);
  result.rec_ext_mc = CalculateSummary(mc_rec_ext, residual_xmin, residual_xmax, true);
  result.rec_true_mc = CalculateSummary(mc_rec_true, kRecTrueXMin, kRecTrueXMax, true);
  result.ext_true_mc = CalculateSummary(mc_ext_true, residual_xmin, residual_xmax, true);
  result.split_data = CalculateSummary(data_split, split_xmin, split_xmax, true);
  result.rec_ext_data = CalculateSummary(data_rec_ext, residual_xmin, residual_xmax, true);

  result.split_chi2 = CalculateChi2Q01Q99(data_split, mc_split, nbins);
  result.rec_ext_chi2 = CalculateChi2Q01Q99(data_rec_ext, mc_rec_ext, nbins);

  if (mc_events.empty()) {
    result.warning = "empty MC sample";
  } else if (data_events.empty()) {
    result.warning = "empty data sample";
  } else if (!result.split_mc.valid || !result.rec_ext_mc.valid ||
             !result.rec_true_mc.valid || !result.ext_true_mc.valid ||
             !result.split_data.valid || !result.rec_ext_data.valid) {
    result.warning = "invalid summary";
  } else {
    result.valid = true;
  }

  return result;
}

const AnalysisResult *FindResult(const std::vector<AnalysisResult> &results,
                                 Axis axis,
                                 int analysis_bin) {
  for (const auto &result : results) {
    if (result.axis == axis && result.analysis_bin == analysis_bin) {
      return &result;
    }
  }
  return nullptr;
}

std::string CsvHeader() {
  return "axis,angle_bin,valid,n_mc,n_data,n_data_split,n_data_rec_ext,"
         "split_data_sigma68_mm,split_mc_sigma68_mm,split_data_minus_mc_mm,"
         "split_relative_diff,split_relative_diff_percent,"
         "rec_ext_data_sigma68_mm,rec_ext_mc_sigma68_mm,rec_ext_data_minus_mc_mm,"
         "rec_ext_relative_diff,rec_ext_relative_diff_percent,"
         "rec_true_mc_sigma68_mm,rec_true_mc_median_mm,rec_true_mc_mean_mm,"
         "rec_true_mc_q16_mm,rec_true_mc_q84_mm,"
         "ext_true_mc_sigma68_mm,ext_true_mc_median_mm,ext_true_mc_mean_mm,"
         "ext_true_mc_q16_mm,ext_true_mc_q84_mm,"
         "split_chi2_q01_mm,split_chi2_q99_mm,split_chi2,split_ndf,split_chi2_ndf,split_p_value,"
         "rec_ext_chi2_q01_mm,rec_ext_chi2_q99_mm,rec_ext_chi2,rec_ext_ndf,rec_ext_chi2_ndf,rec_ext_p_value,"
         "warning";
}

std::string CsvRow(const AnalysisResult &result) {
  const double split_rel = RelativeDifference(result.split_data.sigma68,
                                              result.split_mc.sigma68);
  const double rec_ext_rel = RelativeDifference(result.rec_ext_data.sigma68,
                                                result.rec_ext_mc.sigma68);

  std::ostringstream out;
  out << AxisLabel(result.axis) << ','
      << AnalysisBinLabel(result.analysis_bin) << ','
      << (result.valid ? 1 : 0) << ','
      << result.n_mc << ','
      << result.n_data << ','
      << result.n_data_split << ','
      << result.n_data_rec_ext << ','
      << FormatDouble(result.split_data.sigma68) << ','
      << FormatDouble(result.split_mc.sigma68) << ','
      << FormatDouble(Difference(result.split_data.sigma68, result.split_mc.sigma68)) << ','
      << FormatDouble(split_rel) << ','
      << FormatDouble(100. * split_rel) << ','
      << FormatDouble(result.rec_ext_data.sigma68) << ','
      << FormatDouble(result.rec_ext_mc.sigma68) << ','
      << FormatDouble(Difference(result.rec_ext_data.sigma68, result.rec_ext_mc.sigma68)) << ','
      << FormatDouble(rec_ext_rel) << ','
      << FormatDouble(100. * rec_ext_rel) << ','
      << FormatDouble(result.rec_true_mc.sigma68) << ','
      << FormatDouble(result.rec_true_mc.median) << ','
      << FormatDouble(result.rec_true_mc.mean) << ','
      << FormatDouble(result.rec_true_mc.q16) << ','
      << FormatDouble(result.rec_true_mc.q84) << ','
      << FormatDouble(result.ext_true_mc.sigma68) << ','
      << FormatDouble(result.ext_true_mc.median) << ','
      << FormatDouble(result.ext_true_mc.mean) << ','
      << FormatDouble(result.ext_true_mc.q16) << ','
      << FormatDouble(result.ext_true_mc.q84) << ','
      << FormatDouble(result.split_chi2.q01) << ','
      << FormatDouble(result.split_chi2.q99) << ','
      << FormatDouble(result.split_chi2.chi2) << ','
      << result.split_chi2.ndf << ','
      << FormatDouble(result.split_chi2.chi2_ndf) << ','
      << FormatDouble(result.split_chi2.p_value) << ','
      << FormatDouble(result.rec_ext_chi2.q01) << ','
      << FormatDouble(result.rec_ext_chi2.q99) << ','
      << FormatDouble(result.rec_ext_chi2.chi2) << ','
      << result.rec_ext_chi2.ndf << ','
      << FormatDouble(result.rec_ext_chi2.chi2_ndf) << ','
      << FormatDouble(result.rec_ext_chi2.p_value) << ','
      << result.warning;
  return out.str();
}

TString MakeOutputLogFileName(const char *output_pdf,
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

TString MakeOutputRootLogFileName(const char *output_pdf,
                                  const char *output_log) {
  TString root_file = MakeOutputLogFileName(output_pdf, output_log);
  if (root_file.EndsWith(".log") ||
      root_file.EndsWith(".csv") ||
      root_file.EndsWith(".txt")) {
    root_file.Remove(root_file.Length() - 4);
  }
  root_file += ".root";
  return root_file;
}

void WriteCsvLog(const char *output_log,
                 const std::vector<AnalysisResult> &results) {
  std::ofstream log(output_log);
  if (!log) {
    std::cerr << "Warning: cannot open log file: " << output_log << std::endl;
    return;
  }

  log << CsvHeader() << '\n';
  for (const auto &result : results) {
    log << CsvRow(result) << '\n';
  }
}

void WriteRootLogTree(TDirectory *directory,
                      const std::vector<AnalysisResult> &results) {
  if (!directory) {
    return;
  }

  directory->cd();
  TTree tree("reslog", "FROST position resolution validation log");

  std::string axis;
  std::string angle_bin;
  int valid = 0;
  Long64_t n_mc = 0;
  Long64_t n_data = 0;
  Long64_t n_data_split = 0;
  Long64_t n_data_rec_ext = 0;
  double split_data_sigma68_mm = 0.;
  double split_mc_sigma68_mm = 0.;
  double split_data_minus_mc_mm = 0.;
  double split_relative_diff = 0.;
  double split_relative_diff_percent = 0.;
  double rec_ext_data_sigma68_mm = 0.;
  double rec_ext_mc_sigma68_mm = 0.;
  double rec_ext_data_minus_mc_mm = 0.;
  double rec_ext_relative_diff = 0.;
  double rec_ext_relative_diff_percent = 0.;
  double rec_true_mc_sigma68_mm = 0.;
  double rec_true_mc_median_mm = 0.;
  double rec_true_mc_mean_mm = 0.;
  double rec_true_mc_q16_mm = 0.;
  double rec_true_mc_q84_mm = 0.;
  double ext_true_mc_sigma68_mm = 0.;
  double ext_true_mc_median_mm = 0.;
  double ext_true_mc_mean_mm = 0.;
  double ext_true_mc_q16_mm = 0.;
  double ext_true_mc_q84_mm = 0.;
  double split_chi2_q01_mm = 0.;
  double split_chi2_q99_mm = 0.;
  double split_chi2 = 0.;
  int split_ndf = 0;
  double split_chi2_ndf = 0.;
  double split_p_value = 0.;
  double rec_ext_chi2_q01_mm = 0.;
  double rec_ext_chi2_q99_mm = 0.;
  double rec_ext_chi2 = 0.;
  int rec_ext_ndf = 0;
  double rec_ext_chi2_ndf = 0.;
  double rec_ext_p_value = 0.;
  std::string warning;

  tree.Branch("axis", &axis);
  tree.Branch("angle_bin", &angle_bin);
  tree.Branch("valid", &valid);
  tree.Branch("n_mc", &n_mc);
  tree.Branch("n_data", &n_data);
  tree.Branch("n_data_split", &n_data_split);
  tree.Branch("n_data_rec_ext", &n_data_rec_ext);
  tree.Branch("split_data_sigma68_mm", &split_data_sigma68_mm);
  tree.Branch("split_mc_sigma68_mm", &split_mc_sigma68_mm);
  tree.Branch("split_data_minus_mc_mm", &split_data_minus_mc_mm);
  tree.Branch("split_relative_diff", &split_relative_diff);
  tree.Branch("split_relative_diff_percent", &split_relative_diff_percent);
  tree.Branch("rec_ext_data_sigma68_mm", &rec_ext_data_sigma68_mm);
  tree.Branch("rec_ext_mc_sigma68_mm", &rec_ext_mc_sigma68_mm);
  tree.Branch("rec_ext_data_minus_mc_mm", &rec_ext_data_minus_mc_mm);
  tree.Branch("rec_ext_relative_diff", &rec_ext_relative_diff);
  tree.Branch("rec_ext_relative_diff_percent", &rec_ext_relative_diff_percent);
  tree.Branch("rec_true_mc_sigma68_mm", &rec_true_mc_sigma68_mm);
  tree.Branch("rec_true_mc_median_mm", &rec_true_mc_median_mm);
  tree.Branch("rec_true_mc_mean_mm", &rec_true_mc_mean_mm);
  tree.Branch("rec_true_mc_q16_mm", &rec_true_mc_q16_mm);
  tree.Branch("rec_true_mc_q84_mm", &rec_true_mc_q84_mm);
  tree.Branch("ext_true_mc_sigma68_mm", &ext_true_mc_sigma68_mm);
  tree.Branch("ext_true_mc_median_mm", &ext_true_mc_median_mm);
  tree.Branch("ext_true_mc_mean_mm", &ext_true_mc_mean_mm);
  tree.Branch("ext_true_mc_q16_mm", &ext_true_mc_q16_mm);
  tree.Branch("ext_true_mc_q84_mm", &ext_true_mc_q84_mm);
  tree.Branch("split_chi2_q01_mm", &split_chi2_q01_mm);
  tree.Branch("split_chi2_q99_mm", &split_chi2_q99_mm);
  tree.Branch("split_chi2", &split_chi2);
  tree.Branch("split_ndf", &split_ndf);
  tree.Branch("split_chi2_ndf", &split_chi2_ndf);
  tree.Branch("split_p_value", &split_p_value);
  tree.Branch("rec_ext_chi2_q01_mm", &rec_ext_chi2_q01_mm);
  tree.Branch("rec_ext_chi2_q99_mm", &rec_ext_chi2_q99_mm);
  tree.Branch("rec_ext_chi2", &rec_ext_chi2);
  tree.Branch("rec_ext_ndf", &rec_ext_ndf);
  tree.Branch("rec_ext_chi2_ndf", &rec_ext_chi2_ndf);
  tree.Branch("rec_ext_p_value", &rec_ext_p_value);
  tree.Branch("warning", &warning);

  for (const auto &result : results) {
    const double split_rel = RelativeDifference(result.split_data.sigma68,
                                                result.split_mc.sigma68);
    const double rec_ext_rel = RelativeDifference(result.rec_ext_data.sigma68,
                                                  result.rec_ext_mc.sigma68);

    axis = AxisLabel(result.axis);
    angle_bin = AnalysisBinLabel(result.analysis_bin);
    valid = result.valid ? 1 : 0;
    n_mc = result.n_mc;
    n_data = result.n_data;
    n_data_split = result.n_data_split;
    n_data_rec_ext = result.n_data_rec_ext;
    split_data_sigma68_mm = result.split_data.sigma68;
    split_mc_sigma68_mm = result.split_mc.sigma68;
    split_data_minus_mc_mm = Difference(result.split_data.sigma68, result.split_mc.sigma68);
    split_relative_diff = split_rel;
    split_relative_diff_percent = 100. * split_rel;
    rec_ext_data_sigma68_mm = result.rec_ext_data.sigma68;
    rec_ext_mc_sigma68_mm = result.rec_ext_mc.sigma68;
    rec_ext_data_minus_mc_mm = Difference(result.rec_ext_data.sigma68, result.rec_ext_mc.sigma68);
    rec_ext_relative_diff = rec_ext_rel;
    rec_ext_relative_diff_percent = 100. * rec_ext_rel;
    rec_true_mc_sigma68_mm = result.rec_true_mc.sigma68;
    rec_true_mc_median_mm = result.rec_true_mc.median;
    rec_true_mc_mean_mm = result.rec_true_mc.mean;
    rec_true_mc_q16_mm = result.rec_true_mc.q16;
    rec_true_mc_q84_mm = result.rec_true_mc.q84;
    ext_true_mc_sigma68_mm = result.ext_true_mc.sigma68;
    ext_true_mc_median_mm = result.ext_true_mc.median;
    ext_true_mc_mean_mm = result.ext_true_mc.mean;
    ext_true_mc_q16_mm = result.ext_true_mc.q16;
    ext_true_mc_q84_mm = result.ext_true_mc.q84;
    split_chi2_q01_mm = result.split_chi2.q01;
    split_chi2_q99_mm = result.split_chi2.q99;
    split_chi2 = result.split_chi2.chi2;
    split_ndf = result.split_chi2.ndf;
    split_chi2_ndf = result.split_chi2.chi2_ndf;
    split_p_value = result.split_chi2.p_value;
    rec_ext_chi2_q01_mm = result.rec_ext_chi2.q01;
    rec_ext_chi2_q99_mm = result.rec_ext_chi2.q99;
    rec_ext_chi2 = result.rec_ext_chi2.chi2;
    rec_ext_ndf = result.rec_ext_chi2.ndf;
    rec_ext_chi2_ndf = result.rec_ext_chi2.chi2_ndf;
    rec_ext_p_value = result.rec_ext_chi2.p_value;
    warning = result.warning;
    tree.Fill();
  }

  tree.Write("", TObject::kOverwrite);
}

struct PlotSample {
  std::string label;
  std::vector<WeightedValue> values;
  int color = kBlack;
  int line_style = 1;
  int line_width = 1;
  bool draw_as_points = false;
  int marker_style = 20;
  double marker_size = 0.8;
};

double GetHistogramMaximumWithError(const TH1D *hist) {
  if (!hist) {
    return 0.;
  }
  double maximum = 0.;
  for (int ibin = 1; ibin <= hist->GetNbinsX(); ++ibin) {
    maximum = std::max(maximum, hist->GetBinContent(ibin) + hist->GetBinError(ibin));
  }
  return maximum;
}

void SaveHistogram(TDirectory *hist_dir, TH1D *hist) {
  if (!hist_dir || !hist) {
    return;
  }
  hist_dir->cd();
  hist->Write(hist->GetName(), TObject::kOverwrite);
}

void DrawOverlayHistograms(const std::vector<PlotSample> &samples,
                           const std::vector<std::string> &extra_lines,
                           const char *hist_name_prefix,
                           const char *title,
                           const char *x_axis_title,
                           int nbins,
                           double xmin,
                           double xmax,
                           TDirectory *hist_dir) {
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.13);

  double data_integral = 0.;
  std::vector<TH1D*> histograms;
  histograms.reserve(samples.size());

  const std::string prefix = SanitizeName(hist_name_prefix);
  for (std::size_t isample = 0; isample < samples.size(); ++isample) {
    const auto &sample = samples.at(isample);
    const TString name = Form("%s_%s", prefix.c_str(), SanitizeName(sample.label).c_str());
    auto *hist = new TH1D(name,
                          Form("%s;%s;Number of tracks", title, x_axis_title),
                          nbins,
                          xmin,
                          xmax);
    hist->SetDirectory(nullptr);
    hist->Sumw2();
    hist->SetLineColor(sample.color);
    hist->SetLineStyle(sample.line_style);
    hist->SetLineWidth(sample.line_width);
    hist->SetMarkerColor(sample.color);
    hist->SetMarkerStyle(sample.marker_style);
    hist->SetMarkerSize(sample.marker_size);

    for (const auto &entry : sample.values) {
      if (!std::isfinite(entry.value) ||
          !std::isfinite(entry.weight) ||
          entry.weight <= 0.) {
        continue;
      }
      hist->Fill(entry.value, entry.weight);
    }

    if (sample.draw_as_points && data_integral <= 0.) {
      data_integral = hist->Integral();
    }
    histograms.push_back(hist);
  }

  // Scale MC-like histograms to the Data integral in Data/MC comparison plots.
  if (data_integral > 0.) {
    for (std::size_t i = 0; i < samples.size(); ++i) {
      if (samples.at(i).draw_as_points) {
        continue;
      }
      const double integral = histograms.at(i)->Integral();
      if (integral > 0.) {
        histograms.at(i)->Scale(data_integral / integral);
      }
    }
  }

  double max_y = 0.;
  for (auto *hist : histograms) {
    max_y = std::max(max_y, GetHistogramMaximumWithError(hist));
  }
  for (auto *hist : histograms) {
    hist->SetMaximum(max_y > 0. ? 1.35 * max_y : 1.0);
    hist->GetYaxis()->SetTitleOffset(1.55);
  }

  bool first = true;
  for (std::size_t i = 0; i < samples.size(); ++i) {
    if (samples.at(i).draw_as_points) {
      continue;
    }
    if (first) {
      histograms.at(i)->Draw("hist");
      first = false;
    } else {
      histograms.at(i)->Draw("hist same");
    }
  }
  for (std::size_t i = 0; i < samples.size(); ++i) {
    if (!samples.at(i).draw_as_points) {
      continue;
    }
    if (first) {
      histograms.at(i)->Draw("E1");
      first = false;
    } else {
      histograms.at(i)->Draw("E1 same");
    }
  }

  auto *legend = new TLegend(0.58, 0.72, 0.94, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.033);
  for (std::size_t i = 0; i < histograms.size(); ++i) {
    legend->AddEntry(histograms.at(i),
                     samples.at(i).label.c_str(),
                     samples.at(i).draw_as_points ? "lep" : "l");
  }
  legend->Draw("same");

  if (!extra_lines.empty()) {
    auto *info = new TPaveText(0.18, 0.69, 0.56, 0.88, "NDC");
    info->SetFillColor(0);
    info->SetFillStyle(0);
    info->SetBorderSize(0);
    info->SetTextAlign(12);
    info->SetTextFont(42);
    info->SetTextSize(0.030);
    for (const auto &line : extra_lines) {
      info->AddText(line.c_str());
    }
    info->Draw("same");
  }

  for (auto *hist : histograms) {
    SaveHistogram(hist_dir, hist);
  }
}

std::vector<std::string> MakeDataMcInfoLines(const ResidualSummary &data,
                                             const ResidualSummary &mc,
                                             const Chi2Result &chi2) {
  const double rel = RelativeDifference(data.sigma68, mc.sigma68);
  std::vector<std::string> lines;
  lines.push_back(Form("Data #sigma_{68}=%.3f mm", data.sigma68));
  lines.push_back(Form("MC #sigma_{68}=%.3f mm", mc.sigma68));
  lines.push_back(Form("Data/MC-1 = %.1f%%", 100. * rel));
  if (chi2.valid) {
    lines.push_back(Form("#chi^{2}/ndf = %.1f/%d", chi2.chi2, chi2.ndf));
  } else {
    lines.push_back("#chi^{2}/ndf = n/a");
  }
  return lines;
}

std::vector<std::string> MakeRecTrueInfoLines(const ResidualSummary &summary) {
  std::vector<std::string> lines;
  lines.push_back(Form("MC median = %.3f mm", summary.median));
  lines.push_back(Form("MC #sigma_{68} = %.3f mm", summary.sigma68));
  return lines;
}

void DrawOnePad(const AnalysisResult &result,
                const std::vector<McResidualEvent> &mc_events,
                const std::vector<DataResidualEvent> &data_events,
                PlotKind kind,
                int nbins,
                double residual_xmin,
                double residual_xmax,
                double split_xmin,
                double split_xmax,
                TDirectory *hist_dir) {
  const std::string axis = AxisLabel(result.axis);
  const std::string bin_title = AnalysisBinTitle(result.analysis_bin, result.axis);
  const std::string bin_label = AnalysisBinLabel(result.analysis_bin);

  if (kind == PlotKind::kSplit) {
    std::vector<PlotSample> samples = {
      {"MC", MakeMcSplitValues(mc_events), kBlue + 1, 1, 1, false, 20, 0.8},
      {"Data", MakeDataSplitValues(data_events), kBlack, 1, 1, true, 20, 0.8}
    };
    DrawOverlayHistograms(samples,
                          MakeDataMcInfoLines(result.split_data, result.split_mc, result.split_chi2),
                          Form("split_%s_%s", axis.c_str(), bin_label.c_str()),
                          Form("%s", bin_title.c_str()),
                          Form("%s_{PM-only} - %s_{DWG-only} [mm]", axis.c_str(), axis.c_str()),
                          nbins,
                          split_xmin,
                          split_xmax,
                          hist_dir);
    return;
  }

  if (kind == PlotKind::kRecExt) {
    std::vector<PlotSample> samples = {
      {"MC", MakeMcRecExtValues(mc_events), kBlue + 1, 1, 1, false, 20, 0.8},
      {"Data", MakeDataRecExtValues(data_events), kBlack, 1, 1, true, 20, 0.8}
    };
    DrawOverlayHistograms(samples,
                          MakeDataMcInfoLines(result.rec_ext_data, result.rec_ext_mc, result.rec_ext_chi2),
                          Form("rec_ext_%s_%s", axis.c_str(), bin_label.c_str()),
                          Form("%s", bin_title.c_str()),
                          Form("%s_{rec} - %s_{ext} [mm]", axis.c_str(), axis.c_str()),
                          nbins,
                          residual_xmin,
                          residual_xmax,
                          hist_dir);
    return;
  }

  if (kind == PlotKind::kRecTrue) {
    std::vector<PlotSample> samples = {
      {"MC", MakeMcRecTrueValues(mc_events), kBlue + 1, 1, 2, false, 20, 0.8}
    };
    DrawOverlayHistograms(samples,
                          MakeRecTrueInfoLines(result.rec_true_mc),
                          Form("rec_true_%s_%s", axis.c_str(), bin_label.c_str()),
                          Form("%s", bin_title.c_str()),
                          Form("%s_{rec} - %s_{true} [mm]", axis.c_str(), axis.c_str()),
                          nbins,
                          kRecTrueXMin,
                          kRecTrueXMax,
                          hist_dir);
    return;
  }

  if (kind == PlotKind::kExtTrue) {
    std::vector<PlotSample> samples = {
      {"MC", MakeMcExtTrueValues(mc_events), kBlue + 1, 1, 2, false, 20, 0.8}
    };
    DrawOverlayHistograms(samples,
                          MakeRecTrueInfoLines(result.ext_true_mc),
                          Form("ext_true_%s_%s", axis.c_str(), bin_label.c_str()),
                          Form("%s", bin_title.c_str()),
                          Form("%s_{ext} - %s_{true} [mm]", axis.c_str(), axis.c_str()),
                          nbins,
                          residual_xmin,
                          residual_xmax,
                          hist_dir);
    return;
  }
}

void DrawPlotPage(TCanvas &canvas,
                  const std::vector<AnalysisResult> &results,
                  const ResidualSamples &samples,
                  PlotKind kind,
                  int analysis_bin,
                  const char *output_pdf,
                  int nbins,
                  double residual_xmin,
                  double residual_xmax,
                  double split_xmin,
                  double split_xmax,
                  TDirectory *hist_dir) {
  canvas.Clear();
  canvas.Divide(2, 1);

  const AnalysisResult *x_result = FindResult(results, Axis::kX, analysis_bin);
  const AnalysisResult *y_result = FindResult(results, Axis::kY, analysis_bin);

  canvas.cd(1);
  if (x_result) {
    DrawOnePad(*x_result,
               samples.mc_x.at(analysis_bin),
               samples.data_x.at(analysis_bin),
               kind,
               nbins,
               residual_xmin,
               residual_xmax,
               split_xmin,
               split_xmax,
               hist_dir);
  }

  canvas.cd(2);
  if (y_result) {
    DrawOnePad(*y_result,
               samples.mc_y.at(analysis_bin),
               samples.data_y.at(analysis_bin),
               kind,
               nbins,
               residual_xmin,
               residual_xmax,
               split_xmin,
               split_xmax,
               hist_dir);
  }

  canvas.Print(output_pdf);
}

void DrawCoverPage(TCanvas &canvas,
                   const char *output_pdf,
                   const std::vector<std::string> &mc_files,
                   const std::vector<std::string> &data_files,
                   int nbins,
                   double residual_xmin,
                   double residual_xmax,
                   bool use_position_reweighting) {
  canvas.Clear();
  auto *text = new TPaveText(0.06, 0.08, 0.94, 0.92, "NDC");
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->SetTextFont(42);
  text->SetTextSize(0.030);
  text->AddText("FROST position-resolution validation plots");
  text->AddText("Plots are ordered as: PM-only - DWG-only, rec - ext, MC rec - true, and MC ext - true.");
  text->AddText("Each page shows x on the left and y on the right.");
  text->AddText("Median and sigma68 are calculated from weighted residual values inside the displayed histogram range.");
  text->AddText("Data/MC chi2/ndf values use area-normalized histograms in a common weighted q01-q99 range.");
  text->AddText("");
  text->AddText(Form("MC files: %lu", static_cast<unsigned long>(mc_files.size())));
  text->AddText(Form("Data files: %lu", static_cast<unsigned long>(data_files.size())));
  text->AddText(Form("nbins=%d, residual range=%.1f to %.1f mm", nbins, residual_xmin, residual_xmax));
  text->AddText(Form("split range=%.1f to %.1f mm", 2.0 * residual_xmin, 2.0 * residual_xmax));
  text->AddText(Form("rec-true range=%.1f to %.1f mm", kRecTrueXMin, kRecTrueXMax));
  text->AddText(Form("ext-true range=%.1f to %.1f mm", residual_xmin, residual_xmax));
  text->AddText(Form("position reweighting: %s", use_position_reweighting ? "enabled" : "disabled"));
  text->Draw();
  canvas.Print(output_pdf);
}

void DrawValidationTablePage(TCanvas &canvas,
                             const std::vector<AnalysisResult> &results,
                             const char *output_pdf) {
  canvas.Clear();
  auto *text = new TPaveText(0.03, 0.04, 0.97, 0.96, "NDC");
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->SetTextFont(82);
  text->SetTextSize(kTableTextSize);

  text->AddText("Data/MC validation summary");
  text->AddText("rel. diff. = (Data - MC) / MC.  sigma68 values are in mm.");
  text->AddText("");
  text->AddText("axis bin       nMC     nData  splitD splitMC dSplit   recExtD recExtMC dRecExt  recTrueMC extTrueMC  split chi2/ndf  recExt chi2/ndf");
  text->AddText("-----------------------------------------------------------------------------------------------------------------------------------");

  for (const auto &result : results) {
    const std::string split_chi2 = result.split_chi2.valid
      ? Form("%.1f/%d", result.split_chi2.chi2, result.split_chi2.ndf)
      : "nan";
    const std::string rec_ext_chi2 = result.rec_ext_chi2.valid
      ? Form("%.1f/%d", result.rec_ext_chi2.chi2, result.rec_ext_chi2.ndf)
      : "nan";

    const std::string line = Form("%-4s %-9s %7lld %7lld  %6.2f %6.2f %7s  %7.2f %7.2f %7s  %9.2f %9.2f  %14s %15s",
                                  AxisLabel(result.axis),
                                  AnalysisBinLabel(result.analysis_bin).c_str(),
                                  static_cast<long long>(result.n_mc),
                                  static_cast<long long>(result.n_data),
                                  result.split_data.sigma68,
                                  result.split_mc.sigma68,
                                  FormatPercent(RelativeDifference(result.split_data.sigma68,
                                                                    result.split_mc.sigma68), 1).c_str(),
                                  result.rec_ext_data.sigma68,
                                  result.rec_ext_mc.sigma68,
                                  FormatPercent(RelativeDifference(result.rec_ext_data.sigma68,
                                                                    result.rec_ext_mc.sigma68), 1).c_str(),
                                  result.rec_true_mc.sigma68,
                                  result.ext_true_mc.sigma68,
                                  split_chi2.c_str(),
                                  rec_ext_chi2.c_str());
    text->AddText(line.c_str());
  }

  text->Draw();
  canvas.Print(output_pdf);
}

void DrawRecTrueTablePage(TCanvas &canvas,
                          const std::vector<AnalysisResult> &results,
                          const char *output_pdf,
                          bool close_pdf) {
  canvas.Clear();
  TString pdf_name = output_pdf;
  if (close_pdf) {
    pdf_name += ")";
  }

  auto *text = new TPaveText(0.05, 0.06, 0.95, 0.94, "NDC");
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->SetTextFont(82);
  text->SetTextSize(0.019);

  text->AddText("MC truth summary: rec - true and ext - true");
  text->AddText("rec-true uses the plotted range, -10 mm to 10 mm; ext-true uses the residual range.");
  text->AddText("");
  text->AddText("axis bin       nMC     recTrue_sigma68 recTrue_median  extTrue_sigma68 extTrue_median");
  text->AddText("--------------------------------------------------------------------------------");

  for (const auto &result : results) {
    const std::string line = Form("%-4s %-9s %7lld  %15.3f %14.3f  %15.3f %14.3f",
                                  AxisLabel(result.axis),
                                  AnalysisBinLabel(result.analysis_bin).c_str(),
                                  static_cast<long long>(result.n_mc),
                                  result.rec_true_mc.sigma68,
                                  result.rec_true_mc.median,
                                  result.ext_true_mc.sigma68,
                                  result.ext_true_mc.median);
    text->AddText(line.c_str());
  }

  text->Draw();
  canvas.Print(pdf_name.Data());
}

} // namespace

void CalculatePositionResolution(
    const char *mc_input_dirs="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC_retry/6-TrackMatch_externalfit_PMandDWG,/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/6-TrackMatch_externalfit_PMandDWG",
    // const char *data_input_dirs="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM_f8b4eb/2-rootfile_after_TrackMatch_externalfit_PMandDWG/FHC,/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM_f8b4eb/2-rootfile_after_TrackMatch_externalfit_PMandDWG/RHC",
    // const char *output_pdf="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM_f8b4eb/2-rootfile_after_TrackMatch_externalfit_PMandDWG/analysisplot/FROST_position_resolution_validation_nhitDWG5PM10_onetrack_FHCRHC.pdf",
    // const char *data_input_dirs="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM_f8b4eb/2-rootfile_after_TrackMatch_externalfit_PMandDWG/FHC",
    // const char *output_pdf="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM_f8b4eb/2-rootfile_after_TrackMatch_externalfit_PMandDWG/analysisplot/FROST_position_resolution_validation_nhitDWG5PM10_onetrack_FHC.pdf",
    const char *data_input_dirs="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM_f8b4eb/2-rootfile_after_TrackMatch_externalfit_PMandDWG/RHC",
    const char *output_pdf="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM_f8b4eb/2-rootfile_after_TrackMatch_externalfit_PMandDWG/analysisplot/FROST_position_resolution_validation_nhitDWG5PM10_onetrack_RHC.pdf",
    const char *output_log="",
    int nbins = 200,
    double residual_xmin_mm = -50.,
    double residual_xmax_mm = 50.,
    bool use_position_reweighting = true) {
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);

  const std::vector<std::string> mc_input_dir_list = ParseInputDirectories(mc_input_dirs);
  const std::vector<std::string> data_input_dir_list = ParseInputDirectories(data_input_dirs);
  if (mc_input_dir_list.empty()) {
    std::cerr << "Error: no MC input directory is specified" << std::endl;
    return;
  }
  if (data_input_dir_list.empty()) {
    std::cerr << "Error: no data input directory is specified" << std::endl;
    return;
  }

  const std::vector<std::string> mc_root_files = FindRootFiles(mc_input_dir_list);
  const std::vector<std::string> data_root_files = FindRootFiles(data_input_dir_list);
  if (mc_root_files.empty()) {
    std::cerr << "Error: no MC .root files found" << std::endl;
    return;
  }
  if (data_root_files.empty()) {
    std::cerr << "Error: no data .root files found" << std::endl;
    return;
  }

  std::cout << "MC input files: " << mc_root_files.size() << std::endl;
  std::cout << "Data input files: " << data_root_files.size() << std::endl;

  PositionBinWeights mc_position_weights;
  PositionBinWeights data_position_weights;
  if (use_position_reweighting) {
    for (const auto &file_path : mc_root_files) {
      CountPositionBinsOneFile(file_path, SampleType::kMonteCarlo, mc_position_weights);
    }
    mc_position_weights.CalculateWeights();
    PrintPositionWeightSummary("MC", mc_position_weights);

    for (const auto &file_path : data_root_files) {
      CountPositionBinsOneFile(file_path, SampleType::kData, data_position_weights);
    }
    data_position_weights.CalculateWeights();
    PrintPositionWeightSummary("Data", data_position_weights);
  } else {
    mc_position_weights.Disable();
    data_position_weights.Disable();
    std::cout << "Position reweighting: disabled" << std::endl;
  }

  ResidualSamples samples;
  for (const auto &file_path : mc_root_files) {
    std::cout << "Loading MC: " << file_path << std::endl;
    LoadResidualSamplesOneFile(file_path, SampleType::kMonteCarlo, mc_position_weights, samples);
  }
  for (const auto &file_path : data_root_files) {
    std::cout << "Loading data: " << file_path << std::endl;
    LoadResidualSamplesOneFile(file_path, SampleType::kData, data_position_weights, samples);
  }

  const double split_xmin_mm = 2.0 * residual_xmin_mm;
  const double split_xmax_mm = 2.0 * residual_xmax_mm;

  std::vector<AnalysisResult> results;
  results.reserve(2 * kNAnalysisBins);
  for (int analysis_bin = 0; analysis_bin < kNAnalysisBins; ++analysis_bin) {
    results.push_back(EvaluateAnalysisResult(Axis::kX,
                                             analysis_bin,
                                             samples.mc_x.at(analysis_bin),
                                             samples.data_x.at(analysis_bin),
                                             nbins,
                                             residual_xmin_mm,
                                             residual_xmax_mm,
                                             split_xmin_mm,
                                             split_xmax_mm));
  }
  for (int analysis_bin = 0; analysis_bin < kNAnalysisBins; ++analysis_bin) {
    results.push_back(EvaluateAnalysisResult(Axis::kY,
                                             analysis_bin,
                                             samples.mc_y.at(analysis_bin),
                                             samples.data_y.at(analysis_bin),
                                             nbins,
                                             residual_xmin_mm,
                                             residual_xmax_mm,
                                             split_xmin_mm,
                                             split_xmax_mm));
  }

  const TString output_log_file = MakeOutputLogFileName(output_pdf, output_log);
  WriteCsvLog(output_log_file.Data(), results);
  std::cout << "Wrote log: " << output_log_file.Data() << std::endl;

  const TString output_root_log_file = MakeOutputRootLogFileName(output_pdf, output_log);
  TFile output_root(output_root_log_file.Data(), "RECREATE");
  if (output_root.IsZombie()) {
    std::cerr << "Warning: cannot create ROOT output file: "
              << output_root_log_file.Data() << std::endl;
    return;
  }
  WriteRootLogTree(&output_root, results);
  TDirectory *hist_dir = output_root.mkdir("histograms");

  TCanvas canvas("canvas", "FROST position-resolution validation", 1400, 650);

  TString pdf_open = output_pdf;
  pdf_open += "(";
  DrawCoverPage(canvas,
                pdf_open.Data(),
                mc_root_files,
                data_root_files,
                nbins,
                residual_xmin_mm,
                residual_xmax_mm,
                use_position_reweighting);

  for (int analysis_bin = 0; analysis_bin < kNAnalysisBins; ++analysis_bin) {
    DrawPlotPage(canvas,
                 results,
                 samples,
                 PlotKind::kSplit,
                 analysis_bin,
                 output_pdf,
                 nbins,
                 residual_xmin_mm,
                 residual_xmax_mm,
                 split_xmin_mm,
                 split_xmax_mm,
                 hist_dir);
  }

  for (int analysis_bin = 0; analysis_bin < kNAnalysisBins; ++analysis_bin) {
    DrawPlotPage(canvas,
                 results,
                 samples,
                 PlotKind::kRecExt,
                 analysis_bin,
                 output_pdf,
                 nbins,
                 residual_xmin_mm,
                 residual_xmax_mm,
                 split_xmin_mm,
                 split_xmax_mm,
                 hist_dir);
  }

  for (int analysis_bin = 0; analysis_bin < kNAnalysisBins; ++analysis_bin) {
    DrawPlotPage(canvas,
                 results,
                 samples,
                 PlotKind::kRecTrue,
                 analysis_bin,
                 output_pdf,
                 nbins,
                 residual_xmin_mm,
                 residual_xmax_mm,
                 split_xmin_mm,
                 split_xmax_mm,
                 hist_dir);
  }

  for (int analysis_bin = 0; analysis_bin < kNAnalysisBins; ++analysis_bin) {
    DrawPlotPage(canvas,
                 results,
                 samples,
                 PlotKind::kExtTrue,
                 analysis_bin,
                 output_pdf,
                 nbins,
                 residual_xmin_mm,
                 residual_xmax_mm,
                 split_xmin_mm,
                 split_xmax_mm,
                 hist_dir);
  }

  DrawValidationTablePage(canvas, results, output_pdf);
  DrawRecTrueTablePage(canvas, results, output_pdf, true);

  output_root.cd();
  output_root.Write("", TObject::kOverwrite);
  output_root.Close();

  std::cout << "Wrote ROOT log and histograms: " << output_root_log_file.Data() << std::endl;
  std::cout << "Wrote PDF: " << output_pdf << std::endl;
}
