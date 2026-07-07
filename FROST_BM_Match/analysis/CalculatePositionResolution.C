// CalculatePositionResolution.C
//
// Estimate the nominal FROST position resolution from MC truth and assign
// model systematic uncertainties using residual morphing stress tests.
//
// Usage in ROOT:
//   root -l -b -q 'CalculatePositionResolution.C("mc_dir", "data_dir", "output.pdf")'
//
// Multiple input directories can be given as a comma-separated list:
//   root -l -b -q 'CalculatePositionResolution.C("mc_dir1,mc_dir2", "data_dir1,data_dir2", "output.pdf")'
//
// Output:
//   - One multi-page PDF.
//   - One CSV log file. If output_log is empty, output_pdf with ".log" is used.
//   - One ROOT log file with the same base name as the CSV log file.
//   - The final pages of the PDF contain the compact log table.
//
// Method:
//   For each axis and each angle bin, including the all-angle bin,
//   the nominal resolution is taken from the MC truth residual:
//     R_MC = sigma68(x_rec - x_true).
//
//   Data-observable residuals are used for validation and systematic stress
//   tests:
//     S_data = sigma68(PM-only - DWG-only)
//     O_data = sigma68(rec - ext)
//
//   In MC, define residuals relative to the true FROST position:
//     P = x_PM-only  - x_true
//     D = x_DWG-only - x_true
//     E = x_ext      - x_true
//     R = x_rec      - x_true
//
//   First, k_ext is chosen so that
//     sigma68(P' - D') = S_data,
//   where
//     P' = median(P) + k_ext [P - median(P)]
//     D' = median(D) + k_ext [D - median(D)].
//
//   The propagation of this external discrepancy to the combined external
//   residual is scanned with
//     k_ext_combined(alpha) = 1 + alpha * (k_ext - 1),
//   for alpha = 0.0, 0.1, ..., 1.0.
//
//   For each alpha, k_rec(alpha) is chosen so that
//     sigma68(R'_alpha - E'_alpha) = O_data.
//   The alpha-envelope systematic is
//     delta_alpha = max_alpha |sigma68(R'_alpha) - R_MC|.
//
//   The R-E correlation systematic is evaluated by shuffling R relative to E
//   for each alpha. For a given alpha, the shuffled result is compared with
//   the paired result using the same alpha:
//     delta_corr = max_alpha |R_shuffle(alpha) - R_pair(alpha)|.
//   The final model systematic is
//     delta_model = sqrt(delta_alpha^2 + delta_corr^2).
//
// Notes:
//   - No quadrature subtraction is used.
//   - The residual widths are calculated directly from weighted residual values
//     within the displayed histogram range, not from histogram bin contents.
//   - The overlay plots compare Data with the un-morphed MC prediction.
//   - MC histograms are scaled to the Data integral in each plot.
//   - Data/MC validation chi2/ndf and p-values are calculated in a common
//     weighted q01-q99 range using area-normalized histograms.

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TList.h>
#include <TMath.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TStyle.h>
#include <TText.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int NHIT_DWG_X = 4;
constexpr int NHIT_PM_X = 10;
constexpr int NHIT_DWG_Y = 4;
constexpr int NHIT_PM_Y = 10;

constexpr double ACCEPTANCE_X = 560.0; // mm
constexpr double ACCEPTANCE_Y = 600.0; // mm

constexpr int kNPositionBinsX = 10;
constexpr int kNPositionBinsY = 10;
constexpr int kNPositionBins = kNPositionBinsX * kNPositionBinsY;

constexpr double kNonInitializedThreshold = -9.0e7;

// Angle bins [deg]:
//   [0,10), [10,20), [20,50)
constexpr int kNAngleBins = 3;
const double kAngleBins[kNAngleBins + 1] = {
  0.0, 10.0, 20.0, 50.0
};

constexpr int kNAnalysisBins = kNAngleBins + 1; // all + angle bins

enum class SampleType {
  kMonteCarlo,
  kData
};

enum class Axis {
  kX,
  kY
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
  return Form("%.0f_%.0f_deg",
              kAngleBins[angle_bin],
              kAngleBins[angle_bin + 1]);
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

std::vector<std::string> FindRootFiles(const std::vector<std::string> &input_dirs) {
  std::vector<std::string> root_files;

  for (const auto &input_dir : input_dirs) {
    const std::vector<std::string> files = FindRootFilesInDirectory(input_dir);
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
  double data_integral = 0.;
  double mc_integral = 0.;
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

    if (apply_range &&
        (entry.value < xmin || entry.value >= xmax)) {
      continue;
    }

    selected_values.push_back(entry);
  }

  if (selected_values.empty()) {
    return summary;
  }

  std::sort(selected_values.begin(),
            selected_values.end(),
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
  summary.q16 = WeightedQuantileFromSortedValues(selected_values,
                                                 0.16,
                                                 summary.sum_weight);
  summary.median = WeightedQuantileFromSortedValues(selected_values,
                                                    0.50,
                                                    summary.sum_weight);
  summary.q84 = WeightedQuantileFromSortedValues(selected_values,
                                                 0.84,
                                                 summary.sum_weight);
  summary.sigma68 = 0.5 * (summary.q84 - summary.q16);
  summary.valid = std::isfinite(summary.mean) &&
                  std::isfinite(summary.median) &&
                  std::isfinite(summary.sigma68);

  return summary;
}

double CalculateMedian(const std::vector<WeightedValue> &values) {
  const ResidualSummary summary =
    CalculateSummary(values,
                     -std::numeric_limits<double>::infinity(),
                     std::numeric_limits<double>::infinity(),
                     false);
  return summary.valid ? summary.median : 0.;
}

double CalculateWeightedQuantile(std::vector<WeightedValue> values,
                                 double probability) {
  values.erase(std::remove_if(values.begin(),
                              values.end(),
                              [](const WeightedValue &entry) {
                                return !std::isfinite(entry.value) ||
                                       !std::isfinite(entry.weight) ||
                                       entry.weight <= 0.;
                              }),
               values.end());

  if (values.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::sort(values.begin(),
            values.end(),
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
  std::vector<WeightedValue> combined_values;
  const auto normalized_data = MakeAreaNormalizedValues(data_values);
  const auto normalized_mc = MakeAreaNormalizedValues(mc_values);

  combined_values.reserve(normalized_data.size() + normalized_mc.size());
  combined_values.insert(combined_values.end(),
                         normalized_data.begin(),
                         normalized_data.end());
  combined_values.insert(combined_values.end(),
                         normalized_mc.begin(),
                         normalized_mc.end());
  return combined_values;
}

Chi2Result CalculateChi2Q01Q99(const std::vector<WeightedValue> &data_values,
                               const std::vector<WeightedValue> &mc_values,
                               int nbins) {
  static int counter = 0;

  Chi2Result result;
  if (data_values.empty() || mc_values.empty() || nbins <= 0) {
    return result;
  }

  const auto combined_values =
    MakeCombinedNormalizedValues(data_values, mc_values);
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
        entry.weight <= 0.) {
      continue;
    }
    if (entry.value < result.q01 || entry.value >= result.q99) {
      continue;
    }
    data_hist.Fill(entry.value, entry.weight);
  }

  for (const auto &entry : mc_values) {
    if (!std::isfinite(entry.value) ||
        !std::isfinite(entry.weight) ||
        entry.weight <= 0.) {
      continue;
    }
    if (entry.value < result.q01 || entry.value >= result.q99) {
      continue;
    }
    mc_hist.Fill(entry.value, entry.weight);
  }

  result.data_integral = data_hist.Integral();
  result.mc_integral = mc_hist.Integral();

  if (result.data_integral <= 0. || result.mc_integral <= 0.) {
    return result;
  }

  data_hist.Scale(1. / result.data_integral);
  mc_hist.Scale(1. / result.mc_integral);

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

double EstimatePeakPosition(const std::vector<WeightedValue> &values,
                            int nbins,
                            double xmin,
                            double xmax) {
  if (nbins <= 0 || xmin >= xmax) {
    return CalculateMedian(values);
  }

  std::vector<double> bin_sums(nbins, 0.);
  const double bin_width = (xmax - xmin) / static_cast<double>(nbins);

  for (const auto &entry : values) {
    if (!std::isfinite(entry.value) ||
        !std::isfinite(entry.weight) ||
        entry.weight <= 0. ||
        entry.value < xmin ||
        entry.value >= xmax) {
      continue;
    }

    const int bin = static_cast<int>((entry.value - xmin) / bin_width);
    if (bin >= 0 && bin < nbins) {
      bin_sums.at(bin) += entry.weight;
    }
  }

  const auto max_it = std::max_element(bin_sums.begin(), bin_sums.end());
  if (max_it == bin_sums.end() || *max_it <= 0.) {
    return CalculateMedian(values);
  }

  const int max_bin = static_cast<int>(std::distance(bin_sums.begin(), max_it));
  return xmin + (static_cast<double>(max_bin) + 0.5) * bin_width;
}

struct McResidualEvent {
  double p = 0.;      // PM-only - true
  double d = 0.;      // DWG-only - true
  double e = 0.;      // external combined - true
  double r = 0.;      // reconstructed FROST - true
  double weight = 1.;
};

struct DataResidualEvent {
  double split = 0.;       // PM-only - DWG-only
  double o = 0.;           // reconstructed FROST - external combined
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

void AddDataEvent(
    std::array<std::vector<DataResidualEvent>, kNAnalysisBins> &bins,
    int angle_bin,
    const DataResidualEvent &event) {
  bins.at(0).push_back(event);
  if (angle_bin >= 0) {
    bins.at(angle_bin + 1).push_back(event);
  }
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

std::vector<WeightedValue> MakeMcSplitValues(
    const std::vector<McResidualEvent> &events,
    double median_p,
    double median_d,
    double k_ext) {
  std::vector<WeightedValue> values;
  values.reserve(events.size());
  for (const auto &event : events) {
    const double p_morphed =
      median_p + k_ext * (event.p - median_p);
    const double d_morphed =
      median_d + k_ext * (event.d - median_d);
    values.push_back({p_morphed - d_morphed, event.weight});
  }
  return values;
}

std::vector<WeightedValue> MakeMcRecExtValues(
    const std::vector<McResidualEvent> &events,
    double median_r,
    double median_e,
    double k_rec,
    double k_ext_combined) {
  std::vector<WeightedValue> values;
  values.reserve(events.size());
  for (const auto &event : events) {
    const double r_morphed =
      median_r + k_rec * (event.r - median_r);
    const double e_morphed =
      median_e + k_ext_combined * (event.e - median_e);
    values.push_back({r_morphed - e_morphed, event.weight});
  }
  return values;
}

std::vector<WeightedValue> MakeMcRecTrueValues(
    const std::vector<McResidualEvent> &events,
    double median_r,
    double k_rec) {
  std::vector<WeightedValue> values;
  values.reserve(events.size());
  for (const auto &event : events) {
    const double r_morphed =
      median_r + k_rec * (event.r - median_r);
    values.push_back({r_morphed, event.weight});
  }
  return values;
}

std::vector<WeightedValue> MakeMcComponentValues(
    const std::vector<McResidualEvent> &events,
    char component) {
  std::vector<WeightedValue> values;
  values.reserve(events.size());

  for (const auto &event : events) {
    double value = 0.;
    if (component == 'p') {
      value = event.p;
    } else if (component == 'd') {
      value = event.d;
    } else if (component == 'e') {
      value = event.e;
    } else if (component == 'r') {
      value = event.r;
    }
    values.push_back({value, event.weight});
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
    std::cerr << "Warning: match_info tree is missing in "
              << file_path << std::endl;
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
  ok &= SetVectorBranchAddress(tree, "trackmatch_ninja_track_type",
                               &ninja_track_type);

  ok &= SetVectorBranchAddress(tree, "trackmatch_external_expected_x",
                               &external_expected_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_expected_y",
                               &external_expected_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_pm_only_expected_x",
                               &pm_only_expected_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_pm_only_expected_y",
                               &pm_only_expected_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_dwg_only_expected_x",
                               &dwg_only_expected_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_dwg_only_expected_y",
                               &dwg_only_expected_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_frost_nearest_x",
                               &frost_nearest_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_frost_nearest_y",
                               &frost_nearest_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_tangent_x",
                               &external_tangent_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_tangent_y",
                               &external_tangent_y);

  if (IsMonteCarlo(sample_type)) {
    ok &= SetVectorBranchAddress(tree,
                                 "trackmatch_true_frost_nearest_position_x",
                                 &true_frost_x);
    ok &= SetVectorBranchAddress(tree,
                                 "trackmatch_true_frost_nearest_position_y",
                                 &true_frost_y);
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

      const double theta_x_deg =
        std::atan(std::abs(external_tangent_x->at(itrack))) *
        180.0 / TMath::Pi();
      const double theta_y_deg =
        std::atan(std::abs(external_tangent_y->at(itrack))) *
        180.0 / TMath::Pi();

      const int angle_bin_x = FindAngleBin(theta_x_deg);
      const int angle_bin_y = FindAngleBin(theta_y_deg);

      const bool valid_common_x =
        HasIndex(pm_only_expected_x, itrack) &&
        HasIndex(dwg_only_expected_x, itrack) &&
        HasIndex(external_expected_x, itrack) &&
        HasIndex(frost_nearest_x, itrack) &&
        IsValidValue(pm_only_expected_x->at(itrack)) &&
        IsValidValue(dwg_only_expected_x->at(itrack)) &&
        IsValidValue(external_expected_x->at(itrack)) &&
        IsValidValue(frost_nearest_x->at(itrack));

      const bool valid_common_y =
        HasIndex(pm_only_expected_y, itrack) &&
        HasIndex(dwg_only_expected_y, itrack) &&
        HasIndex(external_expected_y, itrack) &&
        HasIndex(frost_nearest_y, itrack) &&
        IsValidValue(pm_only_expected_y->at(itrack)) &&
        IsValidValue(dwg_only_expected_y->at(itrack)) &&
        IsValidValue(external_expected_y->at(itrack)) &&
        IsValidValue(frost_nearest_y->at(itrack));

      if (IsData(sample_type)) {
        const bool split_x_selected =
          HasIndex(pm_only_expected_x, itrack) &&
          HasIndex(dwg_only_expected_x, itrack) &&
          IsValidValue(pm_only_expected_x->at(itrack)) &&
          IsValidValue(dwg_only_expected_x->at(itrack));

        const bool rec_ext_x_selected =
          HasIndex(external_expected_x, itrack) &&
          HasIndex(frost_nearest_x, itrack) &&
          IsValidValue(external_expected_x->at(itrack)) &&
          IsValidValue(frost_nearest_x->at(itrack));

        if (split_x_selected) {
          const DataResidualEvent event = {
            pm_only_expected_x->at(itrack) - dwg_only_expected_x->at(itrack),
            0.,
            weight,
            true,
            false
          };
          AddDataEvent(samples.data_x, angle_bin_x, event);
        }

        if (rec_ext_x_selected) {
          const DataResidualEvent event = {
            0.,
            frost_nearest_x->at(itrack) - external_expected_x->at(itrack),
            weight,
            false,
            true
          };
          AddDataEvent(samples.data_x, angle_bin_x, event);
        }

        const bool split_y_selected =
          HasIndex(pm_only_expected_y, itrack) &&
          HasIndex(dwg_only_expected_y, itrack) &&
          IsValidValue(pm_only_expected_y->at(itrack)) &&
          IsValidValue(dwg_only_expected_y->at(itrack));

        const bool rec_ext_y_selected =
          HasIndex(external_expected_y, itrack) &&
          HasIndex(frost_nearest_y, itrack) &&
          IsValidValue(external_expected_y->at(itrack)) &&
          IsValidValue(frost_nearest_y->at(itrack));

        if (split_y_selected) {
          const DataResidualEvent event = {
            pm_only_expected_y->at(itrack) - dwg_only_expected_y->at(itrack),
            0.,
            weight,
            true,
            false
          };
          AddDataEvent(samples.data_y, angle_bin_y, event);
        }

        if (rec_ext_y_selected) {
          const DataResidualEvent event = {
            0.,
            frost_nearest_y->at(itrack) - external_expected_y->at(itrack),
            weight,
            false,
            true
          };
          AddDataEvent(samples.data_y, angle_bin_y, event);
        }
      } else {
        if (valid_common_x &&
            HasIndex(true_frost_x, itrack) &&
            IsValidValue(true_frost_x->at(itrack))) {
          const double truth = true_frost_x->at(itrack);
          const McResidualEvent event = {
            pm_only_expected_x->at(itrack) - truth,
            dwg_only_expected_x->at(itrack) - truth,
            external_expected_x->at(itrack) - truth,
            frost_nearest_x->at(itrack) - truth,
            weight
          };
          AddMcEvent(samples.mc_x, angle_bin_x, event);
        }

        if (valid_common_y &&
            HasIndex(true_frost_y, itrack) &&
            IsValidValue(true_frost_y->at(itrack))) {
          const double truth = true_frost_y->at(itrack);
          const McResidualEvent event = {
            pm_only_expected_y->at(itrack) - truth,
            dwg_only_expected_y->at(itrack) - truth,
            external_expected_y->at(itrack) - truth,
            frost_nearest_y->at(itrack) - truth,
            weight
          };
          AddMcEvent(samples.mc_y, angle_bin_y, event);
        }
      }
    }
  }
}

struct ExtModelVariation {
  bool valid = false;
  double alpha = 1.;
  double k_ext_combined = 1.;
  double k_rec = 1.;
  ResidualSummary rec_ext_after;
  ResidualSummary rec_true_after;
  double rec_ext_target_diff = 0.;
  double delta_resolution = 0.;
};

struct RecModelVariation {
  bool valid = false;
  int n_trials = 0;
  double alpha = 1.;
  double k_ext_combined = 1.;
  double reference_rec_true_sigma68 = 0.;
  double k_rec_mean = 0.;
  double k_rec_rms = 0.;
  double rec_ext_after_mean = 0.;
  double rec_true_after_mean = 0.;
  double rec_true_after_rms = 0.;
  double rec_true_after_min = 0.;
  double rec_true_after_max = 0.;
  double delta_resolution = 0.;
  double delta_resolution_max_abs = 0.;
};

struct MorphResult {
  Axis axis = Axis::kX;
  int analysis_bin = 0;

  bool valid = false;

  Long64_t n_mc = 0;
  Long64_t n_data = 0;
  Long64_t n_data_split = 0;
  Long64_t n_data_rec_ext = 0;

  ResidualSummary split_data;
  ResidualSummary split_mc_before;
  ResidualSummary split_mc_after;

  ResidualSummary rec_ext_data;
  ResidualSummary rec_ext_mc_before;
  ResidualSummary rec_ext_mc_after;

  ResidualSummary rec_true_mc_before;
  ResidualSummary rec_true_mc_after;

  double median_p = 0.;
  double median_d = 0.;
  double median_e = 0.;
  double median_r = 0.;

  // k_ext is determined from the PM-DWG split. k_rec is the value for
  // alpha = 1.0 and is kept for diagnostic purposes.
  double k_ext = 1.;
  double k_rec = 1.;
  double rec_ext_target_diff = 0.;

  std::vector<ExtModelVariation> alpha_scan;
  ExtModelVariation alpha_max;
  ExtModelVariation ext_alpha05;
  ExtModelVariation ext_alpha0;

  // delta_ext_model is the alpha-envelope uncertainty.
  double delta_ext_model = 0.;

  RecModelVariation rec_shuffle;
  std::vector<RecModelVariation> rec_shuffle_scan;
  double delta_rec_model = 0.;
  double delta_syst_model = 0.;

  Chi2Result split_chi2;
  Chi2Result rec_ext_chi2;

  std::string warning;
};

double FindBestKExt(const std::vector<McResidualEvent> &mc_events,
                    double median_p,
                    double median_d,
                    double target_sigma68,
                    double xmin,
                    double xmax,
                    double k_min,
                    double k_max) {
  if (mc_events.empty() || !std::isfinite(target_sigma68)) {
    return 1.;
  }

  auto sigma_for_k = [&](double k_ext) -> double {
    const auto values = MakeMcSplitValues(mc_events,
                                          median_p,
                                          median_d,
                                          k_ext);
    const ResidualSummary summary = CalculateSummary(values, xmin, xmax, true);
    return summary.valid ? summary.sigma68
                         : std::numeric_limits<double>::quiet_NaN();
  };

  double best_k = 1.;
  double best_diff = std::numeric_limits<double>::infinity();

  auto scan_range = [&](double min_k, double max_k, double step) {
    for (double k = min_k; k <= max_k + 0.5 * step; k += step) {
      const double sigma = sigma_for_k(k);
      if (!std::isfinite(sigma)) {
        continue;
      }
      const double diff = std::abs(sigma - target_sigma68);
      if (diff < best_diff) {
        best_diff = diff;
        best_k = k;
      }
    }
  };

  // The simple ratio S_data/S_MC would be exact for an affine scaling if all
  // residual values were used. Here sigma68 is evaluated only within the
  // displayed histogram range, so the selected values can change as k_ext
  // changes. Therefore k_ext is determined by a direct scan.
  scan_range(k_min, k_max, 0.05);

  const double refine_min = std::max(k_min, best_k - 0.10);
  const double refine_max = std::min(k_max, best_k + 0.10);
  scan_range(refine_min, refine_max, 0.002);

  return best_k;
}

double FindBestKRec(const std::vector<McResidualEvent> &mc_events,
                    double median_r,
                    double median_e,
                    double k_ext_combined,
                    double target_sigma68,
                    double xmin,
                    double xmax,
                    double k_min,
                    double k_max) {
  if (mc_events.empty() || !std::isfinite(target_sigma68)) {
    return 1.;
  }

  auto sigma_for_k = [&](double k_rec) -> double {
    const auto values = MakeMcRecExtValues(mc_events,
                                           median_r,
                                           median_e,
                                           k_rec,
                                           k_ext_combined);
    const ResidualSummary summary = CalculateSummary(values, xmin, xmax, true);
    return summary.valid ? summary.sigma68
                         : std::numeric_limits<double>::quiet_NaN();
  };

  double best_k = 1.;
  double best_diff = std::numeric_limits<double>::infinity();

  auto scan_range = [&](double min_k, double max_k, double step) {
    for (double k = min_k; k <= max_k + 0.5 * step; k += step) {
      const double sigma = sigma_for_k(k);
      if (!std::isfinite(sigma)) {
        continue;
      }
      const double diff = std::abs(sigma - target_sigma68);
      if (diff < best_diff) {
        best_diff = diff;
        best_k = k;
      }
    }
  };

  // Coarse scan followed by a local fine scan. This avoids assuming an exact
  // quadrature relation or strict Gaussian behavior.
  scan_range(k_min, k_max, 0.05);

  const double refine_min = std::max(k_min, best_k - 0.10);
  const double refine_max = std::min(k_max, best_k + 0.10);
  scan_range(refine_min, refine_max, 0.002);

  return best_k;
}

std::vector<WeightedValue> MakeMcRecExtValuesShuffledR(
    const std::vector<McResidualEvent> &events,
    const std::vector<std::size_t> &shuffled_indices,
    double median_r,
    double median_e,
    double k_rec,
    double k_ext_combined) {
  std::vector<WeightedValue> values;
  values.reserve(events.size());

  if (events.empty() || shuffled_indices.size() != events.size()) {
    return values;
  }

  for (std::size_t i = 0; i < events.size(); ++i) {
    const auto &event_e = events.at(i);
    const auto &event_r = events.at(shuffled_indices.at(i));
    const double r_morphed =
      median_r + k_rec * (event_r.r - median_r);
    const double e_morphed =
      median_e + k_ext_combined * (event_e.e - median_e);
    values.push_back({r_morphed - e_morphed, event_e.weight});
  }

  return values;
}

double FindBestKRecShuffledR(
    const std::vector<McResidualEvent> &mc_events,
    const std::vector<std::size_t> &shuffled_indices,
    double median_r,
    double median_e,
    double k_ext_combined,
    double target_sigma68,
    double xmin,
    double xmax,
    double k_min,
    double k_max) {
  if (mc_events.empty() || !std::isfinite(target_sigma68)) {
    return 1.;
  }

  auto sigma_for_k = [&](double k_rec) -> double {
    const auto values = MakeMcRecExtValuesShuffledR(mc_events,
                                                   shuffled_indices,
                                                   median_r,
                                                   median_e,
                                                   k_rec,
                                                   k_ext_combined);
    const ResidualSummary summary = CalculateSummary(values, xmin, xmax, true);
    return summary.valid ? summary.sigma68
                         : std::numeric_limits<double>::quiet_NaN();
  };

  double best_k = 1.;
  double best_diff = std::numeric_limits<double>::infinity();

  auto scan_range = [&](double min_k, double max_k, double step) {
    for (double k = min_k; k <= max_k + 0.5 * step; k += step) {
      const double sigma = sigma_for_k(k);
      if (!std::isfinite(sigma)) {
        continue;
      }
      const double diff = std::abs(sigma - target_sigma68);
      if (diff < best_diff) {
        best_diff = diff;
        best_k = k;
      }
    }
  };

  scan_range(k_min, k_max, 0.05);

  const double refine_min = std::max(k_min, best_k - 0.10);
  const double refine_max = std::min(k_max, best_k + 0.10);
  scan_range(refine_min, refine_max, 0.002);

  return best_k;
}

ExtModelVariation EvaluateExtModelVariation(
    const std::vector<McResidualEvent> &mc_events,
    double alpha,
    double k_ext,
    double nominal_resolution,
    double median_r,
    double median_e,
    double target_rec_ext_sigma68,
    double residual_xmin,
    double residual_xmax,
    double k_min,
    double k_max) {
  ExtModelVariation variation;
  variation.alpha = alpha;
  variation.k_ext_combined = 1. + alpha * (k_ext - 1.);

  variation.k_rec = FindBestKRec(mc_events,
                                 median_r,
                                 median_e,
                                 variation.k_ext_combined,
                                 target_rec_ext_sigma68,
                                 residual_xmin,
                                 residual_xmax,
                                 k_min,
                                 k_max);

  const auto rec_ext_values = MakeMcRecExtValues(mc_events,
                                                 median_r,
                                                 median_e,
                                                 variation.k_rec,
                                                 variation.k_ext_combined);
  const auto rec_true_values = MakeMcRecTrueValues(mc_events,
                                                   median_r,
                                                   variation.k_rec);

  variation.rec_ext_after = CalculateSummary(rec_ext_values,
                                             residual_xmin,
                                             residual_xmax,
                                             true);
  variation.rec_true_after = CalculateSummary(rec_true_values,
                                              residual_xmin,
                                              residual_xmax,
                                              true);

  if (variation.rec_ext_after.valid && variation.rec_true_after.valid) {
    variation.rec_ext_target_diff =
      variation.rec_ext_after.sigma68 - target_rec_ext_sigma68;
    variation.delta_resolution =
      variation.rec_true_after.sigma68 - nominal_resolution;
    variation.valid = true;
  }

  return variation;
}

RecModelVariation EvaluateRecShuffleVariation(
    const std::vector<McResidualEvent> &mc_events,
    int n_shuffle_trials,
    unsigned int shuffle_seed,
    double alpha,
    double reference_rec_true_sigma68,
    double median_r,
    double median_e,
    double k_ext_combined,
    double target_rec_ext_sigma68,
    double residual_xmin,
    double residual_xmax,
    double k_min,
    double k_max) {
  RecModelVariation variation;
  variation.n_trials = 0;
  variation.alpha = alpha;
  variation.k_ext_combined = k_ext_combined;
  variation.reference_rec_true_sigma68 = reference_rec_true_sigma68;

  if (mc_events.empty() || n_shuffle_trials <= 0) {
    return variation;
  }

  std::vector<std::size_t> indices(mc_events.size());
  std::iota(indices.begin(), indices.end(), 0);

  std::mt19937 rng(shuffle_seed);
  std::vector<double> k_rec_values;
  std::vector<double> rec_ext_values;
  std::vector<double> rec_true_values;
  k_rec_values.reserve(n_shuffle_trials);
  rec_ext_values.reserve(n_shuffle_trials);
  rec_true_values.reserve(n_shuffle_trials);

  for (int itrial = 0; itrial < n_shuffle_trials; ++itrial) {
    std::shuffle(indices.begin(), indices.end(), rng);

    const double k_rec_shuffle = FindBestKRecShuffledR(mc_events,
                                                       indices,
                                                       median_r,
                                                       median_e,
                                                       k_ext_combined,
                                                       target_rec_ext_sigma68,
                                                       residual_xmin,
                                                       residual_xmax,
                                                       k_min,
                                                       k_max);

    const auto rec_ext_shuffled = MakeMcRecExtValuesShuffledR(mc_events,
                                                              indices,
                                                              median_r,
                                                              median_e,
                                                              k_rec_shuffle,
                                                              k_ext_combined);
    const auto rec_true_shuffled = MakeMcRecTrueValues(mc_events,
                                                       median_r,
                                                       k_rec_shuffle);

    const ResidualSummary rec_ext_summary =
      CalculateSummary(rec_ext_shuffled, residual_xmin, residual_xmax, true);
    const ResidualSummary rec_true_summary =
      CalculateSummary(rec_true_shuffled, residual_xmin, residual_xmax, true);

    if (!rec_ext_summary.valid || !rec_true_summary.valid) {
      continue;
    }

    k_rec_values.push_back(k_rec_shuffle);
    rec_ext_values.push_back(rec_ext_summary.sigma68);
    rec_true_values.push_back(rec_true_summary.sigma68);
  }

  if (rec_true_values.empty()) {
    return variation;
  }

  auto mean = [](const std::vector<double> &values) -> double {
    const double sum = std::accumulate(values.begin(), values.end(), 0.);
    return sum / static_cast<double>(values.size());
  };

  auto rms = [&](const std::vector<double> &values, double mean_value) -> double {
    double sum2 = 0.;
    for (const auto value : values) {
      const double diff = value - mean_value;
      sum2 += diff * diff;
    }
    return std::sqrt(sum2 / static_cast<double>(values.size()));
  };

  variation.n_trials = static_cast<int>(rec_true_values.size());
  variation.k_rec_mean = mean(k_rec_values);
  variation.k_rec_rms = rms(k_rec_values, variation.k_rec_mean);
  variation.rec_ext_after_mean = mean(rec_ext_values);
  variation.rec_true_after_mean = mean(rec_true_values);
  variation.rec_true_after_rms = rms(rec_true_values,
                                     variation.rec_true_after_mean);
  variation.rec_true_after_min =
    *std::min_element(rec_true_values.begin(), rec_true_values.end());
  variation.rec_true_after_max =
    *std::max_element(rec_true_values.begin(), rec_true_values.end());
  variation.delta_resolution =
    variation.rec_true_after_mean - reference_rec_true_sigma68;

  variation.delta_resolution_max_abs = 0.;
  for (const auto value : rec_true_values) {
    variation.delta_resolution_max_abs =
      std::max(variation.delta_resolution_max_abs,
               std::abs(value - reference_rec_true_sigma68));
  }

  variation.valid = true;
  return variation;
}

MorphResult EstimateMorphResult(Axis axis,
                                int analysis_bin,
                                const std::vector<McResidualEvent> &mc_events,
                                const std::vector<DataResidualEvent> &data_events,
                                int nbins,
                                double residual_xmin,
                                double residual_xmax,
                                double split_xmin,
                                double split_xmax,
                                double k_min,
                                double k_max,
                                int n_shuffle_trials,
                                unsigned int shuffle_seed) {
  MorphResult result;
  result.axis = axis;
  result.analysis_bin = analysis_bin;
  result.n_mc = static_cast<Long64_t>(mc_events.size());
  result.n_data = static_cast<Long64_t>(data_events.size());

  const auto data_split = MakeDataSplitValues(data_events);
  const auto data_rec_ext = MakeDataRecExtValues(data_events);
  result.n_data_split = static_cast<Long64_t>(data_split.size());
  result.n_data_rec_ext = static_cast<Long64_t>(data_rec_ext.size());

  if (mc_events.empty() || data_split.empty() || data_rec_ext.empty()) {
    result.warning = "empty MC, data split, or data rec-ext sample";
    return result;
  }

  result.split_data = CalculateSummary(data_split, split_xmin, split_xmax, true);
  result.rec_ext_data = CalculateSummary(data_rec_ext,
                                         residual_xmin,
                                         residual_xmax,
                                         true);

  const auto p_values = MakeMcComponentValues(mc_events, 'p');
  const auto d_values = MakeMcComponentValues(mc_events, 'd');
  const auto e_values = MakeMcComponentValues(mc_events, 'e');
  const auto r_values = MakeMcComponentValues(mc_events, 'r');

  result.median_p = CalculateMedian(p_values);
  result.median_d = CalculateMedian(d_values);
  result.median_e = CalculateMedian(e_values);
  result.median_r = CalculateMedian(r_values);

  const auto mc_split_before =
    MakeMcSplitValues(mc_events, result.median_p, result.median_d, 1.0);
  const auto mc_rec_ext_before =
    MakeMcRecExtValues(mc_events, result.median_r, result.median_e, 1.0, 1.0);
  const auto mc_rec_true_before =
    MakeMcRecTrueValues(mc_events, result.median_r, 1.0);

  result.split_mc_before =
    CalculateSummary(mc_split_before, split_xmin, split_xmax, true);
  result.rec_ext_mc_before =
    CalculateSummary(mc_rec_ext_before, residual_xmin, residual_xmax, true);
  result.rec_true_mc_before =
    CalculateSummary(mc_rec_true_before, residual_xmin, residual_xmax, true);

  result.split_chi2 =
    CalculateChi2Q01Q99(data_split, mc_split_before, nbins);
  result.rec_ext_chi2 =
    CalculateChi2Q01Q99(data_rec_ext, mc_rec_ext_before, nbins);

  if (!result.split_data.valid ||
      !result.rec_ext_data.valid ||
      !result.split_mc_before.valid ||
      !result.rec_ext_mc_before.valid ||
      !result.rec_true_mc_before.valid ||
      result.split_mc_before.sigma68 <= 0.) {
    result.warning = "invalid sigma68 input";
    return result;
  }

  result.k_ext =
    FindBestKExt(mc_events,
                 result.median_p,
                 result.median_d,
                 result.split_data.sigma68,
                 split_xmin,
                 split_xmax,
                 k_min,
                 k_max);
  if (!std::isfinite(result.k_ext)) {
    result.warning = "invalid k_ext";
    return result;
  }

  const auto mc_split_after =
    MakeMcSplitValues(mc_events,
                      result.median_p,
                      result.median_d,
                      result.k_ext);
  result.split_mc_after =
    CalculateSummary(mc_split_after, split_xmin, split_xmax, true);

  const double nominal_resolution = result.rec_true_mc_before.sigma68;
  result.delta_ext_model = 0.;
  result.alpha_scan.clear();

  for (int ialpha = 0; ialpha <= 10; ++ialpha) {
    const double alpha = 0.1 * static_cast<double>(ialpha);
    ExtModelVariation variation = EvaluateExtModelVariation(mc_events,
                                                            alpha,
                                                            result.k_ext,
                                                            nominal_resolution,
                                                            result.median_r,
                                                            result.median_e,
                                                            result.rec_ext_data.sigma68,
                                                            residual_xmin,
                                                            residual_xmax,
                                                            k_min,
                                                            k_max);
    result.alpha_scan.push_back(variation);

    if (ialpha == 0) {
      result.ext_alpha0 = variation;
    }
    if (ialpha == 5) {
      result.ext_alpha05 = variation;
    }
    if (ialpha == 10) {
      result.k_rec = variation.k_rec;
      result.rec_ext_mc_after = variation.rec_ext_after;
      result.rec_true_mc_after = variation.rec_true_after;
      result.rec_ext_target_diff = variation.rec_ext_target_diff;
    }

    if (variation.valid) {
      const double abs_delta = std::abs(variation.delta_resolution);
      if (!result.alpha_max.valid || abs_delta > result.delta_ext_model) {
        result.delta_ext_model = abs_delta;
        result.alpha_max = variation;
      }
    }
  }

  const unsigned int seed_offset =
    static_cast<unsigned int>(1000 * (axis == Axis::kY ? 1 : 0) +
                              analysis_bin);
  result.rec_shuffle_scan.clear();
  result.delta_rec_model = 0.;

  for (std::size_t ialpha = 0; ialpha < result.alpha_scan.size(); ++ialpha) {
    const auto &alpha_variation = result.alpha_scan.at(ialpha);
    if (!alpha_variation.valid) {
      continue;
    }

    RecModelVariation shuffle_variation =
      EvaluateRecShuffleVariation(mc_events,
                                  n_shuffle_trials,
                                  shuffle_seed + seed_offset +
                                    static_cast<unsigned int>(100 * ialpha),
                                  alpha_variation.alpha,
                                  alpha_variation.rec_true_after.sigma68,
                                  result.median_r,
                                  result.median_e,
                                  alpha_variation.k_ext_combined,
                                  result.rec_ext_data.sigma68,
                                  residual_xmin,
                                  residual_xmax,
                                  k_min,
                                  k_max);

    result.rec_shuffle_scan.push_back(shuffle_variation);

    if (shuffle_variation.valid) {
      const double abs_delta = std::abs(shuffle_variation.delta_resolution);
      if (!result.rec_shuffle.valid || abs_delta > result.delta_rec_model) {
        result.delta_rec_model = abs_delta;
        result.rec_shuffle = shuffle_variation;
      }
    }
  }

  result.delta_syst_model =
    std::sqrt(result.delta_ext_model * result.delta_ext_model +
              result.delta_rec_model * result.delta_rec_model);

  result.valid = result.rec_true_mc_before.valid;
  return result;
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

std::string FormatDouble(double value, int precision = 6) {
  if (!std::isfinite(value)) {
    return "nan";
  }

  std::ostringstream stream;
  stream << std::fixed << std::setprecision(precision) << value;
  return stream.str();
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

std::string SanitizeCsvText(std::string value) {
  for (char &character : value) {
    if (character == ',') {
      character = ';';
    }
    if (character == '\n' || character == '\r') {
      character = ' ';
    }
  }
  return value;
}

std::string FormatAlphaScanValues(
    const std::vector<ExtModelVariation> &alpha_scan,
    char quantity,
    int precision = 4) {
  std::ostringstream stream;
  bool first = true;
  for (const auto &variation : alpha_scan) {
    if (!variation.valid) {
      continue;
    }

    double value = std::numeric_limits<double>::quiet_NaN();
    if (quantity == 'r') {
      value = variation.rec_true_after.sigma68;
    } else if (quantity == 'd') {
      value = variation.delta_resolution;
    } else if (quantity == 'k') {
      value = variation.k_rec;
    } else if (quantity == 'o') {
      value = variation.rec_ext_after.sigma68;
    }

    if (!std::isfinite(value)) {
      continue;
    }

    if (!first) {
      stream << '|';
    }
    stream << std::fixed << std::setprecision(1) << variation.alpha
           << ':' << std::setprecision(precision) << value;
    first = false;
  }
  return stream.str();
}

std::string FormatRecShuffleScanValues(
    const std::vector<RecModelVariation> &shuffle_scan,
    char quantity,
    int precision = 4) {
  std::ostringstream stream;
  bool first = true;
  for (const auto &variation : shuffle_scan) {
    if (!variation.valid) {
      continue;
    }

    double value = std::numeric_limits<double>::quiet_NaN();
    if (quantity == 'r') {
      value = variation.rec_true_after_mean;
    } else if (quantity == 'd') {
      value = variation.delta_resolution;
    } else if (quantity == 'k') {
      value = variation.k_rec_mean;
    } else if (quantity == 'o') {
      value = variation.rec_ext_after_mean;
    } else if (quantity == 's') {
      value = variation.rec_true_after_rms;
    } else if (quantity == 'b') {
      value = variation.reference_rec_true_sigma68;
    }

    if (!std::isfinite(value)) {
      continue;
    }

    if (!first) {
      stream << '|';
    }
    stream << std::fixed << std::setprecision(1) << variation.alpha
           << ':' << std::setprecision(precision) << value;
    first = false;
  }
  return stream.str();
}

std::string MakeCsvHeader() {
  return "axis,angle_bin,valid,n_mc,n_data,n_data_split,n_data_rec_ext,"
         "nominal_rec_true_sigma68_mm,"
         "nominal_rec_true_mean_mm,"
         "nominal_rec_true_median_mm,"
         "nominal_rec_true_q16_mm,"
         "nominal_rec_true_q84_mm,"
         "delta_alpha_model_mm,alpha_at_max,rec_true_alpha_at_max_mm,"
         "k_rec_alpha_at_max,k_ext_combined_alpha_at_max,"
         "rec_ext_alpha_at_max_sigma68_mm,rec_ext_target_alpha_at_max_diff_mm,"
         "delta_corr_model_mm,delta_model_mm,"
         "rec_shuffle_alpha_at_max,rec_shuffle_ref_rec_true_sigma68_mm,"
         "rec_shuffle_n_trials,rec_shuffle_k_rec_mean,rec_shuffle_k_rec_rms,"
         "rec_shuffle_rec_ext_after_mean_mm,rec_shuffle_rec_true_after_mean_mm,"
         "rec_shuffle_rec_true_after_rms_mm,"
         "split_data_sigma68_mm,split_mc_sigma68_mm,split_data_minus_mc_mm,"
         "split_relative_diff,rec_ext_data_sigma68_mm,rec_ext_mc_sigma68_mm,"
         "rec_ext_data_minus_mc_mm,rec_ext_relative_diff,"
         "split_chi2_q01_mm,split_chi2_q99_mm,split_chi2,split_ndf,"
         "split_chi2_ndf,split_p_value,"
         "rec_ext_chi2_q01_mm,rec_ext_chi2_q99_mm,rec_ext_chi2,"
         "rec_ext_ndf,rec_ext_chi2_ndf,rec_ext_p_value,"
         "alpha_scan_rec_true_sigma68_mm,alpha_scan_delta_resolution_mm,"
         "alpha_scan_k_rec,shuffle_scan_rec_true_sigma68_mm,"
         "shuffle_scan_delta_resolution_mm,shuffle_scan_k_rec,warning";
}

std::string MakeCsvRow(const MorphResult &result) {
  std::ostringstream out;
  out << AxisLabel(result.axis) << ','
      << AnalysisBinLabel(result.analysis_bin) << ','
      << (result.valid ? 1 : 0) << ','
      << result.n_mc << ','
      << result.n_data << ','
      << result.n_data_split << ','
      << result.n_data_rec_ext << ','
      << FormatDouble(result.rec_true_mc_before.sigma68) << ','
      << FormatDouble(result.rec_true_mc_before.mean) << ','
      << FormatDouble(result.rec_true_mc_before.median) << ','
      << FormatDouble(result.rec_true_mc_before.q16) << ','
      << FormatDouble(result.rec_true_mc_before.q84) << ','
      << FormatDouble(result.delta_ext_model) << ','
      << FormatDouble(result.alpha_max.alpha, 1) << ','
      << FormatDouble(result.alpha_max.rec_true_after.sigma68) << ','
      << FormatDouble(result.alpha_max.k_rec) << ','
      << FormatDouble(result.alpha_max.k_ext_combined) << ','
      << FormatDouble(result.alpha_max.rec_ext_after.sigma68) << ','
      << FormatDouble(result.alpha_max.rec_ext_target_diff) << ','
      << FormatDouble(result.delta_rec_model) << ','
      << FormatDouble(result.delta_syst_model) << ','
      << FormatDouble(result.rec_shuffle.alpha, 1) << ','
      << FormatDouble(result.rec_shuffle.reference_rec_true_sigma68) << ','
      << result.rec_shuffle.n_trials << ','
      << FormatDouble(result.rec_shuffle.k_rec_mean) << ','
      << FormatDouble(result.rec_shuffle.k_rec_rms) << ','
      << FormatDouble(result.rec_shuffle.rec_ext_after_mean) << ','
      << FormatDouble(result.rec_shuffle.rec_true_after_mean) << ','
      << FormatDouble(result.rec_shuffle.rec_true_after_rms) << ','
      << FormatDouble(result.split_data.sigma68) << ','
      << FormatDouble(result.split_mc_before.sigma68) << ','
      << FormatDouble(Difference(result.split_data.sigma68,
                                 result.split_mc_before.sigma68)) << ','
      << FormatDouble(RelativeDifference(result.split_data.sigma68,
                                         result.split_mc_before.sigma68)) << ','
      << FormatDouble(result.rec_ext_data.sigma68) << ','
      << FormatDouble(result.rec_ext_mc_before.sigma68) << ','
      << FormatDouble(Difference(result.rec_ext_data.sigma68,
                                 result.rec_ext_mc_before.sigma68)) << ','
      << FormatDouble(RelativeDifference(result.rec_ext_data.sigma68,
                                         result.rec_ext_mc_before.sigma68)) << ','
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
      << FormatAlphaScanValues(result.alpha_scan, 'r') << ','
      << FormatAlphaScanValues(result.alpha_scan, 'd') << ','
      << FormatAlphaScanValues(result.alpha_scan, 'k') << ','
      << FormatRecShuffleScanValues(result.rec_shuffle_scan, 'r') << ','
      << FormatRecShuffleScanValues(result.rec_shuffle_scan, 'd') << ','
      << FormatRecShuffleScanValues(result.rec_shuffle_scan, 'k') << ','
      << SanitizeCsvText(result.warning);
  return out.str();
}

void WriteLogFile(const char *output_log,
                  const std::vector<MorphResult> &results,
                  std::vector<std::string> &log_lines) {
  log_lines.clear();
  log_lines.push_back(MakeCsvHeader());
  for (const auto &result : results) {
    log_lines.push_back(MakeCsvRow(result));
  }

  std::ofstream log(output_log);
  if (!log) {
    std::cerr << "Warning: cannot open log file: "
              << output_log << std::endl;
    return;
  }

  for (const auto &line : log_lines) {
    log << line << '\n';
  }
}

void WriteRootLogFile(const char *output_root,
                      const std::vector<MorphResult> &results) {
  TFile root_file(output_root, "RECREATE");
  if (root_file.IsZombie()) {
    std::cerr << "Warning: cannot open ROOT log file: "
              << output_root << std::endl;
    return;
  }

  TTree tree("reslog", "FROST position resolution log");

  std::string axis;
  std::string angle_bin;
  int valid = 0;
  Long64_t n_mc = 0;
  Long64_t n_data = 0;
  Long64_t n_data_split = 0;
  Long64_t n_data_rec_ext = 0;
  double nominal_rec_true_sigma68_mm = 0.;
  double nominal_rec_true_mean_mm = 0.;
  double nominal_rec_true_median_mm = 0.;
  double nominal_rec_true_q16_mm = 0.;
  double nominal_rec_true_q84_mm = 0.;
  double delta_alpha_model_mm = 0.;
  double alpha_at_max = 0.;
  double rec_true_alpha_at_max_mm = 0.;
  double k_rec_alpha_at_max = 0.;
  double k_ext_combined_alpha_at_max = 0.;
  double rec_ext_alpha_at_max_sigma68_mm = 0.;
  double rec_ext_target_alpha_at_max_diff_mm = 0.;
  double delta_corr_model_mm = 0.;
  double delta_model_mm = 0.;
  double rec_shuffle_alpha_at_max = 0.;
  double rec_shuffle_ref_rec_true_sigma68_mm = 0.;
  int rec_shuffle_n_trials = 0;
  double rec_shuffle_k_rec_mean = 0.;
  double rec_shuffle_k_rec_rms = 0.;
  double rec_shuffle_rec_ext_after_mean_mm = 0.;
  double rec_shuffle_rec_true_after_mean_mm = 0.;
  double rec_shuffle_rec_true_after_rms_mm = 0.;
  double split_data_sigma68_mm = 0.;
  double split_mc_sigma68_mm = 0.;
  double split_data_minus_mc_mm = 0.;
  double split_relative_diff = 0.;
  double rec_ext_data_sigma68_mm = 0.;
  double rec_ext_mc_sigma68_mm = 0.;
  double rec_ext_data_minus_mc_mm = 0.;
  double rec_ext_relative_diff = 0.;
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
  std::string alpha_scan_rec_true_sigma68_mm;
  std::string alpha_scan_delta_resolution_mm;
  std::string alpha_scan_k_rec;
  std::string shuffle_scan_rec_true_sigma68_mm;
  std::string shuffle_scan_delta_resolution_mm;
  std::string shuffle_scan_k_rec;
  std::string warning;

  tree.Branch("axis", &axis);
  tree.Branch("angle_bin", &angle_bin);
  tree.Branch("valid", &valid);
  tree.Branch("n_mc", &n_mc);
  tree.Branch("n_data", &n_data);
  tree.Branch("n_data_split", &n_data_split);
  tree.Branch("n_data_rec_ext", &n_data_rec_ext);
  tree.Branch("nominal_rec_true_sigma68_mm", &nominal_rec_true_sigma68_mm);
  tree.Branch("nominal_rec_true_mean_mm", &nominal_rec_true_mean_mm);
  tree.Branch("nominal_rec_true_median_mm", &nominal_rec_true_median_mm);
  tree.Branch("nominal_rec_true_q16_mm", &nominal_rec_true_q16_mm);
  tree.Branch("nominal_rec_true_q84_mm", &nominal_rec_true_q84_mm);
  tree.Branch("delta_alpha_model_mm", &delta_alpha_model_mm);
  tree.Branch("alpha_at_max", &alpha_at_max);
  tree.Branch("rec_true_alpha_at_max_mm", &rec_true_alpha_at_max_mm);
  tree.Branch("k_rec_alpha_at_max", &k_rec_alpha_at_max);
  tree.Branch("k_ext_combined_alpha_at_max", &k_ext_combined_alpha_at_max);
  tree.Branch("rec_ext_alpha_at_max_sigma68_mm", &rec_ext_alpha_at_max_sigma68_mm);
  tree.Branch("rec_ext_target_alpha_at_max_diff_mm", &rec_ext_target_alpha_at_max_diff_mm);
  tree.Branch("delta_corr_model_mm", &delta_corr_model_mm);
  tree.Branch("delta_model_mm", &delta_model_mm);
  tree.Branch("rec_shuffle_alpha_at_max", &rec_shuffle_alpha_at_max);
  tree.Branch("rec_shuffle_ref_rec_true_sigma68_mm", &rec_shuffle_ref_rec_true_sigma68_mm);
  tree.Branch("rec_shuffle_n_trials", &rec_shuffle_n_trials);
  tree.Branch("rec_shuffle_k_rec_mean", &rec_shuffle_k_rec_mean);
  tree.Branch("rec_shuffle_k_rec_rms", &rec_shuffle_k_rec_rms);
  tree.Branch("rec_shuffle_rec_ext_after_mean_mm", &rec_shuffle_rec_ext_after_mean_mm);
  tree.Branch("rec_shuffle_rec_true_after_mean_mm", &rec_shuffle_rec_true_after_mean_mm);
  tree.Branch("rec_shuffle_rec_true_after_rms_mm", &rec_shuffle_rec_true_after_rms_mm);
  tree.Branch("split_data_sigma68_mm", &split_data_sigma68_mm);
  tree.Branch("split_mc_sigma68_mm", &split_mc_sigma68_mm);
  tree.Branch("split_data_minus_mc_mm", &split_data_minus_mc_mm);
  tree.Branch("split_relative_diff", &split_relative_diff);
  tree.Branch("rec_ext_data_sigma68_mm", &rec_ext_data_sigma68_mm);
  tree.Branch("rec_ext_mc_sigma68_mm", &rec_ext_mc_sigma68_mm);
  tree.Branch("rec_ext_data_minus_mc_mm", &rec_ext_data_minus_mc_mm);
  tree.Branch("rec_ext_relative_diff", &rec_ext_relative_diff);
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
  tree.Branch("alpha_scan_rec_true_sigma68_mm", &alpha_scan_rec_true_sigma68_mm);
  tree.Branch("alpha_scan_delta_resolution_mm", &alpha_scan_delta_resolution_mm);
  tree.Branch("alpha_scan_k_rec", &alpha_scan_k_rec);
  tree.Branch("shuffle_scan_rec_true_sigma68_mm", &shuffle_scan_rec_true_sigma68_mm);
  tree.Branch("shuffle_scan_delta_resolution_mm", &shuffle_scan_delta_resolution_mm);
  tree.Branch("shuffle_scan_k_rec", &shuffle_scan_k_rec);
  tree.Branch("warning", &warning);

  for (const auto &result : results) {
    axis = AxisLabel(result.axis);
    angle_bin = AnalysisBinLabel(result.analysis_bin);
    valid = result.valid ? 1 : 0;
    n_mc = result.n_mc;
    n_data = result.n_data;
    n_data_split = result.n_data_split;
    n_data_rec_ext = result.n_data_rec_ext;
    nominal_rec_true_sigma68_mm = result.rec_true_mc_before.sigma68;
    nominal_rec_true_mean_mm = result.rec_true_mc_before.mean;
    nominal_rec_true_median_mm = result.rec_true_mc_before.median;
    nominal_rec_true_q16_mm = result.rec_true_mc_before.q16;
    nominal_rec_true_q84_mm = result.rec_true_mc_before.q84;
    delta_alpha_model_mm = result.delta_ext_model;
    alpha_at_max = result.alpha_max.alpha;
    rec_true_alpha_at_max_mm = result.alpha_max.rec_true_after.sigma68;
    k_rec_alpha_at_max = result.alpha_max.k_rec;
    k_ext_combined_alpha_at_max = result.alpha_max.k_ext_combined;
    rec_ext_alpha_at_max_sigma68_mm = result.alpha_max.rec_ext_after.sigma68;
    rec_ext_target_alpha_at_max_diff_mm = result.alpha_max.rec_ext_target_diff;
    delta_corr_model_mm = result.delta_rec_model;
    delta_model_mm = result.delta_syst_model;
    rec_shuffle_alpha_at_max = result.rec_shuffle.alpha;
    rec_shuffle_ref_rec_true_sigma68_mm = result.rec_shuffle.reference_rec_true_sigma68;
    rec_shuffle_n_trials = result.rec_shuffle.n_trials;
    rec_shuffle_k_rec_mean = result.rec_shuffle.k_rec_mean;
    rec_shuffle_k_rec_rms = result.rec_shuffle.k_rec_rms;
    rec_shuffle_rec_ext_after_mean_mm = result.rec_shuffle.rec_ext_after_mean;
    rec_shuffle_rec_true_after_mean_mm = result.rec_shuffle.rec_true_after_mean;
    rec_shuffle_rec_true_after_rms_mm = result.rec_shuffle.rec_true_after_rms;
    split_data_sigma68_mm = result.split_data.sigma68;
    split_mc_sigma68_mm = result.split_mc_before.sigma68;
    split_data_minus_mc_mm = Difference(result.split_data.sigma68,
                                        result.split_mc_before.sigma68);
    split_relative_diff = RelativeDifference(result.split_data.sigma68,
                                             result.split_mc_before.sigma68);
    rec_ext_data_sigma68_mm = result.rec_ext_data.sigma68;
    rec_ext_mc_sigma68_mm = result.rec_ext_mc_before.sigma68;
    rec_ext_data_minus_mc_mm = Difference(result.rec_ext_data.sigma68,
                                          result.rec_ext_mc_before.sigma68);
    rec_ext_relative_diff = RelativeDifference(result.rec_ext_data.sigma68,
                                               result.rec_ext_mc_before.sigma68);
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
    alpha_scan_rec_true_sigma68_mm = FormatAlphaScanValues(result.alpha_scan, 'r');
    alpha_scan_delta_resolution_mm = FormatAlphaScanValues(result.alpha_scan, 'd');
    alpha_scan_k_rec = FormatAlphaScanValues(result.alpha_scan, 'k');
    shuffle_scan_rec_true_sigma68_mm = FormatRecShuffleScanValues(result.rec_shuffle_scan, 'r');
    shuffle_scan_delta_resolution_mm = FormatRecShuffleScanValues(result.rec_shuffle_scan, 'd');
    shuffle_scan_k_rec = FormatRecShuffleScanValues(result.rec_shuffle_scan, 'k');
    warning = result.warning;

    tree.Fill();
  }

  root_file.cd();
  tree.Write();
  root_file.Close();
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
    const double value =
      hist->GetBinContent(ibin) + hist->GetBinError(ibin);
    if (value > maximum) {
      maximum = value;
    }
  }
  return maximum;
}

void DrawOverlayHistograms(const std::vector<PlotSample> &samples,
                           const char *title,
                           const char *x_axis_title,
                           int nbins,
                           double xmin,
                           double xmax) {
  static int counter = 0;

  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.13);

  double data_integral = 0.;
  std::vector<TH1D*> histograms;
  histograms.reserve(samples.size());

  for (const auto &sample : samples) {
    TString name = Form("hist_overlay_%d", counter++);
    auto *hist = new TH1D(name,
                          Form("%s;%s;Number of tracks",
                               title,
                               x_axis_title),
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

  // Scale MC-like histograms to the Data integral in this plot.
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
    hist->SetMaximum(max_y > 0. ? 1.30 * max_y : 1.0);
    hist->GetYaxis()->SetTitleOffset(1.55);
  }

  bool first = true;

  // First draw MC-like samples as histograms.
  for (std::size_t i = 0; i < samples.size(); ++i) {
    if (samples.at(i).draw_as_points) {
      continue;
    }

    auto *hist = histograms.at(i);
    if (first) {
      hist->Draw("hist");
      first = false;
    } else {
      hist->Draw("hist same");
    }
  }

  // Then draw data-like samples as points with error bars.
  for (std::size_t i = 0; i < samples.size(); ++i) {
    if (!samples.at(i).draw_as_points) {
      continue;
    }

    auto *hist = histograms.at(i);
    if (first) {
      hist->Draw("E1");
      first = false;
    } else {
      hist->Draw("E1 same");
    }
  }

  auto *legend = new TLegend(0.60, 0.74, 0.94, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.034);
  for (std::size_t i = 0; i < histograms.size(); ++i) {
    legend->AddEntry(histograms.at(i),
                     samples.at(i).label.c_str(),
                     samples.at(i).draw_as_points ? "lep" : "l");
  }
  legend->Draw("same");
}

void DrawResultPage(TCanvas &canvas,
                    const MorphResult &result,
                    const std::vector<McResidualEvent> &mc_events,
                    const std::vector<DataResidualEvent> &data_events,
                    const char *output_pdf,
                    int nbins,
                    double residual_xmin,
                    double residual_xmax,
                    double split_xmin,
                    double split_xmax) {
  canvas.Clear();
  canvas.Divide(2, 1);

  const std::string axis = AxisLabel(result.axis);
  const std::string bin_title = AnalysisBinTitle(result.analysis_bin,
                                                 result.axis);
  const TString page_title =
    Form("%s, %s: R_{MC}=%.3f mm, #delta_{model}=%.3f mm",
         axis.c_str(),
         bin_title.c_str(),
         result.rec_true_mc_before.sigma68,
         result.delta_syst_model);

  canvas.cd(1);
  {
    const auto data_split = MakeDataSplitValues(data_events);
    const auto mc_split =
      MakeMcSplitValues(mc_events,
                        result.median_p,
                        result.median_d,
                        1.0);

    std::vector<PlotSample> samples = {
      {"MC", mc_split, kBlue + 1, 1, 1, false, 20, 0.8},
      {"Data", data_split, kBlack, 1, 1, true, 20, 0.8}
    };

    DrawOverlayHistograms(samples,
                          Form("%s, PM-DWG split", page_title.Data()),
                          Form("%s_{PM-only}-%s_{DWG-only} [mm]",
                               axis.c_str(),
                               axis.c_str()),
                          nbins,
                          split_xmin,
                          split_xmax);
  }

  canvas.cd(2);
  {
    const auto data_rec_ext = MakeDataRecExtValues(data_events);
    const auto mc_rec_ext =
      MakeMcRecExtValues(mc_events,
                         result.median_r,
                         result.median_e,
                         1.0,
                         1.0);

    std::vector<PlotSample> samples = {
      {"MC", mc_rec_ext, kBlue + 1, 1, 1, false, 20, 0.8},
      {"Data", data_rec_ext, kBlack, 1, 1, true, 20, 0.8}
    };

    DrawOverlayHistograms(samples,
                          Form("%s, rec-ext", page_title.Data()),
                          Form("%s_{rec}-%s_{ext} [mm]",
                               axis.c_str(),
                               axis.c_str()),
                          nbins,
                          residual_xmin,
                          residual_xmax);
  }

  canvas.Print(output_pdf);
}

void DrawSummaryPage(TCanvas &canvas,
                     const std::vector<MorphResult> &results,
                     const char *output_pdf) {
  canvas.Clear();

  auto *text = new TPaveText(0.04, 0.04, 0.96, 0.96, "NDC");
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->SetTextFont(82);
  text->SetTextSize(0.020);

  text->AddText("FROST position resolution: MC truth nominal and model systematics");
  text->AddText("delta_alpha = alpha-envelope for rec-ext closure/external propagation");
  text->AddText("delta_corr  = R-E shuffle correlation uncertainty");
  text->AddText("delta_model = sqrt(delta_alpha^2 + delta_corr^2)");
  text->AddText("");
  text->AddText("axis bin        nMC     nData   R_MC   d_alpha alpha* R_alpha  d_corr aC Rpair  Rshuf  d_model valid");
  text->AddText("------------------------------------------------------------------------------------------------");

  for (const auto &result : results) {
    const std::string line =
      Form("%-4s %-9s %7lld %8lld %6.3f %7.3f %5.1f %7.3f %7.3f %3.1f %6.3f %6.3f %7.3f   %d",
           AxisLabel(result.axis),
           AnalysisBinLabel(result.analysis_bin).c_str(),
           static_cast<long long>(result.n_mc),
           static_cast<long long>(result.n_data),
           result.rec_true_mc_before.sigma68,
           result.delta_ext_model,
           result.alpha_max.alpha,
           result.alpha_max.rec_true_after.sigma68,
           result.delta_rec_model,
           result.rec_shuffle.alpha,
           result.rec_shuffle.reference_rec_true_sigma68,
           result.rec_shuffle.rec_true_after_mean,
           result.delta_syst_model,
           result.valid ? 1 : 0);
    text->AddText(line.c_str());
  }

  text->Draw();
  canvas.Print(output_pdf);
}

void DrawValidationSummaryPage(TCanvas &canvas,
                               const std::vector<MorphResult> &results,
                               const char *output_pdf) {
  canvas.Clear();

  auto *text = new TPaveText(0.04, 0.04, 0.96, 0.96, "NDC");
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->SetTextFont(82);
  text->SetTextSize(0.019);

  text->AddText("Data/MC validation using un-morphed MC");
  text->AddText("sigma68 differences are Data - MC.  rel. differences are (Data - MC)/MC.");
  text->AddText("");
  text->AddText("axis bin       S_data S_MC  dS_rel  O_data O_MC  dO_rel  split chi2/ndf  rec-ext chi2/ndf");
  text->AddText("---------------------------------------------------------------------------------------------");

  for (const auto &result : results) {
    const std::string split_chi2_ndf =
      result.split_chi2.valid ? FormatDouble(result.split_chi2.chi2_ndf, 2)
                              : "nan";
    const std::string rec_ext_chi2_ndf =
      result.rec_ext_chi2.valid ? FormatDouble(result.rec_ext_chi2.chi2_ndf, 2)
                                : "nan";

    const std::string line =
      Form("%-4s %-9s %6.2f %6.2f %7.1f%% %6.2f %6.2f %7.1f%% %9s/%-3d %13s/%-3d",
           AxisLabel(result.axis),
           AnalysisBinLabel(result.analysis_bin).c_str(),
           result.split_data.sigma68,
           result.split_mc_before.sigma68,
           100. * RelativeDifference(result.split_data.sigma68,
                                     result.split_mc_before.sigma68),
           result.rec_ext_data.sigma68,
           result.rec_ext_mc_before.sigma68,
           100. * RelativeDifference(result.rec_ext_data.sigma68,
                                     result.rec_ext_mc_before.sigma68),
           split_chi2_ndf.c_str(),
           result.split_chi2.ndf,
           rec_ext_chi2_ndf.c_str(),
           result.rec_ext_chi2.ndf);
    text->AddText(line.c_str());
  }

  text->Draw();
  canvas.Print(output_pdf);
}

void DrawAlphaScanPage(TCanvas &canvas,
                       const std::vector<MorphResult> &results,
                       const char *output_pdf,
                       bool close_pdf) {
  canvas.Clear();

  TString pdf_name = output_pdf;
  if (close_pdf) {
    pdf_name += ")";
  }

  auto *text = new TPaveText(0.04, 0.04, 0.96, 0.96, "NDC");
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->SetTextFont(82);
  text->SetTextSize(0.017);

  text->AddText("Final table: alpha scan stress-test details");
  text->AddText("Each cell is alpha:value.  R_alpha is sigma68(x_rec' - x_true) [mm].");
  text->AddText("");

  for (const auto &result : results) {
    const std::string title =
      Form("%s %-9s  R_MC=%.3f  delta_alpha=%.3f  alpha*=%.1f  delta_corr=%.3f  delta_model=%.3f",
           AxisLabel(result.axis),
           AnalysisBinLabel(result.analysis_bin).c_str(),
           result.rec_true_mc_before.sigma68,
           result.delta_ext_model,
           result.alpha_max.alpha,
           result.delta_rec_model,
           result.delta_syst_model);
    text->AddText(title.c_str());
    text->AddText(Form("  R_alpha: %s",
                       FormatAlphaScanValues(result.alpha_scan, 'r', 3).c_str()));
    text->AddText(Form("  d_alpha: %s",
                       FormatAlphaScanValues(result.alpha_scan, 'd', 3).c_str()));
    text->AddText(Form("  R_shuffle: %s",
                       FormatRecShuffleScanValues(result.rec_shuffle_scan, 'r', 3).c_str()));
    text->AddText(Form("  d_corr:   %s",
                       FormatRecShuffleScanValues(result.rec_shuffle_scan, 'd', 3).c_str()));
    text->AddText("");
  }

  text->Draw();
  canvas.Print(pdf_name.Data());
}

const MorphResult *FindResult(const std::vector<MorphResult> &results,
                              Axis axis,
                              int analysis_bin) {
  for (const auto &result : results) {
    if (result.axis == axis && result.analysis_bin == analysis_bin) {
      return &result;
    }
  }
  return nullptr;
}

const ExtModelVariation *FindAlphaVariation(const MorphResult &result,
                                            double target_alpha) {
  for (const auto &variation : result.alpha_scan) {
    if (variation.valid && std::abs(variation.alpha - target_alpha) < 1.0e-6) {
      return &variation;
    }
  }
  return nullptr;
}

void PrintDiagnosticPage(TCanvas &canvas,
                         const char *output_pdf,
                         bool close_pdf) {
  TString pdf_name = output_pdf;
  if (close_pdf) {
    pdf_name += ")";
  }
  canvas.Print(pdf_name.Data());
}

void DrawMissingDiagnosticPage(TCanvas &canvas,
                               const char *output_pdf,
                               const char *message,
                               bool close_pdf) {
  canvas.Clear();
  auto *text = new TPaveText(0.10, 0.30, 0.90, 0.70, "NDC");
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(22);
  text->SetTextSize(0.035);
  text->AddText(message);
  text->Draw();
  PrintDiagnosticPage(canvas, output_pdf, close_pdf);
}

void DrawAlphaOneAllAngleDiagnosticPages(
    TCanvas &canvas,
    const std::vector<MorphResult> &results,
    const ResidualSamples &samples,
    const char *output_pdf,
    int nbins,
    double residual_xmin,
    double residual_xmax,
    double split_xmin,
    double split_xmax,
    bool close_pdf) {
  constexpr double kRecTrueDiagnosticMin = -10.0;
  constexpr double kRecTrueDiagnosticMax = 10.0;

  int page_index = 0;
  constexpr int kTotalPages = 6;

  for (const auto axis : {Axis::kX, Axis::kY}) {
    const bool is_x = axis == Axis::kX;
    const auto *result = FindResult(results, axis, 0);
    const std::vector<McResidualEvent> &mc_events =
      is_x ? samples.mc_x.at(0) : samples.mc_y.at(0);
    const std::vector<DataResidualEvent> &data_events =
      is_x ? samples.data_x.at(0) : samples.data_y.at(0);

    const ExtModelVariation *alpha_one =
      result ? FindAlphaVariation(*result, 1.0) : nullptr;

    auto should_close = [&]() -> bool {
      ++page_index;
      return close_pdf && page_index == kTotalPages;
    };

    if (!result || !alpha_one) {
      DrawMissingDiagnosticPage(canvas,
                                output_pdf,
                                Form("Missing alpha=1 all-angle result for %s",
                                     AxisLabel(axis)),
                                should_close());
      DrawMissingDiagnosticPage(canvas,
                                output_pdf,
                                Form("Missing alpha=1 all-angle result for %s",
                                     AxisLabel(axis)),
                                should_close());
      DrawMissingDiagnosticPage(canvas,
                                output_pdf,
                                Form("Missing alpha=1 all-angle result for %s",
                                     AxisLabel(axis)),
                                should_close());
      continue;
    }

    const std::string axis_label = AxisLabel(axis);

    canvas.Clear();
    {
      const auto data_split = MakeDataSplitValues(data_events);
      const auto mc_split_before =
        MakeMcSplitValues(mc_events,
                          result->median_p,
                          result->median_d,
                          1.0);
      const auto mc_split_after =
        MakeMcSplitValues(mc_events,
                          result->median_p,
                          result->median_d,
                          result->k_ext);

      std::vector<PlotSample> plot_samples = {
        {"MC before k_{ext}", mc_split_before, kBlue + 1, 1, 1, false, 20, 0.8},
        {"MC after k_{ext}", mc_split_after, kRed + 1, 1, 1, false, 20, 0.8},
        {"Data", data_split, kBlack, 1, 1, true, 20, 0.8}
      };

      DrawOverlayHistograms(plot_samples,
                            Form("%s all, alpha=1 diagnostic: PM-DWG split",
                                 axis_label.c_str()),
                            Form("%s_{PM-only}-%s_{DWG-only} [mm]",
                                 axis_label.c_str(),
                                 axis_label.c_str()),
                            nbins,
                            split_xmin,
                            split_xmax);
      PrintDiagnosticPage(canvas, output_pdf, should_close());
    }

    canvas.Clear();
    {
      const auto data_rec_ext = MakeDataRecExtValues(data_events);
      const auto mc_rec_ext_before =
        MakeMcRecExtValues(mc_events,
                           result->median_r,
                           result->median_e,
                           1.0,
                           1.0);
      const auto mc_rec_ext_ext_only =
        MakeMcRecExtValues(mc_events,
                           result->median_r,
                           result->median_e,
                           1.0,
                           alpha_one->k_ext_combined);
      const auto mc_rec_ext_after =
        MakeMcRecExtValues(mc_events,
                           result->median_r,
                           result->median_e,
                           alpha_one->k_rec,
                           alpha_one->k_ext_combined);

      std::vector<PlotSample> plot_samples = {
        {"MC before corr.", mc_rec_ext_before, kBlue + 1, 1, 1, false, 20, 0.8},
        {"MC ext corr. only", mc_rec_ext_ext_only, kOrange + 7, 1, 1, false, 20, 0.8},
        {"MC ext+rec corr.", mc_rec_ext_after, kRed + 1, 1, 1, false, 20, 0.8},
        {"Data", data_rec_ext, kBlack, 1, 1, true, 20, 0.8}
      };

      DrawOverlayHistograms(plot_samples,
                            Form("%s all, alpha=1 diagnostic: rec-ext",
                                 axis_label.c_str()),
                            Form("%s_{rec}-%s_{ext} [mm]",
                                 axis_label.c_str(),
                                 axis_label.c_str()),
                            nbins,
                            residual_xmin,
                            residual_xmax);
      PrintDiagnosticPage(canvas, output_pdf, should_close());
    }

    canvas.Clear();
    {
      const auto mc_rec_true_before =
        MakeMcRecTrueValues(mc_events,
                            result->median_r,
                            1.0);
      const auto mc_rec_true_after =
        MakeMcRecTrueValues(mc_events,
                            result->median_r,
                            alpha_one->k_rec);

      std::vector<PlotSample> plot_samples = {
        {"MC before k_{rec}", mc_rec_true_before, kBlue + 1, 1, 1, false, 20, 0.8},
        {"MC after k_{rec}", mc_rec_true_after, kRed + 1, 1, 1, false, 20, 0.8}
      };

      DrawOverlayHistograms(plot_samples,
                            Form("%s all, alpha=1 diagnostic: rec-true",
                                 axis_label.c_str()),
                            Form("%s_{rec}-%s_{true} [mm]",
                                 axis_label.c_str(),
                                 axis_label.c_str()),
                            nbins,
                            kRecTrueDiagnosticMin,
                            kRecTrueDiagnosticMax);
      PrintDiagnosticPage(canvas, output_pdf, should_close());
    }
  }
}

}  // namespace

void CalculatePositionResolution(
    const char *mc_input_dirs="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC_retry/6-TrackMatch_externalfit_PMandDWG,/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/6-TrackMatch_externalfit_PMandDWG",
    const char *data_input_dirs="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM/2-rootfile_after_TrackMatch_externalfit_PMandDWG_shift2",
    const char *output_pdf="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM/2-rootfile_after_TrackMatch_externalfit_PMandDWG_shift2/analysisplot/FROST_resolution_model_systematics.pdf",
    const char *output_log="",
    int nbins = 200,
    double residual_xmin_mm = -50.,
    double residual_xmax_mm = 50.,
    bool use_position_reweighting = true,
    double k_min = 0.0,
    double k_max = 3.0,
    int n_shuffle_trials = 50,
    unsigned int shuffle_seed = 12345) {
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);

  const std::vector<std::string> mc_input_dir_list =
    ParseInputDirectories(mc_input_dirs);
  const std::vector<std::string> data_input_dir_list =
    ParseInputDirectories(data_input_dirs);

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
      CountPositionBinsOneFile(file_path,
                               SampleType::kMonteCarlo,
                               mc_position_weights);
    }
    mc_position_weights.CalculateWeights();
    PrintPositionWeightSummary("MC", mc_position_weights);

    for (const auto &file_path : data_root_files) {
      CountPositionBinsOneFile(file_path,
                               SampleType::kData,
                               data_position_weights);
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
    LoadResidualSamplesOneFile(file_path,
                               SampleType::kMonteCarlo,
                               mc_position_weights,
                               samples);
  }

  for (const auto &file_path : data_root_files) {
    std::cout << "Loading data: " << file_path << std::endl;
    LoadResidualSamplesOneFile(file_path,
                               SampleType::kData,
                               data_position_weights,
                               samples);
  }

  const double split_xmin_mm = 2.0 * residual_xmin_mm;
  const double split_xmax_mm = 2.0 * residual_xmax_mm;

  std::vector<MorphResult> results;
  results.reserve(2 * kNAnalysisBins);

  for (int analysis_bin = 0; analysis_bin < kNAnalysisBins; ++analysis_bin) {
    results.push_back(EstimateMorphResult(Axis::kX,
                                          analysis_bin,
                                          samples.mc_x.at(analysis_bin),
                                          samples.data_x.at(analysis_bin),
                                          nbins,
                                          residual_xmin_mm,
                                          residual_xmax_mm,
                                          split_xmin_mm,
                                          split_xmax_mm,
                                          k_min,
                                          k_max,
                                          n_shuffle_trials,
                                          shuffle_seed));
  }

  for (int analysis_bin = 0; analysis_bin < kNAnalysisBins; ++analysis_bin) {
    results.push_back(EstimateMorphResult(Axis::kY,
                                          analysis_bin,
                                          samples.mc_y.at(analysis_bin),
                                          samples.data_y.at(analysis_bin),
                                          nbins,
                                          residual_xmin_mm,
                                          residual_xmax_mm,
                                          split_xmin_mm,
                                          split_xmax_mm,
                                          k_min,
                                          k_max,
                                          n_shuffle_trials,
                                          shuffle_seed));
  }

  const TString output_log_file = MakeOutputLogFileName(output_pdf, output_log);
  std::vector<std::string> log_lines;
  WriteLogFile(output_log_file.Data(), results, log_lines);
  std::cout << "Wrote log: " << output_log_file.Data() << std::endl;

  const TString output_root_log_file =
    MakeOutputRootLogFileName(output_pdf, output_log);
  WriteRootLogFile(output_root_log_file.Data(), results);
  std::cout << "Wrote ROOT log: " << output_root_log_file.Data() << std::endl;

  TCanvas canvas("canvas",
                 "FROST resolution model systematics",
                 1400,
                 650);

  TString pdf_open = output_pdf;
  pdf_open += "(";

  DrawSummaryPage(canvas, results, pdf_open.Data());
  DrawValidationSummaryPage(canvas, results, output_pdf);

  for (const auto &result : results) {
    if (!result.valid) {
      continue;
    }

    if (result.axis == Axis::kX) {
      DrawResultPage(canvas,
                     result,
                     samples.mc_x.at(result.analysis_bin),
                     samples.data_x.at(result.analysis_bin),
                     output_pdf,
                     nbins,
                     residual_xmin_mm,
                     residual_xmax_mm,
                     split_xmin_mm,
                     split_xmax_mm);
    } else {
      DrawResultPage(canvas,
                     result,
                     samples.mc_y.at(result.analysis_bin),
                     samples.data_y.at(result.analysis_bin),
                     output_pdf,
                     nbins,
                     residual_xmin_mm,
                     residual_xmax_mm,
                     split_xmin_mm,
                     split_xmax_mm);
    }
  }

  DrawAlphaScanPage(canvas, results, output_pdf, false);
  DrawAlphaOneAllAngleDiagnosticPages(canvas,
                                       results,
                                       samples,
                                       output_pdf,
                                       nbins,
                                       residual_xmin_mm,
                                       residual_xmax_mm,
                                       split_xmin_mm,
                                       split_xmax_mm,
                                       true);

  std::cout << "Wrote PDF: " << output_pdf << std::endl;
}
