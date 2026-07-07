// DrawPMDWGMinusFROSTComparison.C
//
// Diagnostic macro to compare PM-only - FROST-rec and DWG-only - FROST-rec
// residuals between MC and data, separately for x/y and angle bins.
// It also compares PM-only coordinate phases folded modulo 10 mm and 25 mm.
//
// Usage in ROOT:
//   root -l -b -q 'DrawPMDWGMinusFROSTComparison.C("mc_dir", "data_dir", "out.pdf")'
//
// Multiple input directories can be given as a comma-separated list.
// This macro intentionally keeps the functionality minimal for quick checks.

#include <TCanvas.h>
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
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr int NHIT_DWG_X = 4;
constexpr int NHIT_PM_X = 10;
constexpr int NHIT_DWG_Y = 4;
constexpr int NHIT_PM_Y = 10;

constexpr double ACCEPTANCE_X = 560.0; // mm
constexpr double ACCEPTANCE_Y = 600.0; // mm
constexpr double kNonInitializedThreshold = -9.0e7;

// Angle bins [deg]: [0,5), [5,10), [10,20), [20,50)
constexpr int kNAngleBins = 4;
const double kAngleBins[kNAngleBins + 1] = {
  0.0, 5.0, 10.0, 20.0, 50.0
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
    std::cerr << "Warning: cannot list directory: " << input_dir << std::endl;
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

void AddValue(std::array<std::vector<double>, kNAnalysisBins> &bins,
              int angle_bin,
              double value) {
  if (!std::isfinite(value)) {
    return;
  }
  bins.at(0).push_back(value);
  if (angle_bin >= 0) {
    bins.at(angle_bin + 1).push_back(value);
  }
}

struct ResidualSamples {
  std::array<std::vector<double>, kNAnalysisBins> pm_minus_rec_x;
  std::array<std::vector<double>, kNAnalysisBins> dwg_minus_rec_x;
  std::array<std::vector<double>, kNAnalysisBins> pm_minus_rec_y;
  std::array<std::vector<double>, kNAnalysisBins> dwg_minus_rec_y;

  std::array<std::vector<double>, kNAnalysisBins> pm_only_x;
  std::array<std::vector<double>, kNAnalysisBins> pm_only_y;
};

void LoadResidualSamplesOneFile(const std::string &file_path,
                                SampleType sample_type,
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

  std::vector<int> *has_match = nullptr;
  std::vector<int> *n_dwg_x = nullptr;
  std::vector<int> *n_dwg_y = nullptr;
  std::vector<int> *n_pm_x = nullptr;
  std::vector<int> *n_pm_y = nullptr;
  std::vector<int> *ninja_track_type = nullptr;
  std::vector<double> *external_tangent_x = nullptr;
  std::vector<double> *external_tangent_y = nullptr;
  std::vector<double> *pm_only_expected_x = nullptr;
  std::vector<double> *pm_only_expected_y = nullptr;
  std::vector<double> *dwg_only_expected_x = nullptr;
  std::vector<double> *dwg_only_expected_y = nullptr;
  std::vector<double> *frost_nearest_x = nullptr;
  std::vector<double> *frost_nearest_y = nullptr;

  Int_t bsd_good_spill_flag = 0;
  Int_t detector_flags[8] = {};

  bool ok = true;
  ok &= SetVectorBranchAddress(tree, "trackmatch_has_match", &has_match);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_num_planes_downstream_wagasci_x", &n_dwg_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_num_planes_downstream_wagasci_y", &n_dwg_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_num_planes_proton_module_x", &n_pm_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_num_planes_proton_module_y", &n_pm_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_ninja_track_type", &ninja_track_type);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_tangent_x", &external_tangent_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_tangent_y", &external_tangent_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_pm_only_expected_x", &pm_only_expected_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_pm_only_expected_y", &pm_only_expected_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_dwg_only_expected_x", &dwg_only_expected_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_external_dwg_only_expected_y", &dwg_only_expected_y);
  ok &= SetVectorBranchAddress(tree, "trackmatch_frost_nearest_x", &frost_nearest_x);
  ok &= SetVectorBranchAddress(tree, "trackmatch_frost_nearest_y", &frost_nearest_y);

  if (sample_type == SampleType::kData) {
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

    if (sample_type == SampleType::kData) {
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
      if (!HasIndex(ninja_track_type, itrack) || ninja_track_type->at(itrack) != 1) {
        continue;
      }
      if (!TrackPassesCommonSelection(itrack, n_dwg_x, n_dwg_y, n_pm_x, n_pm_y)) {
        continue;
      }
      if (!HasIndex(frost_nearest_x, itrack) ||
          !HasIndex(frost_nearest_y, itrack) ||
          !HasIndex(external_tangent_x, itrack) ||
          !HasIndex(external_tangent_y, itrack) ||
          !IsValidValue(frost_nearest_x->at(itrack)) ||
          !IsValidValue(frost_nearest_y->at(itrack)) ||
          !IsValidValue(external_tangent_x->at(itrack)) ||
          !IsValidValue(external_tangent_y->at(itrack))) {
        continue;
      }

      const double frost_x = frost_nearest_x->at(itrack);
      const double frost_y = frost_nearest_y->at(itrack);
      if (std::abs(frost_x) >= ACCEPTANCE_X || std::abs(frost_y) >= ACCEPTANCE_Y) {
        continue;
      }

      const double theta_x_deg = std::atan(std::abs(external_tangent_x->at(itrack))) * 180.0 / TMath::Pi();
      const double theta_y_deg = std::atan(std::abs(external_tangent_y->at(itrack))) * 180.0 / TMath::Pi();
      const int angle_bin_x = FindAngleBin(theta_x_deg);
      const int angle_bin_y = FindAngleBin(theta_y_deg);

      if (HasIndex(pm_only_expected_x, itrack) &&
          HasIndex(dwg_only_expected_x, itrack) &&
          IsValidValue(pm_only_expected_x->at(itrack)) &&
          IsValidValue(dwg_only_expected_x->at(itrack))) {
        AddValue(samples.pm_minus_rec_x,
                 angle_bin_x,
                 pm_only_expected_x->at(itrack) - frost_x);
        AddValue(samples.dwg_minus_rec_x,
                 angle_bin_x,
                 dwg_only_expected_x->at(itrack) - frost_x);
        AddValue(samples.pm_only_x,
                 angle_bin_x,
                 pm_only_expected_x->at(itrack));
      }

      if (HasIndex(pm_only_expected_y, itrack) &&
          HasIndex(dwg_only_expected_y, itrack) &&
          IsValidValue(pm_only_expected_y->at(itrack)) &&
          IsValidValue(dwg_only_expected_y->at(itrack))) {
        AddValue(samples.pm_minus_rec_y,
                 angle_bin_y,
                 pm_only_expected_y->at(itrack) - frost_y);
        AddValue(samples.dwg_minus_rec_y,
                 angle_bin_y,
                 dwg_only_expected_y->at(itrack) - frost_y);
        AddValue(samples.pm_only_y,
                 angle_bin_y,
                 pm_only_expected_y->at(itrack));
      }
    }
  }
}

double Mean(const std::vector<double> &values) {
  if (values.empty()) {
    return 0.;
  }
  double sum = 0.;
  for (double value : values) {
    sum += value;
  }
  return sum / static_cast<double>(values.size());
}

double Rms(const std::vector<double> &values, double mean) {
  if (values.empty()) {
    return 0.;
  }
  double sum2 = 0.;
  for (double value : values) {
    const double diff = value - mean;
    sum2 += diff * diff;
  }
  return std::sqrt(sum2 / static_cast<double>(values.size()));
}


double PositiveModulo(double value, double period) {
  if (!std::isfinite(value) || period <= 0.) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double mod_value = std::fmod(value, period);
  if (mod_value < 0.) {
    mod_value += period;
  }
  return mod_value;
}

std::vector<double> MakeModuloValues(const std::vector<double> &values,
                                     double period) {
  std::vector<double> modulo_values;
  modulo_values.reserve(values.size());

  for (double value : values) {
    const double mod_value = PositiveModulo(value, period);
    if (std::isfinite(mod_value)) {
      modulo_values.push_back(mod_value);
    }
  }

  return modulo_values;
}

void DrawDataMcOverlay(const std::vector<double> &data_values,
                       const std::vector<double> &mc_values,
                       const char *title,
                       const char *x_axis_title,
                       int nbins,
                       double xmin,
                       double xmax) {
  static int hist_counter = 0;

  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.13);

  TString h_mc_name = Form("h_mc_%d", hist_counter);
  TString h_data_name = Form("h_data_%d", hist_counter);
  ++hist_counter;

  TH1D *h_mc = new TH1D(h_mc_name,
                        Form("%s;%s;Number of tracks", title, x_axis_title),
                        nbins,
                        xmin,
                        xmax);
  TH1D *h_data = new TH1D(h_data_name,
                          Form("%s;%s;Number of tracks", title, x_axis_title),
                          nbins,
                          xmin,
                          xmax);
  h_mc->SetDirectory(nullptr);
  h_data->SetDirectory(nullptr);
  h_mc->Sumw2();
  h_data->Sumw2();

  for (double value : mc_values) {
    if (std::isfinite(value)) {
      h_mc->Fill(value);
    }
  }
  for (double value : data_values) {
    if (std::isfinite(value)) {
      h_data->Fill(value);
    }
  }

  const double data_integral = h_data->Integral();
  const double mc_integral = h_mc->Integral();
  if (data_integral > 0. && mc_integral > 0.) {
    h_mc->Scale(data_integral / mc_integral);
  }

  h_mc->SetLineColor(kBlue + 1);
  h_mc->SetLineWidth(2);
  h_data->SetMarkerColor(kBlack);
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(0.8);
  h_data->SetLineColor(kBlack);

  const double max_y = std::max(h_mc->GetMaximum(), h_data->GetMaximum());
  h_mc->SetMaximum(max_y > 0. ? 1.35 * max_y : 1.0);
  h_mc->GetYaxis()->SetTitleOffset(1.45);
  h_mc->Draw("hist");
  h_data->Draw("E1 same");

  const double data_mean = Mean(data_values);
  const double mc_mean = Mean(mc_values);
  const double data_rms = Rms(data_values, data_mean);
  const double mc_rms = Rms(mc_values, mc_mean);

  auto *legend = new TLegend(0.55, 0.70, 0.94, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.032);
  legend->AddEntry(h_data,
                   Form("Data: N=%lu, mean=%.2f, RMS=%.2f",
                        static_cast<unsigned long>(data_values.size()),
                        data_mean,
                        data_rms),
                   "lep");
  legend->AddEntry(h_mc,
                   Form("MC: N=%lu, mean=%.2f, RMS=%.2f",
                        static_cast<unsigned long>(mc_values.size()),
                        mc_mean,
                        mc_rms),
                   "l");
  legend->Draw("same");
}

void DrawPage(TCanvas &canvas,
              Axis axis,
              int analysis_bin,
              const ResidualSamples &data,
              const ResidualSamples &mc,
              const char *output_pdf,
              int nbins,
              double xmin,
              double xmax) {
  canvas.Clear();
  canvas.Divide(2, 1);

  const std::string axis_label = AxisLabel(axis);
  const std::string bin_title = AnalysisBinTitle(analysis_bin, axis);

  const auto &data_pm = axis == Axis::kX ? data.pm_minus_rec_x.at(analysis_bin)
                                         : data.pm_minus_rec_y.at(analysis_bin);
  const auto &mc_pm = axis == Axis::kX ? mc.pm_minus_rec_x.at(analysis_bin)
                                       : mc.pm_minus_rec_y.at(analysis_bin);
  const auto &data_dwg = axis == Axis::kX ? data.dwg_minus_rec_x.at(analysis_bin)
                                          : data.dwg_minus_rec_y.at(analysis_bin);
  const auto &mc_dwg = axis == Axis::kX ? mc.dwg_minus_rec_x.at(analysis_bin)
                                        : mc.dwg_minus_rec_y.at(analysis_bin);

  canvas.cd(1);
  DrawDataMcOverlay(data_pm,
                    mc_pm,
                    Form("%s, %s: PM-only - FROST-rec", axis_label.c_str(), bin_title.c_str()),
                    Form("%s_{PM-only} - %s_{rec} [mm]", axis_label.c_str(), axis_label.c_str()),
                    nbins,
                    xmin,
                    xmax);

  canvas.cd(2);
  DrawDataMcOverlay(data_dwg,
                    mc_dwg,
                    Form("%s, %s: DWG-only - FROST-rec", axis_label.c_str(), bin_title.c_str()),
                    Form("%s_{DWG-only} - %s_{rec} [mm]", axis_label.c_str(), axis_label.c_str()),
                    nbins,
                    xmin,
                    xmax);

  canvas.Print(output_pdf);
}

void DrawPmModuloPage(TCanvas &canvas,
                      Axis axis,
                      int analysis_bin,
                      const ResidualSamples &data,
                      const ResidualSamples &mc,
                      const char *output_pdf,
                      int nbins_mod10,
                      int nbins_mod25) {
  canvas.Clear();
  canvas.Divide(2, 1);

  const std::string axis_label = AxisLabel(axis);
  const std::string bin_title = AnalysisBinTitle(analysis_bin, axis);

  const auto &data_pm = axis == Axis::kX ? data.pm_only_x.at(analysis_bin)
                                         : data.pm_only_y.at(analysis_bin);
  const auto &mc_pm = axis == Axis::kX ? mc.pm_only_x.at(analysis_bin)
                                       : mc.pm_only_y.at(analysis_bin);

  const auto data_pm_mod10 = MakeModuloValues(data_pm, 10.0);
  const auto mc_pm_mod10 = MakeModuloValues(mc_pm, 10.0);
  const auto data_pm_mod25 = MakeModuloValues(data_pm, 25.0);
  const auto mc_pm_mod25 = MakeModuloValues(mc_pm, 25.0);

  canvas.cd(1);
  DrawDataMcOverlay(data_pm_mod10,
                    mc_pm_mod10,
                    Form("%s, %s: PM-only mod 10 mm", axis_label.c_str(), bin_title.c_str()),
                    Form("%s_{PM-only} mod 10 mm [mm]", axis_label.c_str()),
                    nbins_mod10,
                    0.0,
                    10.0);

  canvas.cd(2);
  DrawDataMcOverlay(data_pm_mod25,
                    mc_pm_mod25,
                    Form("%s, %s: PM-only mod 25 mm", axis_label.c_str(), bin_title.c_str()),
                    Form("%s_{PM-only} mod 25 mm [mm]", axis_label.c_str()),
                    nbins_mod25,
                    0.0,
                    25.0);

  canvas.Print(output_pdf);
}

void DrawCoverPage(TCanvas &canvas,
                   const char *output_pdf,
                   const std::vector<std::string> &mc_files,
                   const std::vector<std::string> &data_files,
                   int nbins,
                   double xmin,
                   double xmax) {
  canvas.Clear();
  auto *text = new TPaveText(0.06, 0.08, 0.94, 0.92, "NDC");
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->SetTextFont(42);
  text->SetTextSize(0.032);

  text->AddText("PM-only / DWG-only residuals relative to FROST-rec");
  text->AddText("Data and MC are compared with MC scaled to the Data integral in each plot.");
  text->AddText("Angle bins: all, 0-5, 5-10, 10-20, 20-50 deg for each projection.");
  text->AddText("Additional pages show PM-only coordinate phases folded modulo 10 mm and 25 mm.");
  text->AddText("Selection follows the common PM/DWG hit-count selection used in CalculatePositionResolution.C.");
  text->AddText("");
  text->AddText(Form("MC files: %lu", static_cast<unsigned long>(mc_files.size())));
  text->AddText(Form("Data files: %lu", static_cast<unsigned long>(data_files.size())));
  text->AddText(Form("Histogram range: %.1f to %.1f mm, nbins=%d", xmin, xmax, nbins));
  text->Draw();
  canvas.Print(output_pdf);
}

} // namespace

void DrawPMDWGMinusFROSTComparison(
    const char *mc_input_dirs="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC_retry/6-TrackMatch_externalfit_PMandDWG,/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/6-TrackMatch_externalfit_PMandDWG",
    const char *data_input_dirs="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM/2-rootfile_after_TrackMatch_externalfit_PMandDWG_shift2",
    const char *output_pdf="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM/2-rootfile_after_TrackMatch_externalfit_PMandDWG_shift2/analysisplot/PM_DWG_minus_FROST_rec_comparison.pdf",
    int nbins = 200,
    double residual_xmin_mm = -100.,
    double residual_xmax_mm = 100.) {
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);

  const auto mc_input_dir_list = ParseInputDirectories(mc_input_dirs);
  const auto data_input_dir_list = ParseInputDirectories(data_input_dirs);
  const auto mc_root_files = FindRootFiles(mc_input_dir_list);
  const auto data_root_files = FindRootFiles(data_input_dir_list);

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

  ResidualSamples mc_samples;
  ResidualSamples data_samples;

  for (const auto &file_path : mc_root_files) {
    std::cout << "Loading MC: " << file_path << std::endl;
    LoadResidualSamplesOneFile(file_path, SampleType::kMonteCarlo, mc_samples);
  }
  for (const auto &file_path : data_root_files) {
    std::cout << "Loading data: " << file_path << std::endl;
    LoadResidualSamplesOneFile(file_path, SampleType::kData, data_samples);
  }

  TCanvas canvas("canvas", "PM/DWG minus FROST-rec comparison", 1400, 650);

  TString pdf_open = output_pdf;
  pdf_open += "(";
  DrawCoverPage(canvas,
                pdf_open.Data(),
                mc_root_files,
                data_root_files,
                nbins,
                residual_xmin_mm,
                residual_xmax_mm);

  for (int analysis_bin = 0; analysis_bin < kNAnalysisBins; ++analysis_bin) {
    DrawPage(canvas,
             Axis::kX,
             analysis_bin,
             data_samples,
             mc_samples,
             output_pdf,
             nbins,
             residual_xmin_mm,
             residual_xmax_mm);
  }

  for (int analysis_bin = 0; analysis_bin < kNAnalysisBins; ++analysis_bin) {
    DrawPage(canvas,
             Axis::kY,
             analysis_bin,
             data_samples,
             mc_samples,
             output_pdf,
             nbins,
             residual_xmin_mm,
             residual_xmax_mm);
  }

  const int nbins_mod10 = 100;
  const int nbins_mod25 = 125;

  for (int analysis_bin = 0; analysis_bin < kNAnalysisBins; ++analysis_bin) {
    DrawPmModuloPage(canvas,
                     Axis::kX,
                     analysis_bin,
                     data_samples,
                     mc_samples,
                     output_pdf,
                     nbins_mod10,
                     nbins_mod25);
  }

  for (int analysis_bin = 0; analysis_bin < kNAnalysisBins; ++analysis_bin) {
    TString pdf_name = output_pdf;
    if (analysis_bin == kNAnalysisBins - 1) {
      pdf_name += ")";
    }
    DrawPmModuloPage(canvas,
                     Axis::kY,
                     analysis_bin,
                     data_samples,
                     mc_samples,
                     pdf_name.Data(),
                     nbins_mod10,
                     nbins_mod25);
  }

  std::cout << "Wrote PDF: " << output_pdf << std::endl;
}
