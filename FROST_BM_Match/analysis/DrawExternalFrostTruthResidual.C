// DrawExternalFrostTruthResidual.C
//
// Usage in ROOT:
//   root -l -b -q 'DrawExternalFrostTruthResidual.C("input_dir", "output.pdf")'
//
// The macro reads all .root files directly under input_dir, loops over the
// match_info tree, and saves a multi-page PDF.
//
// Residual pages:
//   1. all angles: external expected - true FROST
//   2. all angles: FROST reco - external expected
//   3. all angles: FROST reco - true FROST
//   4-11. angle-binned external expected - true FROST
//   12-19. angle-binned FROST reco - external expected
//   20-27. angle-binned FROST reco - true FROST
//
// Correlation pages for MC:
//   28. all angles:
//       horizontal: external expected - true FROST
//       vertical:   FROST reco - true FROST
//   29-36. same correlation plots in angle bins.
//
// The correlation pages display the ordinary RMS-based covariance and Pearson
// correlation coefficient between
//   E = external expected - true FROST
//   R = FROST reco - true FROST.
//
// Selection:
//   common:
//      trackmatch_has_match == 1
//      trackmatch_ninja_track_type == 1
//      abs(frost_nearest_x) < ACCEPTANCE_X mm
//      abs(frost_nearest_y) < ACCEPTANCE_Y mm
//      trackmatch_external_num_planes_downstream_wagasci_x > NHIT_DWG_X
//      trackmatch_external_num_planes_proton_module_x > NHIT_PM_X
//      trackmatch_external_num_planes_downstream_wagasci_y > NHIT_DWG_Y
//      trackmatch_external_num_planes_proton_module_y > NHIT_PM_Y
//
//   The same common selection is applied to both x and y residuals.

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
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
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
                                       const char *title_prefix,
                                       const char *x_axis_title,
                                       int nbins,
                                       double xmin_mm,
                                       double xmax_mm) {
  std::vector<TH1D*> histograms;
  histograms.reserve(kNAngleBins);

  for (int i = 0; i < kNAngleBins; ++i) {
    histograms.push_back(new TH1D(
      Form("%s_angle_%02d", name_prefix, i),
      Form("%s: %.0f #leq #theta < %.0f deg;%s;Number of tracks",
           title_prefix, kAngleBins[i], kAngleBins[i + 1], x_axis_title),
      nbins,
      xmin_mm,
      xmax_mm));
  }

  return histograms;
}

std::vector<TH2D*> MakeAngleCorrelationHistograms(const char *name_prefix,
                                                  const char *title_prefix,
                                                  int nbins,
                                                  double xmin_mm,
                                                  double xmax_mm) {
  std::vector<TH2D*> histograms;
  histograms.reserve(kNAngleBins);

  for (int i = 0; i < kNAngleBins; ++i) {
    histograms.push_back(new TH2D(
      Form("%s_angle_%02d", name_prefix, i),
      Form("%s: %.0f #leq #theta < %.0f deg;external expected - true FROST [mm];FROST reco - true FROST [mm]",
           title_prefix, kAngleBins[i], kAngleBins[i + 1]),
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
                    TH1D *hist_x,
                    TH1D *hist_y,
                    TH1D *hist_reco_external_x,
                    TH1D *hist_reco_external_y,
                    TH1D *hist_reco_truth_x,
                    TH1D *hist_reco_truth_y,
                    TH2D *hist_correlation_x,
                    TH2D *hist_correlation_y,
                    std::vector<TH1D*> &hist_external_truth_x_by_angle,
                    std::vector<TH1D*> &hist_external_truth_y_by_angle,
                    std::vector<TH1D*> &hist_reco_external_x_by_angle,
                    std::vector<TH1D*> &hist_reco_external_y_by_angle,
                    std::vector<TH1D*> &hist_reco_truth_x_by_angle,
                    std::vector<TH1D*> &hist_reco_truth_y_by_angle,
                    std::vector<TH2D*> &hist_correlation_x_by_angle,
                    std::vector<TH2D*> &hist_correlation_y_by_angle,
                    CorrelationStats &correlation_stats_x,
                    CorrelationStats &correlation_stats_y,
                    std::vector<CorrelationStats> &correlation_stats_x_by_angle,
                    std::vector<CorrelationStats> &correlation_stats_y_by_angle,
                    Long64_t &n_tracks_x,
                    Long64_t &n_tracks_y) {
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
  std::vector<double> *frost_nearest_x = nullptr;
  std::vector<double> *frost_nearest_y = nullptr;
  std::vector<double> *true_frost_x = nullptr;
  std::vector<double> *true_frost_y = nullptr;
  std::vector<double> *baby_mind_tangent_x = nullptr;
  std::vector<double> *baby_mind_tangent_y = nullptr;
  std::vector<int> *ninja_track_type = nullptr;

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
    tree, "trackmatch_frost_nearest_x", &frost_nearest_x);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_frost_nearest_y", &frost_nearest_y);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_true_frost_nearest_position_x", &true_frost_x);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_true_frost_nearest_position_y", &true_frost_y);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_baby_mind_tangent_x", &baby_mind_tangent_x);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_baby_mind_tangent_y", &baby_mind_tangent_y);
  ok &= SetVectorBranchAddress(
    tree, "trackmatch_ninja_track_type", &ninja_track_type);

  if (!ok) {
    std::cerr << "Warning: skip file because required branches are missing: "
              << file_path << std::endl;
    return;
  }

  const Long64_t nentries = tree->GetEntries();
  for (Long64_t entry = 0; entry < nentries; ++entry) {
    tree->GetEntry(entry);

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
      const double angle_deg =
        std::atan(std::sqrt(tangent_x * tangent_x +
                            tangent_y * tangent_y)) * 180.0 / TMath::Pi();

      const int angle_bin = FindAngleBin(angle_deg);
      const bool has_angle_bin = (angle_bin >= 0);

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

      const bool x_base_selected =
        HasIndex(external_expected_x, itrack) &&
        IsValidValue(external_expected_x->at(itrack));

      if (x_base_selected &&
          HasIndex(true_frost_x, itrack) &&
          IsValidValue(true_frost_x->at(itrack))) {
        hist_x->Fill(external_expected_x->at(itrack) -
                     true_frost_x->at(itrack));
        ++n_tracks_x;
        if (has_angle_bin) {
          hist_external_truth_x_by_angle.at(angle_bin)->Fill(
            external_expected_x->at(itrack) - true_frost_x->at(itrack));
        }
      }

      if (x_base_selected &&
          HasIndex(frost_nearest_x, itrack) &&
          IsValidValue(frost_nearest_x->at(itrack))) {
        hist_reco_external_x->Fill(frost_nearest_x->at(itrack) -
                                   external_expected_x->at(itrack));

        if (has_angle_bin) {
          hist_reco_external_x_by_angle.at(angle_bin)->Fill(
            frost_nearest_x->at(itrack) - external_expected_x->at(itrack));
        }
      }

      if (x_base_selected &&
          HasIndex(frost_nearest_x, itrack) &&
          HasIndex(true_frost_x, itrack) &&
          IsValidValue(frost_nearest_x->at(itrack)) &&
          IsValidValue(true_frost_x->at(itrack))) {
        const double external_minus_true_x =
          external_expected_x->at(itrack) - true_frost_x->at(itrack);
        const double reco_minus_true_x =
          frost_nearest_x->at(itrack) - true_frost_x->at(itrack);

        hist_reco_truth_x->Fill(reco_minus_true_x);
        hist_correlation_x->Fill(external_minus_true_x, reco_minus_true_x);
        correlation_stats_x.Fill(external_minus_true_x, reco_minus_true_x);

        if (has_angle_bin) {
          hist_reco_truth_x_by_angle.at(angle_bin)->Fill(reco_minus_true_x);
          hist_correlation_x_by_angle.at(angle_bin)->Fill(
            external_minus_true_x, reco_minus_true_x);
          correlation_stats_x_by_angle.at(angle_bin).Fill(
            external_minus_true_x, reco_minus_true_x);
        }
      }

      const bool y_base_selected =
        HasIndex(external_expected_y, itrack) &&
        IsValidValue(external_expected_y->at(itrack));

      if (y_base_selected &&
          HasIndex(true_frost_y, itrack) &&
          IsValidValue(true_frost_y->at(itrack))) {
        hist_y->Fill(external_expected_y->at(itrack) -
                     true_frost_y->at(itrack));
        ++n_tracks_y;

        if (has_angle_bin) {
          hist_external_truth_y_by_angle.at(angle_bin)->Fill(
            external_expected_y->at(itrack) - true_frost_y->at(itrack));
        }
      }

      if (y_base_selected &&
          HasIndex(frost_nearest_y, itrack) &&
          IsValidValue(frost_nearest_y->at(itrack))) {
        hist_reco_external_y->Fill(frost_nearest_y->at(itrack) -
                                   external_expected_y->at(itrack));

        if (has_angle_bin) {
          hist_reco_external_y_by_angle.at(angle_bin)->Fill(
            frost_nearest_y->at(itrack) - external_expected_y->at(itrack));
        }
      }

      if (y_base_selected &&
          HasIndex(frost_nearest_y, itrack) &&
          HasIndex(true_frost_y, itrack) &&
          IsValidValue(frost_nearest_y->at(itrack)) &&
          IsValidValue(true_frost_y->at(itrack))) {
        const double external_minus_true_y =
          external_expected_y->at(itrack) - true_frost_y->at(itrack);
        const double reco_minus_true_y =
          frost_nearest_y->at(itrack) - true_frost_y->at(itrack);

        hist_reco_truth_y->Fill(reco_minus_true_y);
        hist_correlation_y->Fill(external_minus_true_y, reco_minus_true_y);
        correlation_stats_y.Fill(external_minus_true_y, reco_minus_true_y);

        if (has_angle_bin) {
          hist_reco_truth_y_by_angle.at(angle_bin)->Fill(reco_minus_true_y);
          hist_correlation_y_by_angle.at(angle_bin)->Fill(
            external_minus_true_y, reco_minus_true_y);
          correlation_stats_y_by_angle.at(angle_bin).Fill(
            external_minus_true_y, reco_minus_true_y);
        }
      }
    }
  }
}

bool CalculateCentral68Width(TH1D *hist,
                             double &median,
                             double &sigma68,
                             double &q16,
                             double &q84) {
  if (!hist || hist->GetEntries() <= 0) {
    return false;
  }

  double probabilities[3] = {0.16, 0.50, 0.84};
  double quantiles[3] = {0., 0., 0.};
  hist->GetQuantiles(3, quantiles, probabilities);

  q16 = quantiles[0];
  median = quantiles[1];
  q84 = quantiles[2];
  sigma68 = 0.5 * (q84 - q16);

  return std::isfinite(sigma68);
}

void DrawCentral68Text(TH1D *hist) {
  double median = 0.;
  double sigma68 = 0.;
  double q16 = 0.;
  double q84 = 0.;

  if (!CalculateCentral68Width(hist, median, sigma68, q16, q84)) {
    return;
  }

  auto *text = new TPaveText(0.6, 0.40, 0.85, 0.50, "NDC");
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->SetTextSize(0.035);
  text->AddText(Form("median = %.3f mm", median));
  text->AddText(Form("#sigma_{68} = %.3f mm", sigma68));
  text->Draw("same");
}

void DrawHistogram(TH1D *hist) {
  hist->SetLineWidth(2);
  hist->GetXaxis()->SetTitleOffset(1.15);
  hist->GetYaxis()->SetTitleOffset(1.25);
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

  DrawCentral68Text(hist);
}

void DrawTwoPanelPage(TCanvas &canvas,
                      TH1D *hist_left,
                      TH1D *hist_right) {
  canvas.Clear();
  canvas.Divide(2, 1);

  canvas.cd(1);
  DrawHistogram(hist_left);

  canvas.cd(2);
  DrawHistogram(hist_right);
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
  hist->GetXaxis()->SetTitleOffset(1.15);
  hist->GetYaxis()->SetTitleOffset(1.25);
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
                     bool close_pdf = false) {
  TString pdf_close = output_pdf;
  pdf_close += ")";

  for (int i = 0; i < kNAngleBins; ++i) {
    DrawTwoPanelPage(canvas,
                     hist_x_by_angle.at(i),
                     hist_y_by_angle.at(i));

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

}  // namespace

// void DrawExternalFrostTruthResidual(const char *input_dir="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/6-TrackMatch_externalfit_PMandDWG",
//                                     const char *output_pdf="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/6-TrackMatch_externalfit_PMandDWG/residual_external_true.pdf",
//                                     int nbins = 200,
//                                     double xmin_mm = -50.,
//                                     double xmax_mm = 50.) {
void DrawExternalFrostTruthResidual(const char *input_dir="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM/2-rootfile_after_TrackMatch_externalfit_PMandDWG",
                                    const char *output_pdf="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM/2-rootfile_after_TrackMatch_externalfit_PMandDWG/residual_rec_external.pdf",
                                    int nbins = 200,
                                    double xmin_mm = -50.,
                                    double xmax_mm = 50.) {
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1111);

  const std::vector<std::string> root_files = FindRootFiles(input_dir);
  if (root_files.empty()) {
    std::cerr << "Error: no .root files found in directory: "
              << input_dir << std::endl;
    return;
  }

  auto *hist_x = new TH1D(
    "hist_external_truth_residual_x",
    "External-track prediction minus FROST truth;external expected x - true FROST x [mm];Number of tracks",
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_y = new TH1D(
    "hist_external_truth_residual_y",
    "External-track prediction minus FROST truth;external expected y - true FROST y [mm];Number of tracks",
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_reco_external_x = new TH1D(
    "hist_reco_external_residual_x",
    "FROST reco minus external-track prediction;FROST reco x - external expected x [mm];Number of tracks",
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_reco_external_y = new TH1D(
    "hist_reco_external_residual_y",
    "FROST reco minus external-track prediction;FROST reco y - external expected y [mm];Number of tracks",
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_reco_truth_x = new TH1D(
    "hist_reco_truth_residual_x",
    "FROST reco minus FROST truth;FROST reco x - true FROST x [mm];Number of tracks",
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_reco_truth_y = new TH1D(
    "hist_reco_truth_residual_y",
    "FROST reco minus FROST truth;FROST reco y - true FROST y [mm];Number of tracks",
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_correlation_x = new TH2D(
    "hist_correlation_x",
    "Correlation of residuals, x;external expected x - true FROST x [mm];FROST reco x - true FROST x [mm]",
    nbins,
    xmin_mm,
    xmax_mm,
    nbins,
    xmin_mm,
    xmax_mm);

  auto *hist_correlation_y = new TH2D(
    "hist_correlation_y",
    "Correlation of residuals, y;external expected y - true FROST y [mm];FROST reco y - true FROST y [mm]",
    nbins,
    xmin_mm,
    xmax_mm,
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_correlation_x_by_angle = MakeAngleCorrelationHistograms(
    "hist_correlation_x",
    "Correlation of residuals, x",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_correlation_y_by_angle = MakeAngleCorrelationHistograms(
    "hist_correlation_y",
    "Correlation of residuals, y",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_external_truth_x_by_angle = MakeAngleHistograms(
    "hist_external_truth_residual_x",
    "External-track prediction minus FROST truth, x",
    "external expected x - true FROST x [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_external_truth_y_by_angle = MakeAngleHistograms(
    "hist_external_truth_residual_y",
    "External-track prediction minus FROST truth, y",
    "external expected y - true FROST y [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_reco_external_x_by_angle = MakeAngleHistograms(
    "hist_reco_external_residual_x",
    "FROST reco minus external-track prediction, x",
    "FROST reco x - external expected x [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_reco_external_y_by_angle = MakeAngleHistograms(
    "hist_reco_external_residual_y",
    "FROST reco minus external-track prediction, y",
    "FROST reco y - external expected y [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_reco_truth_x_by_angle = MakeAngleHistograms(
    "hist_reco_truth_residual_x",
    "FROST reco minus FROST truth, x",
    "FROST reco x - true FROST x [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  auto hist_reco_truth_y_by_angle = MakeAngleHistograms(
    "hist_reco_truth_residual_y",
    "FROST reco minus FROST truth, y",
    "FROST reco y - true FROST y [mm]",
    nbins,
    xmin_mm,
    xmax_mm);

  CorrelationStats correlation_stats_x;
  CorrelationStats correlation_stats_y;
  std::vector<CorrelationStats> correlation_stats_x_by_angle(kNAngleBins);
  std::vector<CorrelationStats> correlation_stats_y_by_angle(kNAngleBins);

  Long64_t n_tracks_x = 0;
  Long64_t n_tracks_y = 0;

  for (const auto &file_path : root_files) {
    std::cout << "Processing: " << file_path << std::endl;
    ProcessOneFile(file_path,
                   hist_x,
                   hist_y,
                   hist_reco_external_x,
                   hist_reco_external_y,
                   hist_reco_truth_x,
                   hist_reco_truth_y,
                   hist_correlation_x,
                   hist_correlation_y,
                   hist_external_truth_x_by_angle,
                   hist_external_truth_y_by_angle,
                   hist_reco_external_x_by_angle,
                   hist_reco_external_y_by_angle,
                   hist_reco_truth_x_by_angle,
                   hist_reco_truth_y_by_angle,
                   hist_correlation_x_by_angle,
                   hist_correlation_y_by_angle,
                   correlation_stats_x,
                   correlation_stats_y,
                   correlation_stats_x_by_angle,
                   correlation_stats_y_by_angle,
                   n_tracks_x,
                   n_tracks_y);
  }

  std::cout << "Selected tracks for x: " << n_tracks_x << std::endl;
  std::cout << "Selected tracks for y: " << n_tracks_y << std::endl;
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
      std::cout << Form("Angle %.0f-%.0f deg x correlation: ",
                        kAngleBins[i], kAngleBins[i + 1])
                << "N=" << correlation_stats_x_by_angle.at(i).n
                << ", Cov(E,R)=" << covariance << " mm^2"
                << ", rho=" << correlation << std::endl;
    }
    if (correlation_stats_y_by_angle.at(i).Calculate(
          mean_e, mean_r, rms_e, rms_r, covariance, correlation)) {
      std::cout << Form("Angle %.0f-%.0f deg y correlation: ",
                        kAngleBins[i], kAngleBins[i + 1])
                << "N=" << correlation_stats_y_by_angle.at(i).n
                << ", Cov(E,R)=" << covariance << " mm^2"
                << ", rho=" << correlation << std::endl;
    }
  }

  TCanvas canvas("canvas", "External track residuals", 1200, 600);
  TString pdf_open = output_pdf;
  pdf_open += "(";
  TString pdf_close = output_pdf;
  pdf_close += ")";

  DrawTwoPanelPage(canvas, hist_x, hist_y);
  canvas.Print(pdf_open);

  DrawTwoPanelPage(canvas, hist_reco_external_x, hist_reco_external_y);
  canvas.Print(output_pdf);

  DrawTwoPanelPage(canvas, hist_reco_truth_x, hist_reco_truth_y);
  canvas.Print(output_pdf);

  PrintAnglePages(canvas,
                  output_pdf,
                  hist_external_truth_x_by_angle,
                  hist_external_truth_y_by_angle);

  PrintAnglePages(canvas,
                  output_pdf,
                  hist_reco_external_x_by_angle,
                  hist_reco_external_y_by_angle);

  PrintAnglePages(canvas,
                  output_pdf,
                  hist_reco_truth_x_by_angle,
                  hist_reco_truth_y_by_angle);

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
}
