// File: draw_bm_data_quality_plots.C
//
// Usage:
//
//   root -l -q 'draw_bm_data_quality_plots.C("/path/to/root_dir","dq.pdf","tree_name")'
//
//   root -l -q 'draw_bm_data_quality_plots.C("/path/to/root_dir","dq.pdf","tree_name","2025/11/28/12:00","2025/12/09/12:30")'
//
// This macro reads ROOT files containing the "tree_name" tree and makes
// Baby MIND data quality plots in a multi-page PDF.
//
// The optional time range is specified in local time with the format:
//   yyyy/mm/dd/hh:mm
//
// If start_time_str or end_time_str is empty, no corresponding time cut is applied.

#include <TROOT.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <TChain.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TPaveStats.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <ctime>
#include <cmath>

static std::vector<std::string> GetRootFiles(const char* input_dir)
{
  std::vector<std::string> files;

  TSystemDirectory dir("input_dir", input_dir);
  TList* list = dir.GetListOfFiles();

  if (!list) {
    std::cerr << "Error: cannot open input directory: "
              << input_dir << std::endl;
    return files;
  }

  TIter next(list);
  TSystemFile* file = nullptr;

  while ((file = static_cast<TSystemFile*>(next()))) {
    TString name = file->GetName();

    if (file->IsDirectory()) continue;
    if (!name.EndsWith(".root")) continue;

    TString full_path = TString(input_dir) + "/" + name;
    files.push_back(full_path.Data());
  }

  std::sort(files.begin(), files.end());
  return files;
}

static bool ParseDateTime(const char* str, std::time_t& out_time)
{
  TString s = str;

  if (s.Length() == 0) {
    return false;
  }

  int year = 0;
  int month = 0;
  int day = 0;
  int hour = 0;
  int minute = 0;

  const int nread = std::sscanf(
    s.Data(),
    "%d/%d/%d/%d:%d",
    &year,
    &month,
    &day,
    &hour,
    &minute
  );

  if (nread != 5) {
    std::cerr << "Error: invalid datetime format: " << str << std::endl;
    std::cerr << "Expected format: yyyy/mm/dd/hh:mm" << std::endl;
    return false;
  }

  std::tm tm{};
  tm.tm_year = year - 1900;
  tm.tm_mon  = month - 1;
  tm.tm_mday = day;
  tm.tm_hour = hour;
  tm.tm_min  = minute;
  tm.tm_sec  = 0;
  tm.tm_isdst = -1;

  out_time = std::mktime(&tm);

  if (out_time < 0) {
    std::cerr << "Error: failed to convert datetime: " << str << std::endl;
    return false;
  }

  return true;
}

static std::string UnixtimeToString(std::time_t tt)
{
  std::tm lt{};
#if defined(_WIN32)
  localtime_s(&lt, &tt);
#else
  lt = *std::localtime(&tt);
#endif

  std::ostringstream oss;
  oss << std::setw(4) << std::setfill('0') << lt.tm_year + 1900
      << "/"
      << std::setw(2) << std::setfill('0') << lt.tm_mon + 1
      << "/"
      << std::setw(2) << std::setfill('0') << lt.tm_mday
      << "/"
      << std::setw(2) << std::setfill('0') << lt.tm_hour
      << ":"
      << std::setw(2) << std::setfill('0') << lt.tm_min;

  return oss.str();
}

static void DrawOnePage(TCanvas* c,
                        TH1D* h,
                        const char* output_pdf,
                        bool use_logy = false,
                        bool move_stats = false,
                        double titleoffsety = 1.25)
{
  c->Clear();
  c->SetLogy(use_logy);
  c->SetGrid();
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.12);

  h->SetStats(1);
  h->SetLineColor(kBlue + 1);
  h->SetLineWidth(2);
  h->SetFillStyle(0);

  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(titleoffsety);

  if (use_logy) {
    h->SetMinimum(0.5);
  }

  h->Draw("HIST");

  c->Update();

  if (move_stats) {
    TPaveStats* st = static_cast<TPaveStats*>(
      h->FindObject("stats")
    );

    if (st) {
      // Move the statistics box to the upper-left area.
      st->SetX1NDC(0.63);
      st->SetX2NDC(0.88);
      st->SetY1NDC(0.80);
      st->SetY2NDC(0.92);
    }

    c->Modified();
    c->Update();
  }

  c->SaveAs(output_pdf);
}

// void draw_bm_data_quality_plots(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/latest20260504/rootfile_after_TrackMatch/FHC",
//                                 const char* output_pdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bm_dataquality_plots_FHC.pdf",
//                                 const char* tree_name = "ntbm",
//                                 const char* start_time_str = "",
//                                 const char* end_time_str = "")
// void draw_bm_data_quality_plots(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/latest20260504/rootfile_after_TrackMatch/RHC",
//                                 const char* output_pdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bm_dataquality_plots_RHC.pdf",
//                                 const char* tree_name = "ntbm",
//                                 const char* start_time_str = "",
//                                 const char* end_time_str = "")
// void draw_bm_data_quality_plots(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/latest20260504/rootfile_after_TrackMatch/FHC",
//                                 const char* output_pdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bm_dataquality_plots_FHC_beforeJan24.pdf",
//                                 const char* tree_name = "ntbm",
//                                 const char* start_time_str = "",
//                                 const char* end_time_str = "2026/01/24/15:29")
// void draw_bm_data_quality_plots(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/latest20260504/rootfile_after_TrackMatch/FHC",
//                                 const char* output_pdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bm_dataquality_plots_FHC_afterJan24.pdf",
//                                 const char* tree_name = "ntbm",
//                                 const char* start_time_str = "2026/01/24/15:29",
//                                 const char* end_time_str = "")
// void draw_bm_data_quality_plots(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/E71a_data",
//                                 const char* output_pdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bm_dataquality_plots_E71a.pdf",
//                                 const char* tree_name = "tree",
//                                 const char* start_time_str = "",
//                                 const char* end_time_str = "")
void draw_bm_data_quality_plots(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/rootfile_after_TrackMatch_new/zshift/zshift_0/RHC",
                                const char* output_pdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bm_dataquality_plots_afterJan24_beforedeadplnfix_RHC.pdf",
                                const char* tree_name = "ntbm",
                                const char* start_time_str = "2026/01/24/15:29",
                                const char* end_time_str = "")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);

  std::vector<std::string> root_files = GetRootFiles(input_dir);

  if (root_files.empty()) {
    std::cerr << "Error: no ROOT files found in directory: "
              << input_dir << std::endl;
    return;
  }

  bool use_start_time = false;
  bool use_end_time = false;

  std::time_t start_time = 0;
  std::time_t end_time = 0;

  if (TString(start_time_str).Length() > 0) {
    use_start_time = ParseDateTime(start_time_str, start_time);
    if (!use_start_time) return;
  }

  if (TString(end_time_str).Length() > 0) {
    use_end_time = ParseDateTime(end_time_str, end_time);
    if (!use_end_time) return;
  }

  if (use_start_time && use_end_time && start_time > end_time) {
    std::cerr << "Error: start time is later than end time." << std::endl;
    return;
  }

  TChain chain(tree_name);

  for (const auto& path : root_files) {
    chain.Add(path.c_str());
  }

  std::cout << "Input directory : " << input_dir << std::endl;
  std::cout << "Output PDF      : " << output_pdf << std::endl;
  std::cout << "Number of files : " << root_files.size() << std::endl;
  std::cout << "Number of entries in chain : "
            << chain.GetEntries() << std::endl;

  if (use_start_time) {
    std::cout << "Start time      : "
              << UnixtimeToString(start_time)
              << "  unix = " << start_time << std::endl;
  } else {
    std::cout << "Start time      : none" << std::endl;
  }

  if (use_end_time) {
    std::cout << "End time        : "
              << UnixtimeToString(end_time)
              << "  unix = " << end_time << std::endl;
  } else {
    std::cout << "End time        : none" << std::endl;
  }

  TTreeReader reader(&chain);

  chain.LoadTree(0);
  TTree* first_tree = chain.GetTree();

  if (!first_tree) {
    std::cerr << "Error: failed to load the first tree from TChain." << std::endl;
    return;
  }

  const bool has_nll_minus =
    (first_tree->GetBranch("negative_log_likelihood_minus_") != nullptr);

  const bool has_nll_plus =
    (first_tree->GetBranch("negative_log_likelihood_plus_") != nullptr);

  const bool has_llr_branches = has_nll_minus && has_nll_plus;

  if (has_llr_branches) {
    std::cout << "LLR branches     : found" << std::endl;
  } else {
    std::cout << "LLR branches     : not found. Skip LLR plots." << std::endl;
  }


  // Event-level branch used for the time selection.
  TTreeReaderValue<Double_t> timestamp(reader, "timestamp_");

  // Track-level branches.
  TTreeReaderArray<Double_t> momentum(reader, "momentum_");
  TTreeReaderArray<Int_t> maximum_plane(reader, "baby_mind_maximum_plane_");
  TTreeReaderArray<Int_t> charge(reader, "charge_");

  TTreeReaderArray<Double_t>* nll_minus = nullptr;
  TTreeReaderArray<Double_t>* nll_plus = nullptr;

  if (has_llr_branches) {
    nll_minus = new TTreeReaderArray<Double_t>(
      reader,
      "negative_log_likelihood_minus_"
    );

    nll_plus = new TTreeReaderArray<Double_t>(
      reader,
      "negative_log_likelihood_plus_"
    );
  }

  TH1D* h_momentum = new TH1D(
    "h_momentum",
    ";Momentum by range [MeV/c];Number of tracks",
    1000,
    0.0,
    2000.0
  );

  TH1D* h_maximum_plane = new TH1D(
    "h_maximum_plane",
    ";Maximum plane [plane];Number of tracks",
    18,
    -0.5,
    17.5
  );

  TH1D* h_charge = new TH1D(
    "h_charge",
    ";Charge;Number of tracks",
    3,
    -1.5,
    1.5
  );

  TH1D* h_llr_wide = new TH1D(
    "h_llr_wide",
    ";Log likelihood ratio;Number of tracks",
    100,
    -4000.0,
    4000.0
  );

  TH1D* h_llr_narrow = new TH1D(
    "h_llr_narrow",
    ";Log likelihood ratio;Number of tracks",
    100,
    -10.0,
    10.0
  );

  Long64_t n_entries_total = 0;
  Long64_t n_entries_selected = 0;
  Long64_t n_tracks_filled = 0;
  Long64_t n_tracks_llr_filled = 0;

  while (reader.Next()) {
    ++n_entries_total;

    if (!std::isfinite(*timestamp) || *timestamp <= 0.0) continue;

    const std::time_t event_time =
      static_cast<std::time_t>(std::llround(*timestamp));

    if (use_start_time && event_time < start_time) continue;
    if (use_end_time && event_time > end_time) continue;

    ++n_entries_selected;

    const int n_momentum = momentum.GetSize();
    const int n_maximum_plane = maximum_plane.GetSize();
    const int n_charge = charge.GetSize();

    // The track vectors should have the same size, but each histogram is filled
    // with the available values independently to avoid losing useful entries.
    for (int itrack = 0; itrack < n_momentum; ++itrack) {
      if (std::isfinite(momentum[itrack])) {
        h_momentum->Fill(momentum[itrack]);
      }
    }

    for (int itrack = 0; itrack < n_maximum_plane; ++itrack) {
      h_maximum_plane->Fill(maximum_plane[itrack]);
    }

    for (int itrack = 0; itrack < n_charge; ++itrack) {
      h_charge->Fill(charge[itrack]);
    }

    if (has_llr_branches) {
      const int n_nll_minus = nll_minus->GetSize();
      const int n_nll_plus = nll_plus->GetSize();
      const int n_llr = std::min(n_nll_minus, n_nll_plus);

      for (int itrack = 0; itrack < n_llr; ++itrack) {
        if (!std::isfinite((*nll_minus)[itrack])) continue;
        if (!std::isfinite((*nll_plus)[itrack])) continue;

        const double llr = (*nll_minus)[itrack] - (*nll_plus)[itrack];

        h_llr_wide->Fill(llr);
        h_llr_narrow->Fill(llr);

        ++n_tracks_llr_filled;
      }
    }

    n_tracks_filled += std::max(
      std::max(n_momentum, n_maximum_plane),
      n_charge
    );
  }

  std::cout << "Total entries              : " << n_entries_total << std::endl;
  std::cout << "Selected entries           : " << n_entries_selected << std::endl;
  std::cout << "Track-like values filled   : " << n_tracks_filled << std::endl;
  std::cout << "LLR values filled          : " << n_tracks_llr_filled << std::endl;

  TCanvas* c = new TCanvas(
    "c_bm_dq",
    "Baby MIND data quality plots",
    1000,
    800
  );

  // Open the multi-page PDF.
  c->Print(Form("%s[", output_pdf));

  DrawOnePage(c, h_momentum, output_pdf, false, false, 1.45);
  DrawOnePage(c, h_maximum_plane, output_pdf, false, true);
  DrawOnePage(c, h_charge, output_pdf, false);

  if (has_llr_branches) {
    DrawOnePage(c, h_llr_wide, output_pdf, false);
    DrawOnePage(c, h_llr_wide, output_pdf, true);
    DrawOnePage(c, h_llr_narrow, output_pdf, false,false, 1.45);
  }

  // Close the multi-page PDF.
  c->Print(Form("%s]", output_pdf));

  std::cout << "Saved multi-page PDF to " << output_pdf << std::endl;

  delete c;

  delete h_momentum;
  delete h_maximum_plane;
  delete h_charge;
  delete h_llr_wide;
  delete h_llr_narrow;

  delete nll_minus;
  delete nll_plus;
}
