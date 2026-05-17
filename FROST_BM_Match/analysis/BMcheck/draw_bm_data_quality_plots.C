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
#include <TLegend.h>

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

static bool ParseOptionalTimeRange(const char* start_time_str,
                                   const char* end_time_str,
                                   bool& use_start_time,
                                   bool& use_end_time,
                                   std::time_t& start_time,
                                   std::time_t& end_time)
{
  use_start_time = false;
  use_end_time = false;
  start_time = 0;
  end_time = 0;

  if (TString(start_time_str).Length() > 0) {
    use_start_time = ParseDateTime(start_time_str, start_time);
    if (!use_start_time) return false;
  }

  if (TString(end_time_str).Length() > 0) {
    use_end_time = ParseDateTime(end_time_str, end_time);
    if (!use_end_time) return false;
  }

  if (use_start_time && use_end_time && start_time > end_time) {
    std::cerr << "Error: start time is later than end time." << std::endl;
    return false;
  }

  return true;
}

static TH1D* MakeHistLike(const char* name, const TH1D* ref)
{
  TH1D* h = static_cast<TH1D*>(ref->Clone(name));
  h->Reset("ICES");
  return h;
}

static void SetHistogramStyle(TH1D* h,
                              Color_t color,
                              Style_t line_style = 1)
{
  if (!h) return;

  h->SetStats(1);
  h->SetLineColor(color);
  h->SetLineWidth(2);
  h->SetLineStyle(line_style);
  h->SetFillStyle(0);
}

static void ScaleHistogramToReferenceArea(TH1D* h,
                                          const TH1D* h_ref,
                                          const char* hist_label)
{
  if (!h || !h_ref) return;

  // Normalize the area of h to the area of h_ref.
  // Underflow and overflow bins are not included.
  const double area = h->Integral(1, h->GetNbinsX());
  const double area_ref = h_ref->Integral(1, h_ref->GetNbinsX());

  if (area <= 0.0) {
    std::cerr << "Warning: cannot scale " << hist_label
              << " because the input2 histogram area is zero." << std::endl;
    return;
  }

  if (area_ref <= 0.0) {
    std::cerr << "Warning: cannot scale " << hist_label
              << " because the input1 histogram area is zero." << std::endl;
    return;
  }

  const double scale = area_ref / area;
  h->Scale(scale);

  std::cout << "Scale " << hist_label
            << " input2 histogram by " << scale
            << " so that area(input2) = area(input1)." << std::endl;
}

static bool FillDataQualityHistograms(const char* input_dir,
                                      const char* tree_name,
                                      const char* label,
                                      const char* start_time_str,
                                      const char* end_time_str,
                                      TH1D* h_momentum,
                                      TH1D* h_maximum_plane,
                                      TH1D* h_charge,
                                      TH1D* h_llr_wide,
                                      TH1D* h_llr_narrow,
                                      bool& has_llr_branches)
{
  std::vector<std::string> root_files = GetRootFiles(input_dir);

  if (root_files.empty()) {
    std::cerr << "Error: no ROOT files found in directory: "
              << input_dir << std::endl;
    return false;
  }

  bool use_start_time = false;
  bool use_end_time = false;
  std::time_t start_time = 0;
  std::time_t end_time = 0;

  if (!ParseOptionalTimeRange(start_time_str,
                              end_time_str,
                              use_start_time,
                              use_end_time,
                              start_time,
                              end_time)) {
    return false;
  }

  TChain chain(tree_name);

  for (const auto& path : root_files) {
    chain.Add(path.c_str());
  }

  std::cout << "[" << label << "] Input directory : "
            << input_dir << std::endl;
  std::cout << "[" << label << "] Tree name       : "
            << tree_name << std::endl;
  std::cout << "[" << label << "] Number of files : "
            << root_files.size() << std::endl;
  std::cout << "[" << label << "] Number of entries in chain : "
            << chain.GetEntries() << std::endl;

  if (use_start_time) {
    std::cout << "[" << label << "] Start time      : "
              << UnixtimeToString(start_time)
              << "  unix = " << start_time << std::endl;
  } else {
    std::cout << "[" << label << "] Start time      : none" << std::endl;
  }

  if (use_end_time) {
    std::cout << "[" << label << "] End time        : "
              << UnixtimeToString(end_time)
              << "  unix = " << end_time << std::endl;
  } else {
    std::cout << "[" << label << "] End time        : none" << std::endl;
  }

  chain.LoadTree(0);
  TTree* first_tree = chain.GetTree();

  if (!first_tree) {
    std::cerr << "Error: failed to load the first tree from TChain for "
              << label << "." << std::endl;
    return false;
  }

  const bool has_nll_minus =
    (first_tree->GetBranch("negative_log_likelihood_minus_") != nullptr);

  const bool has_nll_plus =
    (first_tree->GetBranch("negative_log_likelihood_plus_") != nullptr);

  has_llr_branches = has_nll_minus && has_nll_plus;

  if (has_llr_branches) {
    std::cout << "[" << label << "] LLR branches     : found" << std::endl;
  } else {
    std::cout << "[" << label << "] LLR branches     : not found. Skip LLR filling." << std::endl;
  }

  TTreeReader reader(&chain);

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

  std::cout << "[" << label << "] Total entries              : "
            << n_entries_total << std::endl;
  std::cout << "[" << label << "] Selected entries           : "
            << n_entries_selected << std::endl;
  std::cout << "[" << label << "] Track-like values filled   : "
            << n_tracks_filled << std::endl;
  std::cout << "[" << label << "] LLR values filled          : "
            << n_tracks_llr_filled << std::endl;

  delete nll_minus;
  delete nll_plus;

  return true;
}

static void DrawOverlayPage(TCanvas* c,
                            TH1D* h1,
                            TH1D* h2,
                            const char* label1,
                            const char* label2,
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

  h1->GetXaxis()->SetTitleSize(0.045);
  h1->GetYaxis()->SetTitleSize(0.045);
  h1->GetXaxis()->SetLabelSize(0.04);
  h1->GetYaxis()->SetLabelSize(0.04);
  h1->GetXaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetTitleOffset(titleoffsety);

  if (use_logy) {
    h1->SetMinimum(0.5);
    if (h2) h2->SetMinimum(0.5);
  }

  double ymax = h1->GetMaximum();
  if (h2) ymax = std::max(ymax, h2->GetMaximum());

  if (ymax > 0.0) {
    h1->SetMaximum(use_logy ? ymax * 5.0 : ymax * 1.15);
  }

  h1->Draw("HIST");

  if (h2) {
    h2->Draw("HIST SAME");

    TLegend* legend = new TLegend(0.45, 0.85, 0.65, 0.92);
    legend->SetBorderSize(1);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(1001);
    legend->SetTextSize(0.035);
    legend->AddEntry(h1, label1, "l");
    legend->AddEntry(h2, label2, "l");
    legend->Draw();
  }

  c->Update();

  if (move_stats) {
    TPaveStats* st = static_cast<TPaveStats*>(
      h1->FindObject("stats")
    );

    if (st) {
      // Move the statistics box to the upper-right area.
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
// void draw_bm_data_quality_plots(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/rootfile_after_TrackMatch_new/zshift/zshift_0/RHC",
//                                 const char* output_pdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bm_dataquality_plots_afterJan24_beforedeadplnfix_RHC.pdf",
//                                 const char* tree_name = "ntbm",
//                                 const char* start_time_str = "2026/01/24/15:29",
//                                 const char* end_time_str = "",
//                                 const char* label = "sample 1",
//                                 const char* input_dir2 = "",
//                                 const char* label2 = "sample 2",
//                                 const char* tree_name2 = "ntbm",
//                                 const char* start_time_str2 = "",
//                                 const char* end_time_str2 = "")
// void draw_bm_data_quality_plots(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/latest20260504/rootfile_after_TrackMatch/FHC",
//                                 const char* output_pdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bm_dataquality_plots_FHC_compare_E71a_and_E71c.pdf",
//                                 const char* tree_name = "ntbm",
//                                 const char* start_time_str = "",
//                                 const char* end_time_str = "",
//                                 const char* label = "T2K Run15",
//                                 const char* input_dir2 = "/group/nu/ninja/work/otani/FROSTReconData/E71a_data",
//                                 const char* label2 = "T2K Run10",
//                                 const char* tree_name2 = "tree",
//                                 const char* start_time_str2 = "",
//                                 const char* end_time_str2 = "")
void draw_bm_data_quality_plots(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/latest20260504/rootfile_after_TrackMatch/FHC",
                                const char* output_pdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bm_dataquality_plots_FHC_compare_2025_and_2026.pdf",
                                const char* tree_name = "ntbm",
                                const char* start_time_str = "",
                                const char* end_time_str = "2025/12/31/23:59",
                                const char* label = "2025",
                                const char* input_dir2 = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/latest20260504/rootfile_after_TrackMatch/FHC",
                                const char* label2 = "2026",
                                const char* tree_name2 = "ntbm",
                                const char* start_time_str2 = "2026/01/01/00:00",
                                const char* end_time_str2 = "")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const bool has_second_sample = (TString(input_dir2).Length() > 0);

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

  TH1D* h_momentum2 = nullptr;
  TH1D* h_maximum_plane2 = nullptr;
  TH1D* h_charge2 = nullptr;
  TH1D* h_llr_wide2 = nullptr;
  TH1D* h_llr_narrow2 = nullptr;

  if (has_second_sample) {
    h_momentum2 = MakeHistLike("h_momentum2", h_momentum);
    h_maximum_plane2 = MakeHistLike("h_maximum_plane2", h_maximum_plane);
    h_charge2 = MakeHistLike("h_charge2", h_charge);
    h_llr_wide2 = MakeHistLike("h_llr_wide2", h_llr_wide);
    h_llr_narrow2 = MakeHistLike("h_llr_narrow2", h_llr_narrow);
  }

  SetHistogramStyle(h_momentum, kBlue + 1);
  SetHistogramStyle(h_maximum_plane, kBlue + 1);
  SetHistogramStyle(h_charge, kBlue + 1);
  SetHistogramStyle(h_llr_wide, kBlue + 1);
  SetHistogramStyle(h_llr_narrow, kBlue + 1);

  SetHistogramStyle(h_momentum2, kRed + 1, 2);
  SetHistogramStyle(h_maximum_plane2, kRed + 1, 2);
  SetHistogramStyle(h_charge2, kRed + 1, 2);
  SetHistogramStyle(h_llr_wide2, kRed + 1, 2);
  SetHistogramStyle(h_llr_narrow2, kRed + 1, 2);
  bool has_llr_branches = false;
  bool has_llr_branches2 = false;

  if (!FillDataQualityHistograms(input_dir,
                                 tree_name,
                                 label,
                                 start_time_str,
                                 end_time_str,
                                 h_momentum,
                                 h_maximum_plane,
                                 h_charge,
                                 h_llr_wide,
                                 h_llr_narrow,
                                 has_llr_branches)) {
    return;
  }

  if (has_second_sample) {
    if (!FillDataQualityHistograms(input_dir2,
                                   tree_name2,
                                   label2,
                                   start_time_str2,
                                   end_time_str2,
                                   h_momentum2,
                                   h_maximum_plane2,
                                   h_charge2,
                                   h_llr_wide2,
                                   h_llr_narrow2,
                                   has_llr_branches2)) {
      return;
    }
  }

  if (has_second_sample) {
    ScaleHistogramToReferenceArea(h_momentum2,
                                  h_momentum,
                                  "Momentum by range");

    ScaleHistogramToReferenceArea(h_maximum_plane2,
                                  h_maximum_plane,
                                  "Maximum plane");

    ScaleHistogramToReferenceArea(h_charge2,
                                  h_charge,
                                  "Charge");

    if (has_llr_branches && has_llr_branches2) {
      ScaleHistogramToReferenceArea(h_llr_wide2,
                                    h_llr_wide,
                                    "Log-likelihood ratio wide");

      ScaleHistogramToReferenceArea(h_llr_narrow2,
                                    h_llr_narrow,
                                    "Log-likelihood ratio narrow");
    }
  }

  const bool draw_llr_pages =
    has_llr_branches || (has_second_sample && has_llr_branches2);

  TCanvas* c = new TCanvas(
    "c_bm_dq",
    "Baby MIND data quality plots",
    1000,
    800
  );

  // Open the multi-page PDF.
  c->Print(Form("%s[", output_pdf));

  DrawOverlayPage(c, h_momentum, h_momentum2,
                  label, label2, output_pdf,
                  false, false, 1.45);

  DrawOverlayPage(c, h_maximum_plane, h_maximum_plane2,
                  label, label2, output_pdf,
                  false, true);

  DrawOverlayPage(c, h_charge, h_charge2,
                  label, label2, output_pdf,
                  false);

  if (draw_llr_pages) {
    DrawOverlayPage(c,
                    h_llr_wide,
                    has_second_sample ? h_llr_wide2 : nullptr,
                    label,
                    label2,
                    output_pdf,
                    false);

    DrawOverlayPage(c,
                    h_llr_wide,
                    has_second_sample ? h_llr_wide2 : nullptr,
                    label,
                    label2,
                    output_pdf,
                    true);

    DrawOverlayPage(c,
                    h_llr_narrow,
                    has_second_sample ? h_llr_narrow2 : nullptr,
                    label,
                    label2,
                    output_pdf,
                    false,
                    false,
                    1.45);
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

  delete h_momentum2;
  delete h_maximum_plane2;
  delete h_charge2;
  delete h_llr_wide2;
  delete h_llr_narrow2;
}
