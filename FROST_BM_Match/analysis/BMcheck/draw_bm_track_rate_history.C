// File: draw_bm_track_rate_history.C
//
// Usage examples:
//   root -l -q 'draw_bm_track_rate_history.C("/path/to/root_dir","FHC","bm_track_rate_FHC.pdf")'
//   root -l -q 'draw_bm_track_rate_history.C("/path/to/root_dir","RHC","bm_track_rate_RHC.pdf")'
//
// This macro reads ROOT files containing the "ntbm" tree and makes a daily
// history plot of the number of Baby MIND tracks normalized by 10^15 POT.
//
// For FHC mode, only entries with bsd_good_spill_flag_ ==  3 are used.
// For RHC mode, only entries with bsd_good_spill_flag_ == -3 are used.
//
// Daily rate:
//   rate = sum(number_of_tracks_) * 1.0e15 / sum(spill_pot_)

#include <TCanvas.h>
#include <TChain.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TF1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>

#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <ctime>

static long long MakeDayKey(std::time_t tt)
{
  std::tm lt{};
#if defined(_WIN32)
  localtime_s(&lt, &tt);
#else
  lt = *std::localtime(&tt);
#endif

  return
    (static_cast<long long>(lt.tm_year + 1900) * 10000LL) +
    (static_cast<long long>(lt.tm_mon + 1) * 100LL) +
    static_cast<long long>(lt.tm_mday);
}

static std::string MakeDayLabel(long long day_key)
{
  int m = static_cast<int>((day_key / 100LL) % 100LL);
  int d = static_cast<int>(day_key % 100LL);

  std::ostringstream oss;
  oss << std::setw(2) << std::setfill('0') << m
      << "/"
      << std::setw(2) << std::setfill('0') << d;
  return oss.str();
}

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

static bool FillDailyTrackRateMaps(const char* input_dir,
                                   const TString& mode_str,
                                   int target_flag,
                                   const TString& tree_name,
                                   std::unordered_map<long long, double>& pot_per_day,
                                   std::unordered_map<long long, long long>& tracks_per_day,
                                   std::unordered_map<long long, long long>& spills_per_day,
                                   const char* sample_label)
{
  std::vector<std::string> root_files = GetRootFiles(input_dir);

  if (root_files.empty()) {
    std::cerr << "Error: no ROOT files found in directory: "
              << input_dir << std::endl;
    return false;
  }

  TChain chain(tree_name);

  for (const auto& path : root_files) {
    chain.Add(path.c_str());
  }

  std::cout << "[" << sample_label << "] Input directory : "
            << input_dir << std::endl;
  std::cout << "[" << sample_label << "] Tree name       : "
            << tree_name << std::endl;
  std::cout << "[" << sample_label << "] Beam mode       : "
            << mode_str << std::endl;
  std::cout << "[" << sample_label << "] Good spill flag : "
            << target_flag << std::endl;
  std::cout << "[" << sample_label << "] Number of files : "
            << root_files.size() << std::endl;
  std::cout << "[" << sample_label << "] Number of entries in chain : "
            << chain.GetEntries() << std::endl;

  TTreeReader reader(&chain);

  TTreeReaderValue<Double_t> spill_pot(reader, "spill_pot_");
  TTreeReaderValue<Double_t> timestamp(reader, "timestamp_");
  TTreeReaderValue<Int_t> bsd_good_spill_flag(reader, "bsd_good_spill_flag_");
  TTreeReaderValue<Int_t> number_of_tracks(reader, "number_of_tracks_");
  TTreeReaderArray<Int_t> detector_flags(reader, "detector_flags_[8]");

  Long64_t n_accepted = 0;

  while (reader.Next()) {
    if (*bsd_good_spill_flag != target_flag) continue;
    if (detector_flags.GetSize() <= 5) continue;
    if (detector_flags[5] == 0) continue;
    if (!std::isfinite(*spill_pot) || *spill_pot <= 0.0) continue;
    if (!std::isfinite(*timestamp) || *timestamp <= 0.0) continue;
    if (*number_of_tracks < 0) continue;

    std::time_t tt = static_cast<std::time_t>(std::llround(*timestamp));
    long long day_key = MakeDayKey(tt);

    pot_per_day[day_key] += *spill_pot;
    tracks_per_day[day_key] += static_cast<long long>(*number_of_tracks);
    spills_per_day[day_key] += 1;

    ++n_accepted;
  }

  std::cout << "[" << sample_label << "] Accepted entries : "
            << n_accepted << std::endl;

  if (pot_per_day.empty()) {
    std::cerr << "Error: no valid entries for " << sample_label
              << " / " << mode_str << " mode." << std::endl;
    return false;
  }

  return true;
}

static void CollectValidDays(const std::unordered_map<long long, double>& pot_per_day,
                             std::vector<long long>& day_keys)
{
  for (const auto& kv : pot_per_day) {
    const long long day_key = kv.first;
    const double pot = kv.second;

    if (pot <= 0.0) continue;
    day_keys.push_back(day_key);
  }
}

void draw_bm_track_rate_history(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/latest20260504/rootfile_after_TrackMatch/FHC",
                                const char* beam_mode = "FHC",
                                const char* output_pdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bm_track_rate_history_FHC_withdeadpln.pdf",
                                const char* input_dir_with_dead_planes = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/BMdeadpln2025/rootfile_after_TrackMatch")
// void draw_bm_track_rate_history(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/rootfile_after_TrackMatch_afterBMdeadplnfix",
//                                 const char* beam_mode = "FHC",
//                                 const char* output_pdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bm_track_rate_history_FHC.pdf",
//                                 const char* input_dir_with_dead_planes = "")
// void draw_bm_track_rate_history(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/rootfile_after_TrackMatch_new/zshift/zshift_0",
                                // const char* beam_mode = "RHC",
                                // const char* output_pdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bm_track_rate_history_beforedeadplnfix_RHC.pdf")
// void draw_bm_track_rate_history(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/E71a_data",
//                                 const char* beam_mode = "E71A",
//                                 const char* output_pdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bm_track_rate_history_E71a.pdf")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TString mode_str = beam_mode;
  mode_str.ToUpper();

  int target_flag = 0;
  if (mode_str == "FHC") {
    target_flag = 3;
  } else if (mode_str == "RHC") {
    target_flag = -3;
  }else if (mode_str == "E71A") {
    target_flag = 1;
  } else {
    std::cerr << "Error: beam_mode must be FHC, RHC, or E71A. Given: "
              << beam_mode << std::endl;
    return;
  }

  std::vector<std::string> root_files = GetRootFiles(input_dir);

  TString tree_name = "ntbm";

  if (mode_str == "E71A") {
    tree_name = "tree";
  }

  std::unordered_map<long long, double> pot_per_day;
  std::unordered_map<long long, long long> tracks_per_day;
  std::unordered_map<long long, long long> spills_per_day;

  const bool has_dead_plane_sample =
    (TString(input_dir_with_dead_planes).Length() > 0);

  std::unordered_map<long long, double> pot_per_day_dead;
  std::unordered_map<long long, long long> tracks_per_day_dead;
  std::unordered_map<long long, long long> spills_per_day_dead;

  if (!FillDailyTrackRateMaps(input_dir,
                              mode_str,
                              target_flag,
                              tree_name,
                              pot_per_day,
                              tracks_per_day,
                              spills_per_day,
                              "nominal")) {
    return;
  }

  if (has_dead_plane_sample) {
    if (!FillDailyTrackRateMaps(input_dir_with_dead_planes,
                                mode_str,
                                target_flag,
                                tree_name,
                                pot_per_day_dead,
                                tracks_per_day_dead,
                                spills_per_day_dead,
                                "with dead planes")) {
      return;
    }
  }


  std::vector<long long> day_keys;
  day_keys.reserve(pot_per_day.size());

  CollectValidDays(pot_per_day, day_keys);

  if (has_dead_plane_sample) {
    CollectValidDays(pot_per_day_dead, day_keys);
  }

  std::sort(day_keys.begin(), day_keys.end());
  day_keys.erase(std::unique(day_keys.begin(), day_keys.end()),
                 day_keys.end());

  if (day_keys.empty()) {
    std::cerr << "Error: no valid days with positive POT." << std::endl;
    return;
  }

  const int nbins = static_cast<int>(day_keys.size());

  TH1D* h = new TH1D("h_bm_track_rate",
                     Form("%s;Date;Number of Baby MIND tracks / 10^{15} POT",
                          mode_str.Data()),
                     nbins, 0.5, nbins + 0.5);

  TH1D* h_dead = nullptr;
  if (has_dead_plane_sample) {
    h_dead = new TH1D("h_bm_track_rate_with_dead_planes",
                      Form("%s;Date;Number of Baby MIND tracks / 10^{15} POT",
                           mode_str.Data()),
                      nbins, 0.5, nbins + 0.5);
  }

  h->SetStats(0);
  h->SetFillStyle(0);
  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.8);
  h->SetMarkerColor(kBlack);

  if (h_dead) {
    h_dead->SetStats(0);
    h_dead->SetFillStyle(0);
    h_dead->SetLineColor(kBlue + 1);
    h_dead->SetLineWidth(2);
    h_dead->SetMarkerStyle(24);
    h_dead->SetMarkerSize(0.8);
    h_dead->SetMarkerColor(kBlue + 1);
  }

  for (int i = 0; i < nbins; ++i) {
    const long long day_key = day_keys[i];

    h->GetXaxis()->SetBinLabel(i + 1, MakeDayLabel(day_key).c_str());
    if (h_dead) {
      h_dead->GetXaxis()->SetBinLabel(i + 1, MakeDayLabel(day_key).c_str());
    }

    const auto it_pot = pot_per_day.find(day_key);
    const auto it_tracks = tracks_per_day.find(day_key);

    if (it_pot != pot_per_day.end() &&
        it_tracks != tracks_per_day.end() &&
        it_pot->second > 0.0) {
      const double pot = it_pot->second;
      const long long n_tracks = it_tracks->second;
      const double rate = static_cast<double>(n_tracks) * 1.0e15 / pot;

      // Statistical error assuming Poisson uncertainty on the number of tracks.
      const double err_rate =
        (n_tracks > 0)
          ? std::sqrt(static_cast<double>(n_tracks)) * 1.0e15 / pot
          : 0.0;

      h->SetBinContent(i + 1, rate);
      h->SetBinError(i + 1, err_rate);

      std::cout << "[nominal] "
                << MakeDayLabel(day_key)
                << "  POT = " << std::setprecision(8) << pot
                << "  tracks = " << n_tracks
                << "  spills = " << spills_per_day[day_key]
                << "  rate = " << rate
                << " +/- " << err_rate
                << std::endl;
    }

    if (h_dead) {
      const auto it_pot_dead = pot_per_day_dead.find(day_key);
      const auto it_tracks_dead = tracks_per_day_dead.find(day_key);

      if (it_pot_dead != pot_per_day_dead.end() &&
          it_tracks_dead != tracks_per_day_dead.end() &&
          it_pot_dead->second > 0.0) {
        const double pot = it_pot_dead->second;
        const long long n_tracks = it_tracks_dead->second;
        const double rate = static_cast<double>(n_tracks) * 1.0e15 / pot;

        const double err_rate =
          (n_tracks > 0)
            ? std::sqrt(static_cast<double>(n_tracks)) * 1.0e15 / pot
            : 0.0;

        h_dead->SetBinContent(i + 1, rate);
        h_dead->SetBinError(i + 1, err_rate);

        std::cout << "[with dead planes] "
                  << MakeDayLabel(day_key)
                  << "  POT = " << std::setprecision(8) << pot
                  << "  tracks = " << n_tracks
                  << "  spills = " << spills_per_day_dead[day_key]
                  << "  rate = " << rate
                  << " +/- " << err_rate
                  << std::endl;
      }
    }
  }

  TCanvas* c = new TCanvas("c_bm_track_rate",
                           "Baby MIND track rate history",
                           1200, 800);

  c->SetGrid();
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.15);

  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->LabelsOption("v");

  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.1);

  double ymax = h->GetMaximum();
  double ymin = h->GetMinimum();

  if (ymax > ymin) {
    h->SetMinimum(std::max(0.0, ymin * 0.9));
    h->SetMaximum(ymax * 1.1);
  } else {
    h->SetMinimum(0.0);
    h->SetMaximum(ymax > 0.0 ? ymax * 1.5 : 1.0);
  }

  h->Draw("E1");
  if (h_dead) {
    h_dead->Draw("E1 SAME");
  }

  int n_fit_bins = 0;
  for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
    const double y = h->GetBinContent(ibin);
    const double ey = h->GetBinError(ibin);

    if (y > 0.0 && ey > 0.0 &&
        std::isfinite(y) && std::isfinite(ey)) {
      ++n_fit_bins;
    }
  }

  TF1* f_const = nullptr;
  TPaveText* pt = nullptr;

  if (n_fit_bins > 0) {
    f_const = new TF1("f_const", "[0]",
                      h->GetXaxis()->GetXmin(),
                      h->GetXaxis()->GetXmax());

    f_const->SetLineWidth(2);
    f_const->SetLineColor(kRed + 2);

    h->Fit(f_const, "Q0");
    f_const->Draw("SAME");

    const double p0 = f_const->GetParameter(0);
    const double p0_err = f_const->GetParError(0);
    const double chi2 = f_const->GetChisquare();
    const int ndf = f_const->GetNDF();

    pt = new TPaveText(0.70, 0.82, 0.95, 0.95, "NDC");
    pt->SetBorderSize(1);
    pt->SetFillColor(kWhite);
    pt->SetFillStyle(1001);
    pt->SetLineColor(kBlack);
    pt->SetLineWidth(1);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.03);
    pt->AddText(Form("%s mode", mode_str.Data()));
    pt->AddText(Form("#chi^{2} / ndf = %.1f / %d", chi2, ndf));
    pt->AddText(Form("p_{0} = %.3f #pm %.3f", p0, p0_err));
    pt->Draw();
  }

  TLegend* legend = nullptr;
  if (h_dead) {
    legend = new TLegend(0.15, 0.78, 0.45, 0.92);
    legend->SetBorderSize(1);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(1001);
    legend->SetTextSize(0.035);
    legend->AddEntry(h, "nominal", "lep");
    legend->AddEntry(h_dead, "with dead planes", "lep");
    legend->Draw();
  }

  c->SaveAs(output_pdf);

  std::cout << "Saved plot to " << output_pdf << std::endl;

  delete legend;
  delete pt;
  delete f_const;
  delete c;
  delete h_dead;
  delete h;
}
