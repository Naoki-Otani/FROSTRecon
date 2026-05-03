// plot_pot_spill_10min_merged.C
// Usage (ROOT):
//   root -l -q plot_pot_spill_10min_merged.C
// or inside ROOT:
//   .L plot_pot_spill_10min_merged.C
//   plot_pot_spill_10min_merged();

#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGaxis.h>

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <memory>

static inline std::string Trim(const std::string& s)
{
  size_t b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  size_t e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

static bool SplitCSV(const std::string& line, std::vector<std::string>& out)
{
  out.clear();
  std::string token;
  std::stringstream ss(line);
  while (std::getline(ss, token, ',')) out.push_back(Trim(token));
  return !out.empty();
}

// Parse "YYYY/MM/DD/HH:MM" (hour may be 1 or 2 digits)
static bool ParseYMDHM(const std::string& s, std::time_t& out_t)
{
  int Y=0, M=0, D=0, h=0, m=0;
  // Accept e.g. 2026/01/28/17:30
  if (std::sscanf(s.c_str(), "%d/%d/%d/%d:%d", &Y, &M, &D, &h, &m) != 5) return false;

  std::tm tm{};
  tm.tm_year = Y - 1900;
  tm.tm_mon  = M - 1;
  tm.tm_mday = D;
  tm.tm_hour = h;
  tm.tm_min  = m;
  tm.tm_sec  = 0;
  tm.tm_isdst = -1; // let system decide

  std::time_t t = std::mktime(&tm); // local time
  if (t == (std::time_t)(-1)) return false;
  out_t = t;
  return true;
}

// Recorded を「水平: 点線」「縦: 実線」のステップで描画
static void DrawStepRecorded(const std::vector<double>& t,
                             const std::vector<double>& y,
                             Color_t col = kRed,
                             int lw = 2,
                             int ls_horiz = 2, // 点線
                             int ls_vert  = 1) // 実線
{
  const int n = (int)t.size();
  if (n < 2) return;

  for (int i = 0; i < n - 1; ++i) {
    const double t0 = t[i];
    const double t1 = t[i+1];
    const double y0 = y[i];
    const double y1 = y[i+1];

    // 水平区間 (t0,y0)->(t1,y0) を点線
    {
      double xx[2] = {t0, t1};
      double yy[2] = {y0, y0};
      auto g = std::make_unique<TGraph>(2, xx, yy);
      g->SetLineColor(col);
      g->SetLineWidth(lw);
      g->SetLineStyle(ls_horiz);
      g->SetMarkerStyle(0);
      g->Draw("L SAME");
      g.release(); // ROOTに所有させる（マクロ実行ならこれでOK）
    }

    // 縦区間 (t1,y0)->(t1,y1) を実線（変化があるときだけ描く）
    if (y1 != y0) {
      double xx[2] = {t1, t1};
      double yy[2] = {y0, y1};
      auto g = std::make_unique<TGraph>(2, xx, yy);
      g->SetLineColor(col);
      g->SetLineWidth(lw);
      g->SetLineStyle(ls_vert);
      g->SetMarkerStyle(0);
      g->Draw("L SAME");
      g.release();
    }
  }
}

static bool IsFinite(double x) { return std::isfinite(x); }

void accumulatedpot_merge(const char *outputpdf = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/pot/accumulated_pot.pdf", const char* log_file_path = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/pot/accumulated_pot.log")
{
  std::ofstream logofs(log_file_path);
  if (!logofs) {
    std::cerr << "[accumulatedpot_merge] ERROR: cannot open log file: "
              << log_file_path << "\n";
    return;
  }

  std::streambuf* cout_backup = std::cout.rdbuf();
  std::cout.rdbuf(logofs.rdbuf());

  // ---- input CSVs (time order: e71c_data -> e71c_data2 -> e71c_data3) ----
  std::vector<std::string> csvs = {
    "/home/nu/notani/e71c_data/dataquality_withBSD/pot_info/pot_spill_10min.csv",
    "/home/nu/notani/e71c_data2/dataquality_withBSD/pot_info/pot_spill_10min.csv",
    "/home/nu/notani/e71c_data3/dataquality_withBSD/pot_info/pot_spill_10min.csv",
    "/home/nu/notani/e71c_data4/dataquality_withBSD/pot_info/pot_spill_10min.csv"
  };

  std::vector<std::string> csvs_frost_bm = {
    "/home/nu/notani/FROSTReconData/BM_FROST/analysis_plot/pot/RecordedPOT_byFROSTandBM_1.csv",
    "/home/nu/notani/FROSTReconData/BM_FROST/analysis_plot/pot/RecordedPOT_byFROSTandBM_2.csv",
    "/home/nu/notani/FROSTReconData/BM_FROST/analysis_plot/pot/RecordedPOT_byFROSTandBM_3.csv",
    "/home/nu/notani/FROSTReconData/BM_FROST/analysis_plot/pot/RecordedPOT_byFROSTandBM_4.csv"
  };

  std::vector<double> vx;      // unix time [s]
  std::vector<double> vy_del;  // accumulated delivered POT (stitched)
  std::vector<double> vy_rec;  // accumulated recorded  POT (stitched)
  std::vector<double> vx_bm;      // unix time [s]
  std::vector<double> vy_rec_bm;  // accumulated recorded POT by FROST and BM (stitched)

  double offset_del = 0.0;
  double offset_rec = 0.0;
  double offset_rec_bm = 0.0;

  std::time_t last_time_global = 0;
  bool has_last_time = false;

  for (const auto& path : csvs) {
    std::ifstream fin(path);
    if (!fin) {
      std::cerr << "[plot_pot_spill_10min_merged] ERROR: cannot open " << path << "\n";
      continue;
    }

    std::string line;
    double last_global_del_in_this_file = offset_del;
    double last_global_rec_in_this_file = offset_rec;
    bool   saw_valid_row = false;

    while (std::getline(fin, line)) {
      line = Trim(line);
      if (line.empty()) continue;
      if (!line.empty() && line[0] == '#') continue;

      std::vector<std::string> cols;
      if (!SplitCSV(line, cols)) continue;
      // expected:
      // 0: time
      // 1: delivered spill
      // 2: recorded spill
      // 3: delivered POT
      // 4: recorded POT
      // 5: efficiency
      if (cols.size() < 5) continue;

      std::time_t t = 0;
      if (!ParseYMDHM(cols[0], t)) continue;

      // keep monotonic time across files (drop overlaps)
      if (has_last_time && t <= last_time_global) continue;

      double del_pot = 0.0;
      double rec_pot = 0.0;
      try {
        del_pot = std::stod(cols[3]);
        rec_pot = std::stod(cols[4]);
      } catch (...) {
        continue;
      }
      if (!IsFinite(del_pot) || !IsFinite(rec_pot)) continue;

      double gdel = del_pot + offset_del;
      double grec = rec_pot + offset_rec;

      vx.push_back((double)t);
      vy_del.push_back(gdel);
      vy_rec.push_back(grec);

      last_global_del_in_this_file = gdel;
      last_global_rec_in_this_file = grec;
      saw_valid_row = true;

      last_time_global = t;
      has_last_time = true;
    }

    // update offsets for next file
    if (saw_valid_row) {
      offset_del = last_global_del_in_this_file;
      offset_rec = last_global_rec_in_this_file;
    }

    fin.close();
  }

  // ---- merge Recorded POT by FROST and Baby MIND CSVs ----
  std::time_t last_time_global_bm = 0;
  bool has_last_time_bm = false;

  for (const auto& path : csvs_frost_bm) {
    std::ifstream fin(path);
    if (!fin) {
      std::cerr << "[plot_pot_spill_10min_merged] ERROR: cannot open " << path << "\n";
      continue;
    }

    std::string line;
    double last_global_rec_bm_in_this_file = offset_rec_bm;
    bool saw_valid_row = false;

    while (std::getline(fin, line)) {
      line = Trim(line);
      if (line.empty()) continue;
      if (!line.empty() && line[0] == '#') continue;

      std::vector<std::string> cols;
      if (!SplitCSV(line, cols)) continue;
      // expected:
      // 0: time
      // 1: recorded POT by FROST and BM
      if (cols.size() < 2) continue;

      std::time_t t = 0;
      if (!ParseYMDHM(cols[0], t)) continue;

      // keep monotonic time across files (drop overlaps)
      if (has_last_time_bm && t <= last_time_global_bm) continue;

      double rec_bm_pot = 0.0;
      try {
        rec_bm_pot = std::stod(cols[1]);
      } catch (...) {
        continue;
      }
      if (!IsFinite(rec_bm_pot)) continue;

      double grec_bm = rec_bm_pot + offset_rec_bm;

      vx_bm.push_back((double)t);
      vy_rec_bm.push_back(grec_bm);

      last_global_rec_bm_in_this_file = grec_bm;
      saw_valid_row = true;

      last_time_global_bm = t;
      has_last_time_bm = true;
    }

    if (saw_valid_row) {
      offset_rec_bm = last_global_rec_bm_in_this_file;
    }

    fin.close();
  }

  if (vx.empty()) {
    std::cerr << "[plot_pot_spill_10min_merged] ERROR: no valid points.\n";
    return;
  }

  // ---- plot ----
  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas("c_pot10min_merged", "Accumulated POT (10-min CSV merged)", 1600, 800);
  c->SetGrid();
  c->SetLeftMargin(0.08);
  c->SetRightMargin(0.04);
  c->SetBottomMargin(0.17);
  c->SetTopMargin(0.05);

  // draw values in units of 1e20 POT
  std::vector<double> vy_del_1e20(vy_del.size());
  std::vector<double> vy_rec_1e20(vy_rec.size());
  std::vector<double> vy_rec_bm_1e20(vy_rec_bm.size());
  for (size_t i = 0; i < vy_del.size(); ++i) {
    vy_del_1e20[i] = vy_del[i] / 1.0e20;
    vy_rec_1e20[i] = vy_rec[i] / 1.0e20;
  }
  for (size_t i = 0; i < vy_rec_bm.size(); ++i) {
    vy_rec_bm_1e20[i] = vy_rec_bm[i] / 1.0e20;
  }

  TGraph* gr_del = new TGraph((int)vx.size(), vx.data(), vy_del_1e20.data());
  gr_del->SetTitle(";time;Accumulated POT [#times10^{20}]");
  gr_del->SetLineColor(kBlack);
  gr_del->SetLineWidth(2);
  gr_del->SetLineStyle(1);
  gr_del->SetMarkerStyle(0);
  gr_del->SetMinimum(0);

  TAxis* ax = gr_del->GetXaxis();
  ax->SetTimeDisplay(1);
  ax->SetTimeOffset(0, "local");
  ax->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  ax->SetLabelSize(0.05);
  ax->SetTitleSize(0.05);
  ax->SetLabelOffset(0.04);
  ax->SetTitleOffset(1.75);
  // ax->SetNdivisions(520);

  TAxis* ay = gr_del->GetYaxis();
  ay->SetLabelSize(0.05);
  ay->SetTitleSize(0.05);
  ay->SetTitleOffset(0.7);
  ay->SetNoExponent(kTRUE); // 余計な 10^N 表示を消す

  gr_del->Draw("AL");

  // TGraph* gr_rec = new TGraph((int)vx.size(), vx.data(), vy_rec.data());
  // gr_rec->SetLineColor(kRed);
  // gr_rec->SetLineWidth(2);
  // gr_rec->SetLineStyle(2); // <-- Recorded POT を点線に
  // gr_rec->SetMarkerStyle(0);
  // gr_rec->Draw("L SAME");

  // Recorded POT by FROST
  DrawStepRecorded(vx, vy_rec_1e20, kRed, 2, 2, 1);

  // Recorded POT by FROST and Baby MIND
  DrawStepRecorded(vx_bm, vy_rec_bm_1e20, kBlue, 2, 2, 1);

  TGraph* dummy_rec = new TGraph();
  dummy_rec->SetLineColor(kRed);
  dummy_rec->SetLineWidth(2);
  dummy_rec->SetLineStyle(2); // 凡例表示は点線でOK

  TGraph* dummy_rec_bm = new TGraph();
  dummy_rec_bm->SetLineColor(kBlue);
  dummy_rec_bm->SetLineWidth(2);
  dummy_rec_bm->SetLineStyle(2);

  TLegend* leg = new TLegend(0.12, 0.75, 0.45, 0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gr_del, "Delivered POT", "l");
  // leg->AddEntry(gr_rec, "Recorded POT",  "l");
  leg->AddEntry(dummy_rec, "Recorded POT by FROST", "l");
  leg->AddEntry(dummy_rec_bm, "Recorded POT by FROST and Baby MIND", "l");
  leg->SetTextSize(0.05);
  leg->Draw();

  c->SaveAs(outputpdf);

  // ---- last point summary ----
  const double last_del    = vy_del.back();     // Accumulated (Delivered) POT
  const double last_rec    = vy_rec.back();     // Recorded POT by FROST
  const double last_rec_bm = vy_rec_bm.back();  // Recorded POT by FROST and Baby MIND

  const std::time_t last_t    = (std::time_t)vx.back();
  const std::time_t last_t_bm = (std::time_t)vx_bm.back();

  std::tm lt{};
  lt = *std::localtime(&last_t);

  std::tm lt_bm{};
  lt_bm = *std::localtime(&last_t_bm);

  std::cout << std::setprecision(16);
  std::cout << "[LAST POINT]\n";
  std::cout << "  time (Delivered / FROST) = "
            << (lt.tm_year + 1900) << "/"
            << std::setw(2) << std::setfill('0') << (lt.tm_mon + 1) << "/"
            << std::setw(2) << std::setfill('0') << lt.tm_mday << " "
            << std::setw(2) << std::setfill('0') << lt.tm_hour << ":"
            << std::setw(2) << std::setfill('0') << lt.tm_min
            << "\n";
  std::cout << "  time (FROST and BM)      = "
            << (lt_bm.tm_year + 1900) << "/"
            << std::setw(2) << std::setfill('0') << (lt_bm.tm_mon + 1) << "/"
            << std::setw(2) << std::setfill('0') << lt_bm.tm_mday << " "
            << std::setw(2) << std::setfill('0') << lt_bm.tm_hour << ":"
            << std::setw(2) << std::setfill('0') << lt_bm.tm_min
            << "\n";
  std::cout << "  Delivered POT                  = " << last_del << "\n";
  std::cout << "  Recorded POT by FROST          = " << last_rec << "\n";
  std::cout << "  Recorded POT by FROST and BM   = " << last_rec_bm << "\n";
  if (last_del > 0.0) {
    std::cout << "  efficiency by FROST[%]         = "
              << (100.0 * last_rec / last_del) << "\n";
    std::cout << "  efficiency by FROST and BM[%]  = "
              << (100.0 * last_rec_bm / last_del) << "\n";
  }

  std::cout.flush();
  std::cout.rdbuf(cout_backup);
  logofs.close();
}
