#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TParameter.h"

#include "CLI.h"
#include "MultiHitAlgo.h"
#include "ReconstructionUtil.h"
#include "StringUtil.h"
#include "FROSTConstants.h"
#include "TreeReader.h"

namespace {

struct PointWithErr {
  double x = 0.0;
  double y = 0.0;
  double ex = 0.0;
  double ey = 0.0;
};

}  // namespace

int main(int argc, char** argv) {
  try {
    FROST::CLI cli(argc, argv);

    cli.Require("mode");
    cli.Require("alpha");
    cli.Require("out");

    const std::string mode = cli.Get("mode");

    // Keep both the string form (for filename templates) and the numeric form (for calculations).
    const std::string alpha_str = cli.Get("alpha");
    const double alpha = std::stod(alpha_str);
    const std::string out_path = cli.Get("out");
    const int max_events = cli.GetInt("max-events", 0);

    const double threshold1 = cli.GetDouble("threshold1", 10.0);
    const double threshold2 = cli.GetDouble("threshold2", 20.0);

    FROST::MultiHitAlgoConfig mh_cfg;
    mh_cfg.alpha = alpha;
    mh_cfg.threshold1 = threshold1;
    mh_cfg.threshold2 = threshold2;

    if (mode == "single") {
      cli.Require("in-x-pattern");
      cli.Require("in-y-pattern");

      const std::string pat_x = cli.Get("in-x-pattern");
      const std::string pat_y = cli.Get("in-y-pattern");

      // Accumulators at the label positions only.
      std::vector<double> xg_mean(FROST::kNx, 0.0);
      std::vector<double> xg_err(FROST::kNx, 0.0);
      std::vector<std::vector<double>> mcx_avg(FROST::kNfibX, std::vector<double>(FROST::kNx, 0.0));

      for (int q = 0; q < FROST::kNx; ++q) {
        const int label = FROST::kXLabel[q];
        std::unordered_map<std::string, std::string> vars{{"label", std::to_string(label)}};
        const std::string in_path = FROST::ApplyTemplate(pat_x, vars);

        FROST::TreeReader reader(in_path, "wls", /*read_x=*/true, /*read_y=*/false);
        const Long64_t n = (max_events > 0) ? std::min<Long64_t>(reader.Entries(), max_events) : reader.Entries();

        TH1D hx("hx", "xg distribution", 400, label - 60.0, label + 60.0);
        hx.SetDirectory(nullptr);

        std::vector<double> mcx_sum(FROST::kNfibX, 0.0);

        for (Long64_t e = 0; e < n; ++e) {
          reader.GetEntry(e);
          const double* px = reader.PhotonX();

          const double xg = FROST::ComputeXg(px, alpha);
          hx.Fill(xg);

          for (int i = 0; i < FROST::kNfibX; ++i) {
            mcx_sum[i] += px[i];
          }
        }

        xg_mean[q] = hx.GetMean();
        xg_err[q] = hx.GetMeanError();

        if (n > 0) {
          for (int i = 0; i < FROST::kNfibX; ++i) {
            mcx_avg[i][q] = mcx_sum[i] / static_cast<double>(n);
          }
        }

        std::cerr << "[single/x] label=" << label << " entries=" << n << " mean=" << xg_mean[q]
                  << " err=" << xg_err[q] << "\n";
      }

      std::vector<double> yg_mean(FROST::kNy, 0.0);
      std::vector<double> yg_err(FROST::kNy, 0.0);
      std::vector<std::vector<double>> mcy_avg(FROST::kNfibY, std::vector<double>(FROST::kNy, 0.0));

      for (int q = 0; q < FROST::kNy; ++q) {
        const int label = FROST::kYLabel[q];
        std::unordered_map<std::string, std::string> vars{{"label", std::to_string(label)}};
        const std::string in_path = FROST::ApplyTemplate(pat_y, vars);

        FROST::TreeReader reader(in_path, "wls", /*read_x=*/false, /*read_y=*/true);
        const Long64_t n = (max_events > 0) ? std::min<Long64_t>(reader.Entries(), max_events) : reader.Entries();

        TH1D hy("hy", "yg distribution", 400, label - 60.0, label + 60.0);
        hy.SetDirectory(nullptr);

        std::vector<double> mcy_sum(FROST::kNfibY, 0.0);

        for (Long64_t e = 0; e < n; ++e) {
          reader.GetEntry(e);
          const double* py = reader.PhotonY();

          const double yg = FROST::ComputeYg(py, alpha);
          hy.Fill(yg);

          for (int i = 0; i < FROST::kNfibY; ++i) {
            mcy_sum[i] += py[i];
          }
        }

        yg_mean[q] = hy.GetMean();
        yg_err[q] = hy.GetMeanError();

        if (n > 0) {
          for (int i = 0; i < FROST::kNfibY; ++i) {
            mcy_avg[i][q] = mcy_sum[i] / static_cast<double>(n);
          }
        }

        std::cerr << "[single/y] label=" << label << " entries=" << n << " mean=" << yg_mean[q]
                  << " err=" << yg_err[q] << "\n";
      }

      // ----------------------------------------------------------------------
      // Build completed (1 mm grid) mapping graphs and response graphs for
      // normal incidence
      //
      // 1) gx/gy: xg -> true position (mm)
      //    - Fill at label points
      //    - Interpolate the remaining 1-mm points in each 10-mm segment
      //    - Enforce symmetry for negative positions
      //    - Clamp endpoints to +/-1000
      //
      // 2) gmcx_i/gmcy_i: expected light yield per fiber as a function of true position (mm)
      //    - Fill at label points (positive side only is enough for deg=0 branch)
      //    - Interpolate inside 10-mm segments
      //    - Mirror negative side from positive side
      //    - Clamp negative values to small positive (0.01)
      // ----------------------------------------------------------------------

      // ---- (1) gx (xg -> x_mm), gy (yg -> y_mm) ----
      const int xg_size = FROST::kNfibX * 10 + 1;     // 1321, positions [-660..+660]
      const int xg_center = FROST::kNfibX * 5;        // 660
      std::vector<double> posxg(xg_size, 0.0);
      for (int i = 0; i < xg_size; ++i) posxg[i] = static_cast<double>(i - xg_center);

      std::vector<double> xgmc(xg_size, 0.0);
      for (int q = 0; q < FROST::kNx; ++q) {
        const int x = FROST::kXLabel[q];
        const int idx = xg_center + x;
        if (0 <= idx && idx < xg_size) xgmc[idx] = xg_mean[q];
      }
      const int idx_x10 = xg_center + 10;
      for (int i = 10; i <= 630; i += 10) {
        const int idx_i = xg_center + i;
        const int idx_i10 = xg_center + (i + 10);
        for (int j = 1; j <= 9; ++j) {
          const int idx_j = xg_center + j;
          const int idx_ij = xg_center + (i + j);
          xgmc[idx_ij] = xgmc[idx_i] + xgmc[idx_j] * (xgmc[idx_i10] - xgmc[idx_i]) / xgmc[idx_x10];
        }
      }
      for (int i = 0; i < xg_center; ++i) {
        xgmc[xg_center - i] = -xgmc[xg_center + i];
      }
      xgmc[xg_size - 1] = 1000.0;
      xgmc[0] = -1000.0;

      TGraphErrors gx(xg_size);
      gx.SetName("gx");
      for (int i = 0; i < xg_size; ++i) {
        gx.SetPoint(i, xgmc[i], posxg[i]);
      }

      const int yg_size = FROST::kNfibY * 10 + 1;     // 1401, positions [-700..+700]
      const int yg_center = FROST::kNfibY * 5;        // 700
      std::vector<double> posyg(yg_size, 0.0);
      for (int i = 0; i < yg_size; ++i) posyg[i] = static_cast<double>(i - yg_center);

      std::vector<double> ygmc(yg_size, 0.0);
      for (int q = 0; q < FROST::kNy; ++q) {
        const int y = FROST::kYLabel[q];
        const int idx = yg_center + y;
        if (0 <= idx && idx < yg_size) ygmc[idx] = yg_mean[q];
      }
      const int idx_y10 = yg_center + 10;
      for (int i = 10; i <= 670; i += 10) {
        const int idx_i = yg_center + i;
        const int idx_i10 = yg_center + (i + 10);
        for (int j = 1; j <= 9; ++j) {
          const int idx_j = yg_center + j;
          const int idx_ij = yg_center + (i + j);
          ygmc[idx_ij] = ygmc[idx_i] + ygmc[idx_j] * (ygmc[idx_i10] - ygmc[idx_i]) / ygmc[idx_y10];
        }
      }
      for (int i = 0; i < yg_center; ++i) {
        ygmc[yg_center - i] = -ygmc[yg_center + i];
      }
      ygmc[yg_size - 1] = 1000.0;
      ygmc[0] = -1000.0;

      TGraphErrors gy(yg_size);
      gy.SetName("gy");
      for (int i = 0; i < yg_size; ++i) {
        gy.SetPoint(i, ygmc[i], posyg[i]);
      }

      // ---- (2) gmcx_i (x_mm -> expected light yield), gmcy_i (y_mm -> expected light yield) ----
      // mc tables use 1-mm grid WITHOUT the +/- (fiber*10+1) endpoints:
      //   X: [-659..+659] => 1319 points
      //   Y: [-699..+699] => 1399 points
      const int nmcpx = FROST::kNmcPX;     // 1319
      const int nmcpy = FROST::kNmcPY;     // 1399
      const int mcx_center = (nmcpx - 1) / 2;   // 659
      const int mcy_center = (nmcpy - 1) / 2;   // 699

      std::vector<double> posx(nmcpx, 0.0);
      for (int i = 0; i < nmcpx; ++i) posx[i] = static_cast<double>(i - FROST::kNfibX * 5 + 1);
      std::vector<double> posy(nmcpy, 0.0);
      for (int i = 0; i < nmcpy; ++i) posy[i] = static_cast<double>(i - FROST::kNfibY * 5 + 1);

      // Build mcx[fiber][posIndex]
      std::vector<std::vector<double>> mcx(FROST::kNfibX, std::vector<double>(nmcpx, 0.0));
      for (int q = 0; q < FROST::kNx; ++q) {
        const int x = FROST::kXLabel[q];
        const int idx = mcx_center + x;
        if (idx < 0 || idx >= nmcpx) continue;
        for (int k = 0; k < FROST::kNfibX; ++k) {
          mcx[k][idx] = mcx_avg[k][q];
        }
      }
      // Interpolate inside each 10-mm segment (deg<=0 branch in xgsinglehit.C):
      //   if (k - i/10 >= 0) mcx[k][center+i+j] = mcx[k-i/10][center+j]
      //   else average endpoints.
      for (int i = 10; i <= 630; i += 10) {
        for (int j = 1; j <= 9; ++j) {
          for (int k = 0; k < FROST::kNfibX; ++k) {
            const int dst = mcx_center + i + j;
            const int src = mcx_center + j;
            if (dst < 0 || dst >= nmcpx) continue;
            if (k - i / 10 >= 0) {
              mcx[k][dst] = mcx[k - i / 10][src];
            } else {
              mcx[k][dst] = 0.5 * (mcx[k][mcx_center + i] + mcx[k][mcx_center + i + 10]);
            }
          }
        }
      }
      // Mirror negative side:
      for (int i = 0; i < mcx_center; ++i) {
        for (int k = 0; k < FROST::kNfibX; ++k) {
          mcx[k][i] = mcx[FROST::kNfibX - 1 - k][nmcpx - 1 - i];
        }
      }
      for (int k = 0; k < FROST::kNfibX; ++k) {
        for (int i = 0; i < nmcpx; ++i) {
          if (mcx[k][i] < 0.0) mcx[k][i] = 0.01;
        }
      }

      // Build mcy[fiber][posIndex]
      std::vector<std::vector<double>> mcy(FROST::kNfibY, std::vector<double>(nmcpy, 0.0));
      for (int q = 0; q < FROST::kNy; ++q) {
        const int y = FROST::kYLabel[q];
        const int idx = mcy_center + y;
        if (idx < 0 || idx >= nmcpy) continue;
        for (int k = 0; k < FROST::kNfibY; ++k) {
          mcy[k][idx] = mcy_avg[k][q];
        }
      }
      for (int i = 10; i <= 670; i += 10) {
        for (int j = 1; j <= 9; ++j) {
          for (int k = 0; k < FROST::kNfibY; ++k) {
            const int dst = mcy_center + i + j;
            const int src = mcy_center + j;
            if (dst < 0 || dst >= nmcpy) continue;
            if (k - i / 10 >= 0) {
              mcy[k][dst] = mcy[k - i / 10][src];
            } else {
              mcy[k][dst] = 0.5 * (mcy[k][mcy_center + i] + mcy[k][mcy_center + i + 10]);
            }
          }
        }
      }
      for (int i = 0; i < mcy_center; ++i) {
        for (int k = 0; k < FROST::kNfibY; ++k) {
          mcy[k][i] = mcy[FROST::kNfibY - 1 - k][nmcpy - 1 - i];
        }
      }
      for (int k = 0; k < FROST::kNfibY; ++k) {
        for (int i = 0; i < nmcpy; ++i) {
          if (mcy[k][i] < 0.0) mcy[k][i] = 0.01;
        }
      }

      // Convert mc tables to TGraph (x_mm -> expected).
      std::vector<std::unique_ptr<TGraph>> gmcx(FROST::kNfibX);
      for (int k = 0; k < FROST::kNfibX; ++k) {
        auto g = std::make_unique<TGraph>(nmcpx);
        g->SetName((std::string("gmcx_") + std::to_string(k)).c_str());
        for (int i = 0; i < nmcpx; ++i) g->SetPoint(i, posx[i], mcx[k][i]);
        gmcx[k] = std::move(g);
      }
      std::vector<std::unique_ptr<TGraph>> gmcy(FROST::kNfibY);
      for (int k = 0; k < FROST::kNfibY; ++k) {
        auto g = std::make_unique<TGraph>(nmcpy);
        g->SetName((std::string("gmcy_") + std::to_string(k)).c_str());
        for (int i = 0; i < nmcpy; ++i) g->SetPoint(i, posy[i], mcy[k][i]);
        gmcy[k] = std::move(g);
      }

      TFile out(out_path.c_str(), "RECREATE");
      if (out.IsZombie()) throw std::runtime_error("Failed to create output file: " + out_path);

      TParameter<double> p_alpha("alpha", alpha);
      p_alpha.Write();

      gx.Write();
      gy.Write();

      for (auto& g : gmcx) g->Write();
      for (auto& g : gmcy) g->Write();

      out.Close();
      std::cerr << "Wrote single-hit map: " << out_path << "\n";
      return 0;

    } else if (mode == "two") {
      cli.Require("in-x-pattern");
      cli.Require("in-y-pattern");

      const std::string pat_x = cli.Get("in-x-pattern");
      const std::string pat_y = cli.Get("in-y-pattern");

      // Store mean xg for each label and each category.
      std::vector<double> xg_xm(FROST::kNx, 0.0), xg_xm_err(FROST::kNx, 0.0);
      std::vector<double> xg_xp(FROST::kNx, 0.0), xg_xp_err(FROST::kNx, 0.0);

      for (int q = 0; q < FROST::kNx; ++q) {
        const int label = FROST::kXLabel[q];
        std::unordered_map<std::string, std::string> vars{{"label", std::to_string(label)}};
        const std::string in_path = FROST::ApplyTemplate(pat_x, vars);

        FROST::TreeReader reader(in_path, "wls", /*read_x=*/true, /*read_y=*/true);
        const Long64_t n = (max_events > 0) ? std::min<Long64_t>(reader.Entries(), max_events) : reader.Entries();

        TH1D hxm("hxm", "xg (minus candidate)", 400, label - 80.0, label + 80.0);
        TH1D hxp("hxp", "xg (plus candidate)", 400, label - 80.0, label + 80.0);
        hxm.SetDirectory(nullptr);
        hxp.SetDirectory(nullptr);

        TH1D hx("hx", "hx", FROST::kNfibX, -FROST::kNfibX * 5.0, FROST::kNfibX * 5.0);
        hx.SetDirectory(nullptr);

        for (Long64_t e = 0; e < n; ++e) {
          reader.GetEntry(e);
          if (!reader.TrueX().has_value() || !reader.TrueX2().has_value()) {
            throw std::runtime_error("Two-hit X calibration file is missing x_muon1/x_muon2 branches: " + in_path);
          }

          const double x = reader.TrueX().value();
          const double x2 = reader.TrueX2().value();
          const double* px = reader.PhotonX();

          // Fill histogram from current event yields.
          for (int i = 0; i < FROST::kNfibX; ++i) {
            hx.SetBinContent(i + 1, px[i]);
          }

          //  1) DivideGroup1 (low threshold)
          //  2)  DivideGroup2
          //  3) Always run MultiReconstructionX to obtain xgmulti[0], xgmulti[1]
          std::vector<std::vector<int>> groupx;
          FROST::DivideGroup1(&hx, groupx, mh_cfg);
          FROST::DivideGroup2(&hx, groupx, mh_cfg);

          int group2_size = (groupx.size() >= 2) ? static_cast<int>(groupx[1].size()) : 0;

          std::vector<double> xgmulti;
          std::vector<int> rectype;
          FROST::MultiReconstructionX(&hx, xgmulti, rectype, groupx, mh_cfg);
          if (xgmulti.size() < 2) continue;  // safety

          const double xrecmuon = (x < x2) ? xgmulti[0] : xgmulti[1];

          // Original code fills histograms only when group2_size == 0:
          //   if (x > x2) fill 1xp with xrecmuon; else fill 1xm with xrecmuon
          if (group2_size == 0) {
            if (x > x2) {
              hxp.Fill(xrecmuon);
              // if you also keep vectors like xg1xp[q], push_back(xrecmuon) here
            } else {
              hxm.Fill(xrecmuon);
              // if you also keep vectors like xg1xm[q], push_back(xrecmuon) here
            }
          }
        }

        xg_xm[q] = hxm.GetMean();
        xg_xm_err[q] = hxm.GetMeanError();
        xg_xp[q] = hxp.GetMean();
        xg_xp_err[q] = hxp.GetMeanError();

        std::cerr << "[two/x] label=" << label << " entries=" << n << " xm_mean=" << xg_xm[q]
                  << " xp_mean=" << xg_xp[q] << "\n";
      }

      std::vector<double> yg_ym(FROST::kNy, 0.0), yg_ym_err(FROST::kNy, 0.0);
      std::vector<double> yg_yp(FROST::kNy, 0.0), yg_yp_err(FROST::kNy, 0.0);

      for (int q = 0; q < FROST::kNy; ++q) {
        const int label = FROST::kYLabel[q];
        std::unordered_map<std::string, std::string> vars{{"label", std::to_string(label)}};
        const std::string in_path = FROST::ApplyTemplate(pat_y, vars);

        FROST::TreeReader reader(in_path, "wls", /*read_x=*/true, /*read_y=*/true);
        const Long64_t n = (max_events > 0) ? std::min<Long64_t>(reader.Entries(), max_events) : reader.Entries();

        TH1D hym("hym", "yg (minus candidate)", 400, label - 80.0, label + 80.0);
        TH1D hyp("hyp", "yg (plus candidate)", 400, label - 80.0, label + 80.0);
        hym.SetDirectory(nullptr);
        hyp.SetDirectory(nullptr);

        TH1D hy("hy", "hy", FROST::kNfibY, -FROST::kNfibY * 5.0, FROST::kNfibY * 5.0);
        hy.SetDirectory(nullptr);

        for (Long64_t e = 0; e < n; ++e) {
          reader.GetEntry(e);
          if (!reader.TrueY().has_value() || !reader.TrueY2().has_value()) {
            throw std::runtime_error("Two-hit Y calibration file is missing y/y2 branches: " + in_path);
          }

          const double y = reader.TrueY().value();
          const double y2 = reader.TrueY2().value();
          const double* py = reader.PhotonY();

          // Fill histogram from current event yields.
          for (int i = 0; i < FROST::kNfibY; ++i) {
            hy.SetBinContent(i + 1, py[i]);
          }

          //  1) DivideGroup1 (low threshold)
          //  2) DivideGroup2
          //  3) Always run MultiReconstructionY to obtain ygmulti[0], ygmulti[1]
          std::vector<std::vector<int>> groupy;
          FROST::DivideGroup1(&hy, groupy, mh_cfg);
          FROST::DivideGroup2(&hy, groupy, mh_cfg);

          int group2_size = (groupy.size() >= 2) ? static_cast<int>(groupy[1].size()) : 0;

          std::vector<double> ygmulti;
          std::vector<int> rectype;
          FROST::MultiReconstructionY(&hy, ygmulti, rectype, groupy, mh_cfg);
          if (ygmulti.size() < 2) continue;  // safety

          const double yrecmuon = (y < y2) ? ygmulti[0] : ygmulti[1];

          // Original code fills histograms only when group2_size == 0:
          //   if (y > y2) fill 1yp with yrecmuon; else fill 1ym with yrecmuon
          if (group2_size == 0) {
            if (y > y2) {
              hyp.Fill(yrecmuon);
              // if you also keep vectors like yg1yp[q], push_back(yrecmuon) here
            } else {
              hym.Fill(yrecmuon);
              // if you also keep vectors like yg1ym[q], push_back(yrecmuon) here
            }
          }
        }

        yg_ym[q] = hym.GetMean();
        yg_ym_err[q] = hym.GetMeanError();
        yg_yp[q] = hyp.GetMean();
        yg_yp_err[q] = hyp.GetMeanError();

        std::cerr << "[two/y] label=" << label << " entries=" << n << " ym_mean=" << yg_ym[q]
                  << " yp_mean=" << yg_yp[q] << "\n";
      }

      // ----------------------------------------------------------------------
      // Build the completed (1 mm grid) mapping graphs, matching the original
      // FROST_reconstruction.cc logic.
      //
      // The original code first fills the mapping at the calibration labels
      // (0..10, 20, 30, ..., 640, 641..659), then interpolates the remaining
      // 1-mm points inside each 10-mm segment using the shape of 0..10 mm.
      // It also enforces symmetry:
      //   xm(-x) = -xp(+x)
      //   xp(-x) = -xm(+x)
      // and clamps endpoints to +/-1000.
      // ----------------------------------------------------------------------

      const int x_center = FROST::kNfibX * 5;          // 660
      const int x_size = FROST::kNmcPX + 2;            // 1321, positions [-660..+660]
      std::vector<double> posxg(x_size, 0.0);
      for (int i = 0; i < x_size; ++i) posxg[i] = static_cast<double>(i - x_center);

      std::vector<double> xgmc1xm(x_size, 0.0);
      std::vector<double> xgmc1xp(x_size, 0.0);

      // Fill known label points on + side (including 0..10, 20.., 640, 641..659).
      for (int q = 0; q < FROST::kNx; ++q) {
        const int x = FROST::kXLabel[q];
        const int idx = x_center + x;
        if (idx < 0 || idx >= x_size) continue;
        xgmc1xm[idx] = xg_xm[q];
        xgmc1xp[idx] = xg_xp[q];
      }

      // Interpolate 1-mm points inside each 10-mm segment using the 0..10 mm shape.
      //   f(i+j) = f(i) + f(j) * (f(i+10) - f(i)) / f(10)
      // for i=10..630 step 10 and j=1..9.
      const int idx_x10 = x_center + 10;
      for (int i = 10; i <= 630; i += 10) {
        const int idx_i = x_center + i;
        const int idx_i10 = x_center + (i + 10);
        for (int j = 1; j <= 9; ++j) {
          const int idx_j = x_center + j;
          const int idx_ij = x_center + (i + j);
          xgmc1xm[idx_ij] = xgmc1xm[idx_i] + xgmc1xm[idx_j] * (xgmc1xm[idx_i10] - xgmc1xm[idx_i]) / xgmc1xm[idx_x10];
          xgmc1xp[idx_ij] = xgmc1xp[idx_i] + xgmc1xp[idx_j] * (xgmc1xp[idx_i10] - xgmc1xp[idx_i]) / xgmc1xp[idx_x10];
        }
      }

      // Symmetry for negative x.
      for (int i = 0; i < x_center; ++i) {
        xgmc1xm[x_center - i] = -xgmc1xp[x_center + i];
        xgmc1xp[x_center - i] = -xgmc1xm[x_center + i];
      }

      // Clamp endpoints.
      xgmc1xm[x_size - 1] = 1000.0;
      xgmc1xm[0] = -1000.0;
      xgmc1xp[x_size - 1] = 1000.0;
      xgmc1xp[0] = -1000.0;

      // Create TGraphErrors: x-axis = xg, y-axis = true position.
      TGraphErrors g1xm(x_size);
      TGraphErrors g1xp(x_size);
      g1xm.SetName("g1xm");
      g1xp.SetName("g1xp");
      for (int i = 0; i < x_size; ++i) {
        g1xm.SetPoint(i, xgmc1xm[i], posxg[i]);
        g1xp.SetPoint(i, xgmc1xp[i], posxg[i]);
      }

      const int y_center = FROST::kNfibY * 5;          // 700
      const int y_size = FROST::kNmcPY + 2;            // 1401, positions [-700..+700]
      std::vector<double> posyg(y_size, 0.0);
      for (int i = 0; i < y_size; ++i) posyg[i] = static_cast<double>(i - y_center);

      std::vector<double> ygmc1ym(y_size, 0.0);
      std::vector<double> ygmc1yp(y_size, 0.0);

      for (int q = 0; q < FROST::kNy; ++q) {
        const int y = FROST::kYLabel[q];
        const int idx = y_center + y;
        if (idx < 0 || idx >= y_size) continue;
        ygmc1ym[idx] = yg_ym[q];
        ygmc1yp[idx] = yg_yp[q];
      }

      const int idx_y10 = y_center + 10;
      for (int i = 10; i <= 670; i += 10) {
        const int idx_i = y_center + i;
        const int idx_i10 = y_center + (i + 10);
        for (int j = 1; j <= 9; ++j) {
          const int idx_j = y_center + j;
          const int idx_ij = y_center + (i + j);
          ygmc1ym[idx_ij] = ygmc1ym[idx_i] + ygmc1ym[idx_j] * (ygmc1ym[idx_i10] - ygmc1ym[idx_i]) / ygmc1ym[idx_y10];
          ygmc1yp[idx_ij] = ygmc1yp[idx_i] + ygmc1yp[idx_j] * (ygmc1yp[idx_i10] - ygmc1yp[idx_i]) / ygmc1yp[idx_y10];
        }
      }

      for (int i = 0; i < y_center; ++i) {
        ygmc1ym[y_center - i] = -ygmc1yp[y_center + i];
        ygmc1yp[y_center - i] = -ygmc1ym[y_center + i];
      }

      ygmc1ym[y_size - 1] = 1000.0;
      ygmc1ym[0] = -1000.0;
      ygmc1yp[y_size - 1] = 1000.0;
      ygmc1yp[0] = -1000.0;

      TGraphErrors g1ym(y_size);
      TGraphErrors g1yp(y_size);
      g1ym.SetName("g1ym");
      g1yp.SetName("g1yp");
      for (int i = 0; i < y_size; ++i) {
        g1ym.SetPoint(i, ygmc1ym[i], posyg[i]);
        g1yp.SetPoint(i, ygmc1yp[i], posyg[i]);
      }

      std::vector<PointWithErr> g1ym_pts;
      std::vector<PointWithErr> g1yp_pts;
      g1ym_pts.reserve(FROST::kNy * 2);
      g1yp_pts.reserve(FROST::kNy * 2);

      for (int q = 0; q < FROST::kNy; ++q) {
        const int pos = FROST::kYLabel[q];
        if (pos == 0) {
          g1ym_pts.push_back({yg_ym[q], 0.0, yg_ym_err[q], 0.0});
          g1yp_pts.push_back({yg_yp[q], 0.0, yg_yp_err[q], 0.0});
        } else {
          g1ym_pts.push_back({yg_ym[q], static_cast<double>(pos), yg_ym_err[q], 0.0});
          g1yp_pts.push_back({yg_yp[q], static_cast<double>(pos), yg_yp_err[q], 0.0});

          g1ym_pts.push_back({-yg_yp[q], -static_cast<double>(pos), yg_yp_err[q], 0.0});
          g1yp_pts.push_back({-yg_ym[q], -static_cast<double>(pos), yg_ym_err[q], 0.0});
        }
      }

      TFile out(out_path.c_str(), "RECREATE");
      if (out.IsZombie()) throw std::runtime_error("Failed to create output file: " + out_path);

      TParameter<double> p_alpha("alpha", alpha);
      p_alpha.Write();
      TParameter<double> p_thr1("threshold1", threshold1);
      TParameter<double> p_thr2("threshold2", threshold2);
      p_thr1.Write();
      p_thr2.Write();

      g1xm.Write();
      g1xp.Write();
      g1ym.Write();
      g1yp.Write();

      out.Close();
      std::cerr << "Wrote two-hit map: " << out_path << "\n";
      return 0;

    } else {
      throw std::runtime_error("Unknown --mode: " + mode + " (expected 'single' or 'two')");
    }

  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    std::cerr << "\n";
    std::cerr << "Usage examples:\n";
    std::cerr << "  mapfunction_tool --mode single --alpha 1.0 --in-x-pattern <pattern> --in-y-pattern <pattern> --out single_map.root\n";
    std::cerr << "  mapfunction_tool --mode two    --alpha 1.0 --in-x-pattern <pattern> --in-y-pattern <pattern> --out two_map.root\n";
    return 1;
  }
}
