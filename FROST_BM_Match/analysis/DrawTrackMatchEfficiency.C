#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TChain.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <TString.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <streambuf>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <Math/QuantFuncMathCore.h>

// ------------------------------------------------------------
// Draw BM-FROST matching efficiencies and dx/dy distributions.
//
// Usage in ROOT:
//   .L DrawTrackMatchEfficiency.C
//   DrawTrackMatchEfficiency("/path/to/input_dir",
//                            "/path/to/output.pdf",
//                            "/path/to/output.log");
//
//   std::vector<std::string> exclude = {"hoge.root", "hage.root"};
//   DrawTrackMatchEfficiency("/path/to/input_dir",
//                            "/path/to/output.pdf",
//                            "/path/to/output.log",
//                            exclude);
//
// Excluded files are specified by file name relative to inputDir.
//
// Denominator selection (counted per track):
//   (trackmatch_ninja_track_type == 1)
//   && abs(trackmatch_expected_x) < 460
//   && abs(trackmatch_expected_y) < 400
//   && bsd_good_spill_flag != 0
//   && matched == 1
//
// Numerator for BM-FROST matching efficiency:
//   denominator condition && trackmatch_has_match == 1
//
// Numerator for FROST is_hit efficiency:
//   denominator condition && trackmatch_frost_is_hit == 1
//
// Angle definition for efficiency:
//   atan(sqrt(trackmatch_baby_mind_tangent_x^2
//             + trackmatch_baby_mind_tangent_y^2)) [deg]
//
// Angle bins [deg]:
//   [0,5), [5,10), [10,15), [15,20), [20,25),
//   [25,30), [30,35), [35,40), [40,50)
//
// Output pages in the PDF:
//   1. BM-FROST matching efficiency vs angle
//   2. FROST is_hit efficiency vs angle
//   3. dx distribution for all denominator tracks
//   4. dy distribution for all denominator tracks
//   5-13.  dx distributions in bins of atan(abs(tan_x)) [deg]
//   14-22. dx distributions in bins of atan(sqrt(tan_x^2 + tan_y^2)) [deg]
//   23-31. dy distributions in bins of atan(abs(tan_y)) [deg]
//   32-40. dy distributions in bins of atan(sqrt(tan_x^2 + tan_y^2)) [deg]
//
// The terminal output is also written to the specified log file.
// ------------------------------------------------------------

namespace {
  constexpr int kNBins = 9;
  const double kAngleBins[kNBins + 1] = {
    0.0, 5.0, 10.0, 15.0, 20.0,
    25.0, 30.0, 35.0, 40.0, 50.0
  };

  bool HasRootExtension(const std::string &name) {
    return name.size() >= 5 && name.substr(name.size() - 5) == ".root";
  }

  bool IsExcludedFile(const std::string &fileName,
                      const std::vector<std::string> &excludedFiles) {
    return std::find(excludedFiles.begin(), excludedFiles.end(), fileName)
           != excludedFiles.end();
  }

  // Stream buffer that writes to two destinations at the same time.
  class TeeBuf : public std::streambuf {
  public:
    TeeBuf(std::streambuf *sb1, std::streambuf *sb2) : sb1_(sb1), sb2_(sb2) {}

  protected:
    virtual int overflow(int c) override {
      if (c == EOF) return !EOF;
      const int r1 = sb1_ ? sb1_->sputc(c) : c;
      const int r2 = sb2_ ? sb2_->sputc(c) : c;
      return (r1 == EOF || r2 == EOF) ? EOF : c;
    }

    virtual int sync() override {
      const int r1 = sb1_ ? sb1_->pubsync() : 0;
      const int r2 = sb2_ ? sb2_->pubsync() : 0;
      return (r1 == 0 && r2 == 0) ? 0 : -1;
    }

  private:
    std::streambuf *sb1_;
    std::streambuf *sb2_;
  };
}

void DrawTrackMatchEfficiency(
    const char *inputDir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/rootfile_after_TrackMatch",
    const char *outputPdfPath = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/efficiency.pdf",
    const char *logFilePath = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/efficiency.log",
    const std::vector<std::string> &excludedFiles = std::vector<std::string>{"BMPM_track_2025-11-29_13-46-59_Run0_afterTrackMatch.root", "BMPM_track_2025-11-30_13-11-36_Run0_afterTrackMatch.root"}
    ) {

  std::ofstream logFile(logFilePath);
  if (!logFile) {
    std::cerr << "Error: cannot open log file: " << logFilePath << std::endl;
    return;
  }

  TeeBuf teeCoutBuf(std::cout.rdbuf(), logFile.rdbuf());
  TeeBuf teeCerrBuf(std::cerr.rdbuf(), logFile.rdbuf());
  std::streambuf *oldCoutBuf = std::cout.rdbuf(&teeCoutBuf);
  std::streambuf *oldCerrBuf = std::cerr.rdbuf(&teeCerrBuf);

  std::cout << "Log file: " << logFilePath << std::endl;
  // Build a chain from all ROOT files in the directory.
  TChain chain("match_info");

  TSystemDirectory dir("input_dir", inputDir);
  TList *fileList = dir.GetListOfFiles();
  if (!fileList) {
    std::cerr << "Error: cannot open directory: " << inputDir << std::endl;
    std::cout.rdbuf(oldCoutBuf);
    std::cerr.rdbuf(oldCerrBuf);
    return;
  }

  int nFilesAdded = 0;
  int nFilesExcluded = 0;

  TIter next(fileList);
  while (TObject *obj = next()) {
    auto *sysFile = dynamic_cast<TSystemFile *>(obj);
    if (!sysFile) continue;
    if (sysFile->IsDirectory()) continue;

    const std::string fileName = sysFile->GetName();
    if (!HasRootExtension(fileName)) continue;

    if (IsExcludedFile(fileName, excludedFiles)) {
      std::cout << "Excluded file: " << fileName << std::endl;
      ++nFilesExcluded;
      continue;
    }

    const std::string fullPath = std::string(inputDir) + "/" + fileName;
    if (chain.Add(fullPath.c_str(), 0) > 0) {
      ++nFilesAdded;
      std::cout << "Added file: " << fileName << std::endl;
    }
  }

  if (nFilesAdded == 0) {
    std::cerr << "Error: no ROOT files with tree 'match_info' were added from "
              << inputDir << std::endl;
    std::cout.rdbuf(oldCoutBuf);
    std::cerr.rdbuf(oldCerrBuf);
    return;
  }

  std::cout << "Added " << nFilesAdded << " ROOT files." << std::endl;
  std::cout << "Excluded " << nFilesExcluded << " ROOT files." << std::endl;
  std::cout << "Total spills in chain: " << chain.GetEntries() << std::endl;

  // Spill-wise branches.
  Int_t matched = 0;
  Int_t bsd_good_spill_flag = 0;

  // Track-wise branches stored as vectors.
  std::vector<int> *trackmatch_has_match = nullptr;
  std::vector<int> *trackmatch_ninja_track_type = nullptr;
  std::vector<double> *trackmatch_expected_x = nullptr;
  std::vector<double> *trackmatch_expected_y = nullptr;
  std::vector<double> *trackmatch_dx = nullptr;
  std::vector<double> *trackmatch_dy = nullptr;
  std::vector<double> *trackmatch_baby_mind_tangent_x = nullptr;
  std::vector<double> *trackmatch_baby_mind_tangent_y = nullptr;
  std::vector<int> *trackmatch_frost_is_hit = nullptr;

  chain.SetBranchAddress("matched", &matched);
  chain.SetBranchAddress("bsd_good_spill_flag", &bsd_good_spill_flag);
  chain.SetBranchAddress("trackmatch_has_match", &trackmatch_has_match);
  chain.SetBranchAddress("trackmatch_ninja_track_type", &trackmatch_ninja_track_type);
  chain.SetBranchAddress("trackmatch_expected_x", &trackmatch_expected_x);
  chain.SetBranchAddress("trackmatch_expected_y", &trackmatch_expected_y);
  chain.SetBranchAddress("trackmatch_dx", &trackmatch_dx);
  chain.SetBranchAddress("trackmatch_dy", &trackmatch_dy);
  chain.SetBranchAddress("trackmatch_baby_mind_tangent_x", &trackmatch_baby_mind_tangent_x);
  chain.SetBranchAddress("trackmatch_baby_mind_tangent_y", &trackmatch_baby_mind_tangent_y);
  chain.SetBranchAddress("trackmatch_frost_is_hit", &trackmatch_frost_is_hit);

  // Histograms for numerator and denominator counts.
  auto *hDen = new TH1D("hDen", "Denominator;Angle [deg];Tracks",
                        kNBins, kAngleBins);
  auto *hNum = new TH1D("hNum", "Numerator;Angle [deg];Tracks",
                        kNBins, kAngleBins);
  auto *hDenIsHit = new TH1D("hDenIsHit", "Denominator for is_hit;Angle [deg];Tracks",
                             kNBins, kAngleBins);
  auto *hNumIsHit = new TH1D("hNumIsHit", "Numerator for is_hit;Angle [deg];Tracks",
                             kNBins, kAngleBins);

  // dx/dy distributions for all denominator tracks.
  auto *hDxAll = new TH1D("hDxAll", ";dx [mm];Number of events", 200, -500.0, 500.0);
  auto *hDyAll = new TH1D("hDyAll", ";dy [mm];Number of events", 200, -500.0, 500.0);

  // dx distributions binned by atan(abs(tan_x)) [deg].
  std::vector<TH1D*> hDxByAngleX;
  std::vector<TH1D*> hDxByAngleTot;
  std::vector<TH1D*> hDyByAngleY;
  std::vector<TH1D*> hDyByAngleTot;
  for (int i = 0; i < kNBins; ++i) {
    hDxByAngleX.push_back(new TH1D(
        Form("hDx_bin%d", i),
        ";dx [mm];Number of events",
        200, -500.0, 500.0));
    hDxByAngleTot.push_back(new TH1D(
        Form("hDxTot_bin%d", i),
        ";dx [mm];Number of events",
        200, -500.0, 500.0));
    hDyByAngleY.push_back(new TH1D(
        Form("hDy_bin%d", i),
        ";dy [mm];Number of events",
        200, -500.0, 500.0));
    hDyByAngleTot.push_back(new TH1D(
        Form("hDyTot_bin%d", i),
        ";dy [mm];Number of events",
        200, -500.0, 500.0));
  }

  hDen->Sumw2();
  hNum->Sumw2();
  hDenIsHit->Sumw2();
  hNumIsHit->Sumw2();
  hDxAll->Sumw2();
  hDyAll->Sumw2();
  for (int i = 0; i < kNBins; ++i) {
    hDxByAngleX[i]->Sumw2();
    hDxByAngleTot[i]->Sumw2();
    hDyByAngleY[i]->Sumw2();
    hDyByAngleTot[i]->Sumw2();
  }

  const Long64_t nEntries = chain.GetEntries();
  for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
    chain.GetEntry(iEntry);

    // Spill-level conditions.
    if (matched != 1) continue;
    if (bsd_good_spill_flag == 0) continue;

    if (!trackmatch_has_match ||
        !trackmatch_ninja_track_type ||
        !trackmatch_expected_x ||
        !trackmatch_expected_y ||
        !trackmatch_dx ||
        !trackmatch_dy ||
        !trackmatch_baby_mind_tangent_x ||
        !trackmatch_baby_mind_tangent_y ||
        !trackmatch_frost_is_hit) {
      continue;
    }

    const std::size_t nTracks = trackmatch_has_match->size();
    if (trackmatch_ninja_track_type->size() != nTracks ||
        trackmatch_expected_x->size() != nTracks ||
        trackmatch_expected_y->size() != nTracks ||
        trackmatch_dx->size() != nTracks ||
        trackmatch_dy->size() != nTracks ||
        trackmatch_baby_mind_tangent_x->size() != nTracks ||
        trackmatch_baby_mind_tangent_y->size() != nTracks ||
        trackmatch_frost_is_hit->size() != nTracks) {
      std::cerr << "Warning: vector size mismatch at spill entry "
                << iEntry << ". Skip this spill." << std::endl;
      continue;
    }

    // Count each track independently.
    for (std::size_t iTrack = 0; iTrack < nTracks; ++iTrack) {
      const int ninjaTrackType = trackmatch_ninja_track_type->at(iTrack);
      const double expectedX = trackmatch_expected_x->at(iTrack);
      const double expectedY = trackmatch_expected_y->at(iTrack);
      const double dx = trackmatch_dx->at(iTrack);
      const double dy = trackmatch_dy->at(iTrack);
      const double tx = trackmatch_baby_mind_tangent_x->at(iTrack);
      const double ty = trackmatch_baby_mind_tangent_y->at(iTrack);

      // const bool passTrackType =
      //   (ninjaTrackType == 1 || ninjaTrackType == 2);
      const bool passTrackType =
        (ninjaTrackType == 1);
      const bool passPosition =
        (std::abs(expectedX) < 360.0 && std::abs(expectedY) < 500.0);

      if (!passTrackType) continue;
      if (!passPosition) continue;

      const double angleDeg =
        std::atan(std::sqrt(tx * tx + ty * ty)) * 180.0 / TMath::Pi();
      const double angleXDeg =
        std::atan(std::abs(tx)) * 180.0 / TMath::Pi();
      const double angleYDeg =
        std::atan(std::abs(ty)) * 180.0 / TMath::Pi();

      // Fill denominator per track.
      hDen->Fill(angleDeg);
      hDenIsHit->Fill(angleDeg);
      hDxAll->Fill(dx);
      hDyAll->Fill(dy);

      for (int iBin = 0; iBin < kNBins; ++iBin) {
        if (angleXDeg >= kAngleBins[iBin] && angleXDeg < kAngleBins[iBin + 1]) {
          hDxByAngleX[iBin]->Fill(dx);
          break;
        }
      }
      for (int iBin = 0; iBin < kNBins; ++iBin) {
        if (angleDeg >= kAngleBins[iBin] && angleDeg < kAngleBins[iBin + 1]) {
          hDxByAngleTot[iBin]->Fill(dx);
          break;
        }
      }
      for (int iBin = 0; iBin < kNBins; ++iBin) {
        if (angleYDeg >= kAngleBins[iBin] && angleYDeg < kAngleBins[iBin + 1]) {
          hDyByAngleY[iBin]->Fill(dy);
          break;
        }
      }
      for (int iBin = 0; iBin < kNBins; ++iBin) {
        if (angleDeg >= kAngleBins[iBin] && angleDeg < kAngleBins[iBin + 1]) {
          hDyByAngleTot[iBin]->Fill(dy);
          break;
        }
      }

      // Fill numerator per track.
      if (trackmatch_has_match->at(iTrack) == 1) {
        hNum->Fill(angleDeg);
      }
      if (trackmatch_frost_is_hit->at(iTrack) == 1) {
        hNumIsHit->Fill(angleDeg);
      }
    }
  }

  // Build efficiency histogram.
  // The histogram stores the central values only.
  // Asymmetric confidence intervals are calculated separately below.
  auto *hEff = new TH1D("hEff",
                        ";Baby MIND reconstructed angle [deg];Track matching efficiency",
                        kNBins, kAngleBins);
  auto *hEffIsHit = new TH1D("hEffIsHit",
                             ";Baby MIND reconstructed angle [deg];FROST hit efficiency",
                             kNBins, kAngleBins);

  // Arrays for asymmetric efficiency error bars.
  double x[kNBins];
  double y[kNBins];
  double exl[kNBins];
  double exh[kNBins];
  double eyl[kNBins];
  double eyh[kNBins];
  double xIsHit[kNBins];
  double yIsHit[kNBins];
  double exlIsHit[kNBins];
  double exhIsHit[kNBins];
  double eylIsHit[kNBins];
  double eyhIsHit[kNBins];

  // One-sigma equivalent central confidence level.
  const double alpha = 1.0 - 0.682689492137;

  for (int iBin = 1; iBin <= kNBins; ++iBin) {
    const double num = hNum->GetBinContent(iBin);
    const double den = hDen->GetBinContent(iBin);
    const double numIsHit = hNumIsHit->GetBinContent(iBin);
    const double denIsHit = hDenIsHit->GetBinContent(iBin);

    double eff = 0.0;
    double effIsHit = 0.0;
    if (den > 0.0) {
      eff = num / den;
    }
    if (denIsHit > 0.0) {
      effIsHit = numIsHit / denIsHit;
    }

    hEff->SetBinContent(iBin, eff);
    hEffIsHit->SetBinContent(iBin, effIsHit);

    const double low = hEff->GetXaxis()->GetBinLowEdge(iBin);
    const double high = hEff->GetXaxis()->GetBinUpEdge(iBin);
    const double center = 0.5 * (low + high);

    x[iBin - 1] = center;
    y[iBin - 1] = eff;
    exl[iBin - 1] = center - low;
    exh[iBin - 1] = high - center;
    eyl[iBin - 1] = 0.0;
    eyh[iBin - 1] = 0.0;

    xIsHit[iBin - 1] = center;
    yIsHit[iBin - 1] = effIsHit;
    exlIsHit[iBin - 1] = center - low;
    exhIsHit[iBin - 1] = high - center;
    eylIsHit[iBin - 1] = 0.0;
    eyhIsHit[iBin - 1] = 0.0;

    if (den > 0.0) {
      double lower = 0.0;
      double upper = 1.0;

      if (num > 0.0) {
        lower = ROOT::Math::beta_quantile(alpha / 2.0, num, den - num + 1.0);
      } else {
        lower = 0.0;
      }

      if (num < den) {
        upper = ROOT::Math::beta_quantile(1.0 - alpha / 2.0, num + 1.0, den - num);
      } else {
        upper = 1.0;
      }

      eyl[iBin - 1] = eff - lower;
      eyh[iBin - 1] = upper - eff;
    }

    if (denIsHit > 0.0) {
      double lowerIsHit = 0.0;
      double upperIsHit = 1.0;

      if (numIsHit > 0.0) {
        lowerIsHit = ROOT::Math::beta_quantile(alpha / 2.0,
                                               numIsHit,
                                               denIsHit - numIsHit + 1.0);
      } else {
        lowerIsHit = 0.0;
      }

      if (numIsHit < denIsHit) {
        upperIsHit = ROOT::Math::beta_quantile(1.0 - alpha / 2.0,
                                               numIsHit + 1.0,
                                               denIsHit - numIsHit);
      } else {
        upperIsHit = 1.0;
      }

      eylIsHit[iBin - 1] = effIsHit - lowerIsHit;
      eyhIsHit[iBin - 1] = upperIsHit - effIsHit;
    }
  }
  auto *gEff = new TGraphAsymmErrors(kNBins, x, y, exl, exh, eyl, eyh);
  auto *gEffIsHit = new TGraphAsymmErrors(kNBins, xIsHit, yIsHit,
                                          exlIsHit, exhIsHit,
                                          eylIsHit, eyhIsHit);

  // Print summary to stdout.
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Track-based efficiency vs angle" << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  double totalNum = 0.0;
  double totalDen = 0.0;

  for (int iBin = 1; iBin <= kNBins; ++iBin) {
    const double low = hEff->GetXaxis()->GetBinLowEdge(iBin);
    const double high = hEff->GetXaxis()->GetBinUpEdge(iBin);
    const double num = hNum->GetBinContent(iBin);
    const double den = hDen->GetBinContent(iBin);
    const double eff = hEff->GetBinContent(iBin);
    const double errLow = eyl[iBin - 1];
    const double errHigh = eyh[iBin - 1];

    totalNum += num;
    totalDen += den;

    std::cout << "[" << std::setw(2) << low << ", " << std::setw(2) << high << ") deg : "
              << "numerator = " << std::setw(8) << num
              << ", denominator = " << std::setw(8) << den
              << ", efficiency = " << 100.0 * eff
              << " +" << 100.0 * errHigh
              << " -" << 100.0 * errLow << " %" << std::endl;
  }

  std::cout << "----------------------------------------" << std::endl;

  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Track-based is_hit efficiency vs angle" << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  double totalNumIsHit = 0.0;
  double totalDenIsHit = 0.0;

  for (int iBin = 1; iBin <= kNBins; ++iBin) {
    const double low = hEffIsHit->GetXaxis()->GetBinLowEdge(iBin);
    const double high = hEffIsHit->GetXaxis()->GetBinUpEdge(iBin);
    const double num = hNumIsHit->GetBinContent(iBin);
    const double den = hDenIsHit->GetBinContent(iBin);
    const double eff = hEffIsHit->GetBinContent(iBin);
    const double errLow = eylIsHit[iBin - 1];
    const double errHigh = eyhIsHit[iBin - 1];

    totalNumIsHit += num;
    totalDenIsHit += den;

    std::cout << "[" << std::setw(2) << low << ", " << std::setw(2) << high << ") deg : "
              << "numerator = " << std::setw(8) << num
              << ", denominator = " << std::setw(8) << den
              << ", efficiency = " << 100.0 * eff
              << " +" << 100.0 * errHigh
              << " -" << 100.0 * errLow << " %" << std::endl;
  }

  std::cout << "----------------------------------------" << std::endl;
  if (totalDenIsHit > 0.0) {
    const double totalEffIsHit = totalNumIsHit / totalDenIsHit;
    double totalLowerIsHit = 0.0;
    double totalUpperIsHit = 1.0;

    if (totalNumIsHit > 0.0) {
      totalLowerIsHit = ROOT::Math::beta_quantile(alpha / 2.0,
                                                  totalNumIsHit,
                                                  totalDenIsHit - totalNumIsHit + 1.0);
    }
    if (totalNumIsHit < totalDenIsHit) {
      totalUpperIsHit = ROOT::Math::beta_quantile(1.0 - alpha / 2.0,
                                                  totalNumIsHit + 1.0,
                                                  totalDenIsHit - totalNumIsHit);
    }

    std::cout << "Total is_hit efficiency = "
              << 100.0 * totalEffIsHit
              << " +" << 100.0 * (totalUpperIsHit - totalEffIsHit)
              << " -" << 100.0 * (totalEffIsHit - totalLowerIsHit)
              << " %" << std::endl;
  } else {
    std::cout << "Total is_hit efficiency = undefined (denominator = 0)" << std::endl;
  }
  std::cout << "----------------------------------------" << std::endl;

  if (totalDen > 0.0) {
    const double totalEff = totalNum / totalDen;
    double totalLower = 0.0;
    double totalUpper = 1.0;

    if (totalNum > 0.0) {
      totalLower = ROOT::Math::beta_quantile(alpha / 2.0,
                                             totalNum,
                                             totalDen - totalNum + 1.0);
    }
    if (totalNum < totalDen) {
      totalUpper = ROOT::Math::beta_quantile(1.0 - alpha / 2.0,
                                             totalNum + 1.0,
                                             totalDen - totalNum);
    }

    std::cout << "Total efficiency = "
              << 100.0 * totalEff
              << " +" << 100.0 * (totalUpper - totalEff)
              << " -" << 100.0 * (totalEff - totalLower)
              << " %" << std::endl;
  } else {
    std::cout << "Total efficiency = undefined (denominator = 0)" << std::endl;
  }
  std::cout << "----------------------------------------" << std::endl;

  // Draw the efficiency histogram and overlay asymmetric error bars.

  auto *canvas = new TCanvas("c_eff", "c_eff", 900, 700);
  canvas->SetGrid();

  hEff->SetMinimum(0.0);
  hEff->SetMaximum(1.0);
  hEff->SetLineColor(0);
  hEff->SetLineWidth(0);
  hEff->SetFillStyle(0);
  hEff->SetFillColor(0);
  hEff->SetMarkerSize(0);

  gEff->SetMarkerStyle(20);
  gEff->SetMarkerSize(1.2);
  gEff->SetLineWidth(2);
  gEffIsHit->SetMarkerStyle(20);
  gEffIsHit->SetMarkerSize(1.2);
  gEffIsHit->SetLineWidth(2);

  // Write a multi-page PDF.
  canvas->SaveAs((std::string(outputPdfPath) + "[").c_str());

  // Page 1: efficiency.
  gStyle->SetOptStat(0);
  hEff->Draw();
  gEff->Draw("P SAME");
  canvas->SaveAs(outputPdfPath);

  // Page 2: is_hit efficiency.
  canvas->Clear();
  canvas->SetGrid();
  gStyle->SetOptStat(0);
  hEffIsHit->SetMinimum(0.0);
  hEffIsHit->SetMaximum(1.0);
  hEffIsHit->SetLineColor(0);
  hEffIsHit->SetLineWidth(0);
  hEffIsHit->SetFillStyle(0);
  hEffIsHit->SetFillColor(0);
  hEffIsHit->SetMarkerSize(0);
  hEffIsHit->Draw();
  gEffIsHit->Draw("P SAME");
  canvas->SaveAs(outputPdfPath);

  // Page 3: dx for all denominator tracks.
  canvas->Clear();
  canvas->SetGrid();
  gStyle->SetOptStat(1110);
  hDxAll->SetTitle("dx distribution for all angles");
  hDxAll->SetLineWidth(2);
  hDxAll->Draw("HIST");
  canvas->SaveAs(outputPdfPath);

  // Page 4: dy for all denominator tracks.
  canvas->Clear();
  canvas->SetGrid();
  gStyle->SetOptStat(1110);
  hDyAll->SetTitle("dy distribution for all angles");
  hDyAll->SetLineWidth(2);
  hDyAll->Draw("HIST");
  canvas->SaveAs(outputPdfPath);

  // Pages 5+: dx by atan(|tan_x|) bins.
  for (int i = 0; i < kNBins; ++i) {
    canvas->Clear();
    canvas->SetGrid();
    gStyle->SetOptStat(1110);
    hDxByAngleX[i]->SetTitle(
        Form("dx distribution: %.0f #leq #theta_{x} < %.0f deg",
             kAngleBins[i], kAngleBins[i + 1]));
    hDxByAngleX[i]->SetLineWidth(2);
    hDxByAngleX[i]->Draw("HIST");
    canvas->SaveAs(outputPdfPath);
  }

  // Pages after that: dx by atan(sqrt(tan_x^2 + tan_y^2)) bins.
  for (int i = 0; i < kNBins; ++i) {
    canvas->Clear();
    canvas->SetGrid();
    gStyle->SetOptStat(1110);
    hDxByAngleTot[i]->SetTitle(
        Form("dx distribution: %.0f #leq #theta < %.0f deg",
             kAngleBins[i], kAngleBins[i + 1]));
    hDxByAngleTot[i]->SetLineWidth(2);
    hDxByAngleTot[i]->Draw("HIST");
    canvas->SaveAs(outputPdfPath);
  }

  // Pages after that: dy by atan(|tan_y|) bins.
  for (int i = 0; i < kNBins; ++i) {
    canvas->Clear();
    canvas->SetGrid();
    gStyle->SetOptStat(1110);
    hDyByAngleY[i]->SetTitle(
        Form("dy distribution: %.0f #leq #theta_{y} < %.0f deg",
             kAngleBins[i], kAngleBins[i + 1]));
    hDyByAngleY[i]->SetLineWidth(2);
    hDyByAngleY[i]->Draw("HIST");
    canvas->SaveAs(outputPdfPath);
  }

  // Final pages: dy by atan(sqrt(tan_x^2 + tan_y^2)) bins.
  for (int i = 0; i < kNBins; ++i) {
    canvas->Clear();
    canvas->SetGrid();
    gStyle->SetOptStat(1110);
    hDyByAngleTot[i]->SetTitle(
        Form("dy distribution: %.0f #leq #theta < %.0f deg",
             kAngleBins[i], kAngleBins[i + 1]));
    hDyByAngleTot[i]->SetLineWidth(2);
    hDyByAngleTot[i]->Draw("HIST");
    canvas->SaveAs(outputPdfPath);
  }

  canvas->SaveAs((std::string(outputPdfPath) + "]").c_str());

  std::cout << "Saved PDF: " << outputPdfPath << std::endl;
  std::cout << "Saved log: " << logFilePath << std::endl;

  std::cout.rdbuf(oldCoutBuf);
  std::cerr.rdbuf(oldCerrBuf);

  delete canvas;
  delete gEff;
  delete gEffIsHit;
  delete hEff;
  delete hEffIsHit;
  delete hNum;
  delete hDen;
  delete hNumIsHit;
  delete hDenIsHit;
  delete hDxAll;
  delete hDyAll;
  for (int i = 0; i < kNBins; ++i) {
    delete hDxByAngleX[i];
    delete hDxByAngleTot[i];
    delete hDyByAngleY[i];
    delete hDyByAngleTot[i];
  }
}
