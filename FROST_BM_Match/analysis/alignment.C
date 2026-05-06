#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TF1.h>
#include <TMath.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <streambuf>
#include <iomanip>
#include <cmath>
#include <algorithm>

// ------------------------------------------------------------
// Draw dx/dy distributions in 24 FROST detector areas.
//
// Usage in ROOT:
//   .L alignment.C
//   alignment("/path/to/input_dir",
//                     "/path/to/output.pdf",
//                     "/path/to/output.log");
//
//   std::vector<std::string> exclude = {"hoge.root", "hage.root"};
//   alignment("/path/to/input_dir",
//                     "/path/to/output.pdf",
//                     "/path/to/output.log",
//                     exclude);
//
// FROST area definition:
//   -660 <= x <= 660 mm, -700 <= y <= 700 mm
//   x bin edges: -660, -440, -220, 0, 220, 440, 660 mm
//   y bin edges: -700, -350, 0, 350, 700 mm
//
// Area assignment:
//   trackmatch_frost_nearest_x and trackmatch_frost_nearest_y
//
// Fill condition:
//   matched == 1
//   && trackmatch_has_match == 1
//   && bsd_good_spill_flag != 0
//   && atan(sqrt(trackmatch_baby_mind_tangent_x^2
//                + trackmatch_baby_mind_tangent_y^2)) < 10 deg
//
// Output pages in the PDF:
//   1. dx distributions in 24 FROST areas
//   2. dy distributions in 24 FROST areas
//   3. dx and dy distributions for all FROST areas
//
// Each page is divided into 6 x 4 pads.
// The pad layout follows the real FROST geometry:
//   left to right: increasing x
//   top to bottom: decreasing y
// ------------------------------------------------------------

namespace {
  constexpr int kNX = 6;
  constexpr int kNY = 4;
  constexpr int kNAreas = kNX * kNY;

  const double kXEdges[kNX + 1] = {
    -660.0, -440.0, -220.0, 0.0, 220.0, 440.0, 660.0
  };

  const double kYEdges[kNY + 1] = {
    -700.0, -350.0, 0.0, 350.0, 700.0
  };

  bool HasRootExtension(const std::string &name) {
    return name.size() >= 5 && name.substr(name.size() - 5) == ".root";
  }

  bool IsExcludedFile(const std::string &fileName,
                      const std::vector<std::string> &excludedFiles) {
    return std::find(excludedFiles.begin(), excludedFiles.end(), fileName)
           != excludedFiles.end();
  }

  class TeeBuf : public std::streambuf {
  public:
    TeeBuf(std::streambuf *sb1, std::streambuf *sb2) : sb1_(sb1), sb2_(sb2) {}

  protected:
    int overflow(int c) override {
      if (c == EOF) return !EOF;
      const int r1 = sb1_ ? sb1_->sputc(c) : c;
      const int r2 = sb2_ ? sb2_->sputc(c) : c;
      return (r1 == EOF || r2 == EOF) ? EOF : c;
    }

    int sync() override {
      const int r1 = sb1_ ? sb1_->pubsync() : 0;
      const int r2 = sb2_ ? sb2_->pubsync() : 0;
      return (r1 == 0 && r2 == 0) ? 0 : -1;
    }

  private:
    std::streambuf *sb1_;
    std::streambuf *sb2_;
  };

  int FindXBin(double x) {
    for (int ix = 0; ix < kNX; ++ix) {
      if (x >= kXEdges[ix] && x < kXEdges[ix + 1]) return ix;
    }

    // Include the upper edge in the last bin.
    if (x == kXEdges[kNX]) return kNX - 1;
    return -1;
  }

  int FindYBin(double y) {
    for (int iy = 0; iy < kNY; ++iy) {
      if (y >= kYEdges[iy] && y < kYEdges[iy + 1]) return iy;
    }

    // Include the upper edge in the last bin.
    if (y == kYEdges[kNY]) return kNY - 1;
    return -1;
  }

  int AreaIndex(int ix, int iy) {
    return iy * kNX + ix;
  }

  TString AreaTitle(int ix, int iy) {
    return Form("%.0f mm < x < %.0f mm, %.0f mm < y < %.0f mm",
                kXEdges[ix], kXEdges[ix + 1],
                kYEdges[iy], kYEdges[iy + 1]);
  }

  int PadIndexForFrostLayout(int ix, int iy) {
    // ROOT pads are numbered left-to-right, top-to-bottom.
    // y bin index increases from bottom to top, so invert iy for display.
    const int display_row = (kNY - 1) - iy;
    return display_row * kNX + ix + 1;
  }

  void DrawHistogramWithGaussianFit(TH1D *hist) {
    hist->SetLineWidth(2);
    hist->Draw("HIST");

    // Require enough entries and a finite RMS for a meaningful Gaussian fit.
    if (hist->GetEntries() < 5) return;
    if (hist->GetRMS() <= 0.0) return;

    const double x_min = hist->GetXaxis()->GetXmin();
    const double x_max = hist->GetXaxis()->GetXmax();

    auto *fit = new TF1(Form("%s_gaus_fit", hist->GetName()),
                        "gaus",
                        x_min,
                        x_max);
    fit->SetLineWidth(2);

    // R: use the TF1 range.
    // Q: quiet fit.
    // 0: do not let ROOT auto-draw the fit. We draw it explicitly below.
    hist->Fit(fit, "RQ0");
    fit->Draw("SAME");

    gPad->Modified();
    gPad->Update();
  }
}

void alignment(
    const char *inputDir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/latest20260504/rootfile_after_TrackMatch/zshift/zshift_0",
    const char *outputPdfPath = "/group/nu/ninja/work/otani/FROSTReconData//BM_FROST/analysis_plot/alignment/alignment_zshift0.pdf",
    const char *logFilePath = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/alignment/alignment_zshift0.log",
    const std::vector<std::string> &excludedFiles = std::vector<std::string>
    {"BMPM_track_2025-11-29_13-46-59_Run0_afterTrackMatch_zshift0.root",
     "BMPM_track_2025-11-30_13-11-36_Run0_afterTrackMatch_zshift0.root"}) {

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

  Int_t matched = 0;
  Int_t bsd_good_spill_flag = 0;

  std::vector<int> *trackmatch_has_match = nullptr;
  std::vector<double> *trackmatch_frost_nearest_x = nullptr;
  std::vector<double> *trackmatch_frost_nearest_y = nullptr;
  std::vector<double> *trackmatch_dx = nullptr;
  std::vector<double> *trackmatch_dy = nullptr;
  std::vector<double> *trackmatch_baby_mind_tangent_x = nullptr;
  std::vector<double> *trackmatch_baby_mind_tangent_y = nullptr;

  chain.SetBranchAddress("matched", &matched);
  chain.SetBranchAddress("bsd_good_spill_flag", &bsd_good_spill_flag);
  chain.SetBranchAddress("trackmatch_has_match", &trackmatch_has_match);
  chain.SetBranchAddress("trackmatch_frost_nearest_x", &trackmatch_frost_nearest_x);
  chain.SetBranchAddress("trackmatch_frost_nearest_y", &trackmatch_frost_nearest_y);
  chain.SetBranchAddress("trackmatch_dx", &trackmatch_dx);
  chain.SetBranchAddress("trackmatch_dy", &trackmatch_dy);
  chain.SetBranchAddress("trackmatch_baby_mind_tangent_x", &trackmatch_baby_mind_tangent_x);
  chain.SetBranchAddress("trackmatch_baby_mind_tangent_y", &trackmatch_baby_mind_tangent_y);

  std::vector<TH1D*> hDx(kNAreas, nullptr);
  std::vector<TH1D*> hDy(kNAreas, nullptr);
  auto *hDxAll = new TH1D("hDxAll",
                          "dx distribution for all FROST areas;dx [mm];Number of events",
                          100, -500.0, 500.0);
  auto *hDyAll = new TH1D("hDyAll",
                          "dy distribution for all FROST areas;dy [mm];Number of events",
                          100, -500.0, 500.0);
  hDxAll->Sumw2();
  hDyAll->Sumw2();

  for (int iy = 0; iy < kNY; ++iy) {
    for (int ix = 0; ix < kNX; ++ix) {
      const int iarea = AreaIndex(ix, iy);
      const TString title = AreaTitle(ix, iy);

      hDx[iarea] = new TH1D(Form("hDx_area_%d_%d", ix, iy),
                            Form("%s;dx [mm];Number of events", title.Data()),
                            100, -500.0, 500.0);
      hDy[iarea] = new TH1D(Form("hDy_area_%d_%d", ix, iy),
                            Form("%s;dy [mm];Number of events", title.Data()),
                            100, -500.0, 500.0);

      hDx[iarea]->Sumw2();
      hDy[iarea]->Sumw2();
    }
  }

  Long64_t nFilledTracks = 0;
  Long64_t nOutOfFrostArea = 0;
  Long64_t nVectorMismatch = 0;

  const Long64_t nEntries = chain.GetEntries();
  for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
    chain.GetEntry(iEntry);

    if (matched != 1) continue;
    if (bsd_good_spill_flag == 0) continue;

    if (!trackmatch_has_match ||
        !trackmatch_frost_nearest_x ||
        !trackmatch_frost_nearest_y ||
        !trackmatch_dx ||
        !trackmatch_dy ||
        !trackmatch_baby_mind_tangent_x ||
        !trackmatch_baby_mind_tangent_y) {
      continue;
    }

    const std::size_t nTracks = trackmatch_has_match->size();
    if (trackmatch_frost_nearest_x->size() != nTracks ||
        trackmatch_frost_nearest_y->size() != nTracks ||
        trackmatch_dx->size() != nTracks ||
        trackmatch_dy->size() != nTracks ||
        trackmatch_baby_mind_tangent_x->size() != nTracks ||
        trackmatch_baby_mind_tangent_y->size() != nTracks) {
      ++nVectorMismatch;
      std::cerr << "Warning: vector size mismatch at spill entry "
                << iEntry << ". Skip this spill." << std::endl;
      continue;
    }

    for (std::size_t iTrack = 0; iTrack < nTracks; ++iTrack) {
      if (trackmatch_has_match->at(iTrack) != 1) continue;

      const double frostX = trackmatch_frost_nearest_x->at(iTrack);
      const double frostY = trackmatch_frost_nearest_y->at(iTrack);
      const double dx = trackmatch_dx->at(iTrack);
      const double dy = trackmatch_dy->at(iTrack);
      const double tx = trackmatch_baby_mind_tangent_x->at(iTrack);
      const double ty = trackmatch_baby_mind_tangent_y->at(iTrack);

      const double angleDeg =
        std::atan(std::sqrt(tx * tx + ty * ty)) * 180.0 / TMath::Pi();
      if (angleDeg >= 10.0) continue;

      const int ix = FindXBin(frostX);
      const int iy = FindYBin(frostY);

      if (ix < 0 || iy < 0) {
        ++nOutOfFrostArea;
        continue;
      }

      const int iarea = AreaIndex(ix, iy);
      hDx[iarea]->Fill(dx);
      hDy[iarea]->Fill(dy);
      hDxAll->Fill(dx);
      hDyAll->Fill(dy);
      ++nFilledTracks;
    }
  }

  std::cout << "Filled matched tracks: " << nFilledTracks << std::endl;
  std::cout << "Tracks outside FROST area: " << nOutOfFrostArea << std::endl;
  std::cout << "Spills skipped by vector mismatch: " << nVectorMismatch << std::endl;

  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1111);

  auto *canvas = new TCanvas("c_frost_area_dxdy", "c_frost_area_dxdy", 1800, 1200);

  canvas->SaveAs((std::string(outputPdfPath) + "[").c_str());

  // Page 1: dx distributions.
  canvas->Clear();
  canvas->Divide(kNX, kNY, 0.001, 0.001);
  for (int iy = 0; iy < kNY; ++iy) {
    for (int ix = 0; ix < kNX; ++ix) {
      const int iarea = AreaIndex(ix, iy);
      const int ipad = PadIndexForFrostLayout(ix, iy);

      canvas->cd(ipad);
      gPad->SetGrid();
      DrawHistogramWithGaussianFit(hDx[iarea]);
    }
  }
  canvas->SaveAs(outputPdfPath);

  // Page 2: dy distributions.
  canvas->Clear();
  canvas->Divide(kNX, kNY, 0.001, 0.001);
  for (int iy = 0; iy < kNY; ++iy) {
    for (int ix = 0; ix < kNX; ++ix) {
      const int iarea = AreaIndex(ix, iy);
      const int ipad = PadIndexForFrostLayout(ix, iy);

      canvas->cd(ipad);
      gPad->SetGrid();
      DrawHistogramWithGaussianFit(hDy[iarea]);
    }
  }
  canvas->SaveAs(outputPdfPath);

  // Page 3: dx and dy distributions for all FROST areas.
  canvas->Clear();
  canvas->Divide(2, 1);

  canvas->cd(1);
  gPad->SetGrid();
  DrawHistogramWithGaussianFit(hDxAll);

  canvas->cd(2);
  gPad->SetGrid();
  DrawHistogramWithGaussianFit(hDyAll);

  canvas->SaveAs(outputPdfPath);

  canvas->SaveAs((std::string(outputPdfPath) + "]").c_str());

  std::cout << "Saved PDF: " << outputPdfPath << std::endl;
  std::cout << "Saved log: " << logFilePath << std::endl;

  std::cout.rdbuf(oldCoutBuf);
  std::cerr.rdbuf(oldCerrBuf);

  delete canvas;
  delete hDxAll;
  delete hDyAll;
  for (int iarea = 0; iarea < kNAreas; ++iarea) {
    delete hDx[iarea];
    delete hDy[iarea];
  }
}
