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
#include <TLegend.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TText.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <streambuf>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <array>
#include <memory>
#include <Math/QuantFuncMathCore.h>

// ------------------------------------------------------------
// Draw BM-FROST matching efficiencies and dx/dy distributions.
//
// This version reads MC and data at the same time and overlays them.
// Efficiency pages show MC and data with different colors.
// Distribution pages show MC as a histogram and data as points with
// error bars.  By default the MC histogram is scaled to the data integral
// on each distribution page.
//
// Usage in ROOT:
//   .L DrawTrackMatchEfficiency.C
//   DrawTrackMatchEfficiency("/path/to/mc_dir",
//                            "/path/to/data_dir",
//                            "/path/to/output.pdf",
//                            "/path/to/output.log");
//
//   std::vector<std::string> exclude_mc = {"hoge.root"};
//   std::vector<std::string> exclude_data = {"hage.root"};
//   DrawTrackMatchEfficiency("/path/to/mc_dir",
//                            "/path/to/data_dir",
//                            "/path/to/output.pdf",
//                            "/path/to/output.log",
//                            exclude_mc,
//                            exclude_data);
//
// Excluded files are specified by file name relative to each input dir.
// ------------------------------------------------------------

namespace {
  constexpr int kNBins = 8;
  const double kAngleBins[kNBins + 1] = {
    0.0, 5.0, 10.0, 15.0, 20.0,
    25.0, 30.0, 35.0, 40.0
  };

  constexpr double kEffYMin = 0.90;
  constexpr double kEffYMax = 1.00;
  constexpr double kEffYTitleSize = 0.05;
  constexpr double kTruthHitEffYTitleSize = 0.04;
  constexpr double kCorrectedEffYMax = 1.05;

  // Histogram binning for residual distributions.
  // kResidualHistBinsAll is used for all-angle pages.
  // kResidualHistBinsByAngle[i] is used for all residual quantities
  // (dx, dy, dtanx and dtany) in the corresponding angle bin:
  //   [kAngleBins[i], kAngleBins[i + 1]) deg.
  constexpr int kResidualHistBinsAll = 200;
  const int kResidualHistBinsByAngle[kNBins] = {
    200, 200, 200, 200,
    100, 100, 50, 50
  };

  int ResidualHistBinsForAngleBin(int angle_bin) {
    if (angle_bin < 0 || angle_bin >= kNBins) {
      return kResidualHistBinsAll;
    }
    return kResidualHistBinsByAngle[angle_bin];
  }

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

  struct SampleHists {
    std::string label;
    bool is_mc = false;
    int n_files_added = 0;
    int n_files_excluded = 0;
    Long64_t n_spills = 0;

    TH1D *hDen = nullptr;
    TH1D *hNum = nullptr;
    TH1D *hDenIsHit = nullptr;
    TH1D *hNumIsHit = nullptr;
    TH1D *hDenTruthMuon = nullptr;
    TH1D *hNumTruthMuon = nullptr;
    TH1D *hDenIsHitGivenTruth = nullptr;
    TH1D *hNumIsHitGivenTruth = nullptr;
    TH1D *hEff = nullptr;
    TH1D *hEffIsHit = nullptr;
    TH1D *hEffTruthMuon = nullptr;
    TH1D *hEffIsHitGivenTruth = nullptr;
    TGraphAsymmErrors *gEff = nullptr;
    TGraphAsymmErrors *gEffIsHit = nullptr;
    TGraphAsymmErrors *gEffTruthMuon = nullptr;
    TGraphAsymmErrors *gEffIsHitGivenTruth = nullptr;

    TH1D *hDxAll = nullptr;
    TH1D *hDyAll = nullptr;
    TH1D *hDtanxAll = nullptr;
    TH1D *hDtanyAll = nullptr;

    std::vector<TH1D*> hDxByAngleX;
    std::vector<TH1D*> hDxByAngleTot;
    std::vector<TH1D*> hDyByAngleY;
    std::vector<TH1D*> hDyByAngleTot;
    std::vector<TH1D*> hDtanxByAngleX;
    std::vector<TH1D*> hDtanxByAngleTot;
    std::vector<TH1D*> hDtanyByAngleY;
    std::vector<TH1D*> hDtanyByAngleTot;
  };

  TH1D *MakeHist(const TString &name,
                 const TString &title,
                 int nbins,
                 double xmin,
                 double xmax) {
    auto *hist = new TH1D(name, title, nbins, xmin, xmax);
    hist->Sumw2();
    hist->SetDirectory(nullptr);
    return hist;
  }

  void MakeSampleHists(SampleHists &sample,
                       const std::string &prefix,
                       const std::string &label,
                       bool is_mc) {
    sample.label = label;
    sample.is_mc = is_mc;

    sample.hDen = new TH1D(Form("%s_hDen", prefix.c_str()),
                           "Denominator;Angle [deg];Tracks",
                           kNBins, kAngleBins);
    sample.hNum = new TH1D(Form("%s_hNum", prefix.c_str()),
                           "Numerator;Angle [deg];Tracks",
                           kNBins, kAngleBins);
    sample.hDenIsHit = new TH1D(Form("%s_hDenIsHit", prefix.c_str()),
                                "Denominator for is_hit;Angle [deg];Tracks",
                                kNBins, kAngleBins);
    sample.hNumIsHit = new TH1D(Form("%s_hNumIsHit", prefix.c_str()),
                                "Numerator for is_hit;Angle [deg];Tracks",
                                kNBins, kAngleBins);
    sample.hDenTruthMuon = new TH1D(Form("%s_hDenTruthMuon", prefix.c_str()),
                                    "Denominator for MC truth FROST muon;Angle [deg];Tracks",
                                    kNBins, kAngleBins);
    sample.hNumTruthMuon = new TH1D(Form("%s_hNumTruthMuon", prefix.c_str()),
                                    "Numerator for MC truth FROST muon;Angle [deg];Tracks",
                                    kNBins, kAngleBins);
    sample.hDenIsHitGivenTruth = new TH1D(Form("%s_hDenIsHitGivenTruth", prefix.c_str()),
                                          "Denominator for is_hit given MC truth FROST muon;Angle [deg];Tracks",
                                          kNBins, kAngleBins);
    sample.hNumIsHitGivenTruth = new TH1D(Form("%s_hNumIsHitGivenTruth", prefix.c_str()),
                                          "Numerator for is_hit given MC truth FROST muon;Angle [deg];Tracks",
                                          kNBins, kAngleBins);
    sample.hEff = new TH1D(Form("%s_hEff", prefix.c_str()),
                           ";Baby MIND reconstructed angle [deg];Track matching efficiency",
                           kNBins, kAngleBins);
    sample.hEffIsHit = new TH1D(Form("%s_hEffIsHit", prefix.c_str()),
                                ";Baby MIND reconstructed angle [deg];FROST hit efficiency",
                                kNBins, kAngleBins);
    sample.hEffTruthMuon = new TH1D(Form("%s_hEffTruthMuon", prefix.c_str()),
                                    ";Baby MIND reconstructed angle [deg];MC truth FROST-muon fraction",
                                    kNBins, kAngleBins);
    sample.hEffIsHitGivenTruth = new TH1D(Form("%s_hEffIsHitGivenTruth", prefix.c_str()),
                                          ";Baby MIND reconstructed angle [deg];FROST hit efficiency for truth FROST muons",
                                          kNBins, kAngleBins);

    sample.hDen->Sumw2();
    sample.hNum->Sumw2();
    sample.hDenIsHit->Sumw2();
    sample.hNumIsHit->Sumw2();
    sample.hDenTruthMuon->Sumw2();
    sample.hNumTruthMuon->Sumw2();
    sample.hDenIsHitGivenTruth->Sumw2();
    sample.hNumIsHitGivenTruth->Sumw2();
    sample.hEff->Sumw2();
    sample.hEffIsHit->Sumw2();
    sample.hEffTruthMuon->Sumw2();
    sample.hEffIsHitGivenTruth->Sumw2();
    sample.hDen->SetDirectory(nullptr);
    sample.hNum->SetDirectory(nullptr);
    sample.hDenIsHit->SetDirectory(nullptr);
    sample.hNumIsHit->SetDirectory(nullptr);
    sample.hDenTruthMuon->SetDirectory(nullptr);
    sample.hNumTruthMuon->SetDirectory(nullptr);
    sample.hDenIsHitGivenTruth->SetDirectory(nullptr);
    sample.hNumIsHitGivenTruth->SetDirectory(nullptr);
    sample.hEff->SetDirectory(nullptr);
    sample.hEffIsHit->SetDirectory(nullptr);
    sample.hEffTruthMuon->SetDirectory(nullptr);
    sample.hEffIsHitGivenTruth->SetDirectory(nullptr);

    sample.hDxAll = MakeHist(Form("%s_hDxAll", prefix.c_str()),
                             ";dx [mm];Number of tracks", kResidualHistBinsAll, -500.0, 500.0);
    sample.hDyAll = MakeHist(Form("%s_hDyAll", prefix.c_str()),
                             ";dy [mm];Number of tracks", kResidualHistBinsAll, -500.0, 500.0);
    sample.hDtanxAll = MakeHist(Form("%s_hDtanxAll", prefix.c_str()),
                                ";dtanx;Number of tracks", kResidualHistBinsAll, -0.25, 0.25);
    sample.hDtanyAll = MakeHist(Form("%s_hDtanyAll", prefix.c_str()),
                                ";dtany;Number of tracks", kResidualHistBinsAll, -0.25, 0.25);

    for (int i = 0; i < kNBins; ++i) {
      const int nbins = ResidualHistBinsForAngleBin(i);
      sample.hDxByAngleX.push_back(MakeHist(Form("%s_hDx_bin%d", prefix.c_str(), i),
                                            ";dx [mm];Number of tracks", nbins, -500.0, 500.0));
      sample.hDxByAngleTot.push_back(MakeHist(Form("%s_hDxTot_bin%d", prefix.c_str(), i),
                                              ";dx [mm];Number of tracks", nbins, -500.0, 500.0));
      sample.hDyByAngleY.push_back(MakeHist(Form("%s_hDy_bin%d", prefix.c_str(), i),
                                            ";dy [mm];Number of tracks", nbins, -500.0, 500.0));
      sample.hDyByAngleTot.push_back(MakeHist(Form("%s_hDyTot_bin%d", prefix.c_str(), i),
                                              ";dy [mm];Number of tracks", nbins, -500.0, 500.0));
      sample.hDtanxByAngleX.push_back(MakeHist(Form("%s_hDtanx_bin%d", prefix.c_str(), i),
                                               ";dtanx;Number of tracks", nbins, -0.25, 0.25));
      sample.hDtanxByAngleTot.push_back(MakeHist(Form("%s_hDtanxTot_bin%d", prefix.c_str(), i),
                                                 ";dtanx;Number of tracks", nbins, -0.25, 0.25));
      sample.hDtanyByAngleY.push_back(MakeHist(Form("%s_hDtany_bin%d", prefix.c_str(), i),
                                               ";dtany;Number of tracks", nbins, -0.25, 0.25));
      sample.hDtanyByAngleTot.push_back(MakeHist(Form("%s_hDtanyTot_bin%d", prefix.c_str(), i),
                                                 ";dtany;Number of tracks", nbins, -0.25, 0.25));
    }
  }

  void DeleteHistVector(std::vector<TH1D*> &histograms) {
    for (auto *hist : histograms) {
      delete hist;
    }
    histograms.clear();
  }

  void DeleteSampleHists(SampleHists &sample) {
    delete sample.gEff;
    delete sample.gEffIsHit;
    delete sample.gEffTruthMuon;
    delete sample.gEffIsHitGivenTruth;
    delete sample.hEff;
    delete sample.hEffIsHit;
    delete sample.hEffTruthMuon;
    delete sample.hEffIsHitGivenTruth;
    delete sample.hNum;
    delete sample.hDen;
    delete sample.hNumIsHit;
    delete sample.hDenIsHit;
    delete sample.hNumTruthMuon;
    delete sample.hDenTruthMuon;
    delete sample.hNumIsHitGivenTruth;
    delete sample.hDenIsHitGivenTruth;
    delete sample.hDxAll;
    delete sample.hDyAll;
    delete sample.hDtanxAll;
    delete sample.hDtanyAll;
    DeleteHistVector(sample.hDxByAngleX);
    DeleteHistVector(sample.hDxByAngleTot);
    DeleteHistVector(sample.hDyByAngleY);
    DeleteHistVector(sample.hDyByAngleTot);
    DeleteHistVector(sample.hDtanxByAngleX);
    DeleteHistVector(sample.hDtanxByAngleTot);
    DeleteHistVector(sample.hDtanyByAngleY);
    DeleteHistVector(sample.hDtanyByAngleTot);
  }

  bool BuildChain(const char *inputDir,
                  const std::vector<std::string> &excludedFiles,
                  const std::string &label,
                  TChain &chain,
                  int &nFilesAdded,
                  int &nFilesExcluded) {
    TSystemDirectory dir(Form("input_dir_%s", label.c_str()), inputDir);
    std::unique_ptr<TList> fileList(dir.GetListOfFiles());
    if (!fileList) {
      std::cerr << "Error: cannot open " << label << " directory: "
                << inputDir << std::endl;
      return false;
    }

    nFilesAdded = 0;
    nFilesExcluded = 0;

    TIter next(fileList.get());
    while (TObject *obj = next()) {
      auto *sysFile = dynamic_cast<TSystemFile *>(obj);
      if (!sysFile || sysFile->IsDirectory()) {
        continue;
      }

      const std::string fileName = sysFile->GetName();
      if (!HasRootExtension(fileName)) {
        continue;
      }

      if (IsExcludedFile(fileName, excludedFiles)) {
        std::cout << label << " excluded file: " << fileName << std::endl;
        ++nFilesExcluded;
        continue;
      }

      const std::string fullPath = std::string(inputDir) + "/" + fileName;
      if (chain.Add(fullPath.c_str(), 0) > 0) {
        ++nFilesAdded;
        std::cout << label << " added file: " << fileName << std::endl;
      }
    }

    if (nFilesAdded == 0) {
      std::cerr << "Error: no ROOT files with tree 'match_info' were added from "
                << label << " directory: " << inputDir << std::endl;
      return false;
    }

    std::cout << label << " added " << nFilesAdded << " ROOT files." << std::endl;
    std::cout << label << " excluded " << nFilesExcluded << " ROOT files." << std::endl;
    std::cout << label << " total spills in chain: " << chain.GetEntries() << std::endl;
    return true;
  }

  struct EfficiencyResult {
    bool valid = false;
    double numerator = 0.0;
    double denominator = 0.0;
    double value = 0.0;
    double lower = 0.0;
    double upper = 0.0;
    double err_low = 0.0;
    double err_high = 0.0;
  };

  EfficiencyResult CalculateEfficiencyResult(double numerator,
                                             double denominator) {
    EfficiencyResult result;
    result.numerator = numerator;
    result.denominator = denominator;

    if (denominator <= 0.0 ||
        !std::isfinite(numerator) ||
        !std::isfinite(denominator)) {
      return result;
    }

    const double alpha = 1.0 - 0.682689492137;
    result.value = numerator / denominator;
    result.lower = 0.0;
    result.upper = 1.0;

    if (numerator > 0.0) {
      result.lower = ROOT::Math::beta_quantile(alpha / 2.0,
                                               numerator,
                                               denominator - numerator + 1.0);
    }
    if (numerator < denominator) {
      result.upper = ROOT::Math::beta_quantile(1.0 - alpha / 2.0,
                                               numerator + 1.0,
                                               denominator - numerator);
    }

    result.err_low = std::max(0.0, result.value - result.lower);
    result.err_high = std::max(0.0, result.upper - result.value);
    result.valid = std::isfinite(result.value) &&
                   std::isfinite(result.err_low) &&
                   std::isfinite(result.err_high);
    return result;
  }

  EfficiencyResult CalculateHistogramTotalEfficiency(const TH1D *numerator,
                                                     const TH1D *denominator) {
    double total_numerator = 0.0;
    double total_denominator = 0.0;
    if (numerator) {
      for (int iBin = 1; iBin <= numerator->GetNbinsX(); ++iBin) {
        total_numerator += numerator->GetBinContent(iBin);
      }
    }
    if (denominator) {
      for (int iBin = 1; iBin <= denominator->GetNbinsX(); ++iBin) {
        total_denominator += denominator->GetBinContent(iBin);
      }
    }
    return CalculateEfficiencyResult(total_numerator, total_denominator);
  }

  EfficiencyResult CalculateRatioResult(const EfficiencyResult &numerator,
                                        const EfficiencyResult &denominator) {
    EfficiencyResult result;
    result.numerator = numerator.value;
    result.denominator = denominator.value;

    if (!numerator.valid || !denominator.valid || denominator.value <= 0.0) {
      return result;
    }

    result.value = numerator.value / denominator.value;

    const double numerator_low = std::max(0.0, numerator.value - numerator.err_low);
    const double numerator_high = numerator.value + numerator.err_high;
    const double denominator_low = std::max(0.0, denominator.value - denominator.err_low);
    const double denominator_high = denominator.value + denominator.err_high;

    result.lower = denominator_high > 0.0 ? numerator_low / denominator_high : result.value;
    result.upper = denominator_low > 0.0 ? numerator_high / denominator_low
                                         : numerator_high / denominator.value;
    result.err_low = std::max(0.0, result.value - result.lower);
    result.err_high = std::max(0.0, result.upper - result.value);
    result.valid = std::isfinite(result.value) &&
                   std::isfinite(result.err_low) &&
                   std::isfinite(result.err_high);
    return result;
  }

  void PrintEfficiencyValue(const char *value_label,
                            const EfficiencyResult &result) {
    if (!result.valid) {
      std::cout << value_label << " = undefined (denominator = 0)";
      return;
    }

    std::cout << value_label << " = "
              << 100.0 * result.value
              << " +" << 100.0 * result.err_high
              << " -" << 100.0 * result.err_low
              << " %";
  }

  void PrintEfficiencyTable(const std::string &title,
                            const TH1D *numerator,
                            const TH1D *denominator,
                            const char *numerator_label,
                            const char *denominator_label,
                            const char *value_label) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << title << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    if (!numerator || !denominator) {
      std::cout << "No histogram available." << std::endl;
      std::cout << "----------------------------------------" << std::endl;
      return;
    }

    for (int iBin = 1; iBin <= kNBins; ++iBin) {
      const double low = denominator->GetXaxis()->GetBinLowEdge(iBin);
      const double high = denominator->GetXaxis()->GetBinUpEdge(iBin);
      const double num = numerator->GetBinContent(iBin);
      const double den = denominator->GetBinContent(iBin);
      const EfficiencyResult result = CalculateEfficiencyResult(num, den);

      std::cout << "[" << std::setw(2) << low << ", " << std::setw(2) << high << ") deg : "
                << numerator_label << " = " << std::setw(8) << num
                << ", " << denominator_label << " = " << std::setw(8) << den
                << ", ";
      PrintEfficiencyValue(value_label, result);
      std::cout << std::endl;
    }

    const EfficiencyResult total = CalculateHistogramTotalEfficiency(numerator, denominator);
    std::cout << "Total " << value_label << " : "
              << numerator_label << " = " << std::setw(8) << total.numerator
              << ", " << denominator_label << " = " << std::setw(8) << total.denominator
              << ", ";
    PrintEfficiencyValue(value_label, total);
    std::cout << std::endl;
    std::cout << "----------------------------------------" << std::endl;
  }

  TGraphAsymmErrors *BuildEfficiencyGraph(TH1D *numerator,
                                         TH1D *denominator,
                                         TH1D *efficiency_hist,
                                         const TString &graph_name,
                                         int color,
                                         int marker_style) {
    std::array<double, kNBins> x{};
    std::array<double, kNBins> y{};
    std::array<double, kNBins> exl{};
    std::array<double, kNBins> exh{};
    std::array<double, kNBins> eyl{};
    std::array<double, kNBins> eyh{};

    for (int iBin = 1; iBin <= kNBins; ++iBin) {
      const double num = numerator ? numerator->GetBinContent(iBin) : 0.0;
      const double den = denominator ? denominator->GetBinContent(iBin) : 0.0;
      const EfficiencyResult result = CalculateEfficiencyResult(num, den);
      if (efficiency_hist) {
        efficiency_hist->SetBinContent(iBin, result.value);
      }

      const double low = denominator->GetXaxis()->GetBinLowEdge(iBin);
      const double high = denominator->GetXaxis()->GetBinUpEdge(iBin);
      const double center = 0.5 * (low + high);

      x[iBin - 1] = center;
      y[iBin - 1] = result.value;
      exl[iBin - 1] = center - low;
      exh[iBin - 1] = high - center;
      eyl[iBin - 1] = result.valid ? result.err_low : 0.0;
      eyh[iBin - 1] = result.valid ? result.err_high : 0.0;
    }

    auto *graph = new TGraphAsymmErrors(kNBins,
                                        x.data(), y.data(),
                                        exl.data(), exh.data(),
                                        eyl.data(), eyh.data());
    graph->SetName(graph_name);
    graph->SetMarkerStyle(marker_style);
    graph->SetMarkerSize(1.2);
    graph->SetLineWidth(2);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    return graph;
  }

  TGraphAsymmErrors *BuildRatioGraph(const TGraphAsymmErrors *numeratorGraph,
                                      const TGraphAsymmErrors *denominatorGraph,
                                      const TString &graph_name,
                                      int color,
                                      int marker_style) {
    if (!numeratorGraph || !denominatorGraph) {
      return nullptr;
    }

    const int npoints = std::min(numeratorGraph->GetN(), denominatorGraph->GetN());
    if (npoints <= 0) {
      return nullptr;
    }

    std::vector<double> x(npoints, 0.0);
    std::vector<double> y(npoints, 0.0);
    std::vector<double> exl(npoints, 0.0);
    std::vector<double> exh(npoints, 0.0);
    std::vector<double> eyl(npoints, 0.0);
    std::vector<double> eyh(npoints, 0.0);

    for (int ipoint = 0; ipoint < npoints; ++ipoint) {
      double x_num = 0.0;
      double y_num = 0.0;
      double x_den = 0.0;
      double y_den = 0.0;
      numeratorGraph->GetPoint(ipoint, x_num, y_num);
      denominatorGraph->GetPoint(ipoint, x_den, y_den);

      x.at(ipoint) = x_num;
      exl.at(ipoint) = numeratorGraph->GetErrorXlow(ipoint);
      exh.at(ipoint) = numeratorGraph->GetErrorXhigh(ipoint);

      EfficiencyResult numerator;
      numerator.valid = std::isfinite(y_num);
      numerator.value = y_num;
      numerator.err_low = numeratorGraph->GetErrorYlow(ipoint);
      numerator.err_high = numeratorGraph->GetErrorYhigh(ipoint);

      EfficiencyResult denominator;
      denominator.valid = std::isfinite(y_den) && y_den > 0.0;
      denominator.value = y_den;
      denominator.err_low = denominatorGraph->GetErrorYlow(ipoint);
      denominator.err_high = denominatorGraph->GetErrorYhigh(ipoint);

      const EfficiencyResult ratio = CalculateRatioResult(numerator, denominator);
      y.at(ipoint) = ratio.value;
      eyl.at(ipoint) = ratio.valid ? ratio.err_low : 0.0;
      eyh.at(ipoint) = ratio.valid ? ratio.err_high : 0.0;
    }

    auto *graph = new TGraphAsymmErrors(npoints,
                                        x.data(), y.data(),
                                        exl.data(), exh.data(),
                                        eyl.data(), eyh.data());
    graph->SetName(graph_name);
    graph->SetMarkerStyle(marker_style);
    graph->SetMarkerSize(1.2);
    graph->SetLineWidth(2);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    return graph;
  }

  void FillEfficiencyGraphs(SampleHists &sample, int color, int marker_style) {
    sample.gEff = BuildEfficiencyGraph(sample.hNum,
                                       sample.hDen,
                                       sample.hEff,
                                       Form("g_%s_Eff", sample.label.c_str()),
                                       color,
                                       marker_style);
    sample.gEffIsHit = BuildEfficiencyGraph(sample.hNumIsHit,
                                            sample.hDenIsHit,
                                            sample.hEffIsHit,
                                            Form("g_%s_EffIsHit", sample.label.c_str()),
                                            color,
                                            marker_style);

    if (sample.is_mc) {
      sample.gEffTruthMuon = BuildEfficiencyGraph(sample.hNumTruthMuon,
                                                  sample.hDenTruthMuon,
                                                  sample.hEffTruthMuon,
                                                  Form("g_%s_EffTruthMuon", sample.label.c_str()),
                                                  kBlue + 1,
                                                  marker_style);
      sample.gEffIsHitGivenTruth = BuildEfficiencyGraph(sample.hNumIsHitGivenTruth,
                                                        sample.hDenIsHitGivenTruth,
                                                        sample.hEffIsHitGivenTruth,
                                                        Form("g_%s_EffIsHitGivenTruth", sample.label.c_str()),
                                                        kBlue + 2,
                                                        marker_style);
    }
  }

  void PrintEfficiencySummary(const SampleHists &sample) {
    std::cout << std::fixed << std::setprecision(4);
    PrintEfficiencyTable(sample.label + ": Track-based efficiency vs angle",
                         sample.hNum,
                         sample.hDen,
                         "numerator",
                         "denominator",
                         "efficiency");

    PrintEfficiencyTable(sample.label + ": Track-based is_hit efficiency vs angle",
                         sample.hNumIsHit,
                         sample.hDenIsHit,
                         "numerator",
                         "denominator",
                         "is_hit efficiency");

    if (sample.is_mc) {
      PrintEfficiencyTable(sample.label + ": MC truth FROST-muon fraction vs angle",
                           sample.hNumTruthMuon,
                           sample.hDenTruthMuon,
                           "truth muon",
                           "denominator",
                           "truth fraction");

      PrintEfficiencyTable(sample.label + ": FROST is_hit efficiency for MC truth FROST muons vs angle",
                           sample.hNumIsHitGivenTruth,
                           sample.hDenIsHitGivenTruth,
                           "is_hit && truth muon",
                           "truth muon",
                           "truth-muon is_hit efficiency");
    }
  }

  void PrintCorrectedDataHitEfficiencySummary(const TH1D *rawDataNumerator,
                                              const TH1D *rawDataDenominator,
                                              const TH1D *truthNumerator,
                                              const TH1D *truthDenominator) {
    if (!rawDataNumerator || !rawDataDenominator || !truthNumerator || !truthDenominator) {
      return;
    }

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Corrected Data FROST hit efficiency vs angle" << std::endl;
    std::cout << "Correction: Data apparent is_hit efficiency / MC truth FROST-muon fraction" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    for (int iBin = 1; iBin <= kNBins; ++iBin) {
      const double low = rawDataDenominator->GetXaxis()->GetBinLowEdge(iBin);
      const double high = rawDataDenominator->GetXaxis()->GetBinUpEdge(iBin);

      const EfficiencyResult raw = CalculateEfficiencyResult(rawDataNumerator->GetBinContent(iBin),
                                                             rawDataDenominator->GetBinContent(iBin));
      const EfficiencyResult truth = CalculateEfficiencyResult(truthNumerator->GetBinContent(iBin),
                                                               truthDenominator->GetBinContent(iBin));
      const EfficiencyResult corrected = CalculateRatioResult(raw, truth);

      std::cout << "[" << std::setw(2) << low << ", " << std::setw(2) << high << ") deg : ";
      PrintEfficiencyValue("raw Data", raw);
      std::cout << ", ";
      PrintEfficiencyValue("MC truth fraction", truth);
      std::cout << ", ";
      PrintEfficiencyValue("corrected Data", corrected);
      std::cout << std::endl;
    }

    const EfficiencyResult total_raw = CalculateHistogramTotalEfficiency(rawDataNumerator,
                                                                         rawDataDenominator);
    const EfficiencyResult total_truth = CalculateHistogramTotalEfficiency(truthNumerator,
                                                                           truthDenominator);
    const EfficiencyResult total_corrected = CalculateRatioResult(total_raw, total_truth);

    std::cout << "Total corrected Data FROST hit efficiency : ";
    PrintEfficiencyValue("raw Data", total_raw);
    std::cout << ", ";
    PrintEfficiencyValue("MC truth fraction", total_truth);
    std::cout << ", ";
    PrintEfficiencyValue("corrected Data", total_corrected);
    std::cout << std::endl;
    std::cout << "----------------------------------------" << std::endl;
  }

  bool ProcessSample(const char *inputDir,
                     const std::vector<std::string> &excludedFiles,
                     bool isMC,
                     const std::string &label,
                     const std::string &prefix,
                     int color,
                     int marker_style,
                     SampleHists &sample) {
    MakeSampleHists(sample, prefix, label, isMC);

    TChain chain("match_info");
    if (!BuildChain(inputDir,
                    excludedFiles,
                    label,
                    chain,
                    sample.n_files_added,
                    sample.n_files_excluded)) {
      return false;
    }
    sample.n_spills = chain.GetEntries();

    Int_t matched = 0;
    Int_t bsd_good_spill_flag = 0;

    std::vector<int> *trackmatch_has_match = nullptr;
    std::vector<int> *trackmatch_ninja_track_type = nullptr;
    std::vector<double> *trackmatch_expected_x = nullptr;
    std::vector<double> *trackmatch_expected_y = nullptr;
    std::vector<double> *trackmatch_dx = nullptr;
    std::vector<double> *trackmatch_dy = nullptr;
    std::vector<double> *trackmatch_dtanx = nullptr;
    std::vector<double> *trackmatch_dtany = nullptr;
    std::vector<double> *trackmatch_baby_mind_tangent_x = nullptr;
    std::vector<double> *trackmatch_baby_mind_tangent_y = nullptr;
    std::vector<int> *trackmatch_frost_is_hit = nullptr;
    std::vector<int> *trackmatch_true_frost_nearest_particle_id = nullptr;

    if (!isMC) {
      if (!chain.GetBranch("matched")) {
        std::cerr << "Error: required data branch 'matched' was not found." << std::endl;
        return false;
      }
      if (!chain.GetBranch("bsd_good_spill_flag")) {
        std::cerr << "Error: required data branch 'bsd_good_spill_flag' was not found." << std::endl;
        return false;
      }
      chain.SetBranchAddress("matched", &matched);
      chain.SetBranchAddress("bsd_good_spill_flag", &bsd_good_spill_flag);
    }

    if (isMC) {
      if (!chain.GetBranch("trackmatch_true_frost_nearest_particle_id")) {
        std::cerr << "Error: required MC truth branch 'trackmatch_true_frost_nearest_particle_id' was not found." << std::endl;
        return false;
      }
      chain.SetBranchAddress("trackmatch_true_frost_nearest_particle_id",
                             &trackmatch_true_frost_nearest_particle_id);
    }

    const std::vector<std::string> required_vector_branches = {
      "trackmatch_has_match",
      "trackmatch_ninja_track_type",
      "trackmatch_expected_x",
      "trackmatch_expected_y",
      "trackmatch_dx",
      "trackmatch_dy",
      "trackmatch_dtanx",
      "trackmatch_dtany",
      "trackmatch_baby_mind_tangent_x",
      "trackmatch_baby_mind_tangent_y",
      "trackmatch_frost_is_hit"
    };
    for (const auto &branch_name : required_vector_branches) {
      if (!chain.GetBranch(branch_name.c_str())) {
        std::cerr << "Error: required branch '" << branch_name
                  << "' was not found in " << label << " sample." << std::endl;
        return false;
      }
    }

    chain.SetBranchAddress("trackmatch_has_match", &trackmatch_has_match);
    chain.SetBranchAddress("trackmatch_ninja_track_type", &trackmatch_ninja_track_type);
    chain.SetBranchAddress("trackmatch_expected_x", &trackmatch_expected_x);
    chain.SetBranchAddress("trackmatch_expected_y", &trackmatch_expected_y);
    chain.SetBranchAddress("trackmatch_dx", &trackmatch_dx);
    chain.SetBranchAddress("trackmatch_dy", &trackmatch_dy);
    chain.SetBranchAddress("trackmatch_dtanx", &trackmatch_dtanx);
    chain.SetBranchAddress("trackmatch_dtany", &trackmatch_dtany);
    chain.SetBranchAddress("trackmatch_baby_mind_tangent_x", &trackmatch_baby_mind_tangent_x);
    chain.SetBranchAddress("trackmatch_baby_mind_tangent_y", &trackmatch_baby_mind_tangent_y);
    chain.SetBranchAddress("trackmatch_frost_is_hit", &trackmatch_frost_is_hit);

    const Long64_t nEntries = chain.GetEntries();
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
      chain.GetEntry(iEntry);

      if (!isMC) {
        if (matched != 1) continue;
        if (bsd_good_spill_flag == 0) continue;
      }

      if (!trackmatch_has_match ||
          !trackmatch_ninja_track_type ||
          !trackmatch_expected_x ||
          !trackmatch_expected_y ||
          !trackmatch_dx ||
          !trackmatch_dy ||
          !trackmatch_dtanx ||
          !trackmatch_dtany ||
          !trackmatch_baby_mind_tangent_x ||
          !trackmatch_baby_mind_tangent_y ||
          !trackmatch_frost_is_hit ||
          (isMC && !trackmatch_true_frost_nearest_particle_id)) {
        continue;
      }

      const std::size_t nTracks = trackmatch_has_match->size();
      if (trackmatch_ninja_track_type->size() != nTracks ||
          trackmatch_expected_x->size() != nTracks ||
          trackmatch_expected_y->size() != nTracks ||
          trackmatch_dx->size() != nTracks ||
          trackmatch_dy->size() != nTracks ||
          trackmatch_dtanx->size() != nTracks ||
          trackmatch_dtany->size() != nTracks ||
          trackmatch_baby_mind_tangent_x->size() != nTracks ||
          trackmatch_baby_mind_tangent_y->size() != nTracks ||
          trackmatch_frost_is_hit->size() != nTracks ||
          (isMC && trackmatch_true_frost_nearest_particle_id->size() != nTracks)) {
        std::cerr << "Warning: vector size mismatch in " << label
                  << " at spill entry " << iEntry
                  << ". Skip this spill." << std::endl;
        continue;
      }

      for (std::size_t iTrack = 0; iTrack < nTracks; ++iTrack) {
        const int ninjaTrackType = trackmatch_ninja_track_type->at(iTrack);
        const double expectedX = trackmatch_expected_x->at(iTrack);
        const double expectedY = trackmatch_expected_y->at(iTrack);
        const double dx = trackmatch_dx->at(iTrack);
        const double dy = trackmatch_dy->at(iTrack);
        const double dtanx = trackmatch_dtanx->at(iTrack);
        const double dtany = trackmatch_dtany->at(iTrack);
        const double tx = trackmatch_baby_mind_tangent_x->at(iTrack);
        const double ty = trackmatch_baby_mind_tangent_y->at(iTrack);

        const bool passTrackType = (ninjaTrackType == 1);
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

        sample.hDen->Fill(angleDeg);
        sample.hDenIsHit->Fill(angleDeg);

        bool hasTruthFrostMuon = false;
        if (isMC && trackmatch_true_frost_nearest_particle_id) {
          const int truthPid = trackmatch_true_frost_nearest_particle_id->at(iTrack);
          hasTruthFrostMuon = (truthPid == 13 || truthPid == -13);
          sample.hDenTruthMuon->Fill(angleDeg);
          if (hasTruthFrostMuon) {
            sample.hNumTruthMuon->Fill(angleDeg);
            sample.hDenIsHitGivenTruth->Fill(angleDeg);
          }
        }

        sample.hDxAll->Fill(dx);
        sample.hDyAll->Fill(dy);
        sample.hDtanxAll->Fill(dtanx);
        sample.hDtanyAll->Fill(dtany);

        for (int iBin = 0; iBin < kNBins; ++iBin) {
          if (angleXDeg >= kAngleBins[iBin] && angleXDeg < kAngleBins[iBin + 1]) {
            sample.hDxByAngleX[iBin]->Fill(dx);
            sample.hDtanxByAngleX[iBin]->Fill(dtanx);
            break;
          }
        }
        for (int iBin = 0; iBin < kNBins; ++iBin) {
          if (angleDeg >= kAngleBins[iBin] && angleDeg < kAngleBins[iBin + 1]) {
            sample.hDxByAngleTot[iBin]->Fill(dx);
            sample.hDtanxByAngleTot[iBin]->Fill(dtanx);
            break;
          }
        }
        for (int iBin = 0; iBin < kNBins; ++iBin) {
          if (angleYDeg >= kAngleBins[iBin] && angleYDeg < kAngleBins[iBin + 1]) {
            sample.hDyByAngleY[iBin]->Fill(dy);
            sample.hDtanyByAngleY[iBin]->Fill(dtany);
            break;
          }
        }
        for (int iBin = 0; iBin < kNBins; ++iBin) {
          if (angleDeg >= kAngleBins[iBin] && angleDeg < kAngleBins[iBin + 1]) {
            sample.hDyByAngleTot[iBin]->Fill(dy);
            sample.hDtanyByAngleTot[iBin]->Fill(dtany);
            break;
          }
        }

        if (trackmatch_has_match->at(iTrack) == 1) {
          sample.hNum->Fill(angleDeg);
        }
        if (trackmatch_frost_is_hit->at(iTrack) == 1) {
          sample.hNumIsHit->Fill(angleDeg);
          if (isMC && hasTruthFrostMuon) {
            sample.hNumIsHitGivenTruth->Fill(angleDeg);
          }
        }
      }
    }

    FillEfficiencyGraphs(sample, color, marker_style);
    PrintEfficiencySummary(sample);
    return true;
  }

  void StyleEfficiencyFrame(TH1D *frame,
                            double y_min = kEffYMin,
                            double y_max = kEffYMax,
                            double y_title_size = kEffYTitleSize) {
    frame->SetMinimum(y_min);
    frame->SetMaximum(y_max);
    frame->SetLineColor(0);
    frame->SetLineWidth(0);
    frame->SetFillStyle(0);
    frame->SetFillColor(0);
    frame->SetMarkerSize(0);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(y_title_size);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetLabelSize(0.05);
    frame->GetXaxis()->SetTitleOffset(1.0);
    frame->GetYaxis()->SetTitleOffset(1.3);
  }

  void DrawEfficiencyPage(TCanvas *canvas,
                          const char *outputPdfPath,
                          TH1D *frame,
                          TGraphAsymmErrors *mcGraph,
                          TGraphAsymmErrors *dataGraph,
                          const char *title,
                          double y_min = kEffYMin,
                          double y_max = kEffYMax,
                          double y_title_size = kEffYTitleSize) {
    canvas->Clear();
    canvas->SetGrid();
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    frame->SetTitle(title);
    StyleEfficiencyFrame(frame, y_min, y_max, y_title_size);
    frame->Draw();
    if (mcGraph) {
      mcGraph->Draw("P SAME");
    }
    if (dataGraph) {
      dataGraph->Draw("P SAME");
    }

    auto *legend = new TLegend(0.58, 0.20, 0.88, 0.34);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    if (mcGraph) legend->AddEntry(mcGraph, "MC", "lep");
    if (dataGraph) legend->AddEntry(dataGraph, "Data", "lep");
    legend->Draw();
    canvas->SaveAs(outputPdfPath);
  }

  void DrawCorrectedDataHitEfficiencyPage(TCanvas *canvas,
                                          const char *outputPdfPath,
                                          TH1D *frame,
                                          TGraphAsymmErrors *rawDataGraph,
                                          TGraphAsymmErrors *correctedDataGraph) {
    canvas->Clear();
    canvas->SetGrid();
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    frame->SetTitle("Corrected Data FROST hit efficiency vs angle");
    frame->GetYaxis()->SetTitle("FROST hit efficiency");
    StyleEfficiencyFrame(frame, kEffYMin, kCorrectedEffYMax, kEffYTitleSize);
    frame->Draw();

    if (rawDataGraph) {
      rawDataGraph->Draw("P SAME");
    }
    if (correctedDataGraph) {
      correctedDataGraph->Draw("P SAME");
    }

    auto *legend = new TLegend(0.50, 0.20, 0.88, 0.38);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    if (rawDataGraph) {
      legend->AddEntry(rawDataGraph, "Data raw", "lep");
    }
    if (correctedDataGraph) {
      legend->AddEntry(correctedDataGraph,
                       "Data corrected by MC truth fraction",
                       "lep");
    }
    legend->Draw();

    canvas->SaveAs(outputPdfPath);
  }

  double HistMaxWithErrors(const TH1D *hist) {
    if (!hist) {
      return 0.0;
    }
    double maximum = 0.0;
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
      maximum = std::max(maximum,
                         hist->GetBinContent(i) + hist->GetBinError(i));
    }
    return maximum;
  }

  double GraphMaxWithErrors(const TGraphAsymmErrors *graph) {
    if (!graph) {
      return 0.0;
    }

    double maximum = 0.0;
    for (int i = 0; i < graph->GetN(); ++i) {
      double x = 0.0;
      double y = 0.0;
      graph->GetPoint(i, x, y);
      maximum = std::max(maximum, y + graph->GetErrorYhigh(i));
    }
    return maximum;
  }

  void GarwoodInterval(double count,
                       double confidence_level,
                       double &lower,
                       double &upper) {
    lower = 0.0;
    upper = 0.0;

    if (count < 0.0 || !std::isfinite(count)) {
      return;
    }

    const double alpha = 1.0 - confidence_level;

    if (count > 0.0) {
      lower = 0.5 * ROOT::Math::chisquared_quantile(alpha / 2.0,
                                                    2.0 * count);
    } else {
      lower = 0.0;
    }

    upper = 0.5 * ROOT::Math::chisquared_quantile_c(alpha / 2.0,
                                                     2.0 * (count + 1.0));
  }

  TGraphAsymmErrors *MakePoissonDataGraph(const TH1D *dataHist,
                                          const TString &name,
                                          double confidence_level = 0.682689492137) {
    if (!dataHist) {
      return nullptr;
    }

    const int nbins = dataHist->GetNbinsX();
    std::vector<double> x(nbins, 0.0);
    std::vector<double> y(nbins, 0.0);
    std::vector<double> exl(nbins, 0.0);
    std::vector<double> exh(nbins, 0.0);
    std::vector<double> eyl(nbins, 0.0);
    std::vector<double> eyh(nbins, 0.0);

    for (int ibin = 1; ibin <= nbins; ++ibin) {
      const double low = dataHist->GetXaxis()->GetBinLowEdge(ibin);
      const double high = dataHist->GetXaxis()->GetBinUpEdge(ibin);
      const double center = 0.5 * (low + high);
      const double count = dataHist->GetBinContent(ibin);

      double lower = 0.0;
      double upper = 0.0;
      GarwoodInterval(count, confidence_level, lower, upper);

      x.at(ibin - 1) = center;
      y.at(ibin - 1) = count;
      exl.at(ibin - 1) = center - low;
      exh.at(ibin - 1) = high - center;
      eyl.at(ibin - 1) = std::max(0.0, count - lower);
      eyh.at(ibin - 1) = std::max(0.0, upper - count);
    }

    auto *graph = new TGraphAsymmErrors(nbins,
                                        x.data(),
                                        y.data(),
                                        exl.data(),
                                        exh.data(),
                                        eyl.data(),
                                        eyh.data());
    graph->SetName(name);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.8);
    graph->SetMarkerColor(kBlack);
    graph->SetLineColor(kBlack);
    graph->SetLineWidth(1);
    return graph;
  }

  struct ResidualStats {
    bool valid = false;
    double median = 0.0;
    double sigma68 = 0.0;
    double q16 = 0.0;
    double q84 = 0.0;
    double integral = 0.0;
  };

  double HistogramQuantileInVisibleRange(const TH1D *hist,
                                         double probability,
                                         double integral) {
    if (!hist || integral <= 0.0) {
      return 0.0;
    }

    if (probability <= 0.0) {
      return hist->GetXaxis()->GetXmin();
    }
    if (probability >= 1.0) {
      return hist->GetXaxis()->GetXmax();
    }

    const double target = probability * integral;
    double cumulative = 0.0;

    for (int ibin = 1; ibin <= hist->GetNbinsX(); ++ibin) {
      const double content = hist->GetBinContent(ibin);
      if (content <= 0.0 || !std::isfinite(content)) {
        continue;
      }

      const double next_cumulative = cumulative + content;
      if (next_cumulative >= target) {
        const double low = hist->GetXaxis()->GetBinLowEdge(ibin);
        const double high = hist->GetXaxis()->GetBinUpEdge(ibin);
        const double fraction = (target - cumulative) / content;
        return low + std::max(0.0, std::min(1.0, fraction)) * (high - low);
      }
      cumulative = next_cumulative;
    }

    return hist->GetXaxis()->GetXmax();
  }

  ResidualStats CalculateResidualStatsInVisibleRange(const TH1D *hist) {
    ResidualStats stats;
    if (!hist) {
      return stats;
    }

    for (int ibin = 1; ibin <= hist->GetNbinsX(); ++ibin) {
      const double content = hist->GetBinContent(ibin);
      if (std::isfinite(content) && content > 0.0) {
        stats.integral += content;
      }
    }

    if (stats.integral <= 0.0) {
      return stats;
    }

    stats.q16 = HistogramQuantileInVisibleRange(hist, 0.16, stats.integral);
    stats.median = HistogramQuantileInVisibleRange(hist, 0.50, stats.integral);
    stats.q84 = HistogramQuantileInVisibleRange(hist, 0.84, stats.integral);
    stats.sigma68 = 0.5 * (stats.q84 - stats.q16);
    stats.valid = std::isfinite(stats.median) && std::isfinite(stats.sigma68);
    return stats;
  }

  void DrawResidualStatsBox(const TH1D *mcHistOriginal,
                            const TH1D *dataHist) {
    const ResidualStats mcStats = CalculateResidualStatsInVisibleRange(mcHistOriginal);
    const ResidualStats dataStats = CalculateResidualStatsInVisibleRange(dataHist);

    auto *statsBox = new TPaveText(0.18, 0.72, 0.55, 0.88, "NDC");
    statsBox->SetFillColor(0);
    statsBox->SetFillStyle(0);
    statsBox->SetBorderSize(0);
    statsBox->SetTextAlign(12);
    statsBox->SetTextFont(42);
    statsBox->SetTextSize(0.028);

    TText *line = nullptr;
    if (mcStats.valid) {
      line = statsBox->AddText(Form("MC: median = %.3g, #sigma_{68} = %.3g",
                                    mcStats.median,
                                    mcStats.sigma68));
    } else {
      line = statsBox->AddText("MC: median = n/a, #sigma_{68} = n/a");
    }
    if (line) {
      line->SetTextColor(kBlue + 1);
    }

    if (dataStats.valid) {
      line = statsBox->AddText(Form("Data: median = %.3g, #sigma_{68} = %.3g",
                                    dataStats.median,
                                    dataStats.sigma68));
    } else {
      line = statsBox->AddText("Data: median = n/a, #sigma_{68} = n/a");
    }
    if (line) {
      line->SetTextColor(kBlack);
    }

    statsBox->Draw("same");
  }

  void DrawMcDataHistPage(TCanvas *canvas,
                          const char *outputPdfPath,
                          TH1D *mcHistOriginal,
                          TH1D *dataHist,
                          const char *title,
                          bool scaleMCHistToData) {
    canvas->Clear();
    canvas->SetGrid();
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);

    std::unique_ptr<TH1D> mcHist(dynamic_cast<TH1D*>(mcHistOriginal->Clone(
      Form("%s_draw", mcHistOriginal->GetName()))));
    mcHist->SetDirectory(nullptr);

    if (scaleMCHistToData) {
      const double mcIntegral = mcHist->Integral();
      const double dataIntegral = dataHist->Integral();
      if (mcIntegral > 0.0 && dataIntegral > 0.0) {
        mcHist->Scale(dataIntegral / mcIntegral);
      }
    }

    mcHist->SetTitle(title);
    mcHist->SetLineColor(kBlue + 1);
    mcHist->SetLineWidth(2);
    mcHist->SetMarkerSize(0);
    mcHist->GetXaxis()->SetTitleSize(0.045);
    mcHist->GetYaxis()->SetTitleSize(0.045);
    mcHist->GetXaxis()->SetLabelSize(0.045);
    mcHist->GetYaxis()->SetLabelSize(0.045);
    mcHist->GetYaxis()->SetTitleOffset(1.4);

    std::unique_ptr<TGraphAsymmErrors> dataGraph(
      MakePoissonDataGraph(dataHist,
                           Form("%s_poisson_points", dataHist->GetName())));

    const double maxY = std::max(HistMaxWithErrors(mcHist.get()),
                                 GraphMaxWithErrors(dataGraph.get()));
    mcHist->SetMaximum(maxY > 0.0 ? 1.15 * maxY : 1.0);
    mcHist->Draw("HIST");
    if (dataGraph) {
      dataGraph->Draw("P SAME");
    }

    auto *legend = new TLegend(0.58, 0.73, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(mcHist.get(),
                     scaleMCHistToData ? "MC (scaled to Data)" : "MC",
                     "l");
    if (dataGraph) {
      legend->AddEntry(dataGraph.get(), "Data", "lep");
    }
    legend->Draw();

    DrawResidualStatsBox(mcHistOriginal, dataHist);

    canvas->SaveAs(outputPdfPath);
  }

}  // namespace

void DrawTrackMatchEfficiency(
    const char *mcInputDir = "/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/6-TrackMatch_externalfit_PMandDWG",
    const char *dataInputDir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM_f8b4eb/2-rootfile_after_TrackMatch_externalfit_PMandDWG/all",
    const char *outputPdfPath = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM_f8b4eb/2-rootfile_after_TrackMatch_externalfit_PMandDWG/analysisplot/efficiency_mc_data_overlay.pdf",
    const char *logFilePath = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM_f8b4eb/2-rootfile_after_TrackMatch_externalfit_PMandDWG/analysisplot/efficiency_mc_data_overlay.log",
    const std::vector<std::string> &excludedMCFiles = {},
    const std::vector<std::string> &excludedDataFiles = {"b2physics_track_2025-11-29_00-00-00_afterTrackMatch.root","b2physics_track_2025-11-30_00-00-00_afterTrackMatch.root"},
    const bool scaleMCHistToData = true) {

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
  std::cout << "MC input dir: " << mcInputDir << std::endl;
  std::cout << "Data input dir: " << dataInputDir << std::endl;
  std::cout << "Distribution MC scaling: "
            << (scaleMCHistToData ? "MC scaled to Data integral" : "absolute counts")
            << std::endl;
  std::cout << "Data point errors on residual distributions: "
            << "central 68.27% Garwood Poisson intervals" << std::endl;
  std::cout << "MC truth FROST-muon definition: "
            << "trackmatch_true_frost_nearest_particle_id == +/-13" << std::endl;
  std::cout << "Corrected Data hit efficiency: "
            << "Data apparent is_hit efficiency divided by MC truth FROST-muon fraction"
            << std::endl;
  std::cout << "All-angle residual histogram bins: "
            << kResidualHistBinsAll << std::endl;
  std::cout << "Angle-binned residual histogram bins:" << std::endl;
  for (int i = 0; i < kNBins; ++i) {
    std::cout << "  [" << kAngleBins[i] << ", " << kAngleBins[i + 1]
              << ") deg : " << kResidualHistBinsByAngle[i]
              << " bins" << std::endl;
  }

  SampleHists mc;
  SampleHists data;

  const bool mcOk = ProcessSample(mcInputDir,
                                  excludedMCFiles,
                                  true,
                                  "MC",
                                  "mc",
                                  kBlue + 1,
                                  24,
                                  mc);
  const bool dataOk = ProcessSample(dataInputDir,
                                    excludedDataFiles,
                                    false,
                                    "Data",
                                    "data",
                                    kBlack,
                                    20,
                                    data);

  if (!mcOk || !dataOk) {
    std::cerr << "Error: failed to process MC or data sample." << std::endl;
    DeleteSampleHists(mc);
    DeleteSampleHists(data);
    std::cout.rdbuf(oldCoutBuf);
    std::cerr.rdbuf(oldCerrBuf);
    return;
  }

  std::unique_ptr<TGraphAsymmErrors> gDataHitEffCorrected(
    BuildRatioGraph(data.gEffIsHit,
                    mc.gEffTruthMuon,
                    "g_Data_EffIsHitCorrectedByTruthFraction",
                    kRed + 1,
                    21));
  PrintCorrectedDataHitEfficiencySummary(data.hNumIsHit,
                                         data.hDenIsHit,
                                         mc.hNumTruthMuon,
                                         mc.hDenTruthMuon);

  auto *canvas = new TCanvas("c_eff", "c_eff", 900, 700);
  canvas->SaveAs((std::string(outputPdfPath) + "[").c_str());

  DrawEfficiencyPage(canvas,
                     outputPdfPath,
                     data.hEff,
                     mc.gEff,
                     data.gEff,
                     "BM-FROST matching efficiency vs angle");

  DrawEfficiencyPage(canvas,
                     outputPdfPath,
                     data.hEffIsHit,
                     mc.gEffIsHit,
                     data.gEffIsHit,
                     "FROST is_hit efficiency vs angle");

  DrawEfficiencyPage(canvas,
                     outputPdfPath,
                     mc.hEffTruthMuon,
                     mc.gEffTruthMuon,
                     nullptr,
                     "MC truth FROST-muon fraction vs angle",
                     kEffYMin,
                     kEffYMax);

  DrawEfficiencyPage(canvas,
                     outputPdfPath,
                     mc.hEffIsHitGivenTruth,
                     mc.gEffIsHitGivenTruth,
                     nullptr,
                     "FROST is_hit efficiency for MC truth FROST muons",
                     kEffYMin,
                     kEffYMax,
                     kTruthHitEffYTitleSize);

  DrawCorrectedDataHitEfficiencyPage(canvas,
                                      outputPdfPath,
                                      data.hEffIsHit,
                                      data.gEffIsHit,
                                      gDataHitEffCorrected.get());

  DrawMcDataHistPage(canvas,
                     outputPdfPath,
                     mc.hDxAll,
                     data.hDxAll,
                     "dx distribution for all angles",
                     scaleMCHistToData);

  DrawMcDataHistPage(canvas,
                     outputPdfPath,
                     mc.hDyAll,
                     data.hDyAll,
                     "dy distribution for all angles",
                     scaleMCHistToData);

  DrawMcDataHistPage(canvas,
                     outputPdfPath,
                     mc.hDtanxAll,
                     data.hDtanxAll,
                     "dtanx distribution for all angles",
                     scaleMCHistToData);

  DrawMcDataHistPage(canvas,
                     outputPdfPath,
                     mc.hDtanyAll,
                     data.hDtanyAll,
                     "dtany distribution for all angles",
                     scaleMCHistToData);

  for (int i = 0; i < kNBins; ++i) {
    DrawMcDataHistPage(canvas,
                       outputPdfPath,
                       mc.hDxByAngleX[i],
                       data.hDxByAngleX[i],
                       Form("dx distribution: %.0f #leq #theta_{x} < %.0f deg",
                            kAngleBins[i], kAngleBins[i + 1]),
                       scaleMCHistToData);
  }

  for (int i = 0; i < kNBins; ++i) {
    DrawMcDataHistPage(canvas,
                       outputPdfPath,
                       mc.hDxByAngleTot[i],
                       data.hDxByAngleTot[i],
                       Form("dx distribution: %.0f #leq #theta < %.0f deg",
                            kAngleBins[i], kAngleBins[i + 1]),
                       scaleMCHistToData);
  }

  for (int i = 0; i < kNBins; ++i) {
    DrawMcDataHistPage(canvas,
                       outputPdfPath,
                       mc.hDyByAngleY[i],
                       data.hDyByAngleY[i],
                       Form("dy distribution: %.0f #leq #theta_{y} < %.0f deg",
                            kAngleBins[i], kAngleBins[i + 1]),
                       scaleMCHistToData);
  }

  for (int i = 0; i < kNBins; ++i) {
    DrawMcDataHistPage(canvas,
                       outputPdfPath,
                       mc.hDyByAngleTot[i],
                       data.hDyByAngleTot[i],
                       Form("dy distribution: %.0f #leq #theta < %.0f deg",
                            kAngleBins[i], kAngleBins[i + 1]),
                       scaleMCHistToData);
  }

  for (int i = 0; i < kNBins; ++i) {
    DrawMcDataHistPage(canvas,
                       outputPdfPath,
                       mc.hDtanxByAngleX[i],
                       data.hDtanxByAngleX[i],
                       Form("dtanx distribution: %.0f #leq #theta_{x} < %.0f deg",
                            kAngleBins[i], kAngleBins[i + 1]),
                       scaleMCHistToData);
  }

  for (int i = 0; i < kNBins; ++i) {
    DrawMcDataHistPage(canvas,
                       outputPdfPath,
                       mc.hDtanxByAngleTot[i],
                       data.hDtanxByAngleTot[i],
                       Form("dtanx distribution: %.0f #leq #theta < %.0f deg",
                            kAngleBins[i], kAngleBins[i + 1]),
                       scaleMCHistToData);
  }

  for (int i = 0; i < kNBins; ++i) {
    DrawMcDataHistPage(canvas,
                       outputPdfPath,
                       mc.hDtanyByAngleY[i],
                       data.hDtanyByAngleY[i],
                       Form("dtany distribution: %.0f #leq #theta_{y} < %.0f deg",
                            kAngleBins[i], kAngleBins[i + 1]),
                       scaleMCHistToData);
  }

  for (int i = 0; i < kNBins; ++i) {
    DrawMcDataHistPage(canvas,
                       outputPdfPath,
                       mc.hDtanyByAngleTot[i],
                       data.hDtanyByAngleTot[i],
                       Form("dtany distribution: %.0f #leq #theta < %.0f deg",
                            kAngleBins[i], kAngleBins[i + 1]),
                       scaleMCHistToData);
  }

  canvas->SaveAs((std::string(outputPdfPath) + "]").c_str());

  std::cout << "Saved PDF: " << outputPdfPath << std::endl;
  std::cout << "Saved log: " << logFilePath << std::endl;

  delete canvas;
  DeleteSampleHists(mc);
  DeleteSampleHists(data);

  std::cout.rdbuf(oldCoutBuf);
  std::cerr.rdbuf(oldCerrBuf);
}
