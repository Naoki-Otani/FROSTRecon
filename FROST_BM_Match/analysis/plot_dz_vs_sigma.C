#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TPad.h>
#include <TF1.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TString.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace {

std::string Trim(const std::string &s) {
  const char *ws = " \t\r\n";
  const auto begin = s.find_first_not_of(ws);
  if (begin == std::string::npos) return "";
  const auto end = s.find_last_not_of(ws);
  return s.substr(begin, end - begin + 1);
}

bool ReadFitCsv(const std::string &csvPath,
                double &meanX, double &meanXErr,
                double &meanY, double &meanYErr,
                double &sigmaX, double &sigmaXErr,
                double &sigmaY, double &sigmaYErr) {
  std::ifstream fin(csvPath.c_str());
  if (!fin) {
    std::cerr << "Error: cannot open " << csvPath << std::endl;
    return false;
  }

  bool foundDx = false;
  bool foundDy = false;

  std::string line;
  while (std::getline(fin, line)) {
    line = Trim(line);
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::replace(line.begin(), line.end(), ',', ' ');

    std::stringstream ss(line);
    std::string quantity;
    double entries = 0.0;
    int fit_ok = 0;
    double mean = 0.0;
    double mean_error = 0.0;
    double sigma = 0.0;
    double sigma_error = 0.0;

    ss >> quantity >> entries >> fit_ok >> mean >> mean_error >> sigma >> sigma_error;
    if (!ss) {
      std::cerr << "Warning: malformed line in " << csvPath << std::endl;
      continue;
    }

    if (quantity == "dx") {
      meanX = mean;
      meanXErr = mean_error;
      sigmaX = sigma;
      sigmaXErr = sigma_error;
      foundDx = true;

      if (!fit_ok) {
        std::cerr << "Warning: dx fit_ok == 0 in " << csvPath << std::endl;
      }
    } else if (quantity == "dy") {
      meanY = mean;
      meanYErr = mean_error;
      sigmaY = sigma;
      sigmaYErr = sigma_error;
      foundDy = true;

      if (!fit_ok) {
        std::cerr << "Warning: dy fit_ok == 0 in " << csvPath << std::endl;
      }
    }
  }

  if (!foundDx || !foundDy) {
    std::cerr << "Error: dx/dy rows were not found in " << csvPath << std::endl;
    return false;
  }

  return true;
}

void FindYRange(const double *y,
                const double *ye,
                int n,
                double &ymin,
                double &ymax) {
  ymin = std::numeric_limits<double>::max();
  ymax = -std::numeric_limits<double>::max();

  for (int i = 0; i < n; ++i) {
    ymin = std::min(ymin, y[i] - ye[i]);
    ymax = std::max(ymax, y[i] + ye[i]);
  }

  const double margin = 0.1 * (ymax - ymin > 0.0 ? ymax - ymin : 1.0);
  ymin -= margin;
  ymax += margin;
}

struct CommonDzFitResult {
  bool fitOk = false;

  double bestDz = 0.0;
  double bestDzErr = 0.0;

  double curvatureX = 0.0;
  double curvatureXErr = 0.0;
  double minimumX = 0.0;
  double minimumXErr = 0.0;

  double curvatureY = 0.0;
  double curvatureYErr = 0.0;
  double minimumY = 0.0;
  double minimumYErr = 0.0;

  double chi2 = 0.0;
  int ndf = 0;
  int nUsedPoints = 0;
};

int gNFit = 0;
const double *gDzFit = nullptr;
const double *gSigmaXFit = nullptr;
const double *gSigmaXErrFit = nullptr;
const double *gSigmaYFit = nullptr;
const double *gSigmaYErrFit = nullptr;

void CommonDzChi2(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
  const double bestDz = par[0];

  const double curvatureX = par[1];
  const double minimumX = par[2];

  const double curvatureY = par[3];
  const double minimumY = par[4];

  double chi2 = 0.0;

  for (int i = 0; i < gNFit; ++i) {
    const double dz = gDzFit[i];

    if (gSigmaXErrFit[i] > 0.0) {
      const double modelX =
          curvatureX * (dz - bestDz) * (dz - bestDz) + minimumX;
      const double pullX = (gSigmaXFit[i] - modelX) / gSigmaXErrFit[i];
      chi2 += pullX * pullX;
    }

    if (gSigmaYErrFit[i] > 0.0) {
      const double modelY =
          curvatureY * (dz - bestDz) * (dz - bestDz) + minimumY;
      const double pullY = (gSigmaYFit[i] - modelY) / gSigmaYErrFit[i];
      chi2 += pullY * pullY;
    }
  }

  f = chi2;
}

CommonDzFitResult FitCommonBestDz(const int n,
                                  const double *dz,
                                  const double *sigmaX,
                                  const double *sigmaXErr,
                                  const double *sigmaY,
                                  const double *sigmaYErr) {
  gNFit = n;
  gDzFit = dz;
  gSigmaXFit = sigmaX;
  gSigmaXErrFit = sigmaXErr;
  gSigmaYFit = sigmaY;
  gSigmaYErrFit = sigmaYErr;

  CommonDzFitResult result;

  int nUsedPoints = 0;
  for (int i = 0; i < n; ++i) {
    if (sigmaXErr[i] > 0.0) ++nUsedPoints;
    if (sigmaYErr[i] > 0.0) ++nUsedPoints;
  }
  result.nUsedPoints = nUsedPoints;
  result.ndf = nUsedPoints - 5;

  if (result.ndf <= 0) {
    std::cerr << "Error: not enough valid points for common dz fit." << std::endl;
    return result;
  }

  TMinuit minuit(5);
  minuit.SetPrintLevel(0);
  minuit.SetFCN(CommonDzChi2);

  Int_t ierflg = 0;
  Double_t arglist[10];

  arglist[0] = 1.0;
  minuit.mnexcm("SET ERR", arglist, 1, ierflg);

  // Parameters:
  // [0] common best dz
  // [1] x curvature
  // [2] x minimum sigma
  // [3] y curvature
  // [4] y minimum sigma
  //
  // Initial values are hard-coded here.
  minuit.mnparm(0, "best dz",      0.0,    1.0,   -100.0, 100.0, ierflg);
  minuit.mnparm(1, "curvature x",  0.0001,  1e-5,     0.0,   1.0, ierflg);
  minuit.mnparm(2, "minimum x",   44.0,    1.0,      0.0, 500.0, ierflg);
  minuit.mnparm(3, "curvature y",  0.0001,  1e-5,     0.0,   1.0, ierflg);
  minuit.mnparm(4, "minimum y",   20.5,    1.0,      0.0, 500.0, ierflg);

  arglist[0] = 5000;
  arglist[1] = 1.0;
  minuit.mnexcm("MIGRAD", arglist, 2, ierflg);

  result.fitOk = (ierflg == 0);

  minuit.GetParameter(0, result.bestDz, result.bestDzErr);
  minuit.GetParameter(1, result.curvatureX, result.curvatureXErr);
  minuit.GetParameter(2, result.minimumX, result.minimumXErr);
  minuit.GetParameter(3, result.curvatureY, result.curvatureYErr);
  minuit.GetParameter(4, result.minimumY, result.minimumYErr);

  Double_t amin = 0.0;
  Double_t edm = 0.0;
  Double_t errdef = 0.0;
  Int_t nvpar = 0;
  Int_t nparx = 0;
  Int_t icstat = 0;
  minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  result.chi2 = amin;

  return result;
}

TPaveText *MakeCommonFitBox(const CommonDzFitResult &fit,
                            const char *axisLabel,
                            double x1 = 0.30,
                            double y1 = 0.65,
                            double x2 = 0.65,
                            double y2 = 0.88) {
  auto *box = new TPaveText(x1, y1, x2, y2, "NDC");
  box->SetFillColor(0);
  box->SetBorderSize(1);
  box->SetTextAlign(12);

  box->AddText("Common dz fit");
  box->AddText(Form("best dz = %.3f #pm %.3f mm",
                    fit.bestDz,
                    fit.bestDzErr));

  if (std::string(axisLabel) == "x") {
    box->AddText(Form("curv x = %.4g #pm %.4g",
                      fit.curvatureX,
                      fit.curvatureXErr));
    box->AddText(Form("min x = %.3f #pm %.3f mm",
                      fit.minimumX,
                      fit.minimumXErr));
  } else {
    box->AddText(Form("curv y = %.4g #pm %.4g",
                      fit.curvatureY,
                      fit.curvatureYErr));
    box->AddText(Form("min y = %.3f #pm %.3f mm",
                      fit.minimumY,
                      fit.minimumYErr));
  }

  box->AddText(Form("#chi^{2}/ndf = %.2f / %d",
                    fit.chi2,
                    fit.ndf));

  return box;
}

}  // namespace

void plot_dz_vs_sigma(
    const char *inputDir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/alignment",
    const char *outputPdfPath = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/alignment/dz_vs_sigma.pdf") {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  const int N = 21;

  double dz[N];
  double dzErr[N] = {0.0};

  double sigmaX[N];
  double sigmaXErr[N];
  double sigmaY[N];
  double sigmaYErr[N];

  double meanX[N];
  double meanXErr[N];
  double meanY[N];
  double meanYErr[N];

  for (int i = 0; i < N; ++i) {
    dz[i] = -100.0 + 10.0 * i;

    meanX[i] = 0.0;
    meanXErr[i] = 0.0;
    meanY[i] = 0.0;
    meanYErr[i] = 0.0;

    sigmaX[i] = 0.0;
    sigmaXErr[i] = 0.0;
    sigmaY[i] = 0.0;
    sigmaYErr[i] = 0.0;

    const int z = static_cast<int>(dz[i]);

    std::ostringstream oss;
    oss << inputDir << "/alignment_zshift_" << z << "_fit.csv";
    const std::string csvPath = oss.str();

    if (!ReadFitCsv(csvPath,
                    meanX[i], meanXErr[i],
                    meanY[i], meanYErr[i],
                    sigmaX[i], sigmaXErr[i],
                    sigmaY[i], sigmaYErr[i])) {
      std::cerr << "Error while reading: " << csvPath << std::endl;
      return;
    }
  }

  double yminX, ymaxX;
  double yminY, ymaxY;
  FindYRange(sigmaX, sigmaXErr, N, yminX, ymaxX);
  FindYRange(sigmaY, sigmaYErr, N, yminY, ymaxY);

  const CommonDzFitResult commonFit =
      FitCommonBestDz(N,
                      dz,
                      sigmaX,
                      sigmaXErr,
                      sigmaY,
                      sigmaYErr);

  if (!commonFit.fitOk) {
    std::cerr << "Warning: common dz fit may not have converged." << std::endl;
  }

  std::cout << "Common best dz: "
            << commonFit.bestDz << " +/- "
            << commonFit.bestDzErr << " mm" << std::endl;

  std::cout << "Common fit chi2/ndf: "
            << commonFit.chi2 << " / "
            << commonFit.ndf << std::endl;

  TCanvas *c = new TCanvas("c_dz_vs_sigma", "dz vs sigma", 1400, 600);
  c->Divide(2, 1);

  c->cd(1);
  gPad->SetGrid();

  TGraphErrors *grX = new TGraphErrors(N, dz, sigmaX, dzErr, sigmaXErr);
  grX->SetTitle("x;dz [mm];#sigma_{x} [mm]");
  grX->SetMarkerStyle(20);
  grX->SetMarkerSize(1.1);
  grX->SetLineWidth(2);
  grX->GetXaxis()->SetLimits(-105.0, 105.0);
  grX->SetMinimum(yminX);
  grX->SetMaximum(ymaxX);
  grX->GetXaxis()->SetTitleOffset(0.9);
  grX->GetYaxis()->SetTitleOffset(0.9);
  grX->GetXaxis()->SetTitleSize(0.05);
  grX->GetYaxis()->SetTitleSize(0.05);
  grX->GetXaxis()->SetLabelSize(0.05);
  grX->GetYaxis()->SetLabelSize(0.05);
  grX->Draw("AP");

  TF1 *fitX = new TF1("fitX",
                      "[0]*(x-[1])*(x-[1]) + [2]",
                      -100.0,
                      100.0);
  fitX->SetParNames("curvature x", "common best dz", "minimum x");
  fitX->SetParameters(commonFit.curvatureX,
                      commonFit.bestDz,
                      commonFit.minimumX);
  fitX->SetLineWidth(2);
  fitX->Draw("SAME");

  TPaveText *boxX = MakeCommonFitBox(commonFit, "x");
  boxX->Draw();

  c->cd(2);
  gPad->SetGrid();

  TGraphErrors *grY = new TGraphErrors(N, dz, sigmaY, dzErr, sigmaYErr);
  grY->SetTitle("y;dz [mm];#sigma_{y} [mm]");
  grY->SetMarkerStyle(21);
  grY->SetMarkerSize(1.1);
  grY->SetLineWidth(2);
  grY->GetXaxis()->SetLimits(-105.0, 105.0);
  grY->SetMinimum(yminY);
  grY->SetMaximum(ymaxY);
  grY->GetXaxis()->SetTitleOffset(0.9);
  grY->GetYaxis()->SetTitleOffset(0.9);
  grY->GetXaxis()->SetTitleSize(0.05);
  grY->GetYaxis()->SetTitleSize(0.05);
  grY->GetXaxis()->SetLabelSize(0.05);
  grY->GetYaxis()->SetLabelSize(0.05);
  grY->Draw("AP");

  TF1 *fitY = new TF1("fitY",
                      "[0]*(x-[1])*(x-[1]) + [2]",
                      -100.0,
                      100.0);
  fitY->SetParNames("curvature y", "common best dz", "minimum y");
  fitY->SetParameters(commonFit.curvatureY,
                      commonFit.bestDz,
                      commonFit.minimumY);
  fitY->SetLineWidth(2);
  fitY->Draw("SAME");

  TPaveText *boxY = MakeCommonFitBox(commonFit, "y");
  boxY->Draw();

  c->SaveAs(outputPdfPath);
  std::cout << "Saved PDF: " << outputPdfPath << std::endl;
}
