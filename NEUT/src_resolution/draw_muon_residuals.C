// File: draw_muon_residuals.C
//
// Usage example:
//   root -l -q 'draw_muon_residuals.C("input.root","output.pdf","output.root",0)'
//
// This macro reads the TTree "frost" from an input ROOT file and creates
// residual histograms for
//
//   x_rec - x_true
//   y_rec - y_true
//
// for muons (pid = +/-13) contained in the event.
//
// Event / particle selection:
//   A particle is counted as a "hit particle" only if
//
//       vertexposz[i] < -5.25 && energydeposit[i] > 0.1
//
//   The target events are events that contain at least one counted muon and are detected by Baby MIND.
//
// True muon position:
//   The true position is defined at z = 0.
//   It is obtained by extrapolating the true track using
//   (posx, posy, posz, momx, momy, momz):
//
//       x_true = posx + momx / momz * (0 - posz)
//       y_true = posy + momy / momz * (0 - posz)
//
// Reconstructed position:
//   x_rec and y_rec candidates are treated independently.
//   For each muon,
//     - x_rec is chosen as the candidate closest to x_true
//     - y_rec is chosen as the candidate closest to y_true
//
// Histogram categories:
//   For both x and y residuals, histograms are created for
//
//     1) all selected muons
//     2) exactly 1 counted hit particle
//     3) exactly 2 counted hit particles
//     4) exactly 3 counted hit particles
//     5) 4 or more counted hit particles
//
// If an event contains multiple counted muons, each muon contributes
// independently to the histograms.
//
// Output:
//   - A multi-page PDF containing all histograms with Gaussian fits
//   - A ROOT file containing all histograms, fit functions, and canvases
//
// Note:
//   In the current frost tree, x_rec and y_rec are stored as
//   std::vector<std::vector<double>> with bunch index.
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TPad.h>
#include <TH1.h>
#include <TLegend.h>

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm>

TF1* FitHistogramWithGaussian(TH1D* hist, const std::string& fitName, double xmin, double xmax)
{
    if (!hist) return nullptr;
    if (hist->GetEntries() < 5) return nullptr;

    TF1* gaus = new TF1(fitName.c_str(), "gaus", xmin, xmax);
    gaus->SetParent(nullptr);
    hist->Fit(gaus, "RQ");

    return gaus;
}

int FindThetaBin(double tanThetaAbs)
{
    const double edges[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.2};
    const int nBins = 8;

    for (int i = 0; i < nBins; ++i) {
        if (tanThetaAbs >= edges[i] && tanThetaAbs < edges[i + 1]) {
            return i;
        }
    }
    return -1;
}

void SetResolutionBinFromFit(TH1D* hRes, int bin, TF1* fit)
{
    if (!hRes || !fit) return;

    const double mean    = fit->GetParameter(1);
    const double sigma   = fit->GetParameter(2);
    const double meanErr = fit->GetParError(1);
    const double sigErr  = fit->GetParError(2);

    const double res2 = mean * mean + sigma * sigma;
    if (res2 <= 0.0) return;

    const double res = std::sqrt(res2);

    // Error propagation:
    // res = sqrt(mean^2 + sigma^2)
    // dres^2 = (mean/res)^2 * dmean^2 + (sigma/res)^2 * dsigma^2
    const double resErr = std::sqrt(
        (mean  * mean  * meanErr * meanErr +
         sigma * sigma * sigErr  * sigErr) / res2
    );

    hRes->SetBinContent(bin, res);
    hRes->SetBinError(bin, resErr);
}

const double L = 750; //mm //distance b/w tracker and Baby MIND
bool IsBMDetect(double x, double y,double momx, double momy, double momz){
  //x,y at FROST, momx,momy,momz at FROST
  double xbm, ybm; //(x,y) at second layer of Baby MIND
  xbm=x+(momx/momz)*L;
  ybm=y+(momy/momz)*L;
  if(xbm<2155 && xbm>-615 && ybm<965 && ybm>-985){
    return true;
  }else{
    return false;
  }
}

void draw_muon_residuals(const char* inputRootFile="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/FROSTrecon_result/rootfile_afterrecon/frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection_afterrecon_4.0.root",
                         const char* outputPdfFile="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/FROSTrecon_result/resolution/resolution_frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection_afterrecon_4.0.pdf",
                         const char* outputRootFile="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/FROSTrecon_result/resolution/resolution_frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection_afterrecon_4.0.root",
                         int bunch=0)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    gStyle->SetTitleFontSize(0.040);
    TH1::AddDirectory(kFALSE);

    TFile* fin = TFile::Open(inputRootFile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: failed to open input ROOT file: "
                  << inputRootFile << std::endl;
        return;
    }

    TTree* tree = dynamic_cast<TTree*>(fin->Get("frost"));
    if (!tree) {
        std::cerr << "Error: TTree \"frost\" was not found in: "
                  << inputRootFile << std::endl;
        fin->Close();
        return;
    }

    std::vector<int>*    pid            = nullptr;
    std::vector<double>* posx           = nullptr;
    std::vector<double>* posy           = nullptr;
    std::vector<double>* posz           = nullptr;
    std::vector<double>* momx           = nullptr;
    std::vector<double>* momy           = nullptr;
    std::vector<double>* momz           = nullptr;
    std::vector<double>* energydeposit  = nullptr;
    std::vector<double>* vertexposz     = nullptr;
    std::vector<std::vector<double>>* x_rec = nullptr;
    std::vector<std::vector<double>>* y_rec = nullptr;

    tree->SetBranchAddress("pid",           &pid);
    tree->SetBranchAddress("posx",          &posx);
    tree->SetBranchAddress("posy",          &posy);
    tree->SetBranchAddress("posz",          &posz);
    tree->SetBranchAddress("momx",          &momx);
    tree->SetBranchAddress("momy",          &momy);
    tree->SetBranchAddress("momz",          &momz);
    tree->SetBranchAddress("energydeposit", &energydeposit);
    tree->SetBranchAddress("vertexposz",    &vertexposz);
    tree->SetBranchAddress("x_rec",         &x_rec);
    tree->SetBranchAddress("y_rec",         &y_rec);

    const int    nBins = 200;
    const double histMin = -10.0;
    const double histMax =  10.0;

    enum Category {
        kAll = 0,
        kHit1,
        kHit2,
        kHit3,
        kHit4p,
        kNCategory
    };

    const char* catName[kNCategory] = {
        "all",
        "hit1",
        "hit2",
        "hit3",
        "hit4p"
    };

    const char* catTitle[kNCategory] = {
        "All selected muons",
        "Exactly 1 counted hit particle",
        "Exactly 2 counted hit particles",
        "Exactly 3 counted hit particles",
        "4 or more counted hit particles"
    };

    const int nThetaBins = 8;
    const double thetaEdges[nThetaBins + 1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.2};

    const char* thetaLabels[nThetaBins] = {
        "0.0-0.1",
        "0.1-0.2",
        "0.2-0.3",
        "0.3-0.4",
        "0.4-0.5",
        "0.5-0.6",
        "0.6-0.8",
        "0.8-1.2"
    };

    const char* thetaDisplay[nThetaBins] = {
    "[0.0,0.1)",
    "[0.1,0.2)",
    "[0.2,0.3)",
    "[0.3,0.4)",
    "[0.4,0.5)",
    "[0.5,0.6)",
    "[0.6,0.8)",
    "[0.8,1.2)"
    };

    TH1D* hX[kNCategory];
    TH1D* hY[kNCategory];
    TH1D* hXTheta[kNCategory][nThetaBins];
    TH1D* hYTheta[kNCategory][nThetaBins];
    TH1D* hResVsThetaX[kNCategory];
    TH1D* hResVsThetaY[kNCategory];

    for (int i = 0; i < kNCategory; ++i) {
        hX[i] = new TH1D(
            Form("h_xres_%s", catName[i]),
            Form("%s; x_{rec} - x_{true} [mm]; Number of events", catTitle[i]),
            nBins, histMin, histMax
        );
        hX[i]->SetDirectory(nullptr);

        hY[i] = new TH1D(
            Form("h_yres_%s", catName[i]),
            Form("%s; y_{rec} - y_{true} [mm]; Number of events", catTitle[i]),
            nBins, histMin, histMax
        );
        hY[i]->SetDirectory(nullptr);

        hX[i]->SetLineWidth(2);
        hY[i]->SetLineWidth(2);
        hX[i]->SetFillColor(kYellow);
        hY[i]->SetFillColor(kYellow);
    }

    for (int icat = 0; icat < kNCategory; ++icat) {
      for (int i = 0; i < nThetaBins; ++i) {
          hXTheta[icat][i] = new TH1D(
              Form("h_xres_theta_%s_%s", catName[icat], thetaLabels[i]),
              Form("%s, |tan#theta_{x}| in %s; x_{rec} - x_{true} [mm]; Number of events",
                  catTitle[icat], thetaDisplay[i]),
              nBins / 2, histMin, histMax
          );
          hXTheta[icat][i]->SetDirectory(nullptr);
          hXTheta[icat][i]->SetLineWidth(2);
          hXTheta[icat][i]->SetFillColor(kYellow);

          hYTheta[icat][i] = new TH1D(
              Form("h_yres_theta_%s_%s", catName[icat], thetaLabels[i]),
              Form("%s, |tan#theta_{y}| in %s; y_{rec} - y_{true} [mm]; Number of events",
                  catTitle[icat], thetaDisplay[i]),
              nBins / 2, histMin, histMax
          );
          hYTheta[icat][i]->SetDirectory(nullptr);
          hYTheta[icat][i]->SetLineWidth(2);
          hYTheta[icat][i]->SetFillColor(kYellow);
      }
    }

    for (int icat = 0; icat < kNCategory; ++icat) {
      hResVsThetaX[icat] = new TH1D(
          Form("h_res_vs_theta_x_%s", catName[icat]),
          Form("%s; tan#theta; Position resolution [mm]", catTitle[icat]),
          nThetaBins, thetaEdges
      );
      hResVsThetaX[icat]->SetDirectory(nullptr);
      hResVsThetaX[icat]->SetLineColor(kRed);
      hResVsThetaX[icat]->SetMarkerColor(kRed);
      hResVsThetaX[icat]->SetMarkerStyle(20);
      hResVsThetaX[icat]->SetLineWidth(2);
      hResVsThetaX[icat]->SetStats(0);

      hResVsThetaY[icat] = new TH1D(
          Form("h_res_vs_theta_y_%s", catName[icat]),
          Form("%s; tan#theta; Position resolution [mm]", catTitle[icat]),
          nThetaBins, thetaEdges
      );
      hResVsThetaY[icat]->SetDirectory(nullptr);
      hResVsThetaY[icat]->SetLineColor(kBlue);
      hResVsThetaY[icat]->SetMarkerColor(kBlue);
      hResVsThetaY[icat]->SetMarkerStyle(21);
      hResVsThetaY[icat]->SetLineWidth(2);
      hResVsThetaY[icat]->SetStats(0);
  }

    const Long64_t nEntries = tree->GetEntries();

    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        if (!pid || !posx || !posy || !posz || !momx || !momy || !momz ||
            !energydeposit || !vertexposz || !x_rec || !y_rec) {
            continue;
        }

        const size_t nParticles = pid->size();
        if (posx->size()          != nParticles ||
            posy->size()          != nParticles ||
            posz->size()          != nParticles ||
            momx->size()          != nParticles ||
            momy->size()          != nParticles ||
            momz->size()          != nParticles ||
            energydeposit->size() != nParticles ||
            vertexposz->size()    != nParticles) {
            continue;
        }

        if (!x_rec || !y_rec) {
            continue;
        }

        if (bunch < 0 ||
            bunch >= static_cast<int>(x_rec->size()) ||
            bunch >= static_cast<int>(y_rec->size())) {
            continue;
        }

        const std::vector<double>& x_rec_b = (*x_rec)[bunch];
        const std::vector<double>& y_rec_b = (*y_rec)[bunch];

        std::vector<int> countedIndices;
        countedIndices.reserve(nParticles);

        for (size_t i = 0; i < nParticles; ++i) {
            if ((*vertexposz)[i] < -5.25 && (*energydeposit)[i] > 0.1) {
                countedIndices.push_back(static_cast<int>(i));
            }
        }

        const int nCounted = static_cast<int>(countedIndices.size());
        if (nCounted < 1) {
            continue;
        }

        std::vector<int> countedMuonIndices;
        for (int idx : countedIndices) {
            if (std::abs((*pid)[idx]) == 13) {
                countedMuonIndices.push_back(idx);
            }
        }

        if (countedMuonIndices.empty()) {
            continue;
        }

        int multiplicityCategory = -1;
        if (nCounted == 1)      multiplicityCategory = kHit1;
        else if (nCounted == 2) multiplicityCategory = kHit2;
        else if (nCounted == 3) multiplicityCategory = kHit3;
        else                    multiplicityCategory = kHit4p;

        for (int muIdx : countedMuonIndices) {
            const double pz = (*momz)[muIdx];

            if (std::abs(pz) < 1.0e-12) {
                continue;
            }

            const double xTrue = (*posx)[muIdx] + (*momx)[muIdx] / pz * (0.0 - (*posz)[muIdx]);
            const double yTrue = (*posy)[muIdx] + (*momy)[muIdx] / pz * (0.0 - (*posz)[muIdx]);

            if(!IsBMDetect(xTrue, yTrue, (*momx)[muIdx], (*momy)[muIdx], (*momz)[muIdx])){
              continue;
            }

            // Select the x reconstruction candidate independently for the selected bunch
            int bestXIdx = -1;
            double bestAbsDx = std::numeric_limits<double>::max();

            for (size_t ix = 0; ix < x_rec_b.size(); ++ix) {
                const double absDx = std::abs(x_rec_b[ix] - xTrue);
                if (absDx < bestAbsDx) {
                    bestAbsDx = absDx;
                    bestXIdx = static_cast<int>(ix);
                }
            }

            // Select the y reconstruction candidate independently for the selected bunch
            int bestYIdx = -1;
            double bestAbsDy = std::numeric_limits<double>::max();

            for (size_t iy = 0; iy < y_rec_b.size(); ++iy) {
                const double absDy = std::abs(y_rec_b[iy] - yTrue);
                if (absDy < bestAbsDy) {
                    bestAbsDy = absDy;
                    bestYIdx = static_cast<int>(iy);
                }
            }

            if (bestXIdx < 0 || bestYIdx < 0) {
                continue;
            }

            const double xResidual = x_rec_b[bestXIdx] - xTrue;
            const double yResidual = y_rec_b[bestYIdx] - yTrue;

            hX[kAll]->Fill(xResidual);
            hY[kAll]->Fill(yResidual);

            hX[multiplicityCategory]->Fill(xResidual);
            hY[multiplicityCategory]->Fill(yResidual);

            const double tanThetaX = std::abs((*momx)[muIdx] / pz);
            const double tanThetaY = std::abs((*momy)[muIdx] / pz);

            const int thetaBinX = FindThetaBin(tanThetaX);
            const int thetaBinY = FindThetaBin(tanThetaY);

            if (thetaBinX >= 0) {
                hXTheta[kAll][thetaBinX]->Fill(xResidual);
                hXTheta[multiplicityCategory][thetaBinX]->Fill(xResidual);
            }

            if (thetaBinY >= 0) {
                hYTheta[kAll][thetaBinY]->Fill(yResidual);
                hYTheta[multiplicityCategory][thetaBinY]->Fill(yResidual);
            }
        }
    }

    TFile* fout = TFile::Open(outputRootFile, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: failed to create output ROOT file: "
                  << outputRootFile << std::endl;
        fin->Close();
        return;
    }

    TString pdfName(outputPdfFile);
    TCanvas* c = new TCanvas("c_residuals", "Muon residuals", 1400, 600);

    c->Print(pdfName + "[");

    for (int i = 0; i < kNCategory; ++i) {
        c->Clear();
        c->Divide(2, 1);

        c->cd(1);
        TF1* fitX = FitHistogramWithGaussian(hX[i], Form("fit_x_%s", catName[i]), -5, 5);
        hX[i]->Draw("HIST");
        if (fitX) {
            fitX->SetLineWidth(2);
            fitX->Draw("SAME");
        }

        c->cd(2);
        TF1* fitY = FitHistogramWithGaussian(hY[i], Form("fit_y_%s", catName[i]), -5, 5);
        hY[i]->Draw("HIST");
        if (fitY) {
            fitY->SetLineWidth(2);
            fitY->Draw("SAME");
        }

        c->Print(pdfName);

        fout->cd();
        if (fitX) fitX->Write();
        if (fitY) fitY->Write();
    }

    for (int icat = 0; icat < kNCategory; ++icat) {
      for (int i = 0; i < nThetaBins; ++i) {
          c->Clear();
          c->cd();

          TF1* fitX = FitHistogramWithGaussian(
              hXTheta[icat][i],
              Form("fit_x_theta_%s_%d", catName[icat], i),
              -5, 5
          );

          hXTheta[icat][i]->Draw("HIST");
          if (fitX) {
              fitX->SetLineWidth(2);
              fitX->Draw("SAME");
              SetResolutionBinFromFit(hResVsThetaX[icat], i + 1, fitX);
          }

          c->Print(pdfName);

          fout->cd();
          if (fitX) fitX->Write();
      }
    }

    for (int icat = 0; icat < kNCategory; ++icat) {
      for (int i = 0; i < nThetaBins; ++i) {
          c->Clear();
          c->cd();

          TF1* fitY = FitHistogramWithGaussian(
              hYTheta[icat][i],
              Form("fit_y_theta_%s_%d", catName[icat], i),
              -5, 5
          );

          hYTheta[icat][i]->Draw("HIST");
          if (fitY) {
              fitY->SetLineWidth(2);
              fitY->Draw("SAME");
              SetResolutionBinFromFit(hResVsThetaY[icat], i + 1, fitY);
          }

          c->Print(pdfName);

          fout->cd();
          if (fitY) fitY->Write();
      }
    }

    gStyle->SetOptStat(0);

    for (int icat = 0; icat < kNCategory; ++icat) {
        c->Clear();
        c->cd();

        double ymax = 0.0;
        for (int ibin = 1; ibin <= hResVsThetaX[icat]->GetNbinsX(); ++ibin) {
            ymax = std::max(ymax, hResVsThetaX[icat]->GetBinContent(ibin) + hResVsThetaX[icat]->GetBinError(ibin));
            ymax = std::max(ymax, hResVsThetaY[icat]->GetBinContent(ibin) + hResVsThetaY[icat]->GetBinError(ibin));
        }

        hResVsThetaX[icat]->SetMinimum(0.0);
        hResVsThetaX[icat]->SetMaximum(1.2 * ymax);
        hResVsThetaX[icat]->Draw("E1");
        hResVsThetaY[icat]->Draw("E1 SAME");

        auto leg = new TLegend(0.80, 0.20, 0.88, 0.40);
        leg->AddEntry(hResVsThetaX[icat], "x", "lep");
        leg->AddEntry(hResVsThetaY[icat], "y", "lep");
        leg->Draw();

        c->Print(pdfName);
    }

    c->Print(pdfName + "]");

    fout->cd();
    for (int i = 0; i < kNCategory; ++i) {
        hX[i]->Write();
        hY[i]->Write();
    }

    for (int icat = 0; icat < kNCategory; ++icat) {
        for (int i = 0; i < nThetaBins; ++i) {
            hXTheta[icat][i]->Write();
            hYTheta[icat][i]->Write();
        }
        hResVsThetaX[icat]->Write();
        hResVsThetaY[icat]->Write();
    }

    delete c;

    fout->Close();
    fin->Close();

    for (int i = 0; i < kNCategory; ++i) {
        delete hX[i];
        delete hY[i];
    }

    for (int icat = 0; icat < kNCategory; ++icat) {
        for (int i = 0; i < nThetaBins; ++i) {
            delete hXTheta[icat][i];
            delete hYTheta[icat][i];
        }
        delete hResVsThetaX[icat];
        delete hResVsThetaY[icat];
    }


    std::cout << "Done." << std::endl;
    std::cout << "  Output PDF : " << outputPdfFile  << std::endl;
    std::cout << "  Output ROOT: " << outputRootFile << std::endl;
}
