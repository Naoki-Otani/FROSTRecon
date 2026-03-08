// File: MakeLightYieldMaxHists.C
//
// Usage:
//   root -l -b -q 'MakeLightYieldMaxHists.C("/path/to/file.root","out.pdf")'
//
// This macro selects events where:
//  - exactly one particle satisfies (vertexposz < -5.25) AND (energydeposit > 0.1)
//  - and that particle is a muon (abs(pid)==13)
// Then it fills histograms of:
//  - max(lightyieldx[132])
//  - max(lightyieldy[140])
//  - 2nd max(lightyieldx[132])
//  - 2nd max(lightyieldy[140])
//
// Notes:
//  - "2nd max" means the second-largest value (not necessarily distinct).
//  - vertexposz threshold is used as given (units assumed to match your tree).
//  - energydeposit threshold is 0.1 (units assumed to match your tree).

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>

#include <vector>
#include <cmath>
#include <iostream>
#include <string>

static inline void FindMaxAndSecondMax(const double* arr, int n, double& max1, double& max2)
{
  // Find the largest and second-largest values in one pass.
  max1 = -1.0e300;
  max2 = -1.0e300;

  for (int i = 0; i < n; ++i) {
    const double v = arr[i];
    if (v > max1) {
      max2 = max1;
      max1 = v;
    } else if (v > max2) {
      max2 = v;
    }
  }

  // If n==1, max2 stays very negative; clamp to max1 for safety.
  if (n < 2) max2 = max1;
}

void MakeLightyieldMaxHists(
  int threshold_max = 20.,
  int threshold_2ndmax = 5.,
  const char* inFile =
    "/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/FROSToutput/aftermppccorrection/"
    "frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection.root",
  const char* outPdf = "/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/plot/lightyield_max_2ndmax.pdf"
){
  // Open input file
  TFile* fin = TFile::Open(inFile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "[ERROR] Cannot open input file: " << inFile << std::endl;
    return;
  }

  // Get tree
  TTree* t = dynamic_cast<TTree*>(fin->Get("wls"));
  if (!t) {
    std::cerr << "[ERROR] Cannot find TTree 'wls' in: " << inFile << std::endl;
    fin->Close();
    delete fin;
    return;
  }

  // Branches
  Int_t numhitparticle = 0;
  std::vector<int>* pid = nullptr;
  std::vector<double>* energydeposit = nullptr;
  std::vector<double>* vertexposz = nullptr;

  double lightyieldx[132];
  double lightyieldy[140];

  t->SetBranchAddress("numhitparticle", &numhitparticle);
  t->SetBranchAddress("pid", &pid);
  t->SetBranchAddress("energydeposit", &energydeposit);
  t->SetBranchAddress("vertexposz", &vertexposz);
  t->SetBranchAddress("lightyieldx", lightyieldx);
  t->SetBranchAddress("lightyieldy", lightyieldy);

  // Histograms (binning: choose something reasonable; you can tune later)
  auto* hMaxX   = new TH1D("hMaxX",   "Max light yield (x);Max lightyield (x) [p.e.];Number of events", 100, 0.0, 200.0);
  auto* hMaxY   = new TH1D("hMaxY",   "Max light yield (y);Max lightyield (y) [p.e.];Number of events", 100, 0.0, 200.0);
  auto* h2ndMaxX= new TH1D("h2ndMaxX","2nd max lightyield (x);2nd max lightyield (x) [p.e.];Number of events", 100, 0.0, 100.0);
  auto* h2ndMaxY= new TH1D("h2ndMaxY","2nd max lightyield (y);2nd max lightyield (y) [p.e.];Number of events", 100, 0.0, 100.0);

  // Enable proper errors if needed
  hMaxX->Sumw2(); hMaxY->Sumw2(); h2ndMaxX->Sumw2(); h2ndMaxY->Sumw2();

  int countMaxX=0;
  int countMaxY=0;
  int count2ndMaxX=0;
  int count2ndMaxY=0;

  const Long64_t nEntries = t->GetEntries();
  Long64_t nSelected = 0;

  for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {
    t->GetEntry(ievt);

    if (!pid || !energydeposit || !vertexposz) continue;
    if (pid->size() != energydeposit->size() || pid->size() != vertexposz->size()) continue;

    // Count particles satisfying (vertexposz < -5.25) AND (energydeposit > 0.1)
    int count = 0;
    int idxMu = -1;

    for (size_t i = 0; i < pid->size(); ++i) {
      if ((*vertexposz)[i] < -5.25 && (*energydeposit)[i] > 0.1) {
        ++count;
        idxMu = static_cast<int>(i);
      }
    }

    // Require exactly one such particle and it must be a muon
    if (count != 1) continue;
    if (idxMu < 0) continue;
    if (std::abs((*pid)[idxMu]) != 13) continue;

    // Compute max and 2nd max for X/Y arrays
    double maxX, secondX, maxY, secondY;
    FindMaxAndSecondMax(lightyieldx, 132, maxX, secondX);
    FindMaxAndSecondMax(lightyieldy, 140, maxY, secondY);

    if(maxX > threshold_max) ++countMaxX;
    if(maxY > threshold_max) ++countMaxY;
    if(secondX > threshold_2ndmax) ++count2ndMaxX;
    if(secondY > threshold_2ndmax) ++count2ndMaxY;
    hMaxX->Fill(maxX);
    h2ndMaxX->Fill(secondX);
    hMaxY->Fill(maxY);
    h2ndMaxY->Fill(secondY);

    ++nSelected;
  }

  std::cout << "[INFO] Total entries: " << nEntries << "\n";
  std::cout << "[INFO] Selected events: " << nSelected << "\n";
  std::cout << "[INFO] Count maxX > " << threshold_max << ": " << countMaxX << "/" << nSelected << "=" << (nSelected > 0 ? static_cast<double>(countMaxX)*100. / nSelected : 0.0) << "\%\n";
  std::cout << "[INFO] Count maxY > " << threshold_max << ": " << countMaxY << "/" << nSelected << "=" << (nSelected > 0 ? static_cast<double>(countMaxY)*100. / nSelected : 0.0) << "\%\n";
  std::cout << "[INFO] Count 2nd maxX > " << threshold_2ndmax << ": " << count2ndMaxX << "/" << nSelected << "=" << (nSelected > 0 ? static_cast<double>(count2ndMaxX)*100. / nSelected : 0.0) << "\%\n";
  std::cout << "[INFO] Count 2nd maxY > " << threshold_2ndmax << ": " << count2ndMaxY << "/" << nSelected << "=" << (nSelected > 0 ? static_cast<double>(count2ndMaxY)*100. / nSelected : 0.0) << "\%\n";

  // Draw and save PDF
  gStyle->SetOptStat(1110);

  hMaxX->SetFillColor(kYellow);
  h2ndMaxX->SetFillColor(kYellow);
  hMaxY->SetFillColor(kYellow);
  h2ndMaxY->SetFillColor(kYellow);

  TCanvas* c = new TCanvas("c", "lightyield maxima", 1200, 900);
  c->Divide(2,2);

  c->cd(1); hMaxX->Draw("HIST");
  c->cd(2); hMaxY->Draw("HIST");
  c->cd(3); h2ndMaxX->Draw("HIST");
  c->cd(4); h2ndMaxY->Draw("HIST");

  c->SaveAs(outPdf);

  // Cleanup
  fin->Close();
  delete fin;

  // (Histograms/canvas can be kept in memory for interactive ROOT; delete if you prefer)
}
