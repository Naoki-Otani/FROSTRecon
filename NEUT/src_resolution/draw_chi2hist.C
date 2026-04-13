// File: draw_chi2hist.C
//
// Usage example:
//   root -l -q 'draw_chi2hist.C("input.root","/output/hoge.pdf","/output/hoge.txt",0)'
//
// This macro reads the TTree "frost" from an input ROOT file and creates
// three chi2 histograms under the following event selections:
//
// (1) Events with exactly one counted hit particle, and that particle is a muon
//     (pid = +/-13).
//
// (2) Events with exactly two counted hit particles, at least one of which is
//     a muon (pid = +/-13), and the hit-position distance between the two
//     particles is at least 9.2 mm along either the x axis or the y axis.
//     The hit positions are evaluated using posx and posy.
//
// (3) Events with two or more counted hit particles, and at least one of the
//     counted particles is a muon (pid = +/-13).
//
// Important note:
//   The number of hit particles is NOT taken from numhitparticle.
//   Instead, it is defined as the number of particles satisfying
//
//       (*vertexposz)[i] < -5.25 && (*energydeposit)[i] > 0.1
//
//   Only particles satisfying this condition are counted and used in the
//   selections above.

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TLegend.h>

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>

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

void draw_chi2hist(const char* inputRootFile="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/FROSTrecon_result/rootfile_afterrecon/frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection_afterrecon_1.5.root", const char* outputPdfFile="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/FROSTrecon_result/chi2/chi2_distribution_frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection_afterrecon_1.5.pdf", const char* outputTxtFile="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/FROSTrecon_result/chi2/chi2_threshold_results_frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection_afterrecon_1.5.txt", int bunch=0)
{
    // -------------------------------------------------------------------------
    // Basic style setup
    // -------------------------------------------------------------------------
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1110);
    gStyle->SetTitleFontSize(0.04);

    // -------------------------------------------------------------------------
    // Open input ROOT file
    // -------------------------------------------------------------------------
    TFile* fin = TFile::Open(inputRootFile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: failed to open input ROOT file: "
                  << inputRootFile << std::endl;
        return;
    }

    // -------------------------------------------------------------------------
    // Get the TTree "frost"
    // -------------------------------------------------------------------------
    TTree* tree = dynamic_cast<TTree*>(fin->Get("frost"));
    if (!tree) {
        std::cerr << "Error: TTree \"frost\" was not found in file: "
                  << inputRootFile << std::endl;
        fin->Close();
        return;
    }

    // -------------------------------------------------------------------------
    // Branch variables
    // -------------------------------------------------------------------------
    std::vector<int>*    pid           = nullptr;
    std::vector<double>* posx          = nullptr;
    std::vector<double>* posy          = nullptr;
    std::vector<double>* posz          = nullptr;
    std::vector<double>* momx          = nullptr;
    std::vector<double>* momy          = nullptr;
    std::vector<double>* momz          = nullptr;
    std::vector<double>* energydeposit = nullptr;
    std::vector<double>* vertexposz    = nullptr;
    std::vector<double>* chi2          = nullptr;
    int numhitparticle = 0;

    tree->SetBranchAddress("pid",           &pid);
    tree->SetBranchAddress("posx",          &posx);
    tree->SetBranchAddress("posy",          &posy);
    tree->SetBranchAddress("posz",          &posz);
    tree->SetBranchAddress("momx",          &momx);
    tree->SetBranchAddress("momy",          &momy);
    tree->SetBranchAddress("momz",          &momz);
    tree->SetBranchAddress("energydeposit", &energydeposit);
    tree->SetBranchAddress("vertexposz",    &vertexposz);
    tree->SetBranchAddress("chi2",          &chi2);
    tree->SetBranchAddress("numhitparticle", &numhitparticle);

    // -------------------------------------------------------------------------
    // Histogram definitions
    //
    // The binning below is chosen as a reasonable default. Adjust the range
    // if your chi2 distribution is known more precisely.
    // -------------------------------------------------------------------------
    TH1D* h1_small = new TH1D(
        "h1_small",
        "#chi^{2} distribution: exactly 1 hit particle and muon;#chi^{2}/ndf;Number of events",
        200, 0.0, 5.0
    );

    TH1D* h1_large = new TH1D(
        "h1_large",
        "#chi^{2} distribution: exactly 1 hit particle and muon;#chi^{2}/ndf;Number of events",
        200, 0.0, 50.0
    );

    TH1D* h2_small = new TH1D(
        "h2_small",
        "#chi^{2} distribution: exactly 2 hit particles including muon, |#Deltax| or |#Deltay| >= 9.2 mm;#chi^{2}/ndf;Number of events",
        200, 0.0, 5.0
    );

    TH1D* h2_large = new TH1D(
        "h2_large",
        "#chi^{2} distribution: exactly 2 hit particles including muon, |#Deltax| or |#Deltay| >= 9.2 mm;#chi^{2}/ndf;Number of events",
        200, 0.0, 50.0
    );

    TH1D* h3 = new TH1D(
        "h3",
        "#chi^{2} distribution: >=2 hit particles including muon;#chi^{2}/ndf;Number of events",
        200, 0.0, 50.0
    );

    // -------------------------------------------------------------------------
    // chi2 threshold
    // -------------------------------------------------------------------------
    const int nThresholds = 500;
    double chi2_threshold[nThresholds];
    for (int i = 0; i < nThresholds; ++i) {
        chi2_threshold[i] = 1.0 + 0.01 * i;
    }
    int count_under_threshold_singlehit[nThresholds]={0};
    int count_under_threshold_twohit[nThresholds]={0};
    int count_under_threshold_multihit[nThresholds]={0};

    int total_singlehit = 0;
    int total_twohit = 0;
    int total_multihit = 0;

    // -------------------------------------------------------------------------
    // Event loop
    // -------------------------------------------------------------------------
    const Long64_t nEntries = tree->GetEntries();

    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        // Safety check for vector branches
        if (!pid || !posx || !posy || !energydeposit || !vertexposz) {
            continue;
        }
        if (!chi2) {
            continue;
        }
        if (bunch < 0 || bunch >= static_cast<int>(chi2->size())) {
            continue;
        }
        const double chi2_value = (*chi2)[bunch];

        const size_t nParticles = pid->size();

        // Ensure consistent vector sizes before using the data
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

        // ---------------------------------------------------------------------
        // Build a list of particle indices that satisfy the counting condition:
        //
        //   (*vertexposz)[i] < -5.25 && (*energydeposit)[i] > 0.1
        //
        // Only these particles are considered as "counted hit particles".
        // ---------------------------------------------------------------------
        std::vector<int> countedIndices;
        countedIndices.reserve(nParticles);

        for (size_t i = 0; i < nParticles; ++i) {
            if ((*vertexposz)[i] < -5.25 && (*energydeposit)[i] > 0.1) {
                countedIndices.push_back(static_cast<int>(i));
            }
        }

        const int nCounted = static_cast<int>(countedIndices.size());

        // ---------------------------------------------------------------------
        // Check whether at least one counted particle is a muon and is detected by Baby MIND.
        // ---------------------------------------------------------------------
        bool hasBMDetectMuon = false;
        std::vector<int> countedMuonIndices;

        for (int idx : countedIndices) {
            if (std::abs((*pid)[idx]) == 13) {
                countedMuonIndices.push_back(idx);

                const double pz = (*momz)[idx];
                if (std::abs(pz) < 1.0e-12) continue;

                const double xTrue = (*posx)[idx] + (*momx)[idx] / pz * (0.0 - (*posz)[idx]);
                const double yTrue = (*posy)[idx] + (*momy)[idx] / pz * (0.0 - (*posz)[idx]);

                if (IsBMDetect(xTrue, yTrue, (*momx)[idx], (*momy)[idx], (*momz)[idx])) {
                    hasBMDetectMuon = true;
                }
            }
        }
        // ---------------------------------------------------------------------
        // (1) Exactly 1 counted hit particle, and that particle is a muon detected by Baby MIND
        // ---------------------------------------------------------------------
        if (nCounted == 1) {
            const int idx0 = countedIndices[0];
            if (std::abs((*pid)[idx0]) == 13) {
                const double pz = (*momz)[idx0];
                if (std::abs(pz) > 1.0e-12) {
                    const double xTrue = (*posx)[idx0] + (*momx)[idx0] / pz * (0.0 - (*posz)[idx0]);
                    const double yTrue = (*posy)[idx0] + (*momy)[idx0] / pz * (0.0 - (*posz)[idx0]);

                    if (IsBMDetect(xTrue, yTrue, (*momx)[idx0], (*momy)[idx0], (*momz)[idx0])) {
                        h1_small->Fill(chi2_value);
                        h1_large->Fill(chi2_value);
                        total_singlehit++;
                        for (int i = 0; i < nThresholds; ++i) {
                            if (chi2_value < chi2_threshold[i]) {
                                count_under_threshold_singlehit[i]++;
                            }
                        }
                    }
                }
            }
        }

        // ---------------------------------------------------------------------
        // (2) Exactly 2 counted hit particles
        //     - at least one is a muon detected by Baby MIND
        //     - distance along x OR y is at least 9.2 mm
        //
        // The hit positions are evaluated using posx and posy.
        // ---------------------------------------------------------------------
        if (nCounted == 2) {
            const int idx0 = countedIndices[0];
            const int idx1 = countedIndices[1];

            bool muonIncludedBM = false;

            for (int idx : countedIndices) {
                if (std::abs((*pid)[idx]) != 13) continue;

                const double pz = (*momz)[idx];
                if (std::abs(pz) < 1.0e-12) continue;

                const double xTrue = (*posx)[idx] + (*momx)[idx] / pz * (0.0 - (*posz)[idx]);
                const double yTrue = (*posy)[idx] + (*momy)[idx] / pz * (0.0 - (*posz)[idx]);

                if (IsBMDetect(xTrue, yTrue, (*momx)[idx], (*momy)[idx], (*momz)[idx])) {
                    muonIncludedBM = true;
                    break;
                }
            }

            const double dx = std::abs((*posx)[idx0] - (*posx)[idx1]);
            const double dy = std::abs((*posy)[idx0] - (*posy)[idx1]);

            if (muonIncludedBM && (dx >= 9.2 || dy >= 9.2)) {
                total_twohit++;
                h2_small->Fill(chi2_value);
                h2_large->Fill(chi2_value);
                for (int i = 0; i < nThresholds; ++i) {
                    if (chi2_value < chi2_threshold[i]) {
                        count_under_threshold_twohit[i]++;
                    }
                }
            }
        }

        // ---------------------------------------------------------------------
        // (3) Two or more counted hit particles, and at least one is a muon detected by Baby MIND
        // ---------------------------------------------------------------------
        if (nCounted >= 2 && hasBMDetectMuon) {
            h3->Fill(chi2_value);
            total_multihit++;
            for (int i = 0; i < nThresholds; ++i) {
                if (chi2_value < chi2_threshold[i]) {
                    count_under_threshold_multihit[i]++;
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Draw histograms and export them into a multi-page PDF
    // -------------------------------------------------------------------------
    TCanvas* c = new TCanvas("c", "c", 900, 700);

    TString pdfName(outputPdfFile);

    // Open multi-page PDF
    c->Print(pdfName + "[");

    // Page 1
    c->Clear();
    h1_small->SetLineWidth(2);
    h1_small->Draw("HIST");
    c->Print(pdfName);

    // Page 2
    c->Clear();
    h1_large->SetLineWidth(2);
    h1_large->Draw("HIST");
    c->Print(pdfName);

    // Page 3
    c->Clear();
    h2_small->SetLineWidth(2);
    h2_small->Draw("HIST");
    c->Print(pdfName);

    // Page 4
    c->Clear();
    h2_large->SetLineWidth(2);
    h2_large->Draw("HIST");
    c->Print(pdfName);

    // Page 5
    c->Clear();
    h3->SetLineWidth(2);
    h3->Draw("HIST");
    c->Print(pdfName);

    // Close multi-page PDF
    c->Print(pdfName + "]");

    // -------------------------------------------------------------------------
    // Cleanup
    // -------------------------------------------------------------------------
    delete c;
    delete h1_small;
    delete h1_large;
    delete h2_small;
    delete h2_large;
    delete h3;

    fin->Close();
    delete fin;

    std::cout << "Done: output PDF was written to " << outputPdfFile << std::endl;

    // -------------------------------------------------------------------------
    // Write chi2 threshold results to a text file
    // -------------------------------------------------------------------------
    std::ofstream outFile(outputTxtFile);
    outFile << "Threshold\tSingle-hit count[%]\tTwo-hit count[%]\tMulti-hit count[%]\n";
    for (int i = 0; i < nThresholds; ++i) {
        outFile << chi2_threshold[i] << "\t" << (total_singlehit > 0 ? (double)count_under_threshold_singlehit[i] / total_singlehit * 100 : 0) << "\t" << (total_twohit > 0 ? (double)count_under_threshold_twohit[i] / total_twohit * 100 : 0) << "\t" << (total_multihit > 0 ? (double)count_under_threshold_multihit[i] / total_multihit * 100 : 0) << "\n";
    }
    outFile.close();
    std::cout << "Chi2 threshold results were written to " << outputTxtFile << std::endl;
}
