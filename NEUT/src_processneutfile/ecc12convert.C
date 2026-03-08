#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>

#define xoff 2.51805
#define yoff 14.424
#define zoff -1.48418

// Convert vertex coordinates and write a new ROOT file　(9ECC configuration to 12ECC configuration).
// - Keep ALL original branches by cloning the input tree structure.
// - Only modify EvtVtx; other branches are copied as-is.
// - For ECC4-6, duplicate the event with +/- 0.175 shift in Y (same as original code).
void ecc12convert(
    std::string inputfile  = "/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/NEUToutput/rawroot/neut_320kA_1.077e21pot.root",
    std::string outputfile = "/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/NEUToutput/rawroot/neut_320kA_1.077e21pot_ECC12.root"
) {
    // Open input file
    TFile* fin = TFile::Open(inputfile.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: Could not open input file: " << inputfile << std::endl;
        return;
    }

    // Get input tree
    TTree* tin = dynamic_cast<TTree*>(fin->Get("nRooTracker"));
    if (!tin) {
        std::cerr << "Error: Could not find TTree 'nRooTracker' in: " << inputfile << std::endl;
        fin->Close();
        return;
    }

    // We only need to explicitly bind EvtVtx because we modify it.
    // Other branches will be copied as-is when we Fill the cloned tree.
    Double_t EvtVtx[4];
    tin->SetBranchAddress("EvtVtx", EvtVtx);

    // Create output file
    TFile* fout = TFile::Open(outputfile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: Could not create output file: " << outputfile << std::endl;
        fin->Close();
        return;
    }

    // Clone the structure of the input tree (ALL branches), but keep it empty (0 entries).
    // Since we set branch address for EvtVtx on tin before cloning,
    // the cloned tree will use the same memory for branch buffers.
    TTree* tout = tin->CloneTree(0);

    // Helper: determine X bin (0:left, 1:middle, 2:right)
    auto xBin = [](double x) -> int {
        if (x < -2.675) return 0;
        if (x < -2.36)  return 1;  // -2.675 <= x < -2.36
        return 2;                  // x >= -2.36
    };

    // Helper: determine Y region (0:top, 1:middle, 2:bottom)
    // top    : y > -14.25
    // middle : -14.6 < y < -14.25
    // bottom : y < -14.6
    auto yRegion = [](double y) -> int {
        if (y > -14.25) return 0;
        if (y > -14.6)  return 1;  // -14.6 < y <= -14.25 (matches original strict/loose behavior closely)
        return 2;
    };

    // Helper: apply coordinate offsets using the original vertex and fill once
    auto fillFromOriginal = [&](const Double_t vtx0[4], double yShift) {
        EvtVtx[0] = vtx0[0] + xoff;
        EvtVtx[1] = vtx0[1] + yShift + yoff;
        EvtVtx[2] = vtx0[2] + zoff;
        EvtVtx[3] = vtx0[3]; // keep time unchanged
        tout->Fill();
    };

    const Long64_t nEntries = tin->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tin->GetEntry(i);

        // Save the original vertex before modifying it
        const Double_t vtx0[4] = {EvtVtx[0], EvtVtx[1], EvtVtx[2], EvtVtx[3]};

        const int xb = xBin(vtx0[0]);
        const int yr = yRegion(vtx0[1]);

        if (xb < 0 || xb > 2) continue;

        if (yr == 0) {
            // ECC1-3
            fillFromOriginal(vtx0, +0.175);
        } else if (yr == 1) {
            // ECC4-6 (duplicate, symmetric shifts around the original)
            fillFromOriginal(vtx0, +0.175);
            fillFromOriginal(vtx0, -0.175);
        } else {
            // ECC7-9
            fillFromOriginal(vtx0, -0.175);
        }
    }

    fout->cd();
    tout->Write();
    fout->Close();
    fin->Close();

    std::cout << "Wrote: " << outputfile << std::endl;
}
