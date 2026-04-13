// divide_tree_every_1000.C
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <iostream>
#include <string>
#include <sstream>

void divide_tree_every_1000(const char* inputFileName, const char* outputDir)
{
    const Long64_t chunkSize = 1000;

    // Create the output directory
    gSystem->mkdir(outputDir, true);

    // Open the input file
    TFile* inputFile = TFile::Open(inputFileName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: cannot open input file: " << inputFileName << std::endl;
        return;
    }

    // Get the tree
    TTree* tree = dynamic_cast<TTree*>(inputFile->Get("tree"));
    if (!tree) {
        std::cerr << "Error: cannot find TTree named \"tree\" in " << inputFileName << std::endl;
        inputFile->Close();
        return;
    }

    Long64_t nEntries = tree->GetEntries();
    std::cout << "Input file   : " << inputFileName << std::endl;
    std::cout << "Output dir   : " << outputDir << std::endl;
    std::cout << "Total entries: " << nEntries << std::endl;

    // Build the base name from the input file name
    std::string inputPath(inputFileName);
    std::string baseName = gSystem->BaseName(inputPath.c_str());

    // Remove the ".root" extension
    if (baseName.size() >= 5 && baseName.substr(baseName.size() - 5) == ".root") {
        baseName = baseName.substr(0, baseName.size() - 5);
    }

    Long64_t nChunks = (nEntries + chunkSize - 1) / chunkSize;

    for (Long64_t iChunk = 0; iChunk < nChunks; ++iChunk) {
        Long64_t firstEntry = iChunk * chunkSize;
        Long64_t lastEntry  = firstEntry + chunkSize - 1;
        if (lastEntry >= nEntries) lastEntry = nEntries - 1;

        std::ostringstream outName;
        outName << outputDir << "/" << baseName << "_" << iChunk << ".root";

        std::cout << "Writing: " << outName.str()
                  << "  entries " << firstEntry << " - " << lastEntry << std::endl;

        TFile* outFile = TFile::Open(outName.str().c_str(), "RECREATE");
        if (!outFile || outFile->IsZombie()) {
            std::cerr << "Error: cannot create output file: " << outName.str() << std::endl;
            continue;
        }

        outFile->cd();

        // Create an empty clone and fill it with only the required entries
        TTree* outTree = tree->CloneTree(0);
        if (!outTree) {
            std::cerr << "Error: failed to clone tree structure." << std::endl;
            outFile->Close();
            delete outFile;
            continue;
        }

        for (Long64_t i = firstEntry; i <= lastEntry; ++i) {
            tree->GetEntry(i);
            outTree->Fill();
        }

        outTree->Write();
        outFile->Close();
        delete outFile;
    }

    inputFile->Close();
    delete inputFile;

    std::cout << "Done." << std::endl;
}
