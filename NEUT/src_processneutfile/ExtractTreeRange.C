// File: ExtractTreeRange.C
// Usage from ROOT:
//   root -l -q 'ExtractTreeRange.C("input.root","/path/to/out",10,20)'

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TClass.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <algorithm>

// Recursively find the first TTree in a directory (supports arbitrary tree names and nested dirs).
static TTree* FindFirstTreeInDir(TDirectory* dir, TString& outFullPath)
{
  if (!dir) return nullptr;

  TIter nextKey(dir->GetListOfKeys());
  while (TKey* key = (TKey*)nextKey()) {
    const char* className = key->GetClassName();
    TClass* cl = TClass::GetClass(className);
    if (!cl) continue;

    // If it's a TTree (or derived), return it
    if (cl->InheritsFrom(TTree::Class())) {
      TObject* obj = key->ReadObj();
      TTree* tree = dynamic_cast<TTree*>(obj);
      if (tree) {
        outFullPath = TString(dir->GetPath()) + "/" + tree->GetName();
        return tree;
      }
      delete obj;
      continue;
    }

    // If it's a subdirectory, recurse
    if (cl->InheritsFrom(TDirectory::Class())) {
      TObject* obj = key->ReadObj();
      TDirectory* subdir = dynamic_cast<TDirectory*>(obj);
      if (subdir) {
        TString subPath;
        TTree* t = FindFirstTreeInDir(subdir, subPath);
        if (t) {
          outFullPath = subPath;
          return t;
        }
      }
      // Note: subdir is owned by the file; do not delete.
      continue;
    }
  }

  return nullptr;
}

// Extract entries [firstEntry, lastEntry] (inclusive) from the first TTree found in the file
// and write to: <outputDir>/<inputBase>_<first>_<effectiveLast>.root
void ExtractTreeRange(const char* inputRootFile,
                      const char* outputDir,
                      Long64_t firstEntry,
                      Long64_t lastEntry)
{
  // Basic argument checks
  if (!inputRootFile || TString(inputRootFile).IsNull()) {
    std::cerr << "Error: inputRootFile is empty." << std::endl;
    return;
  }
  if (!outputDir || TString(outputDir).IsNull()) {
    std::cerr << "Error: outputDir is empty." << std::endl;
    return;
  }
  if (firstEntry < 0) {
    std::cerr << "Error: firstEntry must be >= 0." << std::endl;
    return;
  }
  if (lastEntry < firstEntry) {
    std::cerr << "Error: lastEntry must be >= firstEntry." << std::endl;
    return;
  }

  // Open input file
  TFile* fin = TFile::Open(inputRootFile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Error: Could not open input file: " << inputRootFile << std::endl;
    return;
  }

  // Find a TTree (any name, anywhere)
  TString treePath;
  TTree* tin = FindFirstTreeInDir(fin, treePath);
  if (!tin) {
    std::cerr << "Error: No TTree found in input file: " << inputRootFile << std::endl;
    fin->Close();
    delete fin;
    return;
  }

  const Long64_t nEntries = tin->GetEntries();
  if (firstEntry >= nEntries) {
    std::cerr << "Error: firstEntry (" << firstEntry
              << ") >= number of entries (" << nEntries << ")." << std::endl;
    fin->Close();
    delete fin;
    return;
  }

  const Long64_t effectiveLast = std::min(lastEntry, nEntries - 1);

  // Ensure output directory exists
  gSystem->mkdir(outputDir, /*recursive=*/kTRUE);

  // Build output filename: <outputDir>/<base>_<first>_<effectiveLast>.root
  TString inPath(inputRootFile);
  TString base = gSystem->BaseName(inPath);   // e.g. "input.root"
  if (base.EndsWith(".root")) base.Resize(base.Length() - 5);  // remove ".root"
  TString outDirStr(outputDir);
  if (!outDirStr.EndsWith("/")) outDirStr += "/";

  TString outFileName;
  outFileName.Form("%s%s_%lld_%lld.root",
                   outDirStr.Data(),
                   base.Data(),
                   (Long64_t)firstEntry,
                   (Long64_t)effectiveLast);

  // Create output file
  TFile* fout = TFile::Open(outFileName, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Error: Could not create output file: " << outFileName << std::endl;
    fin->Close();
    delete fin;
    return;
  }

  // Clone tree structure only, then fill selected entries
  fout->cd();
  TTree* tout = tin->CloneTree(0);

  for (Long64_t i = firstEntry; i <= effectiveLast; ++i) {
    tin->GetEntry(i);
    tout->Fill();
  }

  // Write and close
  tout->Write();
  fout->Close();
  delete fout;

  fin->Close();
  delete fin;

  std::cout << "Wrote: " << outFileName
            << " (entries " << firstEntry << " to " << effectiveLast
            << ", total " << (effectiveLast - firstEntry + 1) << ")" << std::endl;
}
