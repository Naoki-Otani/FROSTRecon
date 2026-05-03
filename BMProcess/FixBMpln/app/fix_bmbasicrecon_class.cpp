#include <iostream>
#include <string>
#include <filesystem>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"

#include "BMBasicRecon.hpp"
#include "BMBeaminfo.hpp"
#include "BMBSD.hpp"

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input_dir> <output_dir>" << std::endl;
    return 1;
  }

  fs::path input_dir = argv[1];
  fs::path output_dir = argv[2];

  if (!fs::exists(input_dir) || !fs::is_directory(input_dir)) {
    std::cerr << "Error: input directory does not exist or is not a directory: "
              << input_dir << std::endl;
    return 1;
  }

  if (!fs::exists(output_dir)) {
    fs::create_directories(output_dir);
  }

  for (const auto& entry : fs::directory_iterator(input_dir)) {
    if (!entry.is_regular_file()) continue;
    if (entry.path().extension() != ".root") continue;

    fs::path input_path = entry.path();
    fs::path output_path = output_dir / (input_path.stem().string() + "_fixed.root");

    std::cout << "Processing: " << input_path << std::endl;

    TFile fin(input_path.string().c_str(), "READ");
    if (fin.IsZombie()) {
      std::cerr << "Error: cannot open " << input_path << std::endl;
      continue;
    }

    TTree* intree = dynamic_cast<TTree*>(fin.Get("tree"));
    if (!intree) {
      std::cerr << "Warning: tree not found in " << input_path << std::endl;
      continue;
    }

    BMBasicRecon* in_bm = nullptr;
    BMBeaminfo* in_beam = nullptr;
    BMBSD* in_bsd = nullptr;

    intree->SetBranchAddress("BMBasicRecon", &in_bm);
    intree->SetBranchAddress("BMBeaminfo", &in_beam);
    intree->SetBranchAddress("BMBSD", &in_bsd);

    TFile fout(output_path.string().c_str(), "RECREATE");
    TTree outtree("tree", "tree");

    auto out_bm = new BMBasicRecon();
    auto out_beam = new BMBeaminfo();
    auto out_bsd = new BMBSD();

    outtree.Branch("BMBasicRecon", "BMBasicRecon", &out_bm);
    outtree.Branch("BMBeaminfo", "BMBeaminfo", &out_beam);
    outtree.Branch("BMBSD", "BMBSD", &out_bsd);

    Long64_t nentries = intree->GetEntries();
    Long64_t fixed_count = 0;

    for (Long64_t i = 0; i < nentries; ++i) {
      intree->GetEntry(i);

      out_bm->Clear();
      *out_bm = *in_bm;
      *out_beam = *in_beam;
      *out_bsd = *in_bsd;

      std::size_t n = std::min({out_bm->mod.size(), out_bm->view.size(), out_bm->pln.size()});
      for (std::size_t j = 0; j < n; ++j) {
        if (out_bm->mod[j] == 5.0 && out_bm->view[j] == 0.0) {
          out_bm->pln[j] -= 1.0;
          ++fixed_count;
        }
      }

      outtree.Fill();
    }

    fout.cd();
    outtree.Write();

    delete out_bm;
    delete out_beam;
    delete out_bsd;

    std::cout << "fixed elements: " << fixed_count << std::endl;
  }

  return 0;
}
