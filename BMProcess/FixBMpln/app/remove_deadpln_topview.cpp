#include <iostream>
#include <string>
#include <filesystem>
#include <vector>
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
    fs::path output_path = output_dir / (input_path.stem().string() + ".root");

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
    Long64_t removed_count = 0;

    for (Long64_t i = 0; i < nentries; ++i) {
      intree->GetEntry(i);

      out_bm->Clear();
      *out_bm = *in_bm;
      *out_beam = *in_beam;
      *out_bsd = *in_bsd;

      std::vector<std::size_t> keep_indices;
      keep_indices.reserve(out_bm->mod.size());

      std::size_t n = std::min({
        out_bm->mod.size(),
        out_bm->view.size(),
        out_bm->pln.size()
      });

      for (std::size_t j = 0; j < n; ++j) {
        bool remove_hit =
          (out_bm->mod[j] == 5.0) &&
          (out_bm->view[j] == 1.0) &&
          (
            out_bm->pln[j] == 9.0 ||
            out_bm->pln[j] == 10.0 ||
            out_bm->pln[j] == 17.0
          );

        if (remove_hit) {
          ++removed_count;
        } else {
          keep_indices.push_back(j);
        }
      }

      auto filter_vector = [&](auto& vec) {
        using T = typename std::decay_t<decltype(vec)>::value_type;
        std::vector<T> filtered;
        filtered.reserve(keep_indices.size());
        for (std::size_t idx : keep_indices) {
          if (idx < vec.size()) filtered.push_back(vec[idx]);
        }
        vec = std::move(filtered);
      };

      filter_vector(out_bm->mod);
      filter_vector(out_bm->view);
      filter_vector(out_bm->pln);
      filter_vector(out_bm->channel);
      filter_vector(out_bm->HG);
      filter_vector(out_bm->LHG);
      filter_vector(out_bm->RHG);
      filter_vector(out_bm->THG);
      filter_vector(out_bm->BHG);
      filter_vector(out_bm->Lgain);
      filter_vector(out_bm->Rgain);
      filter_vector(out_bm->Tgain);
      filter_vector(out_bm->Bgain);
      filter_vector(out_bm->Lpe);
      filter_vector(out_bm->Rpe);
      filter_vector(out_bm->Tpe);
      filter_vector(out_bm->Bpe);
      filter_vector(out_bm->LG);
      filter_vector(out_bm->Ltime);
      filter_vector(out_bm->Ftime);
      filter_vector(out_bm->Htime);
      filter_vector(out_bm->timedif);
      filter_vector(out_bm->bunch);

      outtree.Fill();
    }

    fout.cd();
    outtree.Write();

    delete out_bm;
    delete out_beam;
    delete out_bsd;

    std::cout << "Removed hits: " << removed_count << std::endl;
  }

  return 0;
}
