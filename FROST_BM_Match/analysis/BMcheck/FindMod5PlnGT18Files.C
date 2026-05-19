#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

bool HasMod5PlnGT18(const std::string& root_path) {
  TFile* file = TFile::Open(root_path.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[WARN] Failed to open: " << root_path << std::endl;
    if (file) file->Close();
    return false;
  }

  TTree* tree = nullptr;
  file->GetObject("tree", tree);
  if (!tree) {
    std::cerr << "[WARN] Tree 'tree' was not found: " << root_path << std::endl;
    file->Close();
    return false;
  }

  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<double>> mod(reader, "mod");
  TTreeReaderValue<std::vector<double>> pln(reader, "pln");

  while (reader.Next()) {
    const size_t n_hits = std::min(mod->size(), pln->size());

    for (size_t ihit = 0; ihit < n_hits; ++ihit) {
      // Check whether at least one hit satisfies mod == 5 and pln > 18.
      if (std::fabs(mod->at(ihit) - 5.0) < 1e-6 && pln->at(ihit) > 18.0) {
        file->Close();
        return true;
      }
    }
  }

  file->Close();
  return false;
}

void FindMod5PlnGT18Files(const char* input_dir) {
  TSystemDirectory dir("input_dir", input_dir);
  TList* files = dir.GetListOfFiles();

  if (!files) {
    std::cerr << "[ERROR] Failed to read directory: " << input_dir << std::endl;
    return;
  }

  std::vector<std::string> root_files;

  TIter next(files);
  TSystemFile* file = nullptr;

  while ((file = dynamic_cast<TSystemFile*>(next()))) {
    std::string filename = file->GetName();

    if (file->IsDirectory()) continue;
    if (filename.size() < 5) continue;
    if (filename.substr(filename.size() - 5) != ".root") continue;

    root_files.push_back(filename);
  }

  std::sort(root_files.begin(), root_files.end());

  std::cout << "[INFO] Files with at least one hit satisfying "
            << "mod == 5 && pln > 18:" << std::endl;

  int n_found = 0;

  for (const auto& filename : root_files) {
    std::string root_path = std::string(input_dir) + "/" + filename;

    if (HasMod5PlnGT18(root_path)) {
      std::cout << filename << std::endl;
      ++n_found;
    }
  }

  std::cout << "[INFO] Number of matched files: " << n_found << std::endl;
}
