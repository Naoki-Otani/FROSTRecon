// make_frost_muon_kinematics_sample.C
//
// Create a ROOT file containing the kinematics of FROST-crossing single-muon
// events from WAGASCIMC output files.
//
// Selection:
//   - Use the "frostmc" tree.
//   - Require num_particles == 1.
//   - Require the only particle to be a muon, |particle_pdg_id[0]| == 13.
//
// Output tree:
//   frost_muon_kinematics
//
// Output branches:
//   sample_id
//   source_file
//   source_file_index
//   source_local_entry
//   pid
//   p_mev
//   tanx
//   tany
//
// Usage:
//   root -l -b -q 'make_frost_muon_kinematics_sample.C("input_dir", "output.root")'
//
// Recursive search:
//   root -l -b -q 'make_frost_muon_kinematics_sample.C("input_dir", "output.root", true)'

#include <TChain.h>
#include <TFile.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace {

bool EndsWithRoot(const std::string& filename)
{
  return filename.size() >= 5 &&
         filename.substr(filename.size() - 5) == ".root";
}

std::string JoinPath(const std::string& directory, const std::string& filename)
{
  if (directory.empty()) {
    return filename;
  }

  if (directory.back() == '/') {
    return directory + filename;
  }

  return directory + "/" + filename;
}

void CollectRootFiles(const std::string& directory,
                      bool recursive,
                      std::vector<std::string>& root_files)
{
  void* dir = gSystem->OpenDirectory(directory.c_str());
  if (dir == nullptr) {
    std::cerr << "Error: Could not open directory: "
              << directory << std::endl;
    return;
  }

  const char* entry_name = nullptr;
  while ((entry_name = gSystem->GetDirEntry(dir)) != nullptr) {
    const std::string name(entry_name);

    if (name == "." || name == "..") {
      continue;
    }

    const std::string path = JoinPath(directory, name);

    FileStat_t file_stat;
    if (gSystem->GetPathInfo(path.c_str(), file_stat) != 0) {
      continue;
    }

    if (R_ISDIR(file_stat.fMode)) {
      if (recursive) {
        CollectRootFiles(path, recursive, root_files);
      }
      continue;
    }

    if (EndsWithRoot(name)) {
      root_files.push_back(path);
    }
  }

  gSystem->FreeDirectory(dir);
}

std::vector<std::string> GetSortedRootFiles(const std::string& input_directory,
                                            bool recursive)
{
  std::vector<std::string> root_files;
  CollectRootFiles(input_directory, recursive, root_files);

  std::sort(root_files.begin(), root_files.end());
  return root_files;
}

template <typename T>
bool CheckBranch(TTree* tree, const char* branch_name)
{
  if (tree->GetBranch(branch_name) == nullptr) {
    std::cerr << "Error: Required branch is missing: "
              << branch_name << std::endl;
    return false;
  }

  return true;
}

bool CheckInputBranches(TTree* tree)
{
  if (tree == nullptr) {
    std::cerr << "Error: Input tree is null." << std::endl;
    return false;
  }

  bool ok = true;
  ok &= CheckBranch<int>(tree, "num_particles");
  ok &= CheckBranch<int>(tree, "particle_pdg_id");
  ok &= CheckBranch<double>(tree, "particle_initial_px_mev_c");
  ok &= CheckBranch<double>(tree, "particle_initial_py_mev_c");
  ok &= CheckBranch<double>(tree, "particle_initial_pz_mev_c");

  return ok;
}

}  // namespace

void make_frost_muon_kinematics_sample(
    const char* input_directory,
    const char* output_filename = "frost_muon_kinematics.root",
    bool recursive = false,
    Long64_t max_entries = -1)
{
  const std::vector<std::string> input_files =
      GetSortedRootFiles(input_directory, recursive);

  if (input_files.empty()) {
    std::cerr << "Error: No ROOT files found in directory: "
              << input_directory << std::endl;
    return;
  }

  std::cout << "Input directory: " << input_directory << std::endl;
  std::cout << "Number of input ROOT files: "
            << input_files.size() << std::endl;

  TChain chain("frostmc");

  for (const auto& input_file : input_files) {
    chain.Add(input_file.c_str());
  }

  if (chain.GetEntries() <= 0) {
    std::cerr << "Error: No entries found in the frostmc TChain."
              << std::endl;
    return;
  }

  if (!CheckInputBranches(&chain)) {
    return;
  }

  std::map<std::string, int> file_index_map;
  for (std::size_t i = 0; i < input_files.size(); ++i) {
    file_index_map[input_files[i]] = static_cast<int>(i);
  }

  Int_t num_particles = 0;
  std::vector<Int_t>* particle_pdg_id = nullptr;
  std::vector<Double_t>* particle_initial_px_mev_c = nullptr;
  std::vector<Double_t>* particle_initial_py_mev_c = nullptr;
  std::vector<Double_t>* particle_initial_pz_mev_c = nullptr;

  chain.SetBranchAddress("num_particles", &num_particles);
  chain.SetBranchAddress("particle_pdg_id", &particle_pdg_id);
  chain.SetBranchAddress("particle_initial_px_mev_c",
                         &particle_initial_px_mev_c);
  chain.SetBranchAddress("particle_initial_py_mev_c",
                         &particle_initial_py_mev_c);
  chain.SetBranchAddress("particle_initial_pz_mev_c",
                         &particle_initial_pz_mev_c);

  TFile output_file(output_filename, "RECREATE");
  if (output_file.IsZombie()) {
    std::cerr << "Error: Could not create output file: "
              << output_filename << std::endl;
    return;
  }

  TTree output_tree("frost_muon_kinematics",
                    "FROST-crossing single-muon kinematics");

  Int_t sample_id = 0;
  std::string source_file;
  Int_t source_file_index = -1;
  Long64_t source_local_entry = -1;

  Int_t pid = 0;
  Double_t p_mev = 0.0;
  Double_t tanx = 0.0;
  Double_t tany = 0.0;

  output_tree.Branch("sample_id", &sample_id, "sample_id/I");
  output_tree.Branch("source_file", &source_file);
  output_tree.Branch("source_file_index", &source_file_index,
                     "source_file_index/I");
  output_tree.Branch("source_local_entry", &source_local_entry,
                     "source_local_entry/L");

  output_tree.Branch("pid", &pid, "pid/I");
  output_tree.Branch("p_mev", &p_mev, "p_mev/D");
  output_tree.Branch("tanx", &tanx, "tanx/D");
  output_tree.Branch("tany", &tany, "tany/D");

  const Long64_t total_entries = chain.GetEntries();
  const Long64_t entries_to_process =
      max_entries > 0 ? std::min(max_entries, total_entries) : total_entries;

  Long64_t selected_entries = 0;
  Long64_t skipped_wrong_num_particles = 0;
  Long64_t skipped_non_muon = 0;
  Long64_t skipped_bad_vectors = 0;
  Long64_t skipped_zero_pz = 0;

  for (Long64_t entry = 0; entry < entries_to_process; ++entry) {
    const Long64_t local_entry = chain.LoadTree(entry);
    if (local_entry < 0) {
      continue;
    }

    chain.GetEntry(entry);

    if (entry % 100000 == 0) {
      std::cout << "Processing entry " << entry
                << " / " << entries_to_process << std::endl;
    }

    if (num_particles != 1) {
      ++skipped_wrong_num_particles;
      continue;
    }

    if (particle_pdg_id == nullptr ||
        particle_initial_px_mev_c == nullptr ||
        particle_initial_py_mev_c == nullptr ||
        particle_initial_pz_mev_c == nullptr) {
      ++skipped_bad_vectors;
      continue;
    }

    if (particle_pdg_id->size() != 1 ||
        particle_initial_px_mev_c->size() != 1 ||
        particle_initial_py_mev_c->size() != 1 ||
        particle_initial_pz_mev_c->size() != 1) {
      ++skipped_bad_vectors;
      continue;
    }

    const Int_t input_pid = particle_pdg_id->at(0);
    if (std::abs(input_pid) != 13) {
      ++skipped_non_muon;
      continue;
    }

    const Double_t px = particle_initial_px_mev_c->at(0);
    const Double_t py = particle_initial_py_mev_c->at(0);
    const Double_t pz = particle_initial_pz_mev_c->at(0);

    if (std::abs(pz) < std::numeric_limits<Double_t>::epsilon()) {
      ++skipped_zero_pz;
      continue;
    }

    TFile* current_file = chain.GetFile();

    source_file = current_file != nullptr ? current_file->GetName() : "";
    source_local_entry = local_entry;

    const auto file_index_it = file_index_map.find(source_file);
    source_file_index =
        file_index_it != file_index_map.end() ? file_index_it->second : -1;

    sample_id = static_cast<Int_t>(selected_entries);

    pid = input_pid;
    p_mev = std::sqrt(px * px + py * py + pz * pz);
    tanx = px / pz;
    tany = py / pz;

    output_tree.Fill();

    ++selected_entries;
  }

  output_file.cd();
  output_tree.Write();
  output_file.Close();

  std::cout << "\nSummary\n"
            << "  Processed entries: " << entries_to_process << "\n"
            << "  Selected single-muon entries: " << selected_entries << "\n"
            << "  Skipped by num_particles != 1: "
            << skipped_wrong_num_particles << "\n"
            << "  Skipped by non-muon PID: " << skipped_non_muon << "\n"
            << "  Skipped by malformed particle vectors: "
            << skipped_bad_vectors << "\n"
            << "  Skipped by zero pz: " << skipped_zero_pz << "\n"
            << "  Output file: " << output_filename << std::endl;
}
