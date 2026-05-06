// File: export_bad_detector_flag_timestamps.C
//
// Usage:
//   root -l -q 'export_bad_detector_flag_timestamps.C("/path/to/root_dir","timestamps.csv")'
//
// This macro scans ROOT files in the input directory, reads the "ntbm" tree,
// and writes timestamps for entries satisfying:
//   bsd_good_spill_flag_ != 0 && detector_flags_[5] == 0
//
// The timestamp is converted from Unix time to:
//   YYYY-MM-DD HH:MM:SS

#include <TROOT.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static std::vector<std::string> GetRootFiles(const char* input_dir)
{
  std::vector<std::string> files;

  TSystemDirectory dir("input_dir", input_dir);
  TList* list = dir.GetListOfFiles();

  if (!list) {
    std::cerr << "Error: cannot open input directory: "
              << input_dir << std::endl;
    return files;
  }

  TIter next(list);
  TSystemFile* file = nullptr;

  while ((file = static_cast<TSystemFile*>(next()))) {
    TString name = file->GetName();

    if (file->IsDirectory()) continue;
    if (!name.EndsWith(".root")) continue;

    TString full_path = TString(input_dir) + "/" + name;
    files.push_back(full_path.Data());
  }

  std::sort(files.begin(), files.end());
  return files;
}

static std::string UnixTimeToDateTimeString(std::time_t tt)
{
  std::tm lt{};
#if defined(_WIN32)
  localtime_s(&lt, &tt);
#else
  lt = *std::localtime(&tt);
#endif

  std::ostringstream oss;
  oss << std::setw(4) << std::setfill('0') << lt.tm_year + 1900
      << "-"
      << std::setw(2) << std::setfill('0') << lt.tm_mon + 1
      << "-"
      << std::setw(2) << std::setfill('0') << lt.tm_mday
      << " "
      << std::setw(2) << std::setfill('0') << lt.tm_hour
      << ":"
      << std::setw(2) << std::setfill('0') << lt.tm_min
      << ":"
      << std::setw(2) << std::setfill('0') << lt.tm_sec;

  return oss.str();
}

void export_bad_detector_flag_timestamps(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/latest20260504/rootfile_after_TrackMatch/RHC",
                                          const char* output_csv = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bad_detector_flag_timestamps_RHC.csv")
// void export_bad_detector_flag_timestamps(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/latest20260504/rootfile_after_TrackMatch/FHC",
//                                           const char* output_csv = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/BMcheck/bad_detector_flag_timestamps_FHC.csv")
{
  gROOT->SetBatch(kTRUE);

  std::vector<std::string> root_files = GetRootFiles(input_dir);

  if (root_files.empty()) {
    std::cerr << "Error: no ROOT files found in directory: "
              << input_dir << std::endl;
    return;
  }

  std::ofstream ofs(output_csv);

  if (!ofs.is_open()) {
    std::cerr << "Error: cannot open output CSV file: "
              << output_csv << std::endl;
    return;
  }

  TChain chain("ntbm");

  for (const auto& path : root_files) {
    chain.Add(path.c_str());
  }

  std::cout << "Input directory : " << input_dir << std::endl;
  std::cout << "Output CSV      : " << output_csv << std::endl;
  std::cout << "Number of files : " << root_files.size() << std::endl;
  std::cout << "Number of entries in chain : "
            << chain.GetEntries() << std::endl;

  TTreeReader reader(&chain);

  TTreeReaderValue<Double_t> timestamp(reader, "timestamp_");
  TTreeReaderValue<Int_t> bsd_good_spill_flag(reader, "bsd_good_spill_flag_");

  // This fixed-size array branch appears as detector_flags_[8] in TTree::Print().
  // If your file exposes it as detector_flags_ instead, change the branch name below.
  TTreeReaderArray<Int_t> detector_flags(reader, "detector_flags_[8]");

  ofs << "timestamp\n";

  Long64_t n_entries = 0;
  Long64_t n_selected = 0;

  while (reader.Next()) {
    ++n_entries;

    if (!std::isfinite(*timestamp) || *timestamp <= 0.0) continue;
    if (*bsd_good_spill_flag == 0) continue;
    if (detector_flags.GetSize() <= 5) continue;
    if (detector_flags[5] != 0) continue;

    const std::time_t tt =
      static_cast<std::time_t>(std::llround(*timestamp));

    ofs << UnixTimeToDateTimeString(tt) << "\n";

    ++n_selected;
  }

  ofs.close();

  std::cout << "Scanned entries   : " << n_entries << std::endl;
  std::cout << "Selected entries  : " << n_selected << std::endl;
  std::cout << "Saved CSV to      : " << output_csv << std::endl;
}
