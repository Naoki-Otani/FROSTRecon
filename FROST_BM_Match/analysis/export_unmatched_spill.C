// export_unmatched_spill.C
// Usage example:
// root -l -b -q 'export_unmatched_spill.C("/path/to/root_dir","/path/to/output.csv")'

#include <TFile.h>
#include <TTree.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>

#include <algorithm>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct CsvRow {
  int wagasci_unixtime;
  int wagasci_spill;
};

std::string FormatUnixTime(int unixtime) {
  std::time_t t = static_cast<std::time_t>(unixtime);
  std::tm *lt = std::localtime(&t);
  if (!lt) return "";

  std::ostringstream oss;
  oss << std::setfill('0')
      << (lt->tm_year + 1900) << "/"
      << std::setw(2) << (lt->tm_mon + 1) << "/"
      << std::setw(2) << lt->tm_mday << "/"
      << std::setw(2) << lt->tm_hour << ":"
      << std::setw(2) << lt->tm_min;
  return oss.str();
}

void export_unmatched_spill(const char* input_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/rootfile_after_TrackMatch/zshift/zshift_0", const char* output_csv = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/unmatched_spill.csv") {
  TSystemDirectory dir("input_dir", input_dir);
  TList *fileList = dir.GetListOfFiles();

  if (!fileList) {
    std::cerr << "Error: cannot open input directory: " << input_dir << std::endl;
    return;
  }

  std::vector<std::string> rootFiles;

  TIter next(fileList);
  TSystemFile *file;
  while ((file = (TSystemFile*)next())) {
    std::string fname = file->GetName();
    if (file->IsDirectory()) continue;
    if (fname.size() >= 5 && fname.substr(fname.size() - 5) == ".root") {
      rootFiles.push_back(std::string(input_dir) + "/" + fname);
    }
  }

  if (rootFiles.empty()) {
    std::cerr << "Error: no ROOT files found in " << input_dir << std::endl;
    return;
  }

  std::vector<CsvRow> rows;
  Long64_t totalSelected = 0;

  for (const auto& filepath : rootFiles) {
    std::cout << "Processing: " << filepath << std::endl;

    TFile *fin = TFile::Open(filepath.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
      std::cerr << "Warning: failed to open " << filepath << std::endl;
      if (fin) {
        fin->Close();
        delete fin;
      }
      continue;
    }

    TTree *tree = (TTree*)fin->Get("match_info");
    if (!tree) {
      std::cerr << "Warning: match_info tree not found in " << filepath << std::endl;
      fin->Close();
      delete fin;
      continue;
    }

    int bsd_good_spill_flag = 0;
    int matched = 0;
    int wagasci_unixtime = 0;
    int wagasci_spill = 0;

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("bsd_good_spill_flag", 1);
    tree->SetBranchStatus("matched", 1);
    tree->SetBranchStatus("wagasci_unixtime", 1);
    tree->SetBranchStatus("wagasci_spill", 1);

    tree->SetBranchAddress("bsd_good_spill_flag", &bsd_good_spill_flag);
    tree->SetBranchAddress("matched", &matched);
    tree->SetBranchAddress("wagasci_unixtime", &wagasci_unixtime);
    tree->SetBranchAddress("wagasci_spill", &wagasci_spill);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);

      if (bsd_good_spill_flag != 0 && matched == 0) {
        CsvRow row;
        row.wagasci_unixtime = wagasci_unixtime;
        row.wagasci_spill = wagasci_spill;
        rows.push_back(row);
        ++totalSelected;
      }
    }

    fin->Close();
    delete fin;
  }

  std::sort(rows.begin(), rows.end(),
            [](const CsvRow& a, const CsvRow& b) {
              if (a.wagasci_unixtime != b.wagasci_unixtime) {
                return a.wagasci_unixtime < b.wagasci_unixtime;
              }
              return a.wagasci_spill < b.wagasci_spill;
            });

  std::ofstream ofs(output_csv);
  if (!ofs.is_open()) {
    std::cerr << "Error: cannot open output CSV: " << output_csv << std::endl;
    return;
  }

  ofs << "#time,wagasci_unixtime,wagasci_spillnum\n";

  for (const auto& row : rows) {
    ofs << FormatUnixTime(row.wagasci_unixtime) << ","
        << row.wagasci_unixtime << ","
        << row.wagasci_spill << "\n";
  }

  ofs.close();

  std::cout << "Done." << std::endl;
  std::cout << "Selected entries written to: " << output_csv << std::endl;
  std::cout << "Total selected entries: " << totalSelected << std::endl;
}
