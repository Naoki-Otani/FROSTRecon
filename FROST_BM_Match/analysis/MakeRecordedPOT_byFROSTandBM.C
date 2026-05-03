// MakeRecordedPOT_byFROSTandBM.C
//
// Usage:
// root -l -b -q 'MakeRecordedPOT_byFROSTandBM.C("/path/to/match_info_dir","/path/to/bsd_dir","/path/to/output_dir")'

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

struct PeriodInfo {
  std::string name;
  time_t start;
  time_t end;
};

struct BSDEntry {
  int spillnum_mod;
  int trg_sec2;
  double bunch12_pot;
};

time_t MakeJSTTime(int year, int month, int day, int hour, int min, int sec = 0) {
  std::tm tm = {};
  tm.tm_year = year - 1900;
  tm.tm_mon  = month - 1;
  tm.tm_mday = day;
  tm.tm_hour = hour - 9;  // JST -> UTC
  tm.tm_min  = min;
  tm.tm_sec  = sec;
  return timegm(&tm);
}

std::string FormatJST(time_t unix_time) {
  time_t jst_time = unix_time + 9 * 3600;
  std::tm *tm = gmtime(&jst_time);

  char buf[32];
  std::snprintf(buf, sizeof(buf), "%04d/%02d/%02d/%02d:%02d",
                tm->tm_year + 1900,
                tm->tm_mon + 1,
                tm->tm_mday,
                tm->tm_hour,
                tm->tm_min);
  return std::string(buf);
}

std::vector<std::string> GetRootFiles(const std::string &dir_path) {
  std::vector<std::string> files;

  void *dirp = gSystem->OpenDirectory(dir_path.c_str());
  if (!dirp) {
    std::cerr << "Error: cannot open directory: " << dir_path << std::endl;
    return files;
  }

  const char *entry;
  while ((entry = gSystem->GetDirEntry(dirp))) {
    std::string name = entry;
    if (name == "." || name == "..") continue;
    if (name.size() >= 5 && name.substr(name.size() - 5) == ".root") {
      files.push_back(dir_path + "/" + name);
    }
  }

  gSystem->FreeDirectory(dirp);
  std::sort(files.begin(), files.end());
  return files;
}

int PositiveMod32768(int x) {
  int m = x % 32768;
  if (m < 0) m += 32768;
  return m;
}

std::vector<BSDEntry> LoadBSDEntries(const std::string &bsd_dir) {
  std::vector<BSDEntry> entries;
  std::vector<std::string> root_files = GetRootFiles(bsd_dir);

  if (root_files.empty()) {
    std::cerr << "Warning: no BSD ROOT files found in " << bsd_dir << std::endl;
    return entries;
  }

  for (const auto &file_path : root_files) {
    std::cout << "Loading BSD: " << file_path << std::endl;

    TFile *file = TFile::Open(file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Warning: failed to open BSD file: " << file_path << std::endl;
      if (file) file->Close();
      delete file;
      continue;
    }

    TTree *tree = dynamic_cast<TTree *>(file->Get("bsd"));
    if (!tree) {
      std::cerr << "Warning: bsd tree not found in " << file_path << std::endl;
      file->Close();
      delete file;
      continue;
    }

    Int_t spillnum = 0;
    Int_t trg_sec[3] = {0, 0, 0};
    Double_t ct_np[5][9];

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("spillnum", 1);
    tree->SetBranchStatus("trg_sec", 1);
    tree->SetBranchStatus("ct_np", 1);

    tree->SetBranchAddress("spillnum", &spillnum);
    tree->SetBranchAddress("trg_sec", trg_sec);
    tree->SetBranchAddress("ct_np", ct_np);

    const Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
      tree->GetEntry(i);

      BSDEntry e;
      e.spillnum_mod = PositiveMod32768(spillnum);
      e.trg_sec2 = trg_sec[2];
      e.bunch12_pot = ct_np[4][1] + ct_np[4][2];
      entries.push_back(e);
    }

    file->Close();
    delete file;
  }

  std::sort(entries.begin(), entries.end(),
            [](const BSDEntry &a, const BSDEntry &b) {
              if (a.spillnum_mod != b.spillnum_mod) return a.spillnum_mod < b.spillnum_mod;
              return a.trg_sec2 < b.trg_sec2;
            });

  std::cout << "Loaded " << entries.size() << " BSD entries." << std::endl;
  return entries;
}

bool FindBSDPOT(const std::vector<BSDEntry> &bsd_entries,
                int wagasci_spill,
                int wagasci_unixtime,
                double &pot_out) {
  pot_out = 0.0;

  const int target_spill_mod = PositiveMod32768(wagasci_spill);

  auto lower = std::lower_bound(
      bsd_entries.begin(), bsd_entries.end(), target_spill_mod,
      [](const BSDEntry &e, int value) { return e.spillnum_mod < value; });

  auto upper = std::upper_bound(
      bsd_entries.begin(), bsd_entries.end(), target_spill_mod,
      [](int value, const BSDEntry &e) { return value < e.spillnum_mod; });

  if (lower == upper) return false;

  int best_dt = std::numeric_limits<int>::max();
  bool found = false;

  for (auto it = lower; it != upper; ++it) {
    int dt = std::abs(it->trg_sec2 - wagasci_unixtime);
    if (dt <= 3000 && dt < best_dt) {
      best_dt = dt;
      pot_out = it->bunch12_pot;
      found = true;
    }
  }

  return found;
}

void WriteCSV(const std::string &out_path,
              const PeriodInfo &period,
              const std::vector<double> &bin_sums) {
  std::ofstream ofs(out_path);
  if (!ofs) {
    std::cerr << "Error: cannot open output file: " << out_path << std::endl;
    return;
  }

  ofs << "#time, recorded POT by FROST and BM\n";

  double cumulative = 0.0;
  const int n_bins = static_cast<int>(bin_sums.size());

  for (int i = 0; i <= n_bins; ++i) {
    time_t t = period.start + static_cast<time_t>(i) * 600;
    if (i > 0) cumulative += bin_sums[i - 1];

    ofs << FormatJST(t) << ",";
    ofs << std::setprecision(16) << cumulative << "\n";
  }

  ofs.close();
}

void MakeRecordedPOT_byFROSTandBM(const char *match_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/rootfile_after_TrackMatch_new/zshift/zshift_0",
                                  const char *bsd_dir = "/group/nu/ninja/work/otani/FROST_beamdata/e71c4/bsd",
                                  const char *output_dir = "/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/pot") {
  std::string in_match_dir = match_dir;
  std::string in_bsd_dir   = bsd_dir;
  std::string out_dir      = output_dir;

  gSystem->mkdir(out_dir.c_str(), true);

  std::vector<PeriodInfo> periods = {
    {"RecordedPOT_byFROSTandBM_1.csv",
      MakeJSTTime(2025, 11, 29, 17, 10), MakeJSTTime(2025, 12, 22,  9, 10)},
    {"RecordedPOT_byFROSTandBM_2.csv",
      MakeJSTTime(2026,  1, 17, 19, 40), MakeJSTTime(2026,  1, 28,  5, 40)},
    {"RecordedPOT_byFROSTandBM_3.csv",
      MakeJSTTime(2026,  1, 28, 18,  0), MakeJSTTime(2026,  3, 11,  5, 40)},
    {"RecordedPOT_byFROSTandBM_4.csv",
      MakeJSTTime(2026,  3, 11, 18, 40), MakeJSTTime(2026,  3, 25,  6, 30)}
  };

  std::vector<std::vector<double>> period_bin_sums;
  for (const auto &p : periods) {
    int n_bins = static_cast<int>((p.end - p.start) / 600);
    if ((p.end - p.start) % 600 != 0) {
      std::cerr << "Error: period range is not divisible by 10 minutes: "
                << p.name << std::endl;
      return;
    }
    period_bin_sums.emplace_back(n_bins, 0.0);
  }

  std::vector<BSDEntry> bsd_entries = LoadBSDEntries(in_bsd_dir);

  std::vector<std::string> match_files = GetRootFiles(in_match_dir);
  if (match_files.empty()) {
    std::cerr << "Error: no match_info ROOT files found in " << in_match_dir << std::endl;
    return;
  }

  Long64_t n_run_le_5_total = 0;
  Long64_t n_run_le_5_bsd_found = 0;
  Long64_t n_run_le_5_bsd_missing = 0;

  for (const auto &file_path : match_files) {
    std::cout << "Processing match_info: " << file_path << std::endl;

    TFile *file = TFile::Open(file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Warning: failed to open " << file_path << std::endl;
      if (file) file->Close();
      delete file;
      continue;
    }

    TTree *tree = dynamic_cast<TTree *>(file->Get("match_info"));
    if (!tree) {
      std::cerr << "Warning: match_info tree not found in " << file_path << std::endl;
      file->Close();
      delete file;
      continue;
    }

    Int_t matched = 0;
    Int_t bsd_good_spill_flag = 0;
    Int_t frost_run_number = 0;
    Int_t wagasci_spill = 0;
    Int_t wagasci_unixtime = 0;
    Double_t spill_pot = 0.0;
    Double_t timestamp = 0.0;
    Int_t detector_flags[8] = {0};

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("matched", 1);
    tree->SetBranchStatus("bsd_good_spill_flag", 1);
    tree->SetBranchStatus("frost_run_number", 1);
    tree->SetBranchStatus("wagasci_spill", 1);
    tree->SetBranchStatus("wagasci_unixtime", 1);
    tree->SetBranchStatus("spill_pot", 1);
    tree->SetBranchStatus("timestamp", 1);
    tree->SetBranchStatus("detector_flags", 1);
    tree->SetBranchAddress("matched", &matched);
    tree->SetBranchAddress("bsd_good_spill_flag", &bsd_good_spill_flag);
    tree->SetBranchAddress("frost_run_number", &frost_run_number);
    tree->SetBranchAddress("wagasci_spill", &wagasci_spill);
    tree->SetBranchAddress("wagasci_unixtime", &wagasci_unixtime);
    tree->SetBranchAddress("spill_pot", &spill_pot);
    tree->SetBranchAddress("timestamp", &timestamp);
    tree->SetBranchAddress("detector_flags", detector_flags);

    const Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
      tree->GetEntry(i);

      if (bsd_good_spill_flag == 0) continue;
      if (detector_flags[5] == 0) continue;
      if (matched != 1) continue;

      double pot_to_add = 0.0;

      if (frost_run_number <= 5) {
        ++n_run_le_5_total;
        double bsd_pot = 0.0;
        if (FindBSDPOT(bsd_entries, wagasci_spill, wagasci_unixtime, bsd_pot)) {
          pot_to_add = bsd_pot;
          ++n_run_le_5_bsd_found;
        } else {
          ++n_run_le_5_bsd_missing;
          continue;
        }
      } else {
        pot_to_add = spill_pot;
      }

      time_t ts = static_cast<time_t>(timestamp);

      for (size_t ip = 0; ip < periods.size(); ++ip) {
        const auto &p = periods[ip];
        if (ts < p.start || ts >= p.end) continue;

        int bin = static_cast<int>((ts - p.start) / 600);
        if (bin >= 0 && bin < static_cast<int>(period_bin_sums[ip].size())) {
          period_bin_sums[ip][bin] += pot_to_add;
        }
        break;
      }
    }

    file->Close();
    delete file;
  }

  for (size_t ip = 0; ip < periods.size(); ++ip) {
    const std::string out_path = out_dir + "/" + periods[ip].name;
    WriteCSV(out_path, periods[ip], period_bin_sums[ip]);
    std::cout << "Wrote: " << out_path << std::endl;
  }

  std::cout << "run<=5 entries needing BSD lookup: " << n_run_le_5_total << std::endl;
  std::cout << "run<=5 matched in BSD:           " << n_run_le_5_bsd_found << std::endl;
  std::cout << "run<=5 not found in BSD:         " << n_run_le_5_bsd_missing << std::endl;
  std::cout << "Done." << std::endl;
}
