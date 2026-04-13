#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TEntryList.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {

// Matching window in unixtime.
constexpr int kTimeWindow = 3000;

// Spill number is wrapped in [0, 32767].
constexpr int kSpillModulo = 32768;

struct PMEvent {
  Long64_t global_entry = -1;
  int spill = -1;
  int unixtime = -1;
  bool used = false;
};

struct PMIndexEntry {
  int unixtime = -1;
  std::size_t pm_event_index = 0;
};

bool HasRootExtension(const fs::path& path) {
  return path.extension() == ".root";
}

bool StartsWith(const std::string& s, const std::string& prefix) {
  return s.rfind(prefix, 0) == 0;
}

int NormalizeSpill(int spill) {
  int x = spill % kSpillModulo;
  if (x < 0) x += kSpillModulo;
  return x;
}

std::vector<fs::path> CollectRootFiles(const fs::path& dir,
                                       const std::string& prefix) {
  std::vector<fs::path> files;

  if (!fs::exists(dir) || !fs::is_directory(dir)) {
    throw std::runtime_error("Directory does not exist: " + dir.string());
  }

  for (const auto& entry : fs::directory_iterator(dir)) {
    if (!entry.is_regular_file()) continue;

    const fs::path& path = entry.path();
    const std::string name = path.filename().string();

    if (!HasRootExtension(path)) continue;
    if (!StartsWith(name, prefix)) continue;

    files.push_back(path);
  }

  std::sort(files.begin(), files.end());
  return files;
}

std::string MakeOutputFileNameFromBM(const fs::path& bm_path) {
  const std::string bm_name = bm_path.filename().string();
  const std::string bm_prefix = "BMBSD_";

  if (!StartsWith(bm_name, bm_prefix)) {
    throw std::runtime_error("Unexpected Baby MIND file name: " + bm_name);
  }

  return "PMBSD_" + bm_name.substr(bm_prefix.size());
}

void EnsureOutputDirectory(const fs::path& out_dir) {
  if (fs::exists(out_dir)) {
    if (!fs::is_directory(out_dir)) {
      throw std::runtime_error("Output path exists but is not a directory: " +
                               out_dir.string());
    }
    return;
  }

  fs::create_directories(out_dir);
}

TTree* GetTreeOrThrow(TFile& file, const std::string& file_label) {
  TTree* tree = nullptr;
  file.GetObject("tree", tree);
  if (!tree) {
    throw std::runtime_error("Could not find TTree named 'tree' in " + file_label);
  }
  return tree;
}

void CheckBranchExistsOrThrow(TTree* tree,
                              const std::string& branch_name,
                              const std::string& file_label) {
  if (!tree->GetBranch(branch_name.c_str())) {
    throw std::runtime_error("Missing branch '" + branch_name + "' in " + file_label);
  }
}

std::vector<PMEvent> BuildPMEventTable(
    TChain& pm_chain,
    std::array<std::vector<PMIndexEntry>, kSpillModulo>& pm_index_by_spill) {
  TTreeReader reader(&pm_chain);
  TTreeReaderValue<int> spill(reader, "spill");
  TTreeReaderValue<int> unixtime(reader, "unixtime");

  std::vector<PMEvent> events;
  Long64_t global_entry = 0;

  while (reader.Next()) {
    PMEvent ev;
    ev.global_entry = global_entry;
    ev.spill = NormalizeSpill(*spill);
    ev.unixtime = *unixtime;
    ev.used = false;

    if (global_entry < 10) {
      std::cout << "[PM] entry=" << global_entry
                << " spill=" << *spill
                << " spill_mod=" << ev.spill
                << " unixtime=" << *unixtime
                << std::endl;
    }

    const std::size_t idx = events.size();
    events.push_back(ev);
    pm_index_by_spill[ev.spill].push_back({ev.unixtime, idx});

    ++global_entry;
  }

  for (auto& vec : pm_index_by_spill) {
    std::sort(vec.begin(), vec.end(),
              [](const PMIndexEntry& a, const PMIndexEntry& b) {
                if (a.unixtime != b.unixtime) return a.unixtime < b.unixtime;
                return a.pm_event_index < b.pm_event_index;
              });
  }

  return events;
}

std::size_t FindBestMatchingPMEvent(
    int bm_spill,
    int bm_unixtime,
    std::vector<PMEvent>& pm_events,
    const std::array<std::vector<PMIndexEntry>, kSpillModulo>& pm_index_by_spill) {
  const auto& candidates = pm_index_by_spill[bm_spill];
  if (candidates.empty()) {
    return std::numeric_limits<std::size_t>::max();
  }

  const int lower_time = bm_unixtime - kTimeWindow;
  const int upper_time = bm_unixtime + kTimeWindow;

  auto lower_it = std::lower_bound(
      candidates.begin(), candidates.end(), lower_time,
      [](const PMIndexEntry& entry, int time) {
        return entry.unixtime < time;
      });

  std::size_t best_idx = std::numeric_limits<std::size_t>::max();
  int best_abs_dt = std::numeric_limits<int>::max();
  Long64_t best_global_entry = std::numeric_limits<Long64_t>::max();

  for (auto it = lower_it; it != candidates.end() && it->unixtime <= upper_time; ++it) {
    PMEvent& pm_ev = pm_events[it->pm_event_index];
    if (pm_ev.used) continue;

    const int abs_dt = std::abs(pm_ev.unixtime - bm_unixtime);

    if (abs_dt < best_abs_dt ||
        (abs_dt == best_abs_dt && pm_ev.global_entry < best_global_entry)) {
      best_idx = it->pm_event_index;
      best_abs_dt = abs_dt;
      best_global_entry = pm_ev.global_entry;
    }
  }

  return best_idx;
}

std::unordered_map<std::string, std::vector<Long64_t>> MatchEventsToBMFiles(
    const std::vector<fs::path>& bm_files,
    std::vector<PMEvent>& pm_events,
    const std::array<std::vector<PMIndexEntry>, kSpillModulo>& pm_index_by_spill) {
  std::unordered_map<std::string, std::vector<Long64_t>> matched_entries_by_output;

  for (const auto& bm_path : bm_files) {
    const std::string output_name = MakeOutputFileNameFromBM(bm_path);

    TFile bm_file(bm_path.string().c_str(), "READ");
    if (bm_file.IsZombie()) {
      throw std::runtime_error("Failed to open Baby MIND file: " + bm_path.string());
    }

    TTree* bm_tree = GetTreeOrThrow(bm_file, bm_path.string());
    CheckBranchExistsOrThrow(bm_tree, "spillnum", bm_path.string());
    CheckBranchExistsOrThrow(bm_tree, "unixtime", bm_path.string());

    TTreeReader reader(bm_tree);
    TTreeReaderValue<int> bm_spillnum(reader, "spillnum");
    TTreeReaderValue<int> bm_unixtime(reader, "unixtime");
    auto& matched_entries = matched_entries_by_output[output_name];

    Long64_t i = 0;
    while (reader.Next()) {

      if (i < 10) {
        std::cout << "[BM] file=" << bm_path.filename().string()
                  << " entry=" << i
                  << " spillnum=" << *bm_spillnum
                  << " spill_mod=" << NormalizeSpill(*bm_spillnum)
                  << " unixtime=" << *bm_unixtime
                  << std::endl;
      }

      // Ignore Baby MIND events during beam-off periods.
      if (*bm_unixtime == -1) {
        ++i;
        continue;
      }

      const int bm_spill = NormalizeSpill(*bm_spillnum);

      const std::size_t pm_idx =
          FindBestMatchingPMEvent(bm_spill, *bm_unixtime, pm_events, pm_index_by_spill);

      if (pm_idx == std::numeric_limits<std::size_t>::max()) {
        ++i;
        continue;
      }

      pm_events[pm_idx].used = true;
      matched_entries.push_back(pm_events[pm_idx].global_entry);
      ++i;
    }

    std::sort(matched_entries.begin(), matched_entries.end());
  }

  return matched_entries_by_output;
}

void WriteMatchedPMFiles(
    TChain& pm_chain,
    const fs::path& out_dir,
    const std::unordered_map<std::string, std::vector<Long64_t>>& matched_entries_by_output) {

  for (const auto& [output_name, entry_list] : matched_entries_by_output) {
    if (entry_list.empty()) {
      std::cout << "[INFO] Skip empty output: " << output_name << std::endl;
      continue;
    }

    const fs::path out_path = out_dir / output_name;

    auto elist = std::make_unique<TEntryList>("elist", "elist");
    for (Long64_t entry : entry_list) {
      elist->Enter(entry, &pm_chain);
    }

    pm_chain.SetEntryList(elist.get());


    TFile out_file(out_path.string().c_str(), "RECREATE");
    if (out_file.IsZombie()) {
      throw std::runtime_error("Failed to create output file: " + out_path.string());
    }

    TTree* out_tree = pm_chain.CopyTree("");
    if (!out_tree) {
      pm_chain.SetEntryList(nullptr);
      throw std::runtime_error("Failed to copy PM tree for output: " + out_path.string());
    }

    out_file.cd();
    out_tree->Write();
    out_file.Close();
    pm_chain.SetEntryList(nullptr);

    std::cout << "[INFO] Wrote " << out_path << " with "
              << entry_list.size() << " matched PM events" << std::endl;
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr
        << "Usage: " << argv[0]
        << " <baby_mind_dir> <proton_module_dir> <output_dir>\n"
        << "\n"
        << "Example:\n"
        << "  " << argv[0]
        << " /path/to/babymind/2-BMBSD"
        << " /path/to/protonmodule/4-PMBSD"
        << " /path/to/reformatted_pm\n";
    return 1;
  }

  try {
    const fs::path bm_dir = argv[1];
    const fs::path pm_dir = argv[2];
    const fs::path out_dir = argv[3];

    EnsureOutputDirectory(out_dir);

    const auto bm_files = CollectRootFiles(bm_dir, "BMBSD_");
    const auto pm_files = CollectRootFiles(pm_dir, "PMBSD_");

    if (bm_files.empty()) {
      throw std::runtime_error("No Baby MIND ROOT files found in: " + bm_dir.string());
    }
    if (pm_files.empty()) {
      throw std::runtime_error("No Proton Module ROOT files found in: " + pm_dir.string());
    }

    std::cout << "[INFO] Found " << bm_files.size() << " Baby MIND files" << std::endl;
    std::cout << "[INFO] Found " << pm_files.size() << " Proton Module files" << std::endl;

    TChain pm_chain("tree");
    for (const auto& pm_path : pm_files) {
      pm_chain.Add(pm_path.string().c_str());
    }

    if (!pm_chain.GetBranch("spill") || !pm_chain.GetBranch("unixtime")) {
      throw std::runtime_error(
          "PM TChain does not contain required branches 'spill' and 'unixtime'.");
    }

    std::array<std::vector<PMIndexEntry>, kSpillModulo> pm_index_by_spill;
    std::vector<PMEvent> pm_events = BuildPMEventTable(pm_chain, pm_index_by_spill);

    std::cout << "[INFO] Indexed " << pm_events.size()
              << " Proton Module events" << std::endl;

    auto matched_entries_by_output =
        MatchEventsToBMFiles(bm_files, pm_events, pm_index_by_spill);

    std::size_t total_matched = 0;
    for (const auto& [_, entries] : matched_entries_by_output) {
      total_matched += entries.size();
    }

    std::cout << "[INFO] Total matched PM events: " << total_matched << std::endl;
    std::cout << "[INFO] Writing output files..." << std::endl;

    if (total_matched == 0) {
      std::cout << "[WARNING] No matched PM events were found. No output files will be written." << std::endl;
      return 0;
    }

    WriteMatchedPMFiles(pm_chain, out_dir, matched_entries_by_output);

    std::cout << "[INFO] Done." << std::endl;
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "[ERROR] " << e.what() << std::endl;
    return 1;
  }
}
