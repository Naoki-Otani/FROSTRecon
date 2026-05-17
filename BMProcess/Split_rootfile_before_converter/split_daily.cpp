#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kTreeName = "tree";
constexpr int kSecondsPerDay = 24 * 60 * 60;

struct DayInfo {
  std::time_t start_unixtime = 0;
  std::string date_label;
};

struct SplitConfig {
  std::string label;
  std::string unixtime_branch_name;
};

void PrintUsage(const char* program_name) {
  std::cerr << "Usage:\n"
            << "  " << program_name << " -bm <input_dir> <output_dir>\n"
            << "  " << program_name << " -pm <input_dir> <output_dir>\n\n"
            << "Options:\n"
            << "  -bm   Split Baby MIND ROOT files by day using unixtime.\n"
            << "  -pm   Split PM ROOT files by day using BMBSD.unixtime.\n";
}

SplitConfig GetSplitConfig(const std::string& mode) {
  if (mode == "-bm") {
    return SplitConfig{"BMBSD", "unixtime"};
  }

  if (mode == "-pm") {
    return SplitConfig{"PMBSD", "BMBSD.unixtime"};
  }

  throw std::runtime_error("Unknown mode: " + mode);
}

bool HasRootExtension(const fs::path& path) {
  return path.extension() == ".root";
}

std::vector<fs::path> CollectRootFiles(const fs::path& input_dir) {
  if (!fs::exists(input_dir)) {
    throw std::runtime_error("Input directory does not exist: " + input_dir.string());
  }

  if (!fs::is_directory(input_dir)) {
    throw std::runtime_error("Input path is not a directory: " + input_dir.string());
  }

  std::vector<fs::path> root_files;

  for (const auto& entry : fs::directory_iterator(input_dir)) {
    if (!entry.is_regular_file()) continue;

    const fs::path& path = entry.path();
    if (!HasRootExtension(path)) continue;

    root_files.push_back(path);
  }

  std::sort(root_files.begin(), root_files.end());
  return root_files;
}

void EnsureOutputDirectory(const fs::path& output_dir) {
  if (fs::exists(output_dir)) {
    if (!fs::is_directory(output_dir)) {
      throw std::runtime_error("Output path exists but is not a directory: " +
                               output_dir.string());
    }
    return;
  }

  fs::create_directories(output_dir);
}

std::string FormatDateLabel(const std::tm& tm_value) {
  std::ostringstream oss;
  oss << std::put_time(&tm_value, "%Y-%m-%d");
  return oss.str();
}

DayInfo GetDayInfoFromUnixTime(int unixtime) {
  const std::time_t t = static_cast<std::time_t>(unixtime);

  std::tm local_tm{};
#if defined(_WIN32)
  localtime_s(&local_tm, &t);
#else
  localtime_r(&t, &local_tm);
#endif

  local_tm.tm_hour = 0;
  local_tm.tm_min = 0;
  local_tm.tm_sec = 0;

  const std::time_t day_start = std::mktime(&local_tm);

  DayInfo info;
  info.start_unixtime = day_start;
  info.date_label = FormatDateLabel(local_tm);

  return info;
}

std::string MakeOutputFileName(const std::string& label,
                               const std::string& date_label) {
  return label + "_" + date_label + "_00-00-00.root";
}

void AddFilesToChain(TChain& chain, const std::vector<fs::path>& root_files) {
  for (const auto& path : root_files) {
    const int added = chain.Add(path.string().c_str());
    if (added == 0) {
      throw std::runtime_error("Failed to add ROOT file to chain: " + path.string());
    }
  }
}

std::map<std::time_t, DayInfo> CollectDays(
    TChain& chain,
    const std::string& unixtime_branch_name) {
  TTreeReader reader(&chain);
  TTreeReaderValue<int> unixtime(reader, unixtime_branch_name.c_str());

  std::map<std::time_t, DayInfo> days;

  Long64_t n_invalid_unixtime = 0;
  Long64_t n_read_entries = 0;

  while (reader.Next()) {
    ++n_read_entries;

    // Ignore invalid unix time values.
    if (*unixtime <= 0) {
      ++n_invalid_unixtime;
      continue;
    }

    const DayInfo info = GetDayInfoFromUnixTime(*unixtime);
    days.emplace(info.start_unixtime, info);
  }

  std::cout << "[INFO] Read entries             : "
            << n_read_entries << '\n';
  std::cout << "[INFO] Invalid unixtime entries : "
            << n_invalid_unixtime << '\n';

  return days;
}

void SplitByDay(const fs::path& input_dir,
                const fs::path& output_dir,
                const SplitConfig& config) {
  EnsureOutputDirectory(output_dir);

  const std::vector<fs::path> root_files = CollectRootFiles(input_dir);
  if (root_files.empty()) {
    throw std::runtime_error("No ROOT files found in input directory: " +
                             input_dir.string());
  }

  TChain chain(kTreeName);
  AddFilesToChain(chain, root_files);

  if (!chain.GetBranch(config.unixtime_branch_name.c_str())) {
    throw std::runtime_error("Branch '" + config.unixtime_branch_name +
                             "' was not found in TTree 'tree'.");
  }

  const Long64_t n_entries = chain.GetEntries();
  std::cout << "[INFO] Mode                  : " << config.label << '\n';
  std::cout << "[INFO] Unixtime branch       : "
            << config.unixtime_branch_name << '\n';
  std::cout << "[INFO] Number of input files : " << root_files.size() << '\n';
  std::cout << "[INFO] Number of entries     : " << n_entries << '\n';

  const std::map<std::time_t, DayInfo> days =
      CollectDays(chain, config.unixtime_branch_name);
  std::cout << "[INFO] Number of output days : " << days.size() << '\n';

  for (const auto& [day_start, info] : days) {
    const std::time_t day_end = day_start + kSecondsPerDay;

    const fs::path output_path =
        output_dir / MakeOutputFileName(config.label, info.date_label);

    std::ostringstream selection;
    selection << config.unixtime_branch_name << " > 0"
              << " && " << config.unixtime_branch_name
              << " >= " << static_cast<long long>(day_start)
              << " && " << config.unixtime_branch_name
              << " < " << static_cast<long long>(day_end);

    std::cout << "[INFO] Writing " << output_path
              << " with selection: " << selection.str() << '\n';

    std::unique_ptr<TFile> output_file(
        TFile::Open(output_path.string().c_str(), "RECREATE"));

    if (!output_file || output_file->IsZombie()) {
      throw std::runtime_error("Failed to create output file: " +
                               output_path.string());
    }

    output_file->cd();

    TTree* output_tree = chain.CopyTree(selection.str().c_str());
    if (!output_tree) {
      throw std::runtime_error("Failed to copy tree for day: " + info.date_label);
    }

    output_tree->SetName(kTreeName);
    output_tree->Write();

    std::cout << "       entries: " << output_tree->GetEntries() << '\n';

    output_file->Close();
  }

  std::cout << "[INFO] Done.\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 4) {
    PrintUsage(argv[0]);
    return 1;
  }

  const std::string mode = argv[1];
  const fs::path input_dir = argv[2];
  const fs::path output_dir = argv[3];

  try {
    const SplitConfig config = GetSplitConfig(mode);
    SplitByDay(input_dir, output_dir, config);
  } catch (const std::exception& e) {
    std::cerr << "[ERROR] " << e.what() << '\n';
    return 1;
  }

  return 0;
}
