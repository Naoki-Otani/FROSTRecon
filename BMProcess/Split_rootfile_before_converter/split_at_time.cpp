#include <TChain.h>
#include <TEventList.h>
#include <TFile.h>
#include <TKey.h>
#include <TObject.h>
#include <TTree.h>

#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kTreeName = "tree";
constexpr const char* kWGDirectoryPrefix = "physics_run";

struct SplitConfig {
  std::string label;
  std::string unixtime_branch_name;
};

void PrintUsage(const char* program_name) {
  std::cerr << "Usage:\n"
            << "  " << program_name
            << " -bm <input.root> <split_time> <output_dir>\n"
            << "  " << program_name
            << " -pm <input.root> <split_time> <output_dir>\n"
            << "  " << program_name
            << " -wg <input_dir> <split_time> <output_dir>\n\n"
            << "Examples:\n"
            << "  " << program_name
            << " -bm BMBSD_2026-01-24_00-00-00.root "
            << "2026-01-24_12-00-00 ./output\n"
            << "  " << program_name
            << " -pm PM_2026-01-24_00-00-00.root "
            << "2026-01-24_12-00:00 ./output\n"
            << "  " << program_name
            << " -wg physics_run_2026-01-24_00-00-00_0 "
            << "2026-01-24_12-30-03 ./output\n\n"
            << "Options:\n"
            << "  -bm   Split BM file using unixtime.\n"
            << "  -pm   Split PM file using BMBSD.unixtime.\n"
            << "  -wg   Split WG directory using bsd.timestamp.\n";
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

std::string FormatDateTimeLabel(const std::tm& tm_value) {
  std::ostringstream oss;
  oss << std::put_time(&tm_value, "%Y-%m-%d_%H-%M-%S");
  return oss.str();
}

std::time_t ParseDateTimeToUnixTime(const std::string& datetime_text) {
  // Accept both:
  //   YYYY-MM-DD_HH-MM-SS
  //   YYYY-MM-DD_HH-MM:SS
  const std::regex pattern(
      R"(^([0-9]{4})-([0-9]{2})-([0-9]{2})_([0-9]{2})-([0-9]{2})[:-]([0-9]{2})$)");

  std::smatch match;
  if (!std::regex_match(datetime_text, match, pattern)) {
    throw std::runtime_error(
        "Invalid split_time format. Use YYYY-MM-DD_HH-MM-SS, "
        "for example 2026-01-24_12-00-00.");
  }

  std::tm tm_value{};
  tm_value.tm_year = std::stoi(match[1].str()) - 1900;
  tm_value.tm_mon = std::stoi(match[2].str()) - 1;
  tm_value.tm_mday = std::stoi(match[3].str());
  tm_value.tm_hour = std::stoi(match[4].str());
  tm_value.tm_min = std::stoi(match[5].str());
  tm_value.tm_sec = std::stoi(match[6].str());
  tm_value.tm_isdst = -1;

  const std::time_t unixtime = std::mktime(&tm_value);
  if (unixtime < 0) {
    throw std::runtime_error("Failed to convert split_time to unixtime.");
  }

  return unixtime;
}

std::string MakeOutputFileName(const std::string& label,
                               const std::string& datetime_label) {
  return label + "_" + datetime_label + ".root";
}

std::string ExtractStartDateTimeLabelFromInputName(const fs::path& input_path,
                                                   const std::string& label) {
  // Expected examples:
  //   BMBSD_2026-01-24_00-00-00.root
  //   PM_2026-01-24_00-00-00.root
  const std::string filename = input_path.filename().string();

  const std::regex pattern("^" + label +
                           R"(_([0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{2}-[0-9]{2}-[0-9]{2})\.root$)");

  std::smatch match;
  if (!std::regex_match(filename, match, pattern)) {
    throw std::runtime_error(
        "Input filename does not match expected format: " + label +
        "_YYYY-MM-DD_HH-MM-SS.root");
  }

  return match[1].str();
}

void WriteSelectedTree(TChain& chain,
                       const fs::path& output_path,
                       const std::string& selection) {
  std::unique_ptr<TFile> output_file(
      TFile::Open(output_path.string().c_str(), "RECREATE"));

  if (!output_file || output_file->IsZombie()) {
    throw std::runtime_error("Failed to create output file: " +
                             output_path.string());
  }

  output_file->cd();

  TTree* output_tree = chain.CopyTree(selection.c_str());
  if (!output_tree) {
    throw std::runtime_error("Failed to copy tree with selection: " + selection);
  }

  output_tree->SetName(kTreeName);
  output_tree->Write();

  std::cout << "[INFO] Wrote " << output_path
            << " entries: " << output_tree->GetEntries() << '\n';

  output_file->Close();
}

std::string MakeWGDirectoryName(const std::string& datetime_label) {
  return std::string(kWGDirectoryPrefix) + "_" + datetime_label + "_0";
}

std::string MakeWGFileName(const std::string& datetime_label, int dif_id) {
  std::ostringstream oss;
  oss << MakeWGDirectoryName(datetime_label)
      << "_ecal_dif_" << dif_id << "_tree.root";
  return oss.str();
}

std::string ExtractStartDateTimeLabelFromWGDirectoryName(
    const fs::path& input_dir) {
  // Expected example:
  //   physics_run_2026-01-24_00-00-00_0
  const std::string dirname = input_dir.filename().string();
  const std::regex pattern(
      R"(^physics_run_([0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{2}-[0-9]{2}-[0-9]{2})_0$)");

  std::smatch match;
  if (!std::regex_match(dirname, match, pattern)) {
    throw std::runtime_error(
        "WG input directory name does not match expected format: "
        "physics_run_YYYY-MM-DD_HH-MM-SS_0");
  }

  return match[1].str();
}

std::vector<int> CollectWGDifIds(const fs::path& input_dir,
                                 const std::string& start_datetime_label) {
  std::vector<int> dif_ids;

  for (int dif_id = 0; dif_id <= 7; ++dif_id) {
    const fs::path input_file =
        input_dir / MakeWGFileName(start_datetime_label, dif_id);
    if (fs::exists(input_file) && fs::is_regular_file(input_file)) {
      dif_ids.push_back(dif_id);
    }
  }

  if (dif_ids.empty()) {
    throw std::runtime_error(
        "No WG ROOT files were found. Expected files like " +
        MakeWGFileName(start_datetime_label, 0));
  }

  return dif_ids;
}

TEventList* MakeWGEventList(TFile& input_file,
                             const std::string& selection,
                             const std::string& list_name) {
  TTree* bsd_tree = nullptr;
  input_file.GetObject("bsd", bsd_tree);
  if (!bsd_tree) {
    throw std::runtime_error("Tree 'bsd' was not found in WG ROOT file.");
  }

  input_file.cd();
  const std::string draw_command = ">>" + list_name;
  bsd_tree->Draw(draw_command.c_str(), selection.c_str(), "goff");

  TEventList* event_list =
      dynamic_cast<TEventList*>(input_file.Get(list_name.c_str()));
  if (!event_list) {
    throw std::runtime_error("Failed to create event list for WG selection: " +
                             selection);
  }

  event_list->SetDirectory(nullptr);
  input_file.Remove(event_list);
  return event_list;
}

void WriteSelectedWGFile(const fs::path& input_path,
                         const fs::path& output_path,
                         const std::string& selection,
                         const std::string& list_name) {
  std::unique_ptr<TFile> input_file(
      TFile::Open(input_path.string().c_str(), "READ"));
  if (!input_file || input_file->IsZombie()) {
    throw std::runtime_error("Failed to open WG input file: " +
                             input_path.string());
  }

  std::unique_ptr<TEventList> event_list(
      MakeWGEventList(*input_file, selection, list_name));

  std::unique_ptr<TFile> output_file(
      TFile::Open(output_path.string().c_str(), "RECREATE"));
  if (!output_file || output_file->IsZombie()) {
    throw std::runtime_error("Failed to create WG output file: " +
                             output_path.string());
  }

  TIter next_key(input_file->GetListOfKeys());
  TKey* key = nullptr;
  Long64_t bsd_entries = -1;

  while ((key = dynamic_cast<TKey*>(next_key()))) {
    std::unique_ptr<TObject> object(key->ReadObj());
    TTree* input_tree = dynamic_cast<TTree*>(object.get());
    if (!input_tree) continue;

    output_file->cd();
    input_tree->SetEventList(event_list.get());
    TTree* output_tree = input_tree->CopyTree("");
    input_tree->SetEventList(nullptr);

    if (!output_tree) {
      throw std::runtime_error("Failed to copy WG tree: " +
                               std::string(input_tree->GetName()));
    }

    output_tree->SetName(input_tree->GetName());
    output_tree->Write();

    if (std::string(input_tree->GetName()) == "bsd") {
      bsd_entries = output_tree->GetEntries();
    }
  }

  std::cout << "[INFO] Wrote " << output_path
            << " bsd entries: " << bsd_entries << '\n';

  output_file->Close();
  input_file->Close();
}

void SplitWGAtTime(const fs::path& input_dir,
                   const std::string& split_time_text,
                   const fs::path& output_dir) {
  if (!fs::exists(input_dir)) {
    throw std::runtime_error("WG input directory does not exist: " +
                             input_dir.string());
  }

  if (!fs::is_directory(input_dir)) {
    throw std::runtime_error("WG input path is not a directory: " +
                             input_dir.string());
  }

  EnsureOutputDirectory(output_dir);

  const std::time_t split_unixtime = ParseDateTimeToUnixTime(split_time_text);

  std::tm split_tm{};
#if defined(_WIN32)
  localtime_s(&split_tm, &split_unixtime);
#else
  localtime_r(&split_unixtime, &split_tm);
#endif

  const std::string split_datetime_label = FormatDateTimeLabel(split_tm);
  const std::string start_datetime_label =
      ExtractStartDateTimeLabelFromWGDirectoryName(input_dir);

  const fs::path output_before_dir =
      output_dir / MakeWGDirectoryName(start_datetime_label);
  const fs::path output_after_dir =
      output_dir / MakeWGDirectoryName(split_datetime_label);

  if (fs::exists(output_before_dir) &&
      fs::equivalent(input_dir, output_before_dir)) {
    throw std::runtime_error(
        "WG output directory must be different from the input directory "
        "to avoid overwriting the original files.");
  }

  EnsureOutputDirectory(output_before_dir);
  EnsureOutputDirectory(output_after_dir);

  const std::vector<int> dif_ids =
      CollectWGDifIds(input_dir, start_datetime_label);

  std::ostringstream before_selection;
  before_selection << "timestamp > 0"
                   << " && timestamp < "
                   << static_cast<long long>(split_unixtime);

  std::ostringstream after_selection;
  after_selection << "timestamp > 0"
                  << " && timestamp >= "
                  << static_cast<long long>(split_unixtime);

  std::cout << "[INFO] Mode            : WG\n";
  std::cout << "[INFO] Time branch     : bsd.timestamp\n";
  std::cout << "[INFO] Input directory : " << input_dir << '\n';
  std::cout << "[INFO] Number of DIFs  : " << dif_ids.size() << '\n';
  std::cout << "[INFO] Split time      : " << split_datetime_label
            << " / unixtime = " << static_cast<long long>(split_unixtime)
            << '\n';

  for (const int dif_id : dif_ids) {
    const fs::path input_file =
        input_dir / MakeWGFileName(start_datetime_label, dif_id);
    const fs::path output_before_file =
        output_before_dir / MakeWGFileName(start_datetime_label, dif_id);
    const fs::path output_after_file =
        output_after_dir / MakeWGFileName(split_datetime_label, dif_id);

    WriteSelectedWGFile(input_file,
                        output_before_file,
                        before_selection.str(),
                        "elist_before");
    WriteSelectedWGFile(input_file,
                        output_after_file,
                        after_selection.str(),
                        "elist_after");
  }

  std::cout << "[INFO] Done.\n";
}

void SplitAtTime(const fs::path& input_path,
                 const std::string& split_time_text,
                 const fs::path& output_dir,
                 const SplitConfig& config) {
  if (!fs::exists(input_path)) {
    throw std::runtime_error("Input file does not exist: " + input_path.string());
  }

  if (!fs::is_regular_file(input_path)) {
    throw std::runtime_error("Input path is not a regular file: " +
                             input_path.string());
  }

  EnsureOutputDirectory(output_dir);

  const std::time_t split_unixtime = ParseDateTimeToUnixTime(split_time_text);

  std::tm split_tm{};
#if defined(_WIN32)
  localtime_s(&split_tm, &split_unixtime);
#else
  localtime_r(&split_unixtime, &split_tm);
#endif

  const std::string split_datetime_label = FormatDateTimeLabel(split_tm);
  const std::string start_datetime_label =
      ExtractStartDateTimeLabelFromInputName(input_path, config.label);

  const fs::path output_before =
      output_dir / MakeOutputFileName(config.label, start_datetime_label);

  const fs::path output_after =
      output_dir / MakeOutputFileName(config.label, split_datetime_label);

  if (fs::equivalent(input_path.parent_path(), output_dir)) {
    throw std::runtime_error(
        "Output directory must be different from the input file directory "
        "to avoid overwriting the original file.");
  }

  TChain chain(kTreeName);
  const int added = chain.Add(input_path.string().c_str());
  if (added == 0) {
    throw std::runtime_error("Failed to add input ROOT file to chain: " +
                             input_path.string());
  }

  if (!chain.GetBranch(config.unixtime_branch_name.c_str())) {
    throw std::runtime_error("Branch '" + config.unixtime_branch_name +
                             "' was not found in TTree 'tree'.");
  }

  std::cout << "[INFO] Mode            : " << config.label << '\n';
  std::cout << "[INFO] Unixtime branch : "
            << config.unixtime_branch_name << '\n';
  std::cout << "[INFO] Input file      : " << input_path << '\n';
  std::cout << "[INFO] Entries         : " << chain.GetEntries() << '\n';
  std::cout << "[INFO] Split time      : " << split_datetime_label
            << " / unixtime = " << static_cast<long long>(split_unixtime)
            << '\n';

  std::ostringstream before_selection;
  before_selection << config.unixtime_branch_name << " > 0"
                   << " && " << config.unixtime_branch_name
                   << " < " << static_cast<long long>(split_unixtime);

  std::ostringstream after_selection;
  after_selection << config.unixtime_branch_name << " > 0"
                  << " && " << config.unixtime_branch_name
                  << " >= " << static_cast<long long>(split_unixtime);

  WriteSelectedTree(chain, output_before, before_selection.str());
  WriteSelectedTree(chain, output_after, after_selection.str());

  std::cout << "[INFO] Done.\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 5) {
    PrintUsage(argv[0]);
    return 1;
  }

  const std::string mode = argv[1];
  const fs::path input_path = argv[2];
  const std::string split_time_text = argv[3];
  const fs::path output_dir = argv[4];

  try {
    if (mode == "-wg") {
      SplitWGAtTime(input_path, split_time_text, output_dir);
    } else {
      const SplitConfig config = GetSplitConfig(mode);
      SplitAtTime(input_path, split_time_text, output_dir, config);
    }
  } catch (const std::exception& e) {
    std::cerr << "[ERROR] " << e.what() << '\n';
    PrintUsage(argv[0]);
    return 1;
  }

  return 0;
}
