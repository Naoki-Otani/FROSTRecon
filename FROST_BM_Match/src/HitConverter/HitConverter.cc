// system includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <regex>
#include <cctype>

// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// ROOT includes
#include <TFile.h>
#include <TTree.h>

// B2 includes
#include "B2Reader.hh"
#include "B2Writer.hh"
#include "B2SpillSummary.hh"
#include "B2BeamSummary.hh"

#include "NTBMConst.hh"

namespace logging = boost::log;
namespace sfs = std::filesystem;

namespace {

logging::trivial::severity_level ParseLogLevel(std::string level) {
  std::transform(level.begin(), level.end(), level.begin(),
                  [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

  if (level == "trace")   return logging::trivial::trace;
  if (level == "debug")   return logging::trivial::debug;
  if (level == "info")    return logging::trivial::info;
  if (level == "warning") return logging::trivial::warning;
  if (level == "error")   return logging::trivial::error;
  if (level == "fatal")   return logging::trivial::fatal;

  throw std::invalid_argument("Unknown log level: " + level +
                              " (use trace/debug/info/warning/error/fatal)");
}

struct FrostSource {
  std::string path;
  std::unique_ptr<TFile> file;
  TTree *tree = nullptr;

  // Branch buffers used only for matching
  Double_t unixtime[272] = {0.0};
  Int_t spillnum = 0;

  // Current scan position to preserve BM-driven one-pass behavior
  Long64_t next_entry = 0;
};

struct FrostMatchResult {
  bool found = false;
  int source_index = -1;
  Long64_t entry = -1;
  Int_t frost_spillnum = -1;
  Int_t frost_unixtime = -1;
  Int_t time_diff = NTBM_NON_INITIALIZED_VALUE;
};

struct MatchRecord {
  Int_t bm_index = -1;
  Int_t wagasci_spill = -1;
  Int_t wagasci_unixtime = -1;

  Int_t matched = 0;

  std::string frost_file;
  Int_t frost_run_number = -1;
  Int_t frost_start_event = -1;
  Int_t frost_end_event = -1;
  Long64_t frost_entry = -1;
  Int_t frost_spillnum = -1;
  Int_t frost_unixtime = -1;

  Int_t time_diff = NTBM_NON_INITIALIZED_VALUE;
};

struct FrostFileKey {
  std::string path;
  Int_t run_number = -1;
  Int_t start_event = -1;
  Int_t end_event = -1;
};

bool ParseFrostFileName(const std::string &path,
                        Int_t &run_number,
                        Int_t &start_event,
                        Int_t &end_event) {
  const std::string filename = sfs::path(path).filename().string();

  std::regex pattern(R"(run(\d+)_(\d+)_(\d+)_recon\.root)");
  std::smatch match;
  if (!std::regex_match(filename, match, pattern)) {
    return false;
  }

  run_number = std::stoi(match[1].str());
  start_event = std::stoi(match[2].str());
  end_event = std::stoi(match[3].str());
  return true;
}

std::vector<std::string> CollectFrostFiles(const std::string &directory_path) {
  if (!sfs::exists(directory_path)) {
    throw std::runtime_error("FROST input directory does not exist: " + directory_path);
  }
  if (!sfs::is_directory(directory_path)) {
    throw std::runtime_error("FROST input path is not a directory: " + directory_path);
  }

  std::vector<FrostFileKey> files;
  for (const auto &entry : sfs::directory_iterator(directory_path)) {
    if (!entry.is_regular_file()) continue;
    if (entry.path().extension() != ".root") continue;

    FrostFileKey key;
    key.path = entry.path().string();
    if (!ParseFrostFileName(key.path, key.run_number, key.start_event, key.end_event)) {
      BOOST_LOG_TRIVIAL(warning)
          << "Skipping file with unexpected FROST filename format: " << key.path;
      continue;
    }
    files.push_back(key);
  }

  std::sort(files.begin(), files.end(),
            [](const FrostFileKey &a, const FrostFileKey &b) {
              if (a.run_number != b.run_number) return a.run_number < b.run_number;
              if (a.start_event != b.start_event) return a.start_event < b.start_event;
              return a.end_event < b.end_event;
            });

  std::vector<std::string> root_files;
  root_files.reserve(files.size());
  for (const auto &f : files) root_files.push_back(f.path);
  return root_files;
}

std::vector<FrostSource> OpenFrostSources(const std::string &directory_path) {
  std::vector<std::string> file_paths = CollectFrostFiles(directory_path);
  if (file_paths.empty()) {
    throw std::runtime_error("No .root files found in FROST input directory: " + directory_path);
  }

  std::vector<FrostSource> sources;
  sources.reserve(file_paths.size());

  for (const auto &path : file_paths) {
    auto file = std::make_unique<TFile>(path.c_str(), "READ");
    if (!file || file->IsZombie()) {
      throw std::runtime_error("Failed to open FROST file: " + path);
    }
    TTree *tree = dynamic_cast<TTree *>(file->Get("frost"));
    if (!tree) {
      throw std::runtime_error("Tree 'frost' not found in file: " + path);
    }

    BOOST_LOG_TRIVIAL(info) << "Opening FROST file: " << path;

    sources.emplace_back();
    FrostSource &source = sources.back();
    source.path = path;
    source.file = std::move(file);
    source.tree = tree;

    source.tree->SetBranchAddress("unixtime", source.unixtime);
    source.tree->SetBranchAddress("spillnum", &source.spillnum);
  }

  return sources;
}

int FindSourceIndexByPath(const std::vector<FrostSource> &sources,
                          const std::string &path) {
  for (int i = 0; i < static_cast<int>(sources.size()); ++i) {
    if (sources[i].path == path) return i;
  }
  return -1;
}

Int_t PositiveMod(Int_t value, Int_t mod) {
  Int_t result = value % mod;
  if (result < 0) result += mod;
  return result;
}

FrostMatchResult FindMatchingFrostEvent(const B2SpillSummary &spill,
                                        std::vector<FrostSource> &sources,
                                        int start_source_index) {
  FrostMatchResult best_match;

  const Int_t wagasci_spill = spill.GetBeamSummary().GetWagasciSpillNumber();
  const Int_t wagasci_unixtime =
      static_cast<Int_t>(spill.GetBeamSummary().GetTimestamp());

  const Int_t wagasci_spill_mod = PositiveMod(wagasci_spill, SPILL_MOD);

  BOOST_LOG_TRIVIAL(debug)
      << "Searching FROST match for BM spill=" << wagasci_spill
      << " (mod " << wagasci_spill_mod << "), unixtime=" << wagasci_unixtime;

  for (int isource = std::max(0, start_source_index);
       isource < static_cast<int>(sources.size()); ++isource) {

    auto &source = sources.at(isource);
    const Long64_t nentries = source.tree->GetEntries();

    for (Long64_t ientry = source.next_entry; ientry < nentries; ++ientry) {
      source.tree->GetEntry(ientry);

      const Int_t frost_unixtime = static_cast<Int_t>(source.unixtime[0]);
      const Int_t frost_spill_mod = PositiveMod(source.spillnum, SPILL_MOD);
      const Int_t time_diff = frost_unixtime - wagasci_unixtime;

      if (frost_spill_mod == wagasci_spill_mod &&
          std::abs(time_diff) <= MAX_TIME_DIFF) {
        best_match.found = true;
        best_match.source_index = isource;
        best_match.entry = ientry;
        best_match.frost_spillnum = source.spillnum;
        best_match.frost_unixtime = frost_unixtime;
        best_match.time_diff = time_diff;

        BOOST_LOG_TRIVIAL(debug)
            << "Matched BM spill=" << wagasci_spill
            << " to FROST file=" << source.path
            << " entry=" << ientry
            << " frost_spill=" << source.spillnum
            << " frost_unixtime=" << frost_unixtime
            << " time_diff=" << time_diff;

        return best_match;
      }
    }
  }

  BOOST_LOG_TRIVIAL(debug)
      << "No FROST match found for BM spill=" << wagasci_spill
      << ", unixtime=" << wagasci_unixtime;

  return best_match;
}

void WriteAdditionalTrees(const std::string &output_path,
                          std::vector<FrostSource> &sources,
                          const std::vector<MatchRecord> &records) {
  std::unique_ptr<TFile> outfile(new TFile(output_path.c_str(), "UPDATE"));
  if (!outfile || outfile->IsZombie()) {
    throw std::runtime_error("Failed to reopen output ROOT file for UPDATE: " + output_path);
  }

  // Remove old trees if rerunning on the same file
  if (outfile->Get("frost_match")) outfile->Delete("frost_match;*");
  if (outfile->Get("match_info")) outfile->Delete("match_info;*");

  TTree *frost_match_tree = nullptr;
  int clone_source_index = -1;

  // source tree currently connected to frost_match_tree branch addresses
  int current_fill_source_index = -1;

  for (const auto &record : records) {
    if (record.matched) {
      clone_source_index = FindSourceIndexByPath(sources, record.frost_file);
      break;
    }
  }

  if (clone_source_index >= 0) {
    auto &source = sources.at(clone_source_index);
    frost_match_tree = source.tree->CloneTree(0);
    frost_match_tree->SetName("frost_match");
    frost_match_tree->SetTitle("Matched FROST events");

    // Make sure the cloned tree uses the addresses of the initial source tree
    source.tree->CopyAddresses(frost_match_tree);
    current_fill_source_index = clone_source_index;
  }

  TTree *match_info_tree = new TTree("match_info", "BM-FROST matching information");

  Int_t bm_index = -1;
  Int_t wagasci_spill = -1;
  Int_t wagasci_unixtime = -1;
  Int_t matched = 0;
  Char_t frost_file[1024] = "";
  Int_t frost_run_number = -1;
  Int_t frost_start_event = -1;
  Int_t frost_end_event = -1;
  Long64_t frost_entry = -1;
  Int_t frost_spillnum = -1;
  Int_t frost_unixtime = -1;
  Int_t time_diff = NTBM_NON_INITIALIZED_VALUE;

  match_info_tree->Branch("bm_index", &bm_index, "bm_index/I");
  match_info_tree->Branch("wagasci_spill", &wagasci_spill, "wagasci_spill/I");
  match_info_tree->Branch("wagasci_unixtime", &wagasci_unixtime, "wagasci_unixtime/I");
  match_info_tree->Branch("matched", &matched, "matched/I");
  match_info_tree->Branch("frost_file", frost_file, "frost_file/C");
  match_info_tree->Branch("frost_run_number", &frost_run_number, "frost_run_number/I");
  match_info_tree->Branch("frost_start_event", &frost_start_event, "frost_start_event/I");
  match_info_tree->Branch("frost_end_event", &frost_end_event, "frost_end_event/I");
  match_info_tree->Branch("frost_entry", &frost_entry, "frost_entry/L");
  match_info_tree->Branch("frost_spillnum", &frost_spillnum, "frost_spillnum/I");
  match_info_tree->Branch("frost_unixtime", &frost_unixtime, "frost_unixtime/I");
  match_info_tree->Branch("time_diff", &time_diff, "time_diff/I");

  for (const auto &record : records) {
    bm_index = record.bm_index;
    wagasci_spill = record.wagasci_spill;
    wagasci_unixtime = record.wagasci_unixtime;
    matched = record.matched;
    std::snprintf(frost_file, sizeof(frost_file), "%s", record.frost_file.c_str());
    frost_run_number = record.frost_run_number;
    frost_start_event = record.frost_start_event;
    frost_end_event = record.frost_end_event;
    frost_entry = record.frost_entry;
    frost_spillnum = record.frost_spillnum;
    frost_unixtime = record.frost_unixtime;
    time_diff = record.time_diff;

    match_info_tree->Fill();

    if (matched && frost_match_tree) {
      const int source_index = FindSourceIndexByPath(sources, record.frost_file);
      if (source_index < 0) {
        BOOST_LOG_TRIVIAL(warning)
            << "Could not find FROST source for file: " << record.frost_file;
      } else {
        auto &source = sources.at(source_index);

        // When the matched event comes from a different FROST file,
        // switch the branch addresses used by the cloned output tree.
        if (source_index != current_fill_source_index) {
          BOOST_LOG_TRIVIAL(debug)
              << "Switch frost_match source from index "
              << current_fill_source_index << " to " << source_index
              << " (" << source.path << ")";
          source.tree->CopyAddresses(frost_match_tree);
          current_fill_source_index = source_index;
        }

        source.tree->GetEntry(record.frost_entry);
        frost_match_tree->Fill();
      }
    }
  }

  outfile->cd();
  if (frost_match_tree) frost_match_tree->Write("", TObject::kOverwrite);
  match_info_tree->Write("", TObject::kOverwrite);
  outfile->Write("", TObject::kOverwrite);
}

} // namespace

int main(int argc, char *argv[]) {
  logging::trivial::severity_level log_level = logging::trivial::info;
  if (argc == 5) {
    log_level = ParseLogLevel(argv[4]);
  }

  logging::core::get()->set_filter(
      logging::trivial::severity >= log_level);

  BOOST_LOG_TRIVIAL(info) << "==========FROST Hit Converter Start==========";

  if (argc != 4 && argc != 5) {
    std::cerr << "Usage : " << argv[0]
              << " <input wagasci/BM file path>"
              << " <input FROST file directory>"
              << " <output file path>"
              << " [trace|debug|info|warning|error|fatal]"
              << std::endl;
    std::exit(1);
  }

  try {
    const std::string bm_input_path = argv[1];
    const std::string frost_input_dir = argv[2];
    const std::string output_path = argv[3];

    BOOST_LOG_TRIVIAL(info) << "-----Settings Summary-----";
    BOOST_LOG_TRIVIAL(info) << "Reader file        : " << bm_input_path;
    BOOST_LOG_TRIVIAL(info) << "FROST input dir    : " << frost_input_dir;
    BOOST_LOG_TRIVIAL(info) << "Writer file        : " << output_path;

    auto frost_sources = OpenFrostSources(frost_input_dir);

    std::vector<MatchRecord> match_records;
    int current_source_index = 0;
    int bm_index = 0;

    {
      B2Reader reader(bm_input_path.c_str());
      B2Writer writer(output_path.c_str(), reader);

      while (reader.ReadNextSpill() > 0) {
        auto &output_spill_summary = writer.GetSpillSummary();

        const Int_t wagasci_spill = output_spill_summary.GetBeamSummary().GetWagasciSpillNumber();
        const Int_t wagasci_unixtime =
            static_cast<Int_t>(output_spill_summary.GetBeamSummary().GetTimestamp());

        FrostMatchResult match =
            FindMatchingFrostEvent(output_spill_summary, frost_sources, current_source_index);

        MatchRecord record;
        record.bm_index = bm_index;
        record.wagasci_spill = wagasci_spill;
        record.wagasci_unixtime = wagasci_unixtime;

        if (match.found) {
          record.matched = 1;
          record.frost_file = frost_sources.at(match.source_index).path;
          record.frost_entry = match.entry;
          record.frost_spillnum = match.frost_spillnum;
          record.frost_unixtime = match.frost_unixtime;
          record.time_diff = match.time_diff;

          if (!ParseFrostFileName(record.frost_file,
                                  record.frost_run_number,
                                  record.frost_start_event,
                                  record.frost_end_event)) {
            BOOST_LOG_TRIVIAL(warning)
                << "Failed to parse FROST filename: " << record.frost_file;
          }

          auto &matched_source = frost_sources.at(match.source_index);
          matched_source.next_entry = match.entry + 1;
          current_source_index = match.source_index;
        } else {
          record.matched = 0;
        }

        match_records.push_back(record);
        writer.Fill();
        ++bm_index;
      }
    } // B2Writer closes here

    WriteAdditionalTrees(output_path, frost_sources, match_records);

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  } catch (const std::exception &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Exception : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========FROST Hit Converter Finish==========";
  std::exit(0);
}
