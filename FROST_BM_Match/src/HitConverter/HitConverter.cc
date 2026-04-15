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

constexpr Long64_t SEARCH_WINDOW_ENTRIES = 100;

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

  // Branch buffers
  Int_t evnum = 0;
  Int_t cablenum[272] = {0};
  Int_t rayraw[272] = {0};
  Int_t rayrawch[272] = {0};
  Double_t lightyield[272][8] = {{0.0}};
  Double_t unixtime[272] = {0.0};
  Int_t spillnum = 0;
  Int_t pileup_flag[272] = {0};
  Int_t undershoot_flag[272] = {0};
  Double_t baseline[272] = {0.0};

  std::vector<std::vector<double>> *hit_bunch = nullptr;
  std::vector<std::vector<double>> *undershoot_bunch = nullptr;
  std::vector<std::vector<double>> *overlapped_bunch = nullptr;
  std::vector<std::vector<double>> *leading = nullptr;
  std::vector<std::vector<double>> *trailing = nullptr;
  std::vector<std::vector<double>> *leading_fromadc = nullptr;
  std::vector<std::vector<double>> *trailing_fromadc = nullptr;
  std::vector<int> *is_hit = nullptr;
  std::vector<std::vector<double>> *xg = nullptr;
  std::vector<std::vector<double>> *yg = nullptr;
  std::vector<std::vector<double>> *x_rec = nullptr;
  std::vector<std::vector<double>> *y_rec = nullptr;
  std::vector<double> *chi2 = nullptr;
  std::vector<int> *is_multihit = nullptr;
  std::vector<std::vector<int>> *x_rectype = nullptr;
  std::vector<std::vector<int>> *y_rectype = nullptr;
  std::vector<std::vector<std::vector<int>>> *groupx = nullptr;
  std::vector<std::vector<std::vector<int>>> *groupy = nullptr;

  // Current scan position to preserve BM-driven one-pass behavior
  Long64_t next_entry = 0;
};

struct FrostOutputEvent {
  Int_t evnum = 0;
  Int_t cablenum[272] = {0};
  Int_t rayraw[272] = {0};
  Int_t rayrawch[272] = {0};
  Double_t lightyield[272][8] = {{0.0}};
  Double_t unixtime[272] = {0.0};
  Int_t spillnum = 0;
  Int_t pileup_flag[272] = {0};
  Int_t undershoot_flag[272] = {0};
  Double_t baseline[272] = {0.0};

  std::vector<std::vector<double>> hit_bunch;
  std::vector<std::vector<double>> undershoot_bunch;
  std::vector<std::vector<double>> overlapped_bunch;
  std::vector<std::vector<double>> leading;
  std::vector<std::vector<double>> trailing;
  std::vector<std::vector<double>> leading_fromadc;
  std::vector<std::vector<double>> trailing_fromadc;
  std::vector<int> is_hit;
  std::vector<std::vector<double>> xg;
  std::vector<std::vector<double>> yg;
  std::vector<std::vector<double>> x_rec;
  std::vector<std::vector<double>> y_rec;
  std::vector<double> chi2;
  std::vector<int> is_multihit;
  std::vector<std::vector<int>> x_rectype;
  std::vector<std::vector<int>> y_rectype;
  std::vector<std::vector<std::vector<int>>> groupx;
  std::vector<std::vector<std::vector<int>>> groupy;
};

void ResetFrostOutputEvent(FrostOutputEvent &out) {
  out.evnum = 0;
  std::fill_n(out.cablenum, 272, 0);
  std::fill_n(out.rayraw, 272, 0);
  std::fill_n(out.rayrawch, 272, 0);
  std::fill(&out.lightyield[0][0], &out.lightyield[0][0] + 272 * 8, 0.0);
  std::fill_n(out.unixtime, 272, 0.0);
  out.spillnum = 0;
  std::fill_n(out.pileup_flag, 272, 0);
  std::fill_n(out.undershoot_flag, 272, 0);
  std::fill_n(out.baseline, 272, 0.0);

  // Keep the bunch container shape for unmatched events.
  out.hit_bunch.assign(8, {});
  out.undershoot_bunch.assign(8, {});
  out.overlapped_bunch.assign(8, {});
  out.leading.assign(8, {});
  out.trailing.assign(8, {});
  out.leading_fromadc.assign(8, {});
  out.trailing_fromadc.assign(8, {});
  out.is_hit.clear();
  out.xg.assign(8, {});
  out.yg.assign(8, {});
  out.x_rec.assign(8, {});
  out.y_rec.assign(8, {});
  out.chi2.clear();
  out.is_multihit.clear();
  out.x_rectype.assign(8, {});
  out.y_rectype.assign(8, {});
  out.groupx.assign(8, {});
  out.groupy.assign(8, {});
}

void CopyFrostSourceToOutput(const FrostSource &src, FrostOutputEvent &out) {
  out.evnum = src.evnum;
  std::copy_n(src.cablenum, 272, out.cablenum);
  std::copy_n(src.rayraw, 272, out.rayraw);
  std::copy_n(src.rayrawch, 272, out.rayrawch);
  std::copy(&src.lightyield[0][0], &src.lightyield[0][0] + 272 * 8, &out.lightyield[0][0]);
  std::copy_n(src.unixtime, 272, out.unixtime);
  out.spillnum = src.spillnum;
  std::copy_n(src.pileup_flag, 272, out.pileup_flag);
  std::copy_n(src.undershoot_flag, 272, out.undershoot_flag);
  std::copy_n(src.baseline, 272, out.baseline);

  out.hit_bunch = src.hit_bunch ? *src.hit_bunch : std::vector<std::vector<double>>(8);
  out.undershoot_bunch = src.undershoot_bunch ? *src.undershoot_bunch : std::vector<std::vector<double>>(8);
  out.overlapped_bunch = src.overlapped_bunch ? *src.overlapped_bunch : std::vector<std::vector<double>>(8);
  out.leading = src.leading ? *src.leading : std::vector<std::vector<double>>(8);
  out.trailing = src.trailing ? *src.trailing : std::vector<std::vector<double>>(8);
  out.leading_fromadc = src.leading_fromadc ? *src.leading_fromadc : std::vector<std::vector<double>>(8);
  out.trailing_fromadc = src.trailing_fromadc ? *src.trailing_fromadc : std::vector<std::vector<double>>(8);
  out.is_hit = src.is_hit ? *src.is_hit : std::vector<int>{};
  out.xg = src.xg ? *src.xg : std::vector<std::vector<double>>(8);
  out.yg = src.yg ? *src.yg : std::vector<std::vector<double>>(8);
  out.x_rec = src.x_rec ? *src.x_rec : std::vector<std::vector<double>>(8);
  out.y_rec = src.y_rec ? *src.y_rec : std::vector<std::vector<double>>(8);
  out.chi2 = src.chi2 ? *src.chi2 : std::vector<double>{};
  out.is_multihit = src.is_multihit ? *src.is_multihit : std::vector<int>{};
  out.x_rectype = src.x_rectype ? *src.x_rectype : std::vector<std::vector<int>>(8);
  out.y_rectype = src.y_rectype ? *src.y_rectype : std::vector<std::vector<int>>(8);
  out.groupx = src.groupx ? *src.groupx : std::vector<std::vector<std::vector<int>>>(8);
  out.groupy = src.groupy ? *src.groupy : std::vector<std::vector<std::vector<int>>>(8);
}

struct FrostMatchResult {
  bool found = false;
  int source_index = -1;
  Long64_t entry = -1;
  Int_t frost_spillnum = -1;
  Int_t frost_unixtime = -1;
  Int_t time_diff = NTBM_NON_INITIALIZED_VALUE;
  bool insertion_found = false;
  int insertion_source_index = -1;
  Long64_t insertion_entry = -1;
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
  std::smatch match;

  // Pattern 1: run00006_0_999_recon.root
  static const std::regex split_pattern(R"(run(\d+)_(\d+)_(\d+)_recon\.root)");

  // Pattern 2: run00006_recon.root
  static const std::regex single_pattern(R"(run(\d+)_recon\.root)");

  if (std::regex_match(filename, match, split_pattern)) {
    run_number = std::stoi(match[1].str());
    start_event = std::stoi(match[2].str());
    end_event = std::stoi(match[3].str());
    return true;
  }

  if (std::regex_match(filename, match, single_pattern)) {
    run_number = std::stoi(match[1].str());

    // No event-range information in the filename
    start_event = -1;
    end_event = -1;
    return true;
  }

  return false;
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

    source.tree->SetBranchAddress("evnum", &source.evnum);
    source.tree->SetBranchAddress("cablenum", source.cablenum);
    source.tree->SetBranchAddress("rayraw", source.rayraw);
    source.tree->SetBranchAddress("rayrawch", source.rayrawch);
    source.tree->SetBranchAddress("lightyield", source.lightyield);
    source.tree->SetBranchAddress("unixtime", source.unixtime);
    source.tree->SetBranchAddress("spillnum", &source.spillnum);
    source.tree->SetBranchAddress("pileup_flag", source.pileup_flag);
    source.tree->SetBranchAddress("undershoot_flag", source.undershoot_flag);
    source.tree->SetBranchAddress("baseline", source.baseline);
    source.tree->SetBranchAddress("hit_bunch", &source.hit_bunch);
    source.tree->SetBranchAddress("undershoot_bunch", &source.undershoot_bunch);
    source.tree->SetBranchAddress("overlapped_bunch", &source.overlapped_bunch);
    source.tree->SetBranchAddress("leading", &source.leading);
    source.tree->SetBranchAddress("trailing", &source.trailing);
    source.tree->SetBranchAddress("leading_fromadc", &source.leading_fromadc);
    source.tree->SetBranchAddress("trailing_fromadc", &source.trailing_fromadc);
    source.tree->SetBranchAddress("is_hit", &source.is_hit);
    source.tree->SetBranchAddress("xg", &source.xg);
    source.tree->SetBranchAddress("yg", &source.yg);
    source.tree->SetBranchAddress("x_rec", &source.x_rec);
    source.tree->SetBranchAddress("y_rec", &source.y_rec);
    source.tree->SetBranchAddress("chi2", &source.chi2);
    source.tree->SetBranchAddress("is_multihit", &source.is_multihit);
    source.tree->SetBranchAddress("x_rectype", &source.x_rectype);
    source.tree->SetBranchAddress("y_rectype", &source.y_rectype);
    source.tree->SetBranchAddress("groupx", &source.groupx);
    source.tree->SetBranchAddress("groupy", &source.groupy);
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
                                        int start_source_index,
                                        bool use_limited_window,
                                        Long64_t anchor_entry) {
  FrostMatchResult best_match;

  const Int_t wagasci_spill = spill.GetBeamSummary().GetWagasciSpillNumber();
  const Int_t wagasci_unixtime =
      static_cast<Int_t>(spill.GetBeamSummary().GetTimestamp());

  const Int_t wagasci_spill_mod = PositiveMod(wagasci_spill, SPILL_MOD);

  BOOST_LOG_TRIVIAL(debug)
      << "Searching FROST match for BM spill=" << wagasci_spill
      << " (mod " << wagasci_spill_mod << "), unixtime=" << wagasci_unixtime;

  const int first_source_index = std::max(0, start_source_index);
  int last_source_index = static_cast<int>(sources.size()) - 1;
  if (use_limited_window) {
    last_source_index =
        std::min(first_source_index + 1, static_cast<int>(sources.size()) - 1);
  }

  for (int isource = first_source_index;
       isource <= last_source_index; ++isource) {

    auto &source = sources.at(isource);
    const Long64_t nentries = source.tree->GetEntries();
    if (nentries <= 0) continue;

    Long64_t search_begin = 0;
    if (use_limited_window) {
      if (isource == first_source_index) {
        const Long64_t center = (anchor_entry >= 0 ? anchor_entry : source.next_entry);
        search_begin = std::max<Long64_t>(0, center - SEARCH_WINDOW_ENTRIES);
      } else {
        search_begin = 0;
      }
    } else {
      if (isource == first_source_index) {
        search_begin = source.next_entry;
      } else {
        search_begin = 0;
      }
    }

    if (search_begin >= nentries) continue;

    Long64_t search_end = nentries;
    if (use_limited_window) {
      if (isource == first_source_index) {
        const Long64_t center = (anchor_entry >= 0 ? anchor_entry : source.next_entry);
        search_end = std::min(center + SEARCH_WINDOW_ENTRIES + 1, nentries);
      } else {
        search_end = std::min<Long64_t>(SEARCH_WINDOW_ENTRIES, nentries);
      }
    }

    BOOST_LOG_TRIVIAL(debug)
        << "Scanning FROST file=" << source.path
        << " entries [" << search_begin << ", " << search_end << ")"
        << (use_limited_window ? " with limited window" : " with full scan");

    for (Long64_t ientry = search_begin; ientry < search_end; ++ientry) {
      source.tree->GetEntry(ientry);

      const Int_t frost_unixtime = static_cast<Int_t>(source.unixtime[0]);
      const Int_t frost_spill_mod = PositiveMod(source.spillnum, SPILL_MOD);
      const Int_t time_diff = frost_unixtime - wagasci_unixtime;

      if (!best_match.insertion_found && frost_unixtime > wagasci_unixtime) {
        best_match.insertion_found = true;
        best_match.insertion_source_index = isource;
        best_match.insertion_entry = ientry;
      }

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
      << ", unixtime=" << wagasci_unixtime
      << (use_limited_window ? " within limited search window"
                             : " after full scan");

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

  // Build an aligned frost_match tree with the full FROST schema.
  // We fill one entry per BM spill. When there is no matched FROST event,
  // an empty event is written so that frost_match and match_info have the
  // same number of rows.
  TTree *frost_match_tree = new TTree("frost_match", "Aligned FROST events for BM spills");
  FrostOutputEvent frost_out;
  ResetFrostOutputEvent(frost_out);

  frost_match_tree->Branch("evnum", &frost_out.evnum, "evnum/I");
  frost_match_tree->Branch("cablenum", frost_out.cablenum, "cablenum[272]/I");
  frost_match_tree->Branch("rayraw", frost_out.rayraw, "rayraw[272]/I");
  frost_match_tree->Branch("rayrawch", frost_out.rayrawch, "rayrawch[272]/I");
  frost_match_tree->Branch("lightyield", frost_out.lightyield, "lightyield[272][8]/D");
  frost_match_tree->Branch("unixtime", frost_out.unixtime, "unixtime[272]/D");
  frost_match_tree->Branch("spillnum", &frost_out.spillnum, "spillnum/I");
  frost_match_tree->Branch("pileup_flag", frost_out.pileup_flag, "pileup_flag[272]/I");
  frost_match_tree->Branch("undershoot_flag", frost_out.undershoot_flag, "undershoot_flag[272]/I");
  frost_match_tree->Branch("baseline", frost_out.baseline, "baseline[272]/D");
  frost_match_tree->Branch("hit_bunch", &frost_out.hit_bunch);
  frost_match_tree->Branch("undershoot_bunch", &frost_out.undershoot_bunch);
  frost_match_tree->Branch("overlapped_bunch", &frost_out.overlapped_bunch);
  frost_match_tree->Branch("leading", &frost_out.leading);
  frost_match_tree->Branch("trailing", &frost_out.trailing);
  frost_match_tree->Branch("leading_fromadc", &frost_out.leading_fromadc);
  frost_match_tree->Branch("trailing_fromadc", &frost_out.trailing_fromadc);
  frost_match_tree->Branch("is_hit", &frost_out.is_hit);
  frost_match_tree->Branch("xg", &frost_out.xg);
  frost_match_tree->Branch("yg", &frost_out.yg);
  frost_match_tree->Branch("x_rec", &frost_out.x_rec);
  frost_match_tree->Branch("y_rec", &frost_out.y_rec);
  frost_match_tree->Branch("chi2", &frost_out.chi2);
  frost_match_tree->Branch("is_multihit", &frost_out.is_multihit);
  frost_match_tree->Branch("x_rectype", &frost_out.x_rectype);
  frost_match_tree->Branch("y_rectype", &frost_out.y_rectype);
  frost_match_tree->Branch("groupx", &frost_out.groupx);
  frost_match_tree->Branch("groupy", &frost_out.groupy);

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

    ResetFrostOutputEvent(frost_out);

    if (matched) {
      const int source_index = FindSourceIndexByPath(sources, record.frost_file);
      if (source_index < 0) {
        BOOST_LOG_TRIVIAL(warning)
            << "Could not find FROST source for file: " << record.frost_file;
      } else {
        auto &source = sources.at(source_index);

        source.tree->GetEntry(record.frost_entry);

        CopyFrostSourceToOutput(source, frost_out);
      }
    }

    frost_match_tree->Fill();
  }

  outfile->cd();
  frost_match_tree->Write("", TObject::kOverwrite);
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
    bool found_first_match = false;
    Long64_t current_anchor_entry = -1;
    Int_t previous_wagasci_spill = NTBM_NON_INITIALIZED_VALUE;

    {
      B2Reader reader(bm_input_path.c_str());
      B2Writer writer(output_path.c_str(), reader);

      while (reader.ReadNextSpill() > 0) {
        auto &output_spill_summary = writer.GetSpillSummary();

        const Int_t wagasci_spill = output_spill_summary.GetBeamSummary().GetWagasciSpillNumber();
        const Int_t wagasci_unixtime =
            static_cast<Int_t>(output_spill_summary.GetBeamSummary().GetTimestamp());

        bool bm_spill_jumped = false;
        if (previous_wagasci_spill != NTBM_NON_INITIALIZED_VALUE) {
          bm_spill_jumped = (std::abs(wagasci_spill - previous_wagasci_spill) > 1);
        }

        const bool use_limited_window =
            found_first_match && !bm_spill_jumped;

        BOOST_LOG_TRIVIAL(debug)
            << "BM spill=" << wagasci_spill
            << ", previous BM spill=" << previous_wagasci_spill
            << ", bm_spill_jumped=" << bm_spill_jumped
            << ", found_first_match=" << found_first_match
            << ", use_limited_window=" << use_limited_window
            << ", current_source_index=" << current_source_index
            << ", current_anchor_entry=" << current_anchor_entry;

        FrostMatchResult match =
            FindMatchingFrostEvent(output_spill_summary,
                                   frost_sources,
                                   current_source_index,
                                   use_limited_window,
                                   current_anchor_entry);

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
          found_first_match = true;
          current_anchor_entry = match.entry;
        } else {
          record.matched = 0;
          if (match.insertion_found) {
            current_source_index = match.insertion_source_index;
            current_anchor_entry = match.insertion_entry;
            found_first_match = true;
          }
        }

        match_records.push_back(record);
        writer.Fill();
        previous_wagasci_spill = wagasci_spill;
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
