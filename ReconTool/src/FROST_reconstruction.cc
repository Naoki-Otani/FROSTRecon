#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <memory>
#include <string>
#include <stdexcept>
#include <array>
#include <vector>
#include <algorithm>
#include <optional>

#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TKey.h"
#include "TObject.h"

#include "CLI.h"
#include "Mapping.h"
#include "MultiHitAlgo.h"
#include "PoissonChi2.h"
#include "ReconstructionUtil.h"
#include "FROSTConstants.h"
#include "TreeReader.h"

namespace {

struct Candidate1D {
  std::vector<double> pos;      // candidate positions in mm
  std::vector<double> xg;       // corresponding xg/yg values
  std::vector<int> rectype;     // mapping type (0/1/2)
  std::vector<std::vector<int>> group; // grouping result (vector of groups)
};

inline int Group2Size(const std::vector<std::vector<int>>& g) {
  return (g.size() >= 2) ? static_cast<int>(g[1].size()) : 0;
}

Candidate1D BuildXCandidates(const FROST::SingleHitMap& single_map,
                             const FROST::TwoHitMap& two_map,
                             const TH1D& hx,
                             const FROST::MultiHitAlgoConfig& cfg) {
  Candidate1D out;

  std::vector<std::vector<int>> groupx;
  FROST::DivideGroup1(&hx, groupx, cfg);
  FROST::DivideGroup2(&hx, groupx, cfg);

  std::vector<double> xgmulti;
  std::vector<int> rectype;
  FROST::MultiReconstructionX(&hx, xgmulti, rectype, groupx, cfg);

  out.group = groupx;
  out.xg = xgmulti;
  out.rectype = rectype;
  out.pos.reserve(xgmulti.size());

  for (size_t i = 0; i < xgmulti.size(); ++i) {
    const double xg = xgmulti[i];
    const int rt = rectype[i];

    double x = 0.0;
    if (rt == 0) {
      x = two_map.XFromXgMinus(xg);
    } else if (rt == 1) {
      x = two_map.XFromXgPlus(xg);
    } else {
      x = single_map.XFromXg(xg);
    }
    out.pos.push_back(x);
  }

  return out;
}

Candidate1D BuildYCandidates(const FROST::SingleHitMap& single_map,
                             const FROST::TwoHitMap& two_map,
                             const TH1D& hy,
                             const FROST::MultiHitAlgoConfig& cfg) {
  Candidate1D out;

  std::vector<std::vector<int>> groupy;
  FROST::DivideGroup1(&hy, groupy, cfg);
  FROST::DivideGroup2(&hy, groupy, cfg);

  std::vector<double> ygmulti;
  std::vector<int> rectype;
  FROST::MultiReconstructionY(&hy, ygmulti, rectype, groupy, cfg);

  out.group = groupy;
  out.xg = ygmulti;
  out.rectype = rectype;
  out.pos.reserve(ygmulti.size());

  for (size_t i = 0; i < ygmulti.size(); ++i) {
    const double yg = ygmulti[i];
    const int rt = rectype[i];

    double y = 0.0;
    if (rt == 0) {
      y = two_map.YFromYgMinus(yg);
    } else if (rt == 1) {
      y = two_map.YFromYgPlus(yg);
    } else {
      y = single_map.YFromYg(yg);
    }
    out.pos.push_back(y);
  }

  return out;
}

}  // namespace

namespace {

enum class InputMode {
  kMC,
  kData
};

InputMode ParseInputMode(int argc, char** argv) {
  bool mc = false;
  bool data = false;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "-mc" || a == "--mc") mc = true;
    if (a == "-data" || a == "--data") data = true;
  }

  if (mc && data) {
    throw std::runtime_error("Specify only one of --mc or --data.");
  }
  if (!mc && !data) {
    throw std::runtime_error("Specify either --mc or --data.");
  }
  return data ? InputMode::kData : InputMode::kMC;
}

std::string DefaultTreeName(InputMode mode) {
  return (mode == InputMode::kData) ? "tree" : "wls";
}

constexpr int kNbunch = 8;

// Validate that the map ROOT file contains completed (1-mm grid) graphs.
// This prevents silently using incomplete maps (label-only graphs), which breaks reconstruction.
void ValidateSingleMapFile(const std::string& path) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("Failed to open single-map file: " + path);

  auto* gx = dynamic_cast<TGraph*>(f->Get("gx"));
  auto* gy = dynamic_cast<TGraph*>(f->Get("gy"));
  if (!gx || !gy) throw std::runtime_error("Missing gx/gy in single-map file: " + path);

  const int nx_expected = FROST::kNfibX * 10 + 1;  // [-660..+660] inclusive => 1321
  const int ny_expected = FROST::kNfibY * 10 + 1;  // [-700..+700] inclusive => 1401
  if (gx->GetN() < nx_expected || gy->GetN() < ny_expected) {
    throw std::runtime_error("Single-map gx/gy are not completed (1-mm grid). Regenerate maps with mapfunction_tool. File: " + path);
  }

  // Check at least one response graph to ensure completed 1-mm grid for chi2.
  auto* gmcx0 = dynamic_cast<TGraph*>(f->Get("gmcx_0"));
  auto* gmcy0 = dynamic_cast<TGraph*>(f->Get("gmcy_0"));
  if (!gmcx0 || !gmcy0) throw std::runtime_error("Missing gmcx_*/gmcy_* in single-map file: " + path);

  if (gmcx0->GetN() < FROST::kNmcPX || gmcy0->GetN() < FROST::kNmcPY) {
    throw std::runtime_error("Single-map gmcx/gmcy are not completed (1-mm grid). Regenerate maps with mapfunction_tool. File: " + path);
  }
}

void ValidateTwoMapFile(const std::string& path) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("Failed to open two-map file: " + path);

  auto* g1xm = dynamic_cast<TGraph*>(f->Get("g1xm"));
  auto* g1xp = dynamic_cast<TGraph*>(f->Get("g1xp"));
  auto* g1ym = dynamic_cast<TGraph*>(f->Get("g1ym"));
  auto* g1yp = dynamic_cast<TGraph*>(f->Get("g1yp"));
  if (!g1xm || !g1xp || !g1ym || !g1yp) {
    throw std::runtime_error("Missing g1xm/g1xp/g1ym/g1yp in two-map file: " + path);
  }

  const int nx_expected = FROST::kNmcPX + 2;  // [-660..+660] inclusive => 1321
  const int ny_expected = FROST::kNmcPY + 2;  // [-700..+700] inclusive => 1401
  if (g1xm->GetN() < nx_expected || g1xp->GetN() < nx_expected ||
      g1ym->GetN() < ny_expected || g1yp->GetN() < ny_expected) {
    throw std::runtime_error("Two-map g1* graphs are not completed (1-mm grid). Regenerate maps with mapfunction_tool. File: " + path);
  }
}

// Copy all TTrees except the source tree and the output frost tree.
void CopyAllObjectsExceptSourceTree(TFile& in, TFile& out, const std::string& source_tree_name) {
  TIter nextKey(in.GetListOfKeys());
  while (TObject* keyObj = nextKey()) {
    auto* key = dynamic_cast<TKey*>(keyObj);
    if (!key) continue;

    const std::string name = key->GetName();
    if (name == source_tree_name) continue;  // handled separately
    if (name == "frost") continue;           // avoid name collision

    std::unique_ptr<TObject> obj(key->ReadObj());
    if (!obj) continue;

    out.cd();

    if (auto* tr = dynamic_cast<TTree*>(obj.get())) {
      // Clone all entries.
      TTree* cloned = tr->CloneTree(-1, "fast");
      cloned->Write();
    } else {
      // Copy histograms, TParameters, etc.
      obj->Write();
    }
  }
}

std::optional<double> ReadTParameterDouble(TFile& f, const std::string& name) {
  if (auto* p = dynamic_cast<TParameter<double>*>(f.Get(name.c_str()))) {
    return p->GetVal();
  }
  // Some files may store ints for thresholds; accept that too.
  if (auto* p = dynamic_cast<TParameter<int>*>(f.Get(name.c_str()))) {
    return static_cast<double>(p->GetVal());
  }
  return std::nullopt;
}

struct MapParams {
  double alpha = 0.0;
  double threshold1 = 0.0;
  double threshold2 = 0.0;
};

MapParams LoadParamsFromMapsOrThrow(const std::string& single_map_path,
                                   const std::string& two_map_path) {
  std::unique_ptr<TFile> fs(TFile::Open(single_map_path.c_str(), "READ"));
  if (!fs || fs->IsZombie()) throw std::runtime_error("Failed to open single-map file: " + single_map_path);
  std::unique_ptr<TFile> ft(TFile::Open(two_map_path.c_str(), "READ"));
  if (!ft || ft->IsZombie()) throw std::runtime_error("Failed to open two-map file: " + two_map_path);

  const auto a_single = ReadTParameterDouble(*fs, "alpha");
  const auto a_two    = ReadTParameterDouble(*ft, "alpha");
  if (!a_single.has_value()) throw std::runtime_error("Missing TParameter<double> 'alpha' in single-map: " + single_map_path);
  if (!a_two.has_value())    throw std::runtime_error("Missing TParameter<double> 'alpha' in two-map: " + two_map_path);
  if (std::abs(*a_single - *a_two) > 1e-12) {
    throw std::runtime_error("alpha mismatch between maps: single-map alpha=" + std::to_string(*a_single) +
                             ", two-map alpha=" + std::to_string(*a_two));
  }

  // threshold1/threshold2 are defined by two-map only (as requested).
  const auto t1 = ReadTParameterDouble(*ft, "threshold1");
  const auto t2 = ReadTParameterDouble(*ft, "threshold2");
  if (!t1.has_value()) throw std::runtime_error("Missing TParameter 'threshold1' in two-map: " + two_map_path);
  if (!t2.has_value()) throw std::runtime_error("Missing TParameter 'threshold2' in two-map: " + two_map_path);

  MapParams out;
  out.alpha = *a_single;
  out.threshold1 = *t1;
  out.threshold2 = *t2;
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    FROST::CLI cli(argc, argv);

    cli.Require("in");
    cli.Require("out");
    cli.Require("single-map");
    cli.Require("two-map");

    const InputMode input_mode = ParseInputMode(argc, argv);

    const std::string in_path = cli.Get("in");
    const std::string out_path = cli.Get("out");
    const std::string single_map_path = cli.Get("single-map");
    const std::string two_map_path = cli.Get("two-map");

    const std::string tree_name = cli.Get("tree", DefaultTreeName(input_mode));

    const int max_events = cli.GetInt("max-events", 0);

    const double chi2_threshold = cli.GetDouble("chi2-threshold", 1.26);

    // Validate mapping files (graph completeness checks).
    ValidateSingleMapFile(single_map_path);
    ValidateTwoMapFile(two_map_path);

    // Load alpha/thresholds from map files.
    // alpha must match between single-map and two-map; mismatch => ERROR.
    const MapParams mp = LoadParamsFromMapsOrThrow(single_map_path, two_map_path);

    FROST::MultiHitAlgoConfig cfg;
    cfg.alpha = mp.alpha;
    cfg.threshold1 = mp.threshold1;
    cfg.threshold2 = mp.threshold2;

    // Construct maps (after validation).
    FROST::SingleHitMap single_map(single_map_path);
    FROST::TwoHitMap two_map(two_map_path);

    const double alpha = cfg.alpha;  // use map-defined alpha everywhere

    // Open input ROOT file (needed to clone/copy trees).
    std::unique_ptr<TFile> fin(TFile::Open(in_path.c_str(), "READ"));
    if (!fin || fin->IsZombie()) throw std::runtime_error("Failed to open input file: " + in_path);

    auto* inTree = dynamic_cast<TTree*>(fin->Get(tree_name.c_str()));
    if (!inTree) throw std::runtime_error("Missing TTree '" + tree_name + "' in input file: " + in_path);

    // Reader for photon arrays etc. (uses the same file).
    FROST::TreeReader reader(in_path, tree_name, /*read_x=*/true, /*read_y=*/true);
    const Long64_t nentries_all = reader.Entries();
    const Long64_t nentries = (max_events > 0) ? std::min<Long64_t>(nentries_all, max_events) : nentries_all;

    // Create output file.
    TFile out(out_path.c_str(), "RECREATE");
    if (out.IsZombie()) throw std::runtime_error("Failed to create output file: " + out_path);

    // Copy everything except the source tree.
    CopyAllObjectsExceptSourceTree(*fin, out, tree_name);

    // Configure which input branches are copied into the output frost tree.
    // For Data input, drop the large waveform branch to reduce output size.
    inTree->SetBranchStatus("*", 1);
    if (input_mode == InputMode::kData && inTree->GetBranch("waveform")) {
      inTree->SetBranchStatus("waveform", 0);
    }

    // Clone source tree structure (0 entries), then add new branches.
    out.cd();
    TTree* outFrost = inTree->CloneTree(0);
    outFrost->SetName("frost");

    // New branches to be added to frost.
    // Structure:
    //   xg[bunch][candidate]
    //   x_rec[bunch][candidate]
    //   groupx[bunch][group][bin]
    std::vector<int> out_is_hit;
    std::vector<std::vector<double>> out_xg, out_yg;
    std::vector<std::vector<double>> out_xrec, out_yrec;
    std::vector<double> out_chi2;
    std::vector<int> out_is_multihit;
    std::vector<std::vector<int>> out_x_rectype, out_y_rectype;
    std::vector<std::vector<std::vector<int>>> out_groupx, out_groupy;

    outFrost->Branch("is_hit", &out_is_hit);
    outFrost->Branch("xg", &out_xg);
    outFrost->Branch("yg", &out_yg);
    outFrost->Branch("x_rec", &out_xrec);
    outFrost->Branch("y_rec", &out_yrec);
    outFrost->Branch("chi2", &out_chi2);
    outFrost->Branch("is_multihit", &out_is_multihit);
    outFrost->Branch("x_rectype", &out_x_rectype);
    outFrost->Branch("y_rectype", &out_y_rectype);
    outFrost->Branch("groupx", &out_groupx);
    outFrost->Branch("groupy", &out_groupy);

    // Histograms for grouping logic (not saved).
    TH1D hx("hx_tmp", "hx_tmp", FROST::kNfibX, -FROST::kNfibX * 5.0, FROST::kNfibX * 5.0);
    TH1D hy("hy_tmp", "hy_tmp", FROST::kNfibY, -FROST::kNfibY * 5.0, FROST::kNfibY * 5.0);
    hx.SetDirectory(nullptr);
    hy.SetDirectory(nullptr);

    const int ndf = FROST::kNfibX + FROST::kNfibY - 2;

    // Progress reporting (1% steps; prints to stderr).
    const Long64_t total = nentries;
    Long64_t next_report = 0;
    int last_percent = -1;
    auto report_progress = [&](Long64_t cur) {
      if (total <= 0) return;
      const int pct = static_cast<int>((100.0 * static_cast<double>(cur)) / static_cast<double>(total));
      if (pct != last_percent) {
        last_percent = pct;
        std::cerr << "\rProcessing: " << std::setw(3) << pct << "% (" << cur << "/" << total << ")" << std::flush;
      }
    };

    for (Long64_t e = 0; e < nentries; ++e) {
      // Progress: update roughly every 1% (or at the first entry).
      if (e == next_report || e == 0) {
        report_progress(e);
        // Set next report position (at least +1 to avoid infinite loop for small totals).
        const Long64_t step = std::max<Long64_t>(1, total / 100);
        next_report = e + step;
      }
      // Load all input branches for this entry into memory.
      reader.GetEntry(e);
      inTree->GetEntry(e);

      out_is_hit.assign(kNbunch, 0);
      out_xg.assign(kNbunch, {});
      out_yg.assign(kNbunch, {});
      out_xrec.assign(kNbunch, {});
      out_yrec.assign(kNbunch, {});
      out_chi2.assign(kNbunch, 0.0);
      out_is_multihit.assign(kNbunch, 0);
      out_x_rectype.assign(kNbunch, {});
      out_y_rectype.assign(kNbunch, {});
      out_groupx.assign(kNbunch, {});
      out_groupy.assign(kNbunch, {});

      const int nbunch = (input_mode == InputMode::kData) ? kNbunch : 1;

      for (int b = 0; b < nbunch; ++b) {
        if (input_mode == InputMode::kData) {
          reader.SetBunch(b);
        }

        const double* px = reader.PhotonX();
        const double* py = reader.PhotonY();

        double maxx = 0.0;
        double maxy = 0.0;
        for (int i = 0; i < FROST::kNfibX; ++i) {
          if (px[i] > maxx) maxx = px[i];
        }
        for (int i = 0; i < FROST::kNfibY; ++i) {
          if (py[i] > maxy) maxy = py[i];
        }

        out_is_hit[b] = (maxx >= 10.0 && maxy >= 10.0) ? 1 : 0;

        const double xg_single = FROST::ComputeXg(px, alpha);
        const double yg_single = FROST::ComputeYg(py, alpha);

        const double x_single = single_map.XFromXg(xg_single);
        const double y_single = single_map.YFromYg(yg_single);

        const double chi2_single = FROST::PoissonChi2(single_map, px, py, x_single, y_single);
        const double chi2red_single = chi2_single / static_cast<double>(ndf);

        out_chi2[b] = chi2red_single;
        out_is_multihit[b] = (chi2red_single > chi2_threshold) ? 1 : 0;

        if (!out_is_multihit[b]) {
          // Single-hit: store as a single candidate. groupx/groupy are empty by design.
          out_xg[b].push_back(xg_single);
          out_xrec[b].push_back(x_single);
          out_x_rectype[b].push_back(2);

          out_yg[b].push_back(yg_single);
          out_yrec[b].push_back(y_single);
          out_y_rectype[b].push_back(2);
        } else {
          for (int i = 0; i < FROST::kNfibX; ++i) hx.SetBinContent(i + 1, px[i]);
          for (int i = 0; i < FROST::kNfibY; ++i) hy.SetBinContent(i + 1, py[i]);

          const Candidate1D xcand = BuildXCandidates(single_map, two_map, hx, cfg);
          const Candidate1D ycand = BuildYCandidates(single_map, two_map, hy, cfg);

          // Store all candidates (aligned by index).
          out_xg[b] = xcand.xg;
          out_xrec[b] = xcand.pos;
          out_x_rectype[b] = xcand.rectype;
          out_groupx[b] = xcand.group;

          out_yg[b] = ycand.xg;
          out_yrec[b] = ycand.pos;
          out_y_rectype[b] = ycand.rectype;
          out_groupy[b] = ycand.group;

          // Fallback safety: ensure at least one candidate exists.
          if (out_xg[b].empty()) {
            out_xg[b].push_back(xg_single);
            out_xrec[b].push_back(x_single);
            out_x_rectype[b].push_back(2);
          }
          if (out_yg[b].empty()) {
            out_yg[b].push_back(yg_single);
            out_yrec[b].push_back(y_single);
            out_y_rectype[b].push_back(2);
          }
        }
      }

      // Fill output frost (contains all original branches + new branches).
      outFrost->Fill();
    }

    // Final progress line.
    report_progress(total);
    std::cerr << "\n";

    out.cd();
    outFrost->Write();

    // Write parameters to output file as TParameter.
    {
      TParameter<double> p_alpha("alpha", cfg.alpha);
      TParameter<double> p_thr1("threshold1", cfg.threshold1);
      TParameter<double> p_thr2("threshold2", cfg.threshold2);
      TParameter<double> p_chi2_threshold("chi2_threshold", chi2_threshold);
      p_alpha.Write();
      p_thr1.Write();
      p_thr2.Write();
      p_chi2_threshold.Write();
    }

    out.Close();

    std::cerr << "Wrote augmented ROOT file: " << out_path << " (entries=" << nentries << ")\n";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n\n";
    std::cerr << "Usage example:\n";
    std::cerr << "  FROST_reconstruction --mc|--data --in input.root --out out.root --single-map single_map.root --two-map two_map.root [--chi2-threshold <chi2_threshold>] [--tree wls|tree]\n";
    return 1;
  }
}
