// system includes
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <memory>
#include <map>
#include <string>
#include <limits>
#include <cctype>
#include <stdexcept>
#include <tuple>

// boost includes
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// root includes
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TString.h>
#include <TCanvas.h>
#include <TRandom.h>

// B2 includes
#include <B2Reader.hh>
#include <B2Writer.hh>
#include <B2Enum.hh>
#include <B2Dimension.hh>
#include <B2SpillSummary.hh>
#include <B2BeamSummary.hh>
#include <B2HitSummary.hh>
#include <B2VertexSummary.hh>
#include <B2ClusterSummary.hh>
#include <B2TrackSummary.hh>
#include <B2EventSummary.hh>
#include <B2EmulsionSummary.hh>
#include <B2Pdg.hh>
#include "NTBMSummary.hh"
#include "NTBMConst.hh"

#include "TrackMatch.hpp"

namespace logging = boost::log;

enum class BeamMode {
  kAuto,
  kFhc,
  kRhc
};

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

  throw std::invalid_argument(
      "Unknown log level: " + level +
      " (use trace/debug/info/warning/error/fatal)");
}

void PrintUsage(const char *program_name) {
  std::cerr
    << "Usage : " << program_name
    << " <input B2 ROOT file>"
    << " <output ROOT file>"
    << " <z shift>"
    << " <MC(0)/data(1)>"
    << " [fhc|rhc]"
    << " [trace|debug|info|warning|error|fatal]"
    << std::endl
    << std::endl
    << "Examples:" << std::endl
    << "  Data:" << std::endl
    << "    " << program_name
    << " afterHitConverter.root afterTrackMatch.root 0.0 1 info"
    << std::endl
    << std::endl
    << "  MC FHC:" << std::endl
    << "    " << program_name
    << " frost_recon_mc.root afterTrackMatchMC.root 0.0 0 fhc info"
    << std::endl
    << std::endl
    << "  MC RHC:" << std::endl
    << "    " << program_name
    << " frost_recon_mc.root afterTrackMatchMC.root 0.0 0 rhc info"
    << std::endl
    << std::endl
    << "Notes:" << std::endl
    << "  - Use 0 for MC and 1 for real data." << std::endl
    << "  - For MC input, fhc or rhc is required explicitly." << std::endl
    << "  - For real data, fhc/rhc is ignored and BsdGoodSpillFlag is used."
    << std::endl;
}

bool IsBeamModeToken(std::string mode) {
  std::transform(mode.begin(), mode.end(), mode.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return mode == "auto" || mode == "fhc" || mode == "rhc";
}

BeamMode ParseBeamMode(std::string mode) {
  std::transform(mode.begin(), mode.end(), mode.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

  if (mode == "auto") return BeamMode::kAuto;
  if (mode == "fhc")  return BeamMode::kFhc;
  if (mode == "rhc")  return BeamMode::kRhc;

  throw std::invalid_argument(
      "Unknown beam mode: " + mode + " (use fhc/rhc)");
}

int BeamModeToBsdGoodSpillFlag(BeamMode beam_mode) {
  if (beam_mode == BeamMode::kFhc) return 1;
  if (beam_mode == BeamMode::kRhc) return -1;

  throw std::invalid_argument(
      "MC beam mode must be explicitly set to fhc or rhc");
}

} // namespace

struct TrackMatchRow {
  Int_t has_match = 0;
  Int_t bm_track_id = -1;
  Int_t frost_match_bunch = -1;
  Int_t ninja_track_type = -1;
  Double_t baby_mind_tangent_x = B2_NON_INITIALIZED_VALUE;
  Double_t baby_mind_tangent_y = B2_NON_INITIALIZED_VALUE;
  Double_t expected_x = B2_NON_INITIALIZED_VALUE;
  Double_t expected_y = B2_NON_INITIALIZED_VALUE;
  Double_t frost_x = B2_NON_INITIALIZED_VALUE;
  Double_t frost_y = B2_NON_INITIALIZED_VALUE;
  Double_t dx = B2_NON_INITIALIZED_VALUE;
  Double_t dy = B2_NON_INITIALIZED_VALUE;
  Double_t tangent_x = B2_NON_INITIALIZED_VALUE;
  Double_t tangent_y = B2_NON_INITIALIZED_VALUE;
  Double_t dtanx = B2_NON_INITIALIZED_VALUE;
  Double_t dtany = B2_NON_INITIALIZED_VALUE;
  Int_t frost_is_hit = -1;
  Double_t external_expected_x = B2_NON_INITIALIZED_VALUE;
  Double_t external_expected_y = B2_NON_INITIALIZED_VALUE;
  Double_t external_tangent_x = B2_NON_INITIALIZED_VALUE;
  Double_t external_tangent_y = B2_NON_INITIALIZED_VALUE;
  Double_t external_chi2_x = B2_NON_INITIALIZED_VALUE;
  Double_t external_chi2_y = B2_NON_INITIALIZED_VALUE;
  Int_t external_ndof_x = B2_NON_INITIALIZED_VALUE;
  Int_t external_ndof_y = B2_NON_INITIALIZED_VALUE;
  Int_t external_num_planes_proton_module_x = 0;
  Int_t external_num_planes_proton_module_y = 0;
  Int_t external_num_planes_downstream_wagasci_x = 0;
  Int_t external_num_planes_downstream_wagasci_y = 0;
  Double_t external_pm_only_expected_x = B2_NON_INITIALIZED_VALUE;
  Double_t external_pm_only_expected_y = B2_NON_INITIALIZED_VALUE;
  Double_t external_pm_only_tangent_x = B2_NON_INITIALIZED_VALUE;
  Double_t external_pm_only_tangent_y = B2_NON_INITIALIZED_VALUE;
  Double_t external_pm_only_chi2_x = B2_NON_INITIALIZED_VALUE;
  Double_t external_pm_only_chi2_y = B2_NON_INITIALIZED_VALUE;
  Int_t external_pm_only_ndof_x = B2_NON_INITIALIZED_VALUE;
  Int_t external_pm_only_ndof_y = B2_NON_INITIALIZED_VALUE;
  Double_t external_dwg_only_expected_x = B2_NON_INITIALIZED_VALUE;
  Double_t external_dwg_only_expected_y = B2_NON_INITIALIZED_VALUE;
  Double_t external_dwg_only_tangent_x = B2_NON_INITIALIZED_VALUE;
  Double_t external_dwg_only_tangent_y = B2_NON_INITIALIZED_VALUE;
  Double_t external_dwg_only_chi2_x = B2_NON_INITIALIZED_VALUE;
  Double_t external_dwg_only_chi2_y = B2_NON_INITIALIZED_VALUE;
  Int_t external_dwg_only_ndof_x = B2_NON_INITIALIZED_VALUE;
  Int_t external_dwg_only_ndof_y = B2_NON_INITIALIZED_VALUE;
  Int_t true_frost_nearest_particle_id = B2_NON_INITIALIZED_VALUE;
  Double_t true_frost_nearest_position_x = B2_NON_INITIALIZED_VALUE;
  Double_t true_frost_nearest_position_y = B2_NON_INITIALIZED_VALUE;
  Double_t true_frost_nearest_tangent_x = B2_NON_INITIALIZED_VALUE;
  Double_t true_frost_nearest_tangent_y = B2_NON_INITIALIZED_VALUE;
};

struct NearestFrostPositionResult {
  bool found = false;
  double frost_position = B2_NON_INITIALIZED_VALUE;
  double diff = B2_NON_INITIALIZED_VALUE;
};

struct ExternalFitGroupKey {
  int detector_id = -1;
  int view = -1;
  int plane = -1;

  bool operator<(const ExternalFitGroupKey &rhs) const {
    return std::tie(detector_id, view, plane) <
           std::tie(rhs.detector_id, rhs.view, rhs.plane);
  }
};

struct ExternalFitRawHit {
  int detector_id = -1;
  int view = -1;
  int plane = -1;
  double position = B2_NON_INITIALIZED_VALUE;
  double z = B2_NON_INITIALIZED_VALUE;
  double position_error = B2_NON_INITIALIZED_VALUE;
  double z_error = B2_NON_INITIALIZED_VALUE;
};

struct ExternalFitPoint {
  int detector_id = -1;
  int view = -1;
  int plane = -1;
  double position = B2_NON_INITIALIZED_VALUE;
  double z = B2_NON_INITIALIZED_VALUE;
  double position_error_low = B2_NON_INITIALIZED_VALUE;
  double position_error_high = B2_NON_INITIALIZED_VALUE;
  double z_error_low = B2_NON_INITIALIZED_VALUE;
  double z_error_high = B2_NON_INITIALIZED_VALUE;
};

struct ExternalFitOneViewResult {
  bool success = false;
  double tangent = B2_NON_INITIALIZED_VALUE;
  double intercept = B2_NON_INITIALIZED_VALUE;
  double chi2 = B2_NON_INITIALIZED_VALUE;
  int ndof = B2_NON_INITIALIZED_VALUE;
};

struct ExternalTrackFitResult {
  ExternalFitOneViewResult x_fit;
  ExternalFitOneViewResult y_fit;
  double expected_x = B2_NON_INITIALIZED_VALUE;
  double expected_y = B2_NON_INITIALIZED_VALUE;
  int num_planes_proton_module_x = 0;
  int num_planes_proton_module_y = 0;
  int num_planes_downstream_wagasci_x = 0;
  int num_planes_downstream_wagasci_y = 0;
};

struct ExternalTrackFitSet {
  ExternalTrackFitResult combined;
  ExternalTrackFitResult proton_module_only;
  ExternalTrackFitResult downstream_wagasci_only;
};

enum class ExternalFitDetectorSelection {
  kProtonModuleAndDownstreamWagasci,
  kProtonModuleOnly,
  kDownstreamWagasciOnly
};

struct TrueFrostParticleInfo {
  std::vector<int> particle_id;
  std::vector<double> position_x;
  std::vector<double> position_y;
  std::vector<double> tangent_x;
  std::vector<double> tangent_y;
};

struct NearestTrueFrostParticleResult {
  bool found = false;
  Int_t particle_id = B2_NON_INITIALIZED_VALUE;
  Double_t position_x = B2_NON_INITIALIZED_VALUE;
  Double_t position_y = B2_NON_INITIALIZED_VALUE;
  Double_t tangent_x = B2_NON_INITIALIZED_VALUE;
  Double_t tangent_y = B2_NON_INITIALIZED_VALUE;
};

bool IsMuonPdg(int particle_id) {
  return particle_id == 13 || particle_id == -13;
}

bool IsFrostHitAtBunch(const std::vector<int> *is_hit, int bm_bunch) {
  // BM bunch is 1..8, while FROST is_hit index is 0..7.
  const int frost_bunch_index = bm_bunch - 1;
  if (!is_hit) return false;
  if (frost_bunch_index < 0) return false;
  if (frost_bunch_index >= static_cast<int>(is_hit->size())) return false;
  return is_hit->at(frost_bunch_index) == 1;
}

int GetFrostMatchIndex(int bm_bunch, int datatype) {
  if (datatype == B2DataType::kMonteCarlo) {
    // FROSTRecon MC output has no bunch structure. All meaningful
    // reconstructed quantities are stored in index 0.
    return 0;
  }

  // Data: BM bunch is 1..8, FROST vector index is 0..7.
  return bm_bunch - 1;
}

bool IsFrostHitForTrack(const std::vector<int> *is_hit,
                        int bm_bunch,
                        int datatype) {
  if (!is_hit) {
    return false;
  }

  const int frost_index = GetFrostMatchIndex(bm_bunch, datatype);
  if (frost_index < 0) {
    return false;
  }
  if (frost_index >= static_cast<int>(is_hit->size())) {
    return false;
  }

  return is_hit->at(frost_index) == 1;
}

std::vector<int> GetFrostMatchIndices(const FrostEntryData &frost,
                                      const std::vector<int> *is_hit,
                                      int bm_bunch,
                                      int datatype) {
  std::vector<int> indices;

  const int frost_index = GetFrostMatchIndex(bm_bunch, datatype);
  if (frost_index < 0) {
    return indices;
  }

  // Data: require is_hit[bm_bunch - 1].
  // MC: require is_hit[0]. FROSTRecon MC output stores all meaningful
  // reconstructed quantities in index 0.
  if (!IsFrostHitForTrack(is_hit, bm_bunch, datatype)) {
    return indices;
  }

  bool has_position = false;
  if (frost.x_rec &&
      frost_index < static_cast<int>(frost.x_rec->size()) &&
      !frost.x_rec->at(frost_index).empty()) {
    has_position = true;
  }
  if (frost.y_rec &&
      frost_index < static_cast<int>(frost.y_rec->size()) &&
      !frost.y_rec->at(frost_index).empty()) {
    has_position = true;
  }

  if (has_position) {
    indices.push_back(frost_index);
  }

  return indices;
}

Int_t GetFrostHitValueForTrack(const FrostEntryData &frost,
                               const std::vector<int> *is_hit,
                               int bm_bunch,
                               int datatype) {
  return IsFrostHitForTrack(is_hit, bm_bunch, datatype) ? 1 : 0;
}

double GetFrostPositionScale(int view, int datatype) {
  if (datatype != B2DataType::kRealData) {
    return 1.0;
  }

  if (view == B2View::kTopView) {
    return FROST_DATA_X_POSITION_SCALE;
  }
  if (view == B2View::kSideView) {
    return FROST_DATA_Y_POSITION_SCALE;
  }

  return 1.0;
}

double CorrectFrostPosition(double frost_position, int view, int datatype) {
  return frost_position * GetFrostPositionScale(view, datatype);
}

FrostInputTrees OpenFrostInputTrees(TFile *file, int datatype) {
  FrostInputTrees trees;
  trees.is_mc = (datatype == B2DataType::kMonteCarlo);

  if (trees.is_mc) {
    trees.frost_input = dynamic_cast<TTree*>(file->Get("frost"));
    trees.frostmc = dynamic_cast<TTree*>(file->Get("frostmc"));
    trees.nRooTracker = dynamic_cast<TTree*>(file->Get("nRooTracker"));

    if (!trees.frost_input) {
      throw std::runtime_error("MC input file must contain frost tree");
    }
    return trees;
  }

  trees.frost_input = dynamic_cast<TTree*>(file->Get("frost_match"));
  trees.match_info = dynamic_cast<TTree*>(file->Get("match_info"));
  if (!trees.frost_input || !trees.match_info) {
    throw std::runtime_error(
      "Data input file must contain frost_match and match_info trees. "
      "Run HitConverter before TrackMatch for data."
    );
  }
  return trees;
}

template <typename BranchType>
void SetRequiredBranchAddress(TTree *tree, const char *branch_name,
                              BranchType **branch_address) {
  if (!tree) {
    throw std::runtime_error("Cannot set branch address on a null TTree");
  }
  if (!tree->GetBranch(branch_name)) {
    throw std::runtime_error(
      std::string("Required branch ") + branch_name +
      " is not found in tree " + tree->GetName());
  }
  tree->SetBranchAddress(branch_name, branch_address);
}

bool FrostMcEntryHasAllBranches(const FrostMcEntryData &frostmc) {
  return frostmc.particle_pdg_id &&
         frostmc.particle_local_first_x_mm &&
         frostmc.particle_local_first_y_mm &&
         frostmc.particle_local_first_z_mm &&
         frostmc.particle_local_last_x_mm &&
         frostmc.particle_local_last_y_mm &&
         frostmc.particle_local_last_z_mm &&
         frostmc.particle_initial_px_mev_c &&
         frostmc.particle_initial_py_mev_c &&
         frostmc.particle_initial_pz_mev_c;
}

std::size_t GetFrostMcParticleCount(const FrostMcEntryData &frostmc) {
  if (!FrostMcEntryHasAllBranches(frostmc)) {
    return 0;
  }

  std::size_t n = frostmc.particle_pdg_id->size();
  n = std::min(n, frostmc.particle_local_first_x_mm->size());
  n = std::min(n, frostmc.particle_local_first_y_mm->size());
  n = std::min(n, frostmc.particle_local_first_z_mm->size());
  n = std::min(n, frostmc.particle_local_last_x_mm->size());
  n = std::min(n, frostmc.particle_local_last_y_mm->size());
  n = std::min(n, frostmc.particle_local_last_z_mm->size());
  n = std::min(n, frostmc.particle_initial_px_mev_c->size());
  n = std::min(n, frostmc.particle_initial_py_mev_c->size());
  n = std::min(n, frostmc.particle_initial_pz_mev_c->size());
  return n;
}

bool ExtrapolateFrostMcParticleToZ0(const FrostMcEntryData &frostmc,
                                    std::size_t iparticle,
                                    double &x_at_z0,
                                    double &y_at_z0) {
  const double x_first = frostmc.particle_local_first_x_mm->at(iparticle);
  const double y_first = frostmc.particle_local_first_y_mm->at(iparticle);
  const double z_first = frostmc.particle_local_first_z_mm->at(iparticle);
  const double x_last = frostmc.particle_local_last_x_mm->at(iparticle);
  const double y_last = frostmc.particle_local_last_y_mm->at(iparticle);
  const double z_last = frostmc.particle_local_last_z_mm->at(iparticle);

  const double dz = z_last - z_first;
  constexpr double kEpsilon = 1.e-9;

  if (std::fabs(dz) < kEpsilon) {
    if (std::fabs(z_first) < kEpsilon) {
      x_at_z0 = x_first;
      y_at_z0 = y_first;
      return true;
    }
    return false;
  }

  const double scale_to_z0 = -z_first / dz;
  x_at_z0 = x_first + scale_to_z0 * (x_last - x_first);
  y_at_z0 = y_first + scale_to_z0 * (y_last - y_first);
  return true;
}

std::size_t GetTrueFrostParticleCount(const TrueFrostParticleInfo &truth) {
  std::size_t n = truth.particle_id.size();
  n = std::min(n, truth.position_x.size());
  n = std::min(n, truth.position_y.size());
  n = std::min(n, truth.tangent_x.size());
  n = std::min(n, truth.tangent_y.size());
  return n;
}

TrueFrostParticleInfo CollectTrueFrostParticleInfo(
    const FrostMcEntryData &frostmc) {
  TrueFrostParticleInfo truth;

  const std::size_t nparticles = GetFrostMcParticleCount(frostmc);
  truth.particle_id.reserve(nparticles);
  truth.position_x.reserve(nparticles);
  truth.position_y.reserve(nparticles);
  truth.tangent_x.reserve(nparticles);
  truth.tangent_y.reserve(nparticles);

  constexpr double kEpsilon = 1.e-9;
  for (std::size_t iparticle = 0; iparticle < nparticles; ++iparticle) {
    const double pz = frostmc.particle_initial_pz_mev_c->at(iparticle);
    if (std::fabs(pz) < kEpsilon) {
      continue;
    }

    double x_at_z0 = B2_NON_INITIALIZED_VALUE;
    double y_at_z0 = B2_NON_INITIALIZED_VALUE;
    if (!ExtrapolateFrostMcParticleToZ0(
          frostmc, iparticle, x_at_z0, y_at_z0)) {
      continue;
    }

    truth.particle_id.push_back(frostmc.particle_pdg_id->at(iparticle));
    truth.position_x.push_back(x_at_z0);
    truth.position_y.push_back(y_at_z0);
    truth.tangent_x.push_back(
      frostmc.particle_initial_px_mev_c->at(iparticle) / pz);
    truth.tangent_y.push_back(
      frostmc.particle_initial_py_mev_c->at(iparticle) / pz);
  }

  return truth;
}

void FillTrueFrostParticleInfo(const TrueFrostParticleInfo &truth,
                               NTBMSummary *ntbm_summary) {
  ntbm_summary->SetTrueFrostParticleInfo(
    truth.particle_id,
    truth.position_x,
    truth.position_y,
    truth.tangent_x,
    truth.tangent_y);
}

NearestTrueFrostParticleResult FindNearestTrueFrostParticle(
    const TrueFrostParticleInfo &truth,
    double expected_x,
    double expected_y) {
  NearestTrueFrostParticleResult result;

  const std::size_t nparticles = GetTrueFrostParticleCount(truth);
  if (nparticles == 0) {
    return result;
  }

  bool has_muon = false;
  for (std::size_t iparticle = 0; iparticle < nparticles; ++iparticle) {
    if (IsMuonPdg(truth.particle_id.at(iparticle))) {
      has_muon = true;
      break;
    }
  }

  double best_distance2 = std::numeric_limits<double>::max();
  for (std::size_t iparticle = 0; iparticle < nparticles; ++iparticle) {
    const int particle_id = truth.particle_id.at(iparticle);
    if (has_muon && !IsMuonPdg(particle_id)) {
      continue;
    }

    const double dx = truth.position_x.at(iparticle) - expected_x;
    const double dy = truth.position_y.at(iparticle) - expected_y;
    const double distance2 = dx * dx + dy * dy;

    if (!result.found || distance2 < best_distance2) {
      result.found = true;
      result.particle_id = particle_id;
      result.position_x = truth.position_x.at(iparticle);
      result.position_y = truth.position_y.at(iparticle);
      result.tangent_x = truth.tangent_x.at(iparticle);
      result.tangent_y = truth.tangent_y.at(iparticle);
      best_distance2 = distance2;
    }
  }

  return result;
}

void CloneTreeIfExists(TFile *input_file, TFile *output_file, const char *tree_name) {
  TTree *input_tree = dynamic_cast<TTree*>(input_file->Get(tree_name));
  if (!input_tree) {
    BOOST_LOG_TRIVIAL(warning)
      << tree_name << " tree not found in input file. Skip cloning.";
    return;
  }

  output_file->cd();
  TTree *output_tree = input_tree->CloneTree(-1);
  output_tree->SetName(tree_name);
  output_tree->SetDirectory(output_file);
  output_tree->Write("", TObject::kOverwrite);
}

NearestFrostPositionResult FindNearestFrostPosition(const FrostEntryData &frost,
                                                    int bm_bunch,
                                                    int view,
                                                    double expected_position,
                                                    const std::vector<int> *is_hit,
                                                    int datatype) {
  NearestFrostPositionResult result;

  const std::vector<std::vector<double>> *source = nullptr;

  if (view == B2View::kTopView) {
    source = frost.x_rec;
  } else if (view == B2View::kSideView) {
    source = frost.y_rec;
  } else {
    return result;
  }

  if (!source) return result;

  const std::vector<int> frost_indices =
    GetFrostMatchIndices(frost, is_hit, bm_bunch, datatype);
  if (frost_indices.empty()) return result;

  double best_abs_diff = std::numeric_limits<double>::max();
  for (const int frost_index : frost_indices) {
    if (frost_index < 0 ||
        frost_index >= static_cast<int>(source->size())) {
      continue;
    }

    const std::vector<double> &positions = source->at(frost_index);
    for (std::size_t i = 0; i < positions.size(); ++i) {
      const double frost_position =
        CorrectFrostPosition(positions.at(i), view, datatype);
      const double diff = frost_position - expected_position;
      const double abs_diff = std::fabs(diff);
      if (!result.found || abs_diff < best_abs_diff) {
        result.found = true;
        result.frost_position = frost_position;
        result.diff = diff;
        best_abs_diff = abs_diff;
      }

    }
  }
  return result;
}

bool CompareBabyMindHitInOneTrack(const B2HitSummary* lhs, const B2HitSummary *rhs) {
  if (lhs->GetView() != rhs->GetView())
    return lhs->GetView() < rhs->GetView();
  if (lhs->GetPlane() != rhs->GetPlane())
    return lhs->GetPlane() < rhs->GetPlane();

  const TVector3 &lpos = lhs->GetScintillatorPosition().GetValue();
  const TVector3 &rpos = rhs->GetScintillatorPosition().GetValue();

  if (lhs->GetView() == B2View::kSideView) {
    if (lpos.Z() != rpos.Z()) return lpos.Z() < rpos.Z();
    return lpos.Y() < rpos.Y();
  }
  if (lhs->GetView() == B2View::kTopView) {
    if (lpos.Z() != rpos.Z()) return lpos.Z() < rpos.Z();
    return lpos.X() < rpos.X();
  }

  return lhs->GetHitId() < rhs->GetHitId();
}

std::vector<std::vector<double> > CalcMergedOnePlanePositionAndError(std::vector<std::vector<double> > position, int view) {

  std::vector<double> xy_position = position.at(0);
  std::sort(xy_position.begin(), xy_position.end());
  std::vector<double> z_position = position.at(1);
  std::sort(z_position.begin(), z_position.end());

  std::vector<std::vector<double> > position_and_error(3);
  for (int i = 0; i < 3; i++) position_and_error.at(i).resize(2);
  // position_and_error.at(pos/higherr/lowerr).at(xy/z)
  const std::size_t number_of_hits = xy_position.size();

  // Calculate position
  // X/Y
  position_and_error.at(0).at(0) = std::accumulate(xy_position.begin(), xy_position.end(), 0.);
  position_and_error.at(0).at(0) /= (double) number_of_hits;
  // Z
  position_and_error.at(0).at(1) = ( z_position.front() + z_position.back() ) / 2.;

  // Calculate error
  switch (view) {
    double xy_area_max, xy_area_min;
    double z_area_max, z_area_min;
  case B2View::kSideView :
    xy_area_max = xy_position.back()  + 0.5 * BM_HORIZONTAL_SCINTI_LARGE / 3.;
    xy_area_min = xy_position.front() - 0.5 * BM_HORIZONTAL_SCINTI_LARGE / 3.;
    z_area_max = z_position.back()  + 0.5 * BM_HORIZONTAL_SCINTI_THICK;
    z_area_min = z_position.front() - 0.5 * BM_HORIZONTAL_SCINTI_THICK;
    // y errors
    position_and_error.at(1).at(0) = xy_area_max - position_and_error.at(0).at(0);
    position_and_error.at(2).at(0) = position_and_error.at(0).at(0) - xy_area_min;
    // z errors
    position_and_error.at(1).at(1) = z_area_max - position_and_error.at(0).at(1);
    position_and_error.at(2).at(1) = position_and_error.at(0).at(1) - z_area_min;
    break;
  case B2View::kTopView :
    if ( xy_position.size() == 2 &&
	 std::fabs( xy_position.front() - xy_position.back() ) < BM_VERTICAL_SCINTI_LARGE ) {
      double overlap = BM_VERTICAL_SCINTI_LARGE - std::fabs( xy_position.front() - xy_position.back() );
      xy_area_max = position_and_error.at(0).at(0) + 0.5 * overlap;
      xy_area_min = position_and_error.at(0).at(0) - 0.5 * overlap;
    } else {
      xy_area_max = xy_position.back()  + 0.5 * BM_VERTICAL_SCINTI_LARGE;
      xy_area_min = xy_position.front() - 0.5 * BM_VERTICAL_SCINTI_LARGE;
    }
    z_area_max = z_position.back()  + 0.5 * BM_VERTICAL_SCINTI_THICK;
    z_area_min = z_position.front() - 0.5 * BM_VERTICAL_SCINTI_THICK;
    // x errors
    position_and_error.at(1).at(0) = xy_area_max - position_and_error.at(0).at(0);
    position_and_error.at(2).at(0) = position_and_error.at(0).at(0) - xy_area_min;
    // z errors
    position_and_error.at(1).at(1) = z_area_max - position_and_error.at(0).at(1);
    position_and_error.at(2).at(1) = position_and_error.at(0).at(1) - z_area_min;
    break;
  default :
    BOOST_LOG_TRIVIAL(error) << "View is not correctly assigned : " << view;
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(trace) << "Position (XY) : "   << position_and_error.at(0).at(0) << ", "
			   << "Position (Z) : "    << position_and_error.at(0).at(1) << ", "
			   << "Error (XY high) : " << position_and_error.at(1).at(0) << ", "
			   << "Error (XY low) : "  << position_and_error.at(2).at(0) << ", "
			   << "Error (Z) : "       << position_and_error.at(1).at(1);

  return position_and_error;

}


std::vector<std::vector<std::vector<std::vector<double> > > > GenerateMergedPositionAndErrors(std::vector<const B2HitSummary* > hits, int datatype){

  std::sort(hits.begin(), hits.end(), CompareBabyMindHitInOneTrack);

  BOOST_LOG_TRIVIAL(trace) << "New track information with hits";
  BOOST_LOG_TRIVIAL(trace) << "Number of Baby MIND hits used for fitting : " << hits.size();

  std::vector<std::vector<std::vector<std::vector<double> > > > merged_position_and_error(2);
  // merged_position_and_error.at(view).at(pos/higherr/lowerr).at(xy/z).at(plane)
  for ( int iview = 0; iview < 2; iview++ ) {
    merged_position_and_error.at(iview).resize(3);
    for ( int iposerr = 0; iposerr < 3; iposerr++ ) {
      merged_position_and_error.at(iview).at(iposerr).resize(2);
    }
  }

  std::vector<std::vector<double> > position_tmp(2);
  // position_tmp.at(xy/z).at(hits)

  for ( int ihit = 0; ihit < hits.size(); ihit++ ) {
    const auto hit = hits.at(ihit);

    int view = hit->GetView();
    int plane = hit->GetPlane();
    int channel = hit->GetSlot().GetValue(hit->GetReadout1());
    BOOST_LOG_TRIVIAL(trace) << "Detector : " << DETECTOR_NAMES.at(hit->GetDetectorId()) << ", "
			     << "View : "     << VIEW_NAMES.at(hit->GetView()) << ", "
			     << "Plane : "    << hit->GetPlane() << ", "
			     << "Channel : "  << hit->GetSlot().GetValue(hit->GetReadout1());

    const TVector3 &pos = hit->GetScintillatorPosition().GetValue();

    if ( datatype == B2DataType::kRealData && plane >= 2 )
      position_tmp.at(1).push_back(pos.Z() + BM_SCI_CORRECTION);
    else
      position_tmp.at(1).push_back(pos.Z());

    switch (view) {
    case B2View::kSideView :
      position_tmp.at(0).push_back(pos.Y());
      break;
    case B2View::kTopView :
      position_tmp.at(0).push_back(pos.X());
      break;
    default :
      BOOST_LOG_TRIVIAL(error) << "View is not correctly assigned";
      std::exit(1);
    }

    if ( ( hit != hits.back() &&
	   ( view != hits.at(ihit+1)->GetView() ||
	     plane != hits.at(ihit+1)->GetPlane() ) ) ||
	 hit == hits.back() ) {
      std::vector<std::vector<double> > one_plane_pos_and_err = CalcMergedOnePlanePositionAndError(position_tmp, view);
      for ( int iposerr = 0; iposerr < 3; iposerr++ )
	for ( int ixyz = 0; ixyz < 2; ixyz++ )
      merged_position_and_error.at(view).at(iposerr).at(ixyz).push_back(one_plane_pos_and_err.at(iposerr).at(ixyz));
      position_tmp.clear(); position_tmp.resize(2);
    }

  }

  return merged_position_and_error;

}


std::vector<std::vector<double> > FitBabyMind(const B2TrackSummary *track, int datatype) {
  std::vector<std::vector<double> > param(2);
  param.at(0).resize(2); param.at(1).resize(2);
  for ( int iview = 0; iview < 2; iview++ )
    for ( int iparam = 0; iparam < 2; iparam++ )
      param.at(iview).at(iparam) = -9999; // If fitting cannot be done correctly, ignore the track

  TF1 *linear[2];
  TGraphAsymmErrors *hit_graph[2];
  for ( int iview = 0; iview < 2; iview++ ) {
    linear[iview] = new TF1(Form("linear %d", iview), "[0] * x + [1]", -2000., 2000.);
    linear[iview]->SetParameter(0, 0.);
    linear[iview]->SetParameter(1, 0.);
  }

  std::vector<const B2HitSummary* > hits;

  auto it_cluster = track->BeginCluster();
  while ( const auto *cluster = it_cluster.Next() ) {
    auto it_hit = cluster->BeginHit();
    while ( const auto *hit = it_hit.Next() ) {
      if ( hit->GetDetectorId() != B2Detector::kBabyMind ) continue;
      if ( hit->GetView() == B2View::kSideView &&
	   hit->GetPlane() > 2) continue; // We only use upstream three planes for sideview
      hits.push_back(hit);
    } // hit
  } // cluster

  std::vector<std::vector<std::vector<std::vector<double> > > > position_and_errors = GenerateMergedPositionAndErrors(hits, datatype);
  // position_and_errors.at(view).at(pos/higherr/lowerr).at(xy/z).at(plane)
  // sideview vectors
  std::vector<Double_t> position_side_y = position_and_errors.at(0).at(0).at(0);
  std::vector<Double_t> position_side_z = position_and_errors.at(0).at(0).at(1);
  std::vector<Double_t> higherr_side_y = position_and_errors.at(0).at(1).at(0);
  std::vector<Double_t> higherr_side_z = position_and_errors.at(0).at(1).at(1);
  std::vector<Double_t> lowerr_side_y = position_and_errors.at(0).at(2).at(0);
  std::vector<Double_t> lowerr_side_z = position_and_errors.at(0).at(2).at(1);
  // topview vectors
  std::vector<Double_t> position_top_x = position_and_errors.at(1).at(0).at(0);
  std::vector<Double_t> position_top_z = position_and_errors.at(1).at(0).at(1);
  std::vector<Double_t> higherr_top_x = position_and_errors.at(1).at(1).at(0);
  std::vector<Double_t> higherr_top_z = position_and_errors.at(1).at(1).at(1);
  std::vector<Double_t> lowerr_top_x = position_and_errors.at(1).at(2).at(0);
  std::vector<Double_t> lowerr_top_z = position_and_errors.at(1).at(2).at(1);


  for ( int iview = 0; iview < 2; iview++ ) {
    if ( iview == B2View::kSideView ) {
      hit_graph[iview] = new TGraphAsymmErrors(position_side_z.size(),
					       &position_side_z[0],
					       &position_side_y[0],
					       &lowerr_side_z[0],
					       &higherr_side_z[0],
					       &lowerr_side_y[0],
					       &higherr_side_y[0]);
    } else if ( iview == B2View::kTopView ) {
      hit_graph[iview] = new TGraphAsymmErrors(position_top_z.size(),
					       &position_top_z[0],
					       &position_top_x[0],
					       &lowerr_top_z[0],
					       &higherr_top_z[0],
					       &lowerr_top_x[0],
					       &higherr_top_x[0]);
    }

    hit_graph[iview]->Fit(linear[iview],"Q","");
    param.at(iview).at(0) = linear[iview]->GetParameter(0);
    param.at(iview).at(1) = linear[iview]->GetParameter(1);
  }

  delete linear[0];
  delete linear[1];

  return param;

}

std::vector<double> GetBabyMindInitialDirectionAndPosition(const B2TrackSummary *track, int datatype) {

  std::vector<double> initial_direction_and_position(4);
  // 0 : tan Y, 1 : tan X, 2 : pos Y, 3 : pos X

  std::vector<std::vector<double> > param = FitBabyMind(track, datatype);
  for (int iview = 0; iview < 2; iview++) {
    initial_direction_and_position.at(iview) = param.at(iview).at(0);
    initial_direction_and_position.at(iview + 2)
      = param.at(iview).at(1) + param.at(iview).at(0) * BM_SECOND_LAYER_POS;
  }

  return initial_direction_and_position;

}

std::vector<double> CalculateBabyMindSecondLayerPositionInFrostCoordinate(
    NTBMSummary *ntbm, int itrack) {

  const std::vector<double> bm2_position = ntbm->GetBabyMindPosition(itrack);

  std::vector<double> position(2);
  const std::vector<double> baby_mind_position = {BABYMIND_POS_Y, BABYMIND_POS_X};
  const std::vector<double> ninja_overall_position = {NINJA_POS_Y, NINJA_POS_X};
  const std::vector<double> ninja_frost_position = {NINJA_FROST_POS_Y, NINJA_FROST_POS_X};

  for (int iview = 0; iview < 2; ++iview) {
    position.at(iview) = bm2_position.at(iview)
      + baby_mind_position.at(iview)
      - ninja_overall_position.at(iview)
      - ninja_frost_position.at(iview);
  }

  return position;
}

std::vector<double> CalculateExpectedPosition(NTBMSummary *ntbm, int itrack, double z_shift) {

  // Pre reconstructed position/direction in BM coordinate
  std::vector<double> baby_mind_pre_direction = ntbm->GetBabyMindTangent(itrack);
  std::vector<double> baby_mind_pre_position = ntbm->GetBabyMindPosition(itrack);

  std::vector<double> position(2);
  std::vector<double> distance(2);

  std::vector<double> baby_mind_position = {BABYMIND_POS_Y, BABYMIND_POS_X};
  std::vector<double> ninja_overall_position = {NINJA_POS_Y, NINJA_POS_X};
  std::vector<double> ninja_frost_position = {NINJA_FROST_POS_Y, NINJA_FROST_POS_X};

  for ( int iview = 0; iview < 2; iview++ ) {
    // extrapolate Baby MIND track to the tracker position
    distance.at(iview) = BABYMIND_POS_Z + BM_SECOND_LAYER_POS
      - NINJA_POS_Z_FROST - NINJA_FROST_POS_Z - (2 * iview - 1) * 10. + z_shift;
    position.at(iview) = baby_mind_pre_position.at(iview) - baby_mind_pre_direction.at(iview) * distance.at(iview);
    // convert coordinate from BM to the tracker
    position.at(iview) = position.at(iview) + baby_mind_position.at(iview)
      - ninja_overall_position.at(iview) - ninja_frost_position.at(iview);
  }

  return position;

}

bool IsDetectorSelectedForExternalFit(
    int detector_id,
    ExternalFitDetectorSelection detector_selection) {
  switch (detector_selection) {
    case ExternalFitDetectorSelection::kProtonModuleAndDownstreamWagasci:
      return detector_id == B2Detector::kProtonModule ||
             detector_id == B2Detector::kWagasciDownstream;
    case ExternalFitDetectorSelection::kProtonModuleOnly:
      return detector_id == B2Detector::kProtonModule;
    case ExternalFitDetectorSelection::kDownstreamWagasciOnly:
      return detector_id == B2Detector::kWagasciDownstream;
    default:
      return false;
  }
}

TVector3 GetDetectorCenterPosition(int detector_id) {
  switch (detector_id) {
    case B2Detector::kWagasciUpstream:
      return TVector3(WAGASCI_UPSTREAM_POS_X,
                      WAGASCI_UPSTREAM_POS_Y,
                      WAGASCI_UPSTREAM_POS_Z);
    case B2Detector::kProtonModule:
      return TVector3(PROTON_MODULE_POS_X,
                      PROTON_MODULE_POS_Y,
                      PROTON_MODULE_POS_Z);
    case B2Detector::kWagasciDownstream:
      return TVector3(WAGASCI_DOWNSTREAM_POS_X,
                      WAGASCI_DOWNSTREAM_POS_Y,
                      WAGASCI_DOWNSTREAM_POS_Z);
    case B2Detector::kBabyMind:
      return TVector3(BABYMIND_POS_X,
                      BABYMIND_POS_Y,
                      BABYMIND_POS_Z);
    default:
      return TVector3(0., 0., 0.);
  }
}

TVector3 ConvertDetectorLocalToFrostCoordinate(int detector_id,
                                               const TVector3 &local_position) {
  const TVector3 detector_center = GetDetectorCenterPosition(detector_id);
  TVector3 position = local_position + detector_center;

  position.SetX(position.X() - NINJA_POS_X - NINJA_FROST_POS_X);
  position.SetY(position.Y() - NINJA_POS_Y - NINJA_FROST_POS_Y);
  position.SetZ(position.Z() - NINJA_POS_Z_FROST - NINJA_FROST_POS_Z);

  return position;
}

bool UseHitForExternalTrackFit(
    int detector_id,
    int view,
    ExternalFitDetectorSelection detector_selection) {
  if (!IsDetectorSelectedForExternalFit(detector_id, detector_selection)) {
    return false;
  }
  if (view != B2View::kTopView && view != B2View::kSideView) {
    return false;
  }

  return true;
}

void IncrementExternalPlaneCount(ExternalTrackFitResult &result,
                                 int detector_id,
                                 int view) {
  const bool is_x_view = (view == B2View::kTopView);
  const bool is_y_view = (view == B2View::kSideView);

  if (!is_x_view && !is_y_view) {
    return;
  }

  switch (detector_id) {
    case B2Detector::kProtonModule:
      if (is_x_view) {
        result.num_planes_proton_module_x++;
      } else {
        result.num_planes_proton_module_y++;
      }
      break;
    case B2Detector::kWagasciDownstream:
      if (is_x_view) {
        result.num_planes_downstream_wagasci_x++;
      } else {
        result.num_planes_downstream_wagasci_y++;
      }
      break;
    default:
      break;
  }
}

std::vector<ExternalFitRawHit> CollectExternalFitRawHits(
    const B2TrackSummary *track,
    int datatype,
    B2Dimension &dimension,
    ExternalFitDetectorSelection detector_selection) {
  std::vector<ExternalFitRawHit> raw_hits;

  auto it_cluster = track->BeginCluster();
  while (const auto *cluster = it_cluster.Next()) {
    auto it_hit = cluster->BeginHit();
    while (const auto *hit = it_hit.Next()) {
      const int detector_id = hit->GetDetectorId();
      const int view = hit->GetView();
      const int plane = hit->GetPlane();

      if (!UseHitForExternalTrackFit(detector_id, view, detector_selection)) {
        continue;
      }

      const unsigned int slot =
        hit->GetSlot().GetValue(hit->GetReadout1());

      TVector3 local_position;
      if (!dimension.GetPosition(
            static_cast<B2Detector>(detector_id),
            static_cast<B2View>(view),
            static_cast<unsigned int>(plane),
            slot,
            local_position)) {
        BOOST_LOG_TRIVIAL(debug)
          << "Skip external-fit hit because position lookup failed"
          << " : detector=" << detector_id
          << ", view=" << view
          << ", plane=" << plane
          << ", slot=" << slot;
        continue;
      }

      TVector3 local_error;
      if (!B2Dimension::GetError(
            static_cast<B2Detector>(detector_id),
            static_cast<B2View>(view),
            static_cast<unsigned int>(plane),
            slot,
            local_error)) {
        BOOST_LOG_TRIVIAL(debug)
          << "Skip external-fit hit because error lookup failed"
          << " : detector=" << detector_id
          << ", view=" << view
          << ", plane=" << plane
          << ", slot=" << slot;
        continue;
      }

      const TVector3 frost_position =
        ConvertDetectorLocalToFrostCoordinate(detector_id, local_position);

      ExternalFitRawHit raw_hit;
      raw_hit.detector_id = detector_id;
      raw_hit.view = view;
      raw_hit.plane = plane;
      raw_hit.z = frost_position.Z();
      raw_hit.z_error = local_error.Z();

      if (view == B2View::kTopView) {
        raw_hit.position = frost_position.X();
        raw_hit.position_error = local_error.X();
      } else if (view == B2View::kSideView) {
        raw_hit.position = frost_position.Y();
        raw_hit.position_error = local_error.Y();
      }

      raw_hits.push_back(raw_hit);
    }
  }

  return raw_hits;
}

std::vector<ExternalFitPoint> BuildExternalFitPoints(
    const std::vector<ExternalFitRawHit> &raw_hits) {
  std::map<ExternalFitGroupKey, std::vector<ExternalFitRawHit>> grouped_hits;

  for (const auto &raw_hit : raw_hits) {
    ExternalFitGroupKey key;
    key.detector_id = raw_hit.detector_id;
    key.view = raw_hit.view;
    key.plane = raw_hit.plane;
    grouped_hits[key].push_back(raw_hit);
  }

  std::vector<ExternalFitPoint> fit_points;
  for (const auto &group : grouped_hits) {
    const auto &key = group.first;
    const auto &hits = group.second;
    if (hits.empty()) {
      continue;
    }

    double position_sum = 0.;
    double position_min = std::numeric_limits<double>::max();
    double position_max = -std::numeric_limits<double>::max();
    double z_min = std::numeric_limits<double>::max();
    double z_max = -std::numeric_limits<double>::max();

    for (const auto &hit : hits) {
      position_sum += hit.position;
      position_min = std::min(position_min,
                              hit.position - hit.position_error);
      position_max = std::max(position_max,
                              hit.position + hit.position_error);
      z_min = std::min(z_min, hit.z - hit.z_error);
      z_max = std::max(z_max, hit.z + hit.z_error);
    }

    ExternalFitPoint point;
    point.detector_id = key.detector_id;
    point.view = key.view;
    point.plane = key.plane;
    point.position = position_sum / static_cast<double>(hits.size());
    point.z = 0.5 * (z_min + z_max);
    point.position_error_low = point.position - position_min;
    point.position_error_high = position_max - point.position;
    point.z_error_low = point.z - z_min;
    point.z_error_high = z_max - point.z;
    fit_points.push_back(point);
  }

  return fit_points;
}

ExternalFitOneViewResult FitExternalOneView(
    const std::vector<ExternalFitPoint> &fit_points,
    int view) {
  ExternalFitOneViewResult result;

  std::vector<Double_t> z_values;
  std::vector<Double_t> position_values;
  std::vector<Double_t> z_error_low;
  std::vector<Double_t> z_error_high;
  std::vector<Double_t> position_error_low;
  std::vector<Double_t> position_error_high;

  for (const auto &point : fit_points) {
    if (point.view != view) {
      continue;
    }
    z_values.push_back(point.z);
    position_values.push_back(point.position);
    z_error_low.push_back(point.z_error_low);
    z_error_high.push_back(point.z_error_high);
    position_error_low.push_back(point.position_error_low);
    position_error_high.push_back(point.position_error_high);
  }

  if (z_values.size() < 2) {
    return result;
  }

  const auto minmax_z =
    std::minmax_element(z_values.begin(), z_values.end());
  if (std::fabs(*minmax_z.second - *minmax_z.first) < 1.e-9) {
    return result;
  }

  TGraphAsymmErrors graph(
    static_cast<int>(z_values.size()),
    z_values.data(),
    position_values.data(),
    z_error_low.data(),
    z_error_high.data(),
    position_error_low.data(),
    position_error_high.data());

  static int function_id = 0;
  TF1 linear(Form("external_linear_%d_%d", view, function_id++),
             "[0] * x + [1]", -5000., 5000.);
  linear.SetParameter(0, 0.);
  linear.SetParameter(1, 0.);

  const int fit_status = graph.Fit(&linear, "Q");
  if (fit_status != 0) {
    return result;
  }

  result.success = true;
  result.tangent = linear.GetParameter(0);
  result.intercept = linear.GetParameter(1);
  result.chi2 = linear.GetChisquare();
  result.ndof = linear.GetNDF();

  return result;
}

ExternalTrackFitResult FitExternalTrackToFrost(
    const B2TrackSummary *track,
    int datatype,
    B2Dimension &dimension,
    ExternalFitDetectorSelection detector_selection) {
  ExternalTrackFitResult result;

  const std::vector<ExternalFitRawHit> raw_hits =
    CollectExternalFitRawHits(track, datatype, dimension, detector_selection);
  const std::vector<ExternalFitPoint> fit_points =
    BuildExternalFitPoints(raw_hits);

  for (const auto &point : fit_points) {
    IncrementExternalPlaneCount(result, point.detector_id, point.view);
  }

  result.x_fit = FitExternalOneView(fit_points, B2View::kTopView);
  result.y_fit = FitExternalOneView(fit_points, B2View::kSideView);

  // The external fit is performed in the FROST local coordinate system.
  // The FROST scintillator center is z = 0, so the expected position at
  // FROST is the intercept of the fitted line.
  if (result.x_fit.success) {
    result.expected_x = result.x_fit.intercept;
  }
  if (result.y_fit.success) {
    result.expected_y = result.y_fit.intercept;
  }

  return result;
}

void FillRowExternalFitResult(
    const ExternalTrackFitResult &external_fit,
    Double_t &expected_x,
    Double_t &expected_y,
    Double_t &tangent_x,
    Double_t &tangent_y,
    Double_t &chi2_x,
    Double_t &chi2_y,
    Int_t &ndof_x,
    Int_t &ndof_y) {
  expected_x = external_fit.expected_x;
  expected_y = external_fit.expected_y;
  tangent_x = external_fit.x_fit.tangent;
  tangent_y = external_fit.y_fit.tangent;
  chi2_x = external_fit.x_fit.chi2;
  chi2_y = external_fit.y_fit.chi2;
  ndof_x = external_fit.x_fit.ndof;
  ndof_y = external_fit.y_fit.ndof;
}

bool NinjaHitExpected(NTBMSummary *ntbm, int itrack, double z_shift) {

  //top view: x, side view: y
  std::vector<double> hit_expected_position = CalculateExpectedPosition(ntbm, itrack, z_shift);
  // Extrapolated position inside tracker area TODO
  if ( std::fabs(hit_expected_position.at(B2View::kTopView))
          > NINJA_FROST_SCI_WIDTH/2. - TEMPORAL_ALLOWANCE[B2View::kTopView] ||
       std::fabs(hit_expected_position.at(B2View::kSideView))
          > NINJA_FROST_SCI_HEIGHT/2. - TEMPORAL_ALLOWANCE[B2View::kSideView] )
    return false;

  // Downstream WAGASCI interaction
  if ( ntbm->GetNinjaTrackType(itrack) == -1 )
    return false;

  return true;

}

FrostTrackCandidates CollectFrostTrackCandidates(NTBMSummary* ntbm, int itrack,
                                                 const FrostEntryData &frost,
                                                 double z_shift,
                                                 const std::vector<int> *is_hit,
                                                 int datatype) {
  FrostTrackCandidates candidates;

  const std::vector<double> hit_expected_position =
    CalculateExpectedPosition(ntbm, itrack, z_shift);

  const int bm_bunch = ntbm->GetBunch(itrack);
  const std::vector<int> frost_indices =
    GetFrostMatchIndices(frost, is_hit, bm_bunch, datatype);

  if (frost_indices.empty()) {
    BOOST_LOG_TRIVIAL(debug)
      << "Skip BM track " << itrack
      << " because no usable FROST candidate was found"
      << " : bm_bunch=" << bm_bunch
      << ", datatype=" << datatype;
    return candidates;
  }

  const double dz_y =
    BABYMIND_POS_Z + BM_SECOND_LAYER_POS - NINJA_POS_Z_FROST - NINJA_FROST_POS_Z
    - (2 * B2View::kSideView - 1) * 10. + z_shift;
  const double dz_x =
    BABYMIND_POS_Z + BM_SECOND_LAYER_POS - NINJA_POS_Z_FROST - NINJA_FROST_POS_Z
    - (2 * B2View::kTopView - 1) * 10. + z_shift;

  if (frost.y_rec) {
    for (const int frost_index : frost_indices) {
      if (frost_index < 0 ||
          frost_index >= static_cast<int>(frost.y_rec->size())) {
        continue;
      }
      const std::vector<double> &ybunch = frost.y_rec->at(frost_index);
      for (std::size_t j = 0; j < ybunch.size(); ++j) {
        const double frost_y =
          CorrectFrostPosition(ybunch.at(j), B2View::kSideView, datatype);
        const double expected_y = hit_expected_position.at(B2View::kSideView);
        const std::vector<double> bm2_position_frost =
          CalculateBabyMindSecondLayerPositionInFrostCoordinate(ntbm, itrack);
        const double bm2_y = bm2_position_frost.at(B2View::kSideView);
        const double dy_tangent = bm2_y - frost_y;

        const double dy = frost_y - expected_y;
        candidates.frost_y_candidates.push_back(frost_y);
        candidates.expected_y_candidates.push_back(expected_y);
        candidates.difference_y_candidates.push_back(dy);
        candidates.tangent_y_candidates.push_back(dy_tangent / dz_y);
      }
    }
  }

  if (frost.x_rec) {
    for (const int frost_index : frost_indices) {
      if (frost_index < 0 ||
          frost_index >= static_cast<int>(frost.x_rec->size())) {
        continue;
      }
      const std::vector<double> &xbunch = frost.x_rec->at(frost_index);
      for (std::size_t j = 0; j < xbunch.size(); ++j) {
        const double frost_x =
          CorrectFrostPosition(xbunch.at(j), B2View::kTopView, datatype);
        const double expected_x = hit_expected_position.at(B2View::kTopView);
        const std::vector<double> bm2_position_frost =
          CalculateBabyMindSecondLayerPositionInFrostCoordinate(ntbm, itrack);
        const double bm2_x = bm2_position_frost.at(B2View::kTopView);
        const double dx_tangent = bm2_x - frost_x;

        const double dx = frost_x - expected_x;
        candidates.frost_x_candidates.push_back(frost_x);
        candidates.expected_x_candidates.push_back(expected_x);
        candidates.difference_x_candidates.push_back(dx);
        candidates.tangent_x_candidates.push_back(dx_tangent / dz_x);
      }
    }
  }

  return candidates;
}

FrostMatchCount CountAcceptedFrostMatches(const FrostTrackCandidates &candidates) {
  FrostMatchCount count;
  count.x_match_count = static_cast<int>(candidates.frost_x_candidates.size());
  count.y_match_count = static_cast<int>(candidates.frost_y_candidates.size());
  return count;
}

void AppendFrostTrackCandidates(NTBMSummary *ntbm, int bm_track_id, int bunch,
                                const FrostTrackCandidates &candidates) {
  if (candidates.frost_x_candidates.empty() && candidates.frost_y_candidates.empty()) {
    return;
  }

  const int nentry = ntbm->GetNumberOfFrostMatchEntries();
  ntbm->SetNumberOfFrostMatchEntries(nentry + 1);

  ntbm->SetBabyMindTrackId(nentry, bm_track_id);
  ntbm->SetFrostMatchBunch(nentry, bunch);
  ntbm->SetNumberOfFrostYCandidates(
    nentry, static_cast<int>(candidates.frost_y_candidates.size()));
  ntbm->SetNumberOfFrostXCandidates(
    nentry, static_cast<int>(candidates.frost_x_candidates.size()));
  ntbm->SetFrostYCandidates(nentry, candidates.frost_y_candidates);
  ntbm->SetFrostXCandidates(nentry, candidates.frost_x_candidates);
  ntbm->SetExpectedYCandidates(nentry, candidates.expected_y_candidates);
  ntbm->SetExpectedXCandidates(nentry, candidates.expected_x_candidates);
  ntbm->SetDifferenceYCandidates(nentry, candidates.difference_y_candidates);
  ntbm->SetDifferenceXCandidates(nentry, candidates.difference_x_candidates);
  ntbm->SetTangentYCandidates(nentry, candidates.tangent_y_candidates);
  ntbm->SetTangentXCandidates(nentry, candidates.tangent_x_candidates);
}

// Transfer B2Summary information

void TransferBeamInfo(const B2SpillSummary &spill_summary, NTBMSummary *ntbm_summary) {
  auto beam_summary = spill_summary.GetBeamSummary();
  ntbm_summary->SetSpillPot(beam_summary.GetSpillPot());
  for (std::size_t bunch = 0; bunch < 8; bunch++){
    ntbm_summary->SetBunchPot(bunch, beam_summary.GetBunchPot(bunch));
  }
    ntbm_summary->SetBsdSpillNumber(beam_summary.GetBsdSpillNumber());
    ntbm_summary->SetTimestamp(beam_summary.GetTimestamp());
    ntbm_summary->SetBsdGoodSpillFlag(beam_summary.GetBsdGoodSpillFlag());
    ntbm_summary->SetWagasciGoodSpillFlag(beam_summary.GetWagasciGoodSpillFlag());
  for (int i = 0; i < 8; i++) {
    ntbm_summary->SetDetectorFlags(i, beam_summary.GetDetectorFlags().at(i));
  }
}

double MyFuncCalculateTrackLength(const B2TrackSummary* track, double ax, double ay,
				  B2Dimension &dimension) {

  double track_length_ = 0;

  const double num_magnet[19] = {0,3,0,1,1,1,2,2,2,2,4,0,3,0,4,4,4,0,0};

  double posx[18] = {};
  double posy[18] = {};
  double posz[18] = {};
  int sizex[18] = {};
  int sizey[18] = {};

  std::vector<UInt_t > used_hit;

  auto it_cluster = track->BeginCluster();
  int debug_cluster_index = 0;
  while ( auto *cluster = it_cluster.Next() ) {
    BOOST_LOG_TRIVIAL(debug)
      << "MyFuncCalculateTrackLength: enter cluster "
      << debug_cluster_index
      << ", cluster ptr=" << static_cast<const void*>(cluster);

    auto it_hit = cluster->BeginHit();
    int debug_hit_index = 0;
    while ( auto *hit = it_hit.Next() ) {
      BOOST_LOG_TRIVIAL(debug)
        << "MyFuncCalculateTrackLength: enter hit "
        << debug_hit_index
        << " in cluster " << debug_cluster_index
        << ", hit ptr=" << static_cast<const void*>(hit);

      BOOST_LOG_TRIVIAL(debug)
        << "MyFuncCalculateTrackLength: about to call GetDetectorId()";
      if ( hit->GetDetectorId() != B2Detector::kBabyMind ) continue;

      BOOST_LOG_TRIVIAL(debug)
        << "MyFuncCalculateTrackLength: passed GetDetectorId(), about to call GetHitId()";
      if ( std::find(used_hit.begin(), used_hit.end(), hit->GetHitId())
	   != used_hit.end() ) continue;

      BOOST_LOG_TRIVIAL(debug)
        << "MyFuncCalculateTrackLength: passed GetHitId(), about to push used_hit";
      used_hit.push_back(hit->GetHitId());

      BOOST_LOG_TRIVIAL(debug)
        << "MyFuncCalculateTrackLength: about to call GetPlane()";
      int plane = hit->GetPlane();

      BOOST_LOG_TRIVIAL(debug)
        << "MyFuncCalculateTrackLength: plane=" << plane
        << ", about to call GetView()";

      const auto view = hit->GetView();

      BOOST_LOG_TRIVIAL(debug)
        << "MyFuncCalculateTrackLength: view=" << view
        << ", about to call GetSlot().GetValue(...)";
      TVector3 position;

      dimension.GetPosBm(view,
			  hit->GetPlane(),
			  hit->GetSlot().GetValue(hit->GetReadout1()),
			  position);
      BOOST_LOG_TRIVIAL(debug)
        << "MyFuncCalculateTrackLength: got position"
        << " x=" << position.X()
        << ", y=" << position.Y()
        << ", z=" << position.Z();

      if ( view == B2View::kSideView ) {
        posy[plane] += position.Y();
        posz[plane] += position.Z();
        sizey[plane]++;
      }
      else if ( view == B2View::kTopView ) {
        posx[plane] += position.X();
        posz[plane] += position.Z();
        sizex[plane]++;
      }
      ++debug_hit_index;
    }
    ++debug_cluster_index;
  }

  std::vector<double > posx_mod;
  std::vector<double > posy_mod;
  std::vector<double > posz_mod;
  std::vector<int > plane_vec;

  for ( int iplane = 0; iplane < 18; iplane++ ) {
    if ( sizex[iplane] < 1 ||
	 sizey[iplane] < 1 ) continue;
    posx[iplane] /= sizex[iplane];
    posy[iplane] /= sizey[iplane];
    posz[iplane] /= (sizex[iplane] + sizey[iplane]);

    posx_mod.push_back(posx[iplane]);
    posy_mod.push_back(posy[iplane]);
    posz_mod.push_back(posz[iplane]);
    plane_vec.push_back(iplane);
  }

  for ( int i = 0; i < plane_vec.size(); i++ ) {

    auto plane = plane_vec.at(i);

    double ax_ = ax;
    double ay_;

    int num_iron = 0;

    if ( plane == plane_vec.front() ) {

      ay_ = ay;

      for ( int iplane = 0; iplane <= plane; iplane++ ) {
	num_iron += num_magnet[iplane];
      }

      track_length_ += num_iron * 3. * 7.86 * std::sqrt(ax_ * ax_ + ay_ * ay_ + 1);
      track_length_ += (sizex[plane] + sizey[plane]) * 0.75 * 1. * std::sqrt(ax_ * ax_ + ay_ * ay_ + 1.) * 1.289;

    }
    else if ( plane != plane_vec.back() ) {

      ay_ = (posy_mod.at(i) - posy_mod.at(i-1)) / (posz_mod.at(i) - posz_mod.at(i-1));

      int startplane = plane_vec.at(i-1);
      for ( int iplane = startplane+1; iplane <= plane; iplane++ ) {
	num_iron += num_magnet[iplane];
      }

      track_length_ += num_iron * 3. * 7.86 * std::sqrt(ax_ * ax_ + ay_ * ay_ + 1);
      track_length_ += (sizex[plane] + sizey[plane]) * 0.75 * 1. * std::sqrt(ax_ * ax_ + ay_ * ay_ + 1.) * 1.289;

    }
    else {

      ay_ = (posy_mod.at(i) - posy_mod.at(i-1)) / (posz_mod.at(i) - posz_mod.at(i-1));

      int startplane = plane_vec.at(i-1);
      for ( int iplane = startplane+1; iplane <= plane; iplane++ ) {
	num_iron += num_magnet[iplane];
      }
      //std::cout << "Num iron at the last plane : " << num_iron + num_magnet[plane+1] * 0.5 << std::endl;
      track_length_ += (num_iron + num_magnet[plane+1] * 0.5) * 3. * 7.86 * std::sqrt(ax_ * ax_ + ay_ * ay_ + 1.);
      track_length_ += (sizex[plane] + sizey[plane]) * 0.75 * 1. * std::sqrt(ax_ * ax_ + ay_ * ay_ + 1) * 1.289;

    }

    //std::cout << "Plane : " << plane << ", # of iron" << num_iron << std::endl;
    //std::cout << "Direction : " << "(" << ax_ << ", " << ay_ << ")" << std::endl;
    //std::cout << "track length : " << track_length_ << std::endl;
  }

  return track_length_;

}

void TransferBabyMindTrackInfo(const B2SpillSummary &spill_summary,
			       NTBMSummary *ntbm_summary,
			       int datatype,
			       B2Dimension &dimension,
			       BeamMode mc_beam_mode,
			       std::vector<ExternalTrackFitSet> *external_fit_results) {

  int itrack = 0;
  int bsd_good_spill_flag = B2_NON_INITIALIZED_VALUE;
  if (datatype == B2DataType::kMonteCarlo) {
    // MC input does not carry a meaningful BSD FHC/RHC flag.
    // Use the explicit command-line beam mode instead.
    bsd_good_spill_flag = BeamModeToBsdGoodSpillFlag(mc_beam_mode);
  } else {
    bsd_good_spill_flag =
      spill_summary.GetBeamSummary().GetBsdGoodSpillFlag();
  }

  auto it_recon_vertex = spill_summary.BeginReconVertex();
  while ( auto *vertex = it_recon_vertex.Next() ) {
    auto it_outgoing_track = vertex->BeginTrack();
    while ( auto *track = it_outgoing_track.Next() ) {
      if ( track->GetTrackType() == B2TrackType::kPrimaryTrack ) {
	if ( track->GetPrimaryTrackType() == B2PrimaryTrackType::kBabyMind3DTrack ) {
	  ntbm_summary->SetNinjaTrackType(itrack, 0); // ECC interaction candidate (or sand muon)
	} else if ( track->GetPrimaryTrackType() == B2PrimaryTrackType::kMatchingTrack &&
		    track->HasDetector(B2Detector::kBabyMind) ) {
	  auto vertex = track->GetParentVertex();
	  auto position = vertex.GetRelativePosition().GetValue();
	  if ( vertex.GetDetector() == B2Detector::kProtonModule ) {
	    if ( std::fabs(position.X()) < 500. &&
		 std::fabs(position.Y()) < 500. &&
		 std::fabs(position.Z()) < 300. ) {
	      ntbm_summary->SetNinjaTrackType(itrack, 2);
	    }
	    else {
	      ntbm_summary->SetNinjaTrackType(itrack, 1);
	    }
	  }
	  else if ( vertex.GetDetector() == B2Detector::kWagasciUpstream ) {
	    if ( std::fabs(position.X()) < 425. &&
		 std::fabs(position.Y()) < 425. &&
		 -155. < position.Z() && position.Z() < 60. ) {
	      ntbm_summary->SetNinjaTrackType(itrack, 2);
	    }
	    else {
	      ntbm_summary->SetNinjaTrackType(itrack, 1);
	    }
	  }
	  else {
	    ntbm_summary->SetNinjaTrackType(itrack, 0);
	  }
	}
	else continue;
      } else { // not primary track
	continue;
      }
      BOOST_LOG_TRIVIAL(debug)
        << "Processing BM track: itrack=" << itrack
        << ", track type=" << track->GetTrackType()
        << ", primary type=" << track->GetPrimaryTrackType()
        << ", bunch=" << track->GetBunch();

      const auto &downstream_hit = track->GetDownstreamHit();
      BOOST_LOG_TRIVIAL(debug)
        << "About to access downstream hit plane";
      // ntbm_summary->SetBabyMindMaximumPlane(itrack, track->GetDownstreamHit().GetPlane());
      int baby_mind_max_plane = -1;

      auto it_cluster_for_plane = track->BeginCluster();
      while (auto *cluster_for_plane = it_cluster_for_plane.Next()) {
        auto it_hit_for_plane = cluster_for_plane->BeginHit();
        while (auto *hit_for_plane = it_hit_for_plane.Next()) {
          if (hit_for_plane->GetDetectorId() != B2Detector::kBabyMind) continue;
          baby_mind_max_plane = std::max(baby_mind_max_plane, hit_for_plane->GetPlane());
        }
      }

      if (baby_mind_max_plane < 0) {
        BOOST_LOG_TRIVIAL(warning)
          << "Skip BM track because no valid Baby MIND hit was found"
          << " : itrack=" << itrack
          << ", primary type=" << track->GetPrimaryTrackType()
          << ", bunch=" << track->GetBunch();
        continue;
      }

      ntbm_summary->SetBabyMindMaximumPlane(itrack, baby_mind_max_plane);

      ntbm_summary->SetTrackLengthTotal(itrack, track->GetTrackLengthTotal());
      double nll_plus = track->GetNegativeLogLikelihoodPlus();
      double nll_minus = track->GetNegativeLogLikelihoodMinus();
      ntbm_summary->SetNegativeLogLikelihoodPlus(itrack, nll_plus);
      ntbm_summary->SetNegativeLogLikelihoodMinus(itrack, nll_minus);

      if (bsd_good_spill_flag > 0) {
        // FHC spill.
        if ( nll_minus - nll_plus >= 4 )
          ntbm_summary->SetCharge(itrack, 1);
        else
          ntbm_summary->SetCharge(itrack, -1);
      } else if (bsd_good_spill_flag < 0) {
        // RHC spill.
        if ( nll_plus - nll_minus >= 4 ) //temporal threshold, to be optimized using MC!!!
          ntbm_summary->SetCharge(itrack, -1);
        else
          ntbm_summary->SetCharge(itrack, 1);
      } else {
        // Beam mode is unknown for bad/undefined spill.
        ntbm_summary->SetCharge(itrack, B2_NON_INITIALIZED_VALUE);
      }

      ntbm_summary->SetBunch(itrack, track->GetBunch());

      TVector3 final_position = track->GetFinalPosition().GetValue();
      if ( std::fabs(final_position.X()) < 1100. &&
	   std::fabs(final_position.Y()) < 900. &&
	   final_position.Z() < 1500. ) {
	//      if ( track->GetIsStopping() )
	ntbm_summary->SetMomentumType(itrack, 0); // Baby MIND range method
      }
      else {
	ntbm_summary->SetMomentumType(itrack, 1); // should be curvature type but not yet implemented
      }
      ntbm_summary->SetMomentum(itrack, track->GetReconMomByRange());
      ntbm_summary->SetMomentumError(itrack, track->GetReconMomByCurve());
      std::vector<Double_t> direction_and_position = GetBabyMindInitialDirectionAndPosition(track, datatype);
      for (int view = 0; view < 2; view++) {
	ntbm_summary->SetBabyMindPosition(itrack, view, direction_and_position.at(view+2));
	ntbm_summary->SetBabyMindTangent(itrack, view, direction_and_position.at(view));
      }

      double track_length = MyFuncCalculateTrackLength(track,
						     direction_and_position.at(1),
						     direction_and_position.at(0),
						     dimension);

      //std::cout << track->GetTrackLengthTotal() << ", " << track_length << std::endl;

      ntbm_summary->SetTrackLengthTotal(itrack, track_length);
      // Use the original B2 track length here.
      // Some tracks contain invalid Baby MIND hit summaries that can crash
      // MyFuncCalculateTrackLength() when accessing hit->GetPlane().
      // ntbm_summary->SetTrackLengthTotal(itrack, track->GetTrackLengthTotal());

      if (external_fit_results) {
        ExternalTrackFitSet external_fit_set;
        external_fit_set.combined = FitExternalTrackToFrost(
          track,
          datatype,
          dimension,
          ExternalFitDetectorSelection::kProtonModuleAndDownstreamWagasci);
        external_fit_set.proton_module_only = FitExternalTrackToFrost(
          track,
          datatype,
          dimension,
          ExternalFitDetectorSelection::kProtonModuleOnly);
        external_fit_set.downstream_wagasci_only = FitExternalTrackToFrost(
          track,
          datatype,
          dimension,
          ExternalFitDetectorSelection::kDownstreamWagasciOnly);
        external_fit_results->push_back(external_fit_set);
      }

      itrack++;

    } // while track
  } // while vertex

}

void TransferMCInfo(const B2SpillSummary &spill_summary, NTBMSummary *ntbm_summary) {
  auto it_event = spill_summary.BeginTrueEvent();
  const auto *event = it_event.Next();

  if (!event) {
    BOOST_LOG_TRIVIAL(warning)
      << "No true event found in this spill. "
      << "Skip MC normalization and total cross section transfer.";
    ntbm_summary->SetNormalization(B2_NON_INITIALIZED_VALUE);
    ntbm_summary->SetTotalCrossSection(B2_NON_INITIALIZED_VALUE);
    return;
  }

  ntbm_summary->SetNormalization(event->GetNormalization());
  ntbm_summary->SetTotalCrossSection(event->GetTotalCrossSection());
}

// main

int main(int argc, char *argv[]) {

  gErrorIgnoreLevel = kError;
  gRandom->SetSeed(1);

  B2Dimension dimension_((std::string)"/opt/wagasci_mc/WagasciMC/etc/wagasci/b2/geometry");

  logging::trivial::severity_level log_level = logging::trivial::info;

  BOOST_LOG_TRIVIAL(info) << "==========FROST-Baby MIND Track Matching Start==========";

  if (argc < 5 || argc > 7) {
    PrintUsage(argv[0]);
    std::exit(1);
  }

  try {
    double z_shift = std::stof(argv[3]);
    int datatype = std::stoi(argv[4]);
    BeamMode mc_beam_mode = BeamMode::kAuto;

    if (datatype != B2DataType::kMonteCarlo && datatype != B2DataType::kRealData) {
      throw std::invalid_argument("Invalid data type. Use 0 for MC and 1 for data.");
    }

    for (int iarg = 5; iarg < argc; ++iarg) {
      const std::string token = argv[iarg];
      if (IsBeamModeToken(token)) {
        if (mc_beam_mode != BeamMode::kAuto) {
          throw std::invalid_argument("Beam mode was specified more than once");
        }
        mc_beam_mode = ParseBeamMode(token);
      } else {
        log_level = ParseLogLevel(token);
      }
    }

    if (datatype == B2DataType::kMonteCarlo &&
        mc_beam_mode == BeamMode::kAuto) {
      PrintUsage(argv[0]);
      throw std::invalid_argument(
        "MC input requires an explicit beam mode argument: fhc or rhc");
    }

    if (datatype == B2DataType::kRealData &&
        mc_beam_mode != BeamMode::kAuto) {
      BOOST_LOG_TRIVIAL(warning)
        << "Beam mode argument is ignored for real data. "
        << "BsdGoodSpillFlag from B2BeamSummary is used instead.";
    }

    logging::core::get()->set_filter(
        logging::trivial::severity >= log_level);

    B2Reader reader(argv[1]);
    std::unique_ptr<TFile> input_root(new TFile(argv[1], "READ"));
    if (!input_root || input_root->IsZombie()) {
      throw std::runtime_error("Cannot open input ROOT file");
    }
    FrostInputTrees frost_trees = OpenFrostInputTrees(input_root.get(), datatype);

    std::vector<std::vector<double>> *x_rec = nullptr;
    std::vector<std::vector<double>> *y_rec = nullptr;
    std::vector<int> *is_hit = nullptr;
    SetRequiredBranchAddress(frost_trees.frost_input, "x_rec", &x_rec);
    SetRequiredBranchAddress(frost_trees.frost_input, "y_rec", &y_rec);
    SetRequiredBranchAddress(frost_trees.frost_input, "is_hit", &is_hit);

    FrostMcEntryData frostmc_entry;
    if (datatype == B2DataType::kMonteCarlo) {
      if (!frost_trees.frostmc) {
        throw std::runtime_error(
          "MC input file must contain frostmc tree to fill true FROST particles");
      }

      SetRequiredBranchAddress(
        frost_trees.frostmc, "particle_pdg_id",
        &frostmc_entry.particle_pdg_id);
      SetRequiredBranchAddress(
        frost_trees.frostmc, "particle_local_first_x_mm",
        &frostmc_entry.particle_local_first_x_mm);
      SetRequiredBranchAddress(
        frost_trees.frostmc, "particle_local_first_y_mm",
        &frostmc_entry.particle_local_first_y_mm);
      SetRequiredBranchAddress(
        frost_trees.frostmc, "particle_local_first_z_mm",
        &frostmc_entry.particle_local_first_z_mm);
      SetRequiredBranchAddress(
        frost_trees.frostmc, "particle_local_last_x_mm",
        &frostmc_entry.particle_local_last_x_mm);
      SetRequiredBranchAddress(
        frost_trees.frostmc, "particle_local_last_y_mm",
        &frostmc_entry.particle_local_last_y_mm);
      SetRequiredBranchAddress(
        frost_trees.frostmc, "particle_local_last_z_mm",
        &frostmc_entry.particle_local_last_z_mm);
      SetRequiredBranchAddress(
        frost_trees.frostmc, "particle_initial_px_mev_c",
        &frostmc_entry.particle_initial_px_mev_c);
      SetRequiredBranchAddress(
        frost_trees.frostmc, "particle_initial_py_mev_c",
        &frostmc_entry.particle_initial_py_mev_c);
      SetRequiredBranchAddress(
        frost_trees.frostmc, "particle_initial_pz_mev_c",
        &frostmc_entry.particle_initial_pz_mev_c);
    }

    TTree *ntbm_tree = new TTree("ntbm", "NINJA BabyMIND Original Summary");
    ntbm_tree->SetDirectory(nullptr);
    NTBMSummary* my_ntbm = new NTBMSummary();
    ntbm_tree->Branch("NTBMSummary", &my_ntbm);

    int nspill = 0;
    std::vector<std::vector<TrackMatchRow>> spill_match_rows;

    //count matched tracks
    int total_tracks = 0;
    int matched_tracks = 0;

    {
      B2Writer writer(argv[2], reader);

      while ( reader.ReadNextSpill() > 0) {

      my_ntbm->SetEntryInDailyFile(reader.GetEntryNumber());

      auto &input_spill_summary = writer.GetSpillSummary();
      BOOST_LOG_TRIVIAL(debug) << "entry : " << reader.GetEntryNumber();
      BOOST_LOG_TRIVIAL(debug) << "timestamp : "
                               << input_spill_summary.GetBeamSummary().GetTimestamp();

      TransferBeamInfo(input_spill_summary, my_ntbm);
      if (datatype == B2DataType::kMonteCarlo) {
        TransferMCInfo(input_spill_summary, my_ntbm);
      }

      FrostEntryData frost_entry;
      if (nspill < frost_trees.frost_input->GetEntries()) {
        frost_trees.frost_input->GetEntry(nspill);
        frost_entry.x_rec = x_rec;
        frost_entry.y_rec = y_rec;
      }

      bool has_frostmc_entry = false;
      if (datatype == B2DataType::kMonteCarlo) {
        if (nspill < frost_trees.frostmc->GetEntries()) {
          frost_trees.frostmc->GetEntry(nspill);
          has_frostmc_entry = true;
        } else {
          BOOST_LOG_TRIVIAL(warning)
            << "No frostmc entry for spill " << nspill;
        }
      }

      TrueFrostParticleInfo true_frost_particles;
      if (datatype == B2DataType::kMonteCarlo &&
          has_frostmc_entry) {
        true_frost_particles = CollectTrueFrostParticleInfo(frostmc_entry);
        FillTrueFrostParticleInfo(true_frost_particles, my_ntbm);
      }

      // Collect all BM 3d tracks
      int number_of_tracks = 0;

      auto it_recon_vertex = input_spill_summary.BeginReconVertex();
      while ( auto *vertex = it_recon_vertex.Next() ) {
	      auto it_outgoing_track = vertex->BeginTrack();
        while ( auto *track = it_outgoing_track.Next() ) {
          if ( track->GetTrackType() == B2TrackType::kPrimaryTrack ) {
            // not start from the other WAGASCI modules
            if ( track->GetPrimaryTrackType() == B2PrimaryTrackType::kBabyMind3DTrack ) {
              number_of_tracks++;
            // start from the other modules and have hits in Baby MIND
            } else if ( track->GetPrimaryTrackType() == B2PrimaryTrackType::kMatchingTrack &&
              track->HasDetector(B2Detector::kBabyMind) ) {
                if ( track->HasDetector(B2Detector::kProtonModule) ||
                  track->HasDetector(B2Detector::kWagasciUpstream) ||
                  track->HasDetector(B2Detector::kWagasciDownstream) ) {
                number_of_tracks++;
                }
            }
          }
        } // track
      } // vertex

      my_ntbm->SetNumberOfTracks(number_of_tracks);

      // Extrapolate BabyMIND tracks to the NINJA FROST position
      // and get the positions to match each BabyMIND track
      if ( number_of_tracks > 0 ) {
        std::vector<ExternalTrackFitSet> external_fit_results;

	      TransferBabyMindTrackInfo(
          input_spill_summary, my_ntbm, datatype, dimension_, mc_beam_mode, &external_fit_results);

        std::vector<TrackMatchRow> rows;

        for ( int ibmtrack = 0; ibmtrack < my_ntbm->GetNumberOfTracks(); ibmtrack++ ) {
          TrackMatchRow row;
          row.bm_track_id = ibmtrack;
          row.ninja_track_type = my_ntbm->GetNinjaTrackType(ibmtrack);

          const int bm_bunch = my_ntbm->GetBunch(ibmtrack);
          row.frost_match_bunch =
            (datatype == B2DataType::kMonteCarlo) ? 0 : bm_bunch;
          row.frost_is_hit =
            GetFrostHitValueForTrack(frost_entry, is_hit, bm_bunch, datatype);

          row.baby_mind_tangent_y = my_ntbm->GetBabyMindTangent(ibmtrack, B2View::kSideView);
          row.baby_mind_tangent_x = my_ntbm->GetBabyMindTangent(ibmtrack, B2View::kTopView);


          const std::vector<double> expected_position =
            CalculateExpectedPosition(my_ntbm, ibmtrack, z_shift);
          my_ntbm->SetExtrapolatedPosition(ibmtrack, expected_position);
          row.expected_y = expected_position.at(B2View::kSideView);
          row.expected_x = expected_position.at(B2View::kTopView);

          if (ibmtrack < static_cast<int>(external_fit_results.size())) {
            const ExternalTrackFitSet &external_fit_set =
              external_fit_results.at(ibmtrack);
            const ExternalTrackFitResult &external_fit =
              external_fit_set.combined;

            FillRowExternalFitResult(
              external_fit,
              row.external_expected_x,
              row.external_expected_y,
              row.external_tangent_x,
              row.external_tangent_y,
              row.external_chi2_x,
              row.external_chi2_y,
              row.external_ndof_x,
              row.external_ndof_y);

            row.external_num_planes_proton_module_x =
              external_fit.num_planes_proton_module_x;
            row.external_num_planes_proton_module_y =
              external_fit.num_planes_proton_module_y;
            row.external_num_planes_downstream_wagasci_x =
              external_fit.num_planes_downstream_wagasci_x;
            row.external_num_planes_downstream_wagasci_y =
              external_fit.num_planes_downstream_wagasci_y;

            FillRowExternalFitResult(
              external_fit_set.proton_module_only,
              row.external_pm_only_expected_x,
              row.external_pm_only_expected_y,
              row.external_pm_only_tangent_x,
              row.external_pm_only_tangent_y,
              row.external_pm_only_chi2_x,
              row.external_pm_only_chi2_y,
              row.external_pm_only_ndof_x,
              row.external_pm_only_ndof_y);

            FillRowExternalFitResult(
              external_fit_set.downstream_wagasci_only,
              row.external_dwg_only_expected_x,
              row.external_dwg_only_expected_y,
              row.external_dwg_only_tangent_x,
              row.external_dwg_only_tangent_y,
              row.external_dwg_only_chi2_x,
              row.external_dwg_only_chi2_y,
              row.external_dwg_only_ndof_x,
              row.external_dwg_only_ndof_y);
          }

          if (datatype == B2DataType::kMonteCarlo &&
              has_frostmc_entry) {
            const NearestTrueFrostParticleResult nearest_true_frost =
              FindNearestTrueFrostParticle(
                true_frost_particles,
                row.expected_x,
                row.expected_y);

            if (nearest_true_frost.found) {
              row.true_frost_nearest_particle_id =
                nearest_true_frost.particle_id;
              row.true_frost_nearest_position_x =
                nearest_true_frost.position_x;
              row.true_frost_nearest_position_y =
                nearest_true_frost.position_y;
              row.true_frost_nearest_tangent_x =
                nearest_true_frost.tangent_x;
              row.true_frost_nearest_tangent_y =
                nearest_true_frost.tangent_y;
            }
          }

          const NearestFrostPositionResult nearest_y =
            FindNearestFrostPosition(frost_entry, bm_bunch, B2View::kSideView,
                                     row.expected_y, is_hit, datatype);
          const NearestFrostPositionResult nearest_x =
            FindNearestFrostPosition(frost_entry, bm_bunch, B2View::kTopView,
                                     row.expected_x, is_hit, datatype);

          const double dz_y =
            BABYMIND_POS_Z + BM_SECOND_LAYER_POS - NINJA_POS_Z_FROST - NINJA_FROST_POS_Z
            - (2 * B2View::kSideView - 1) * 10. + z_shift;
          const double dz_x =
            BABYMIND_POS_Z + BM_SECOND_LAYER_POS - NINJA_POS_Z_FROST - NINJA_FROST_POS_Z
            - (2 * B2View::kTopView - 1) * 10. + z_shift;

          const std::vector<double> bm2_position_frost =
            CalculateBabyMindSecondLayerPositionInFrostCoordinate(my_ntbm, ibmtrack);

          if (nearest_x.found) {
            row.frost_x = nearest_x.frost_position;
            row.dx = nearest_x.diff;
            row.tangent_x =
              (bm2_position_frost.at(B2View::kTopView) - row.frost_x) / dz_x;
            row.dtanx = row.tangent_x - row.baby_mind_tangent_x;
          }

          if (nearest_y.found) {
            row.frost_y = nearest_y.frost_position;
            row.dy = nearest_y.diff;
            row.tangent_y =
              (bm2_position_frost.at(B2View::kSideView) - row.frost_y) / dz_y;
            row.dtany = row.tangent_y - row.baby_mind_tangent_y;
          }

          const bool hit_expected = NinjaHitExpected(my_ntbm, ibmtrack, z_shift);
          if (hit_expected) {
            total_tracks++;
          }

          const FrostTrackCandidates candidates =
            CollectFrostTrackCandidates(my_ntbm, ibmtrack, frost_entry, z_shift, is_hit, datatype);
          const FrostMatchCount match_count =
            CountAcceptedFrostMatches(candidates);

          AppendFrostTrackCandidates(
            my_ntbm, ibmtrack, row.frost_match_bunch, candidates);

          if (match_count.x_match_count > 0 && match_count.y_match_count > 0) {
            row.has_match = 1;
            if (hit_expected) {
              matched_tracks++;
            }
          } else {
            row.has_match = 0;
          }

          rows.push_back(row);
        } // ibmtrack

        spill_match_rows.push_back(rows);
      } else {
        spill_match_rows.push_back(std::vector<TrackMatchRow>{});
      }

      // Create output tree
      BOOST_LOG_TRIVIAL(debug) << *my_ntbm;
      ntbm_tree->Fill();

      my_ntbm->Clear("C");
        writer.Fill();
        nspill++;
      }
    } // writer is closed here

    TFile *output_b2_file = new TFile(argv[2], "UPDATE");
    if (!output_b2_file || output_b2_file->IsZombie()) {
      throw std::runtime_error(std::string("Failed to reopen output ROOT file for UPDATE: ") + argv[2]);
    }
    output_b2_file->cd();
    ntbm_tree->SetDirectory(output_b2_file);
    ntbm_tree->Write("", TObject::kOverwrite);

    TTree *frost_match_out = frost_trees.frost_input->CloneTree(-1);
    // For MC, clone input frost as frost_match to keep the output schema
    // compatible with the data path and downstream jobs.
    frost_match_out->SetName("frost_match");
    frost_match_out->SetDirectory(output_b2_file);

    TTree *match_info_out = nullptr;
    if (frost_trees.match_info) {
      match_info_out = frost_trees.match_info->CloneTree(0);
    } else {
      match_info_out = new TTree("match_info", "FROST-BabyMIND track match info");
    }
    match_info_out->SetName("match_info");
    match_info_out->SetDirectory(output_b2_file);

    std::vector<Int_t> trackmatch_has_match;
    std::vector<Int_t> trackmatch_bm_track_id;
    std::vector<Int_t> trackmatch_frost_match_bunch;
    std::vector<Int_t> trackmatch_ninja_track_type;
    std::vector<Double_t> trackmatch_baby_mind_tangent_x;
    std::vector<Double_t> trackmatch_baby_mind_tangent_y;
    std::vector<Double_t> trackmatch_expected_x;
    std::vector<Double_t> trackmatch_expected_y;
    std::vector<Double_t> trackmatch_frost_nearest_x;
    std::vector<Double_t> trackmatch_frost_nearest_y;
    std::vector<Double_t> trackmatch_dx;
    std::vector<Double_t> trackmatch_dy;
    std::vector<Double_t> trackmatch_tangent_x;
    std::vector<Double_t> trackmatch_tangent_y;
    std::vector<Double_t> trackmatch_dtanx;
    std::vector<Double_t> trackmatch_dtany;
    std::vector<Int_t> trackmatch_frost_is_hit;
    std::vector<Double_t> trackmatch_external_expected_x;
    std::vector<Double_t> trackmatch_external_expected_y;
    std::vector<Double_t> trackmatch_external_tangent_x;
    std::vector<Double_t> trackmatch_external_tangent_y;
    std::vector<Double_t> trackmatch_external_chi2_x;
    std::vector<Double_t> trackmatch_external_chi2_y;
    std::vector<Int_t> trackmatch_external_ndof_x;
    std::vector<Int_t> trackmatch_external_ndof_y;
    std::vector<Int_t> trackmatch_external_num_planes_proton_module_x;
    std::vector<Int_t> trackmatch_external_num_planes_proton_module_y;
    std::vector<Int_t> trackmatch_external_num_planes_downstream_wagasci_x;
    std::vector<Int_t> trackmatch_external_num_planes_downstream_wagasci_y;
    std::vector<Double_t> trackmatch_external_pm_only_expected_x;
    std::vector<Double_t> trackmatch_external_pm_only_expected_y;
    std::vector<Double_t> trackmatch_external_pm_only_tangent_x;
    std::vector<Double_t> trackmatch_external_pm_only_tangent_y;
    std::vector<Double_t> trackmatch_external_pm_only_chi2_x;
    std::vector<Double_t> trackmatch_external_pm_only_chi2_y;
    std::vector<Int_t> trackmatch_external_pm_only_ndof_x;
    std::vector<Int_t> trackmatch_external_pm_only_ndof_y;
    std::vector<Double_t> trackmatch_external_dwg_only_expected_x;
    std::vector<Double_t> trackmatch_external_dwg_only_expected_y;
    std::vector<Double_t> trackmatch_external_dwg_only_tangent_x;
    std::vector<Double_t> trackmatch_external_dwg_only_tangent_y;
    std::vector<Double_t> trackmatch_external_dwg_only_chi2_x;
    std::vector<Double_t> trackmatch_external_dwg_only_chi2_y;
    std::vector<Int_t> trackmatch_external_dwg_only_ndof_x;
    std::vector<Int_t> trackmatch_external_dwg_only_ndof_y;
    Double_t matchinfo_spill_pot = B2_NON_INITIALIZED_VALUE;
    Double_t matchinfo_bunch_pot[NUMBER_OF_BUNCHES];
    Int_t matchinfo_bsd_spill_number = B2_NON_INITIALIZED_VALUE;
    Double_t matchinfo_timestamp = B2_NON_INITIALIZED_VALUE;
    Int_t matchinfo_bsd_good_spill_flag = B2_NON_INITIALIZED_VALUE;
    Int_t matchinfo_wagasci_good_spill_flag = B2_NON_INITIALIZED_VALUE;
    Int_t matchinfo_detector_flags[8];
    std::vector<Int_t> trackmatch_true_frost_nearest_particle_id;
    std::vector<Double_t> trackmatch_true_frost_nearest_position_x;
    std::vector<Double_t> trackmatch_true_frost_nearest_position_y;
    std::vector<Double_t> trackmatch_true_frost_nearest_tangent_x;
    std::vector<Double_t> trackmatch_true_frost_nearest_tangent_y;

    for (int i = 0; i < NUMBER_OF_BUNCHES; ++i) {
      matchinfo_bunch_pot[i] = B2_NON_INITIALIZED_VALUE;
    }
    for (int i = 0; i < 8; ++i) {
      matchinfo_detector_flags[i] = B2_NON_INITIALIZED_VALUE;
    }

    match_info_out->Branch("trackmatch_has_match", &trackmatch_has_match);
    match_info_out->Branch("trackmatch_bm_track_id", &trackmatch_bm_track_id);
    match_info_out->Branch("trackmatch_frost_match_bunch", &trackmatch_frost_match_bunch);
    match_info_out->Branch("trackmatch_ninja_track_type", &trackmatch_ninja_track_type);
    match_info_out->Branch("trackmatch_baby_mind_tangent_x", &trackmatch_baby_mind_tangent_x);
    match_info_out->Branch("trackmatch_baby_mind_tangent_y", &trackmatch_baby_mind_tangent_y);
    match_info_out->Branch("trackmatch_expected_x", &trackmatch_expected_x);
    match_info_out->Branch("trackmatch_expected_y", &trackmatch_expected_y);
    match_info_out->Branch("trackmatch_frost_nearest_x", &trackmatch_frost_nearest_x);
    match_info_out->Branch("trackmatch_frost_nearest_y", &trackmatch_frost_nearest_y);
    match_info_out->Branch("trackmatch_dx", &trackmatch_dx);
    match_info_out->Branch("trackmatch_dy", &trackmatch_dy);
    match_info_out->Branch("trackmatch_tangent_x", &trackmatch_tangent_x);
    match_info_out->Branch("trackmatch_tangent_y", &trackmatch_tangent_y);
    match_info_out->Branch("trackmatch_dtanx", &trackmatch_dtanx);
    match_info_out->Branch("trackmatch_dtany", &trackmatch_dtany);
    match_info_out->Branch("trackmatch_frost_is_hit", &trackmatch_frost_is_hit);
    match_info_out->Branch("trackmatch_external_expected_x",
                           &trackmatch_external_expected_x);
    match_info_out->Branch("trackmatch_external_expected_y",
                           &trackmatch_external_expected_y);
    match_info_out->Branch("trackmatch_external_tangent_x",
                           &trackmatch_external_tangent_x);
    match_info_out->Branch("trackmatch_external_tangent_y",
                           &trackmatch_external_tangent_y);
    match_info_out->Branch("trackmatch_external_chi2_x",
                           &trackmatch_external_chi2_x);
    match_info_out->Branch("trackmatch_external_chi2_y",
                           &trackmatch_external_chi2_y);
    match_info_out->Branch("trackmatch_external_ndof_x",
                           &trackmatch_external_ndof_x);
    match_info_out->Branch("trackmatch_external_ndof_y",
                           &trackmatch_external_ndof_y);
    match_info_out->Branch("trackmatch_external_num_planes_proton_module_x",
                           &trackmatch_external_num_planes_proton_module_x);
    match_info_out->Branch("trackmatch_external_num_planes_proton_module_y",
                           &trackmatch_external_num_planes_proton_module_y);
    match_info_out->Branch("trackmatch_external_num_planes_downstream_wagasci_x",
                           &trackmatch_external_num_planes_downstream_wagasci_x);
    match_info_out->Branch("trackmatch_external_num_planes_downstream_wagasci_y",
                           &trackmatch_external_num_planes_downstream_wagasci_y);
    match_info_out->Branch("trackmatch_external_pm_only_expected_x",
                           &trackmatch_external_pm_only_expected_x);
    match_info_out->Branch("trackmatch_external_pm_only_expected_y",
                           &trackmatch_external_pm_only_expected_y);
    match_info_out->Branch("trackmatch_external_pm_only_tangent_x",
                           &trackmatch_external_pm_only_tangent_x);
    match_info_out->Branch("trackmatch_external_pm_only_tangent_y",
                           &trackmatch_external_pm_only_tangent_y);
    match_info_out->Branch("trackmatch_external_pm_only_chi2_x",
                           &trackmatch_external_pm_only_chi2_x);
    match_info_out->Branch("trackmatch_external_pm_only_chi2_y",
                           &trackmatch_external_pm_only_chi2_y);
    match_info_out->Branch("trackmatch_external_pm_only_ndof_x",
                           &trackmatch_external_pm_only_ndof_x);
    match_info_out->Branch("trackmatch_external_pm_only_ndof_y",
                           &trackmatch_external_pm_only_ndof_y);
    match_info_out->Branch("trackmatch_external_dwg_only_expected_x",
                           &trackmatch_external_dwg_only_expected_x);
    match_info_out->Branch("trackmatch_external_dwg_only_expected_y",
                           &trackmatch_external_dwg_only_expected_y);
    match_info_out->Branch("trackmatch_external_dwg_only_tangent_x",
                           &trackmatch_external_dwg_only_tangent_x);
    match_info_out->Branch("trackmatch_external_dwg_only_tangent_y",
                           &trackmatch_external_dwg_only_tangent_y);
    match_info_out->Branch("trackmatch_external_dwg_only_chi2_x",
                           &trackmatch_external_dwg_only_chi2_x);
    match_info_out->Branch("trackmatch_external_dwg_only_chi2_y",
                           &trackmatch_external_dwg_only_chi2_y);
    match_info_out->Branch("trackmatch_external_dwg_only_ndof_x",
                           &trackmatch_external_dwg_only_ndof_x);
    match_info_out->Branch("trackmatch_external_dwg_only_ndof_y",
                           &trackmatch_external_dwg_only_ndof_y);
    match_info_out->Branch("spill_pot", &matchinfo_spill_pot, "spill_pot/D");
    match_info_out->Branch("bunch_pot", matchinfo_bunch_pot,
                           Form("bunch_pot[%d]/D", NUMBER_OF_BUNCHES));
    match_info_out->Branch("bsd_spill_number", &matchinfo_bsd_spill_number,
                           "bsd_spill_number/I");
    match_info_out->Branch("timestamp", &matchinfo_timestamp, "timestamp/D");
    match_info_out->Branch("bsd_good_spill_flag", &matchinfo_bsd_good_spill_flag,
                           "bsd_good_spill_flag/I");
    match_info_out->Branch("wagasci_good_spill_flag", &matchinfo_wagasci_good_spill_flag,
                           "wagasci_good_spill_flag/I");
    match_info_out->Branch("detector_flags", matchinfo_detector_flags,
                           "detector_flags[8]/I");
    match_info_out->Branch("trackmatch_true_frost_nearest_particle_id",
                           &trackmatch_true_frost_nearest_particle_id);
    match_info_out->Branch("trackmatch_true_frost_nearest_position_x",
                           &trackmatch_true_frost_nearest_position_x);
    match_info_out->Branch("trackmatch_true_frost_nearest_position_y",
                           &trackmatch_true_frost_nearest_position_y);
    match_info_out->Branch("trackmatch_true_frost_nearest_tangent_x",
                           &trackmatch_true_frost_nearest_tangent_x);
    match_info_out->Branch("trackmatch_true_frost_nearest_tangent_y",
                           &trackmatch_true_frost_nearest_tangent_y);

    const Long64_t ninfo = frost_trees.match_info
      ? frost_trees.match_info->GetEntries()
      : ntbm_tree->GetEntries();
    const std::vector<TrackMatchRow> empty_rows;
    for (Long64_t i = 0; i < ninfo; ++i) {
      if (frost_trees.match_info) {
        frost_trees.match_info->GetEntry(i);
      }
      if (i < frost_trees.frost_input->GetEntries()) {
        frost_trees.frost_input->GetEntry(i);
      } else {
        is_hit = nullptr;
      }

      const std::vector<TrackMatchRow> *rows_ptr = &empty_rows;
      if (i < static_cast<Long64_t>(spill_match_rows.size())) {
        rows_ptr = &spill_match_rows.at(i);
      }
      const std::vector<TrackMatchRow> &rows = *rows_ptr;

      trackmatch_has_match.clear();
      trackmatch_bm_track_id.clear();
      trackmatch_frost_match_bunch.clear();
      trackmatch_ninja_track_type.clear();
      trackmatch_baby_mind_tangent_x.clear();
      trackmatch_baby_mind_tangent_y.clear();
      trackmatch_expected_x.clear();
      trackmatch_expected_y.clear();
      trackmatch_frost_nearest_x.clear();
      trackmatch_frost_nearest_y.clear();
      trackmatch_dx.clear();
      trackmatch_dy.clear();
      trackmatch_tangent_x.clear();
      trackmatch_tangent_y.clear();
      trackmatch_dtanx.clear();
      trackmatch_dtany.clear();
      trackmatch_frost_is_hit.clear();
      trackmatch_external_expected_x.clear();
      trackmatch_external_expected_y.clear();
      trackmatch_external_tangent_x.clear();
      trackmatch_external_tangent_y.clear();
      trackmatch_external_chi2_x.clear();
      trackmatch_external_chi2_y.clear();
      trackmatch_external_ndof_x.clear();
      trackmatch_external_ndof_y.clear();
      trackmatch_external_num_planes_proton_module_x.clear();
      trackmatch_external_num_planes_proton_module_y.clear();
      trackmatch_external_num_planes_downstream_wagasci_x.clear();
      trackmatch_external_num_planes_downstream_wagasci_y.clear();
      trackmatch_external_pm_only_expected_x.clear();
      trackmatch_external_pm_only_expected_y.clear();
      trackmatch_external_pm_only_tangent_x.clear();
      trackmatch_external_pm_only_tangent_y.clear();
      trackmatch_external_pm_only_chi2_x.clear();
      trackmatch_external_pm_only_chi2_y.clear();
      trackmatch_external_pm_only_ndof_x.clear();
      trackmatch_external_pm_only_ndof_y.clear();
      trackmatch_external_dwg_only_expected_x.clear();
      trackmatch_external_dwg_only_expected_y.clear();
      trackmatch_external_dwg_only_tangent_x.clear();
      trackmatch_external_dwg_only_tangent_y.clear();
      trackmatch_external_dwg_only_chi2_x.clear();
      trackmatch_external_dwg_only_chi2_y.clear();
      trackmatch_external_dwg_only_ndof_x.clear();
      trackmatch_external_dwg_only_ndof_y.clear();
      trackmatch_true_frost_nearest_particle_id.clear();
      trackmatch_true_frost_nearest_position_x.clear();
      trackmatch_true_frost_nearest_position_y.clear();
      trackmatch_true_frost_nearest_tangent_x.clear();
      trackmatch_true_frost_nearest_tangent_y.clear();

      if (i < ntbm_tree->GetEntries()) {
        ntbm_tree->GetEntry(i);
        matchinfo_spill_pot = my_ntbm->GetSpillPot();
        for (int ibunch = 0; ibunch < NUMBER_OF_BUNCHES; ++ibunch) {
          matchinfo_bunch_pot[ibunch] = my_ntbm->GetBunchPot(ibunch);
        }
        matchinfo_bsd_spill_number = my_ntbm->GetBsdSpillNumber();
        matchinfo_timestamp = my_ntbm->GetTimestamp();
        matchinfo_bsd_good_spill_flag = my_ntbm->GetBsdGoodSpillFlag();
        matchinfo_wagasci_good_spill_flag = my_ntbm->GetWagasciGoodSpillFlag();
        for (int idet = 0; idet < 8; ++idet) {
          matchinfo_detector_flags[idet] = my_ntbm->GetDetectorFlags(idet);
        }
      } else {
        matchinfo_spill_pot = B2_NON_INITIALIZED_VALUE;
        for (int ibunch = 0; ibunch < NUMBER_OF_BUNCHES; ++ibunch) {
          matchinfo_bunch_pot[ibunch] = B2_NON_INITIALIZED_VALUE;
        }
        matchinfo_bsd_spill_number = B2_NON_INITIALIZED_VALUE;
        matchinfo_timestamp = B2_NON_INITIALIZED_VALUE;
        matchinfo_bsd_good_spill_flag = B2_NON_INITIALIZED_VALUE;
        matchinfo_wagasci_good_spill_flag = B2_NON_INITIALIZED_VALUE;
        for (int idet = 0; idet < 8; ++idet) {
          matchinfo_detector_flags[idet] = B2_NON_INITIALIZED_VALUE;
        }
      }


      trackmatch_has_match.reserve(rows.size());
      trackmatch_bm_track_id.reserve(rows.size());
      trackmatch_frost_match_bunch.reserve(rows.size());
      trackmatch_ninja_track_type.reserve(rows.size());
      trackmatch_baby_mind_tangent_x.reserve(rows.size());
      trackmatch_baby_mind_tangent_y.reserve(rows.size());
      trackmatch_expected_x.reserve(rows.size());
      trackmatch_expected_y.reserve(rows.size());
      trackmatch_frost_nearest_x.reserve(rows.size());
      trackmatch_frost_nearest_y.reserve(rows.size());
      trackmatch_dx.reserve(rows.size());
      trackmatch_dy.reserve(rows.size());
      trackmatch_tangent_x.reserve(rows.size());
      trackmatch_tangent_y.reserve(rows.size());
      trackmatch_dtanx.reserve(rows.size());
      trackmatch_dtany.reserve(rows.size());
      trackmatch_frost_is_hit.reserve(rows.size());
      trackmatch_external_expected_x.reserve(rows.size());
      trackmatch_external_expected_y.reserve(rows.size());
      trackmatch_external_tangent_x.reserve(rows.size());
      trackmatch_external_tangent_y.reserve(rows.size());
      trackmatch_external_chi2_x.reserve(rows.size());
      trackmatch_external_chi2_y.reserve(rows.size());
      trackmatch_external_ndof_x.reserve(rows.size());
      trackmatch_external_ndof_y.reserve(rows.size());
      trackmatch_external_num_planes_proton_module_x.reserve(rows.size());
      trackmatch_external_num_planes_proton_module_y.reserve(rows.size());
      trackmatch_external_num_planes_downstream_wagasci_x.reserve(rows.size());
      trackmatch_external_num_planes_downstream_wagasci_y.reserve(rows.size());
      trackmatch_external_pm_only_expected_x.reserve(rows.size());
      trackmatch_external_pm_only_expected_y.reserve(rows.size());
      trackmatch_external_pm_only_tangent_x.reserve(rows.size());
      trackmatch_external_pm_only_tangent_y.reserve(rows.size());
      trackmatch_external_pm_only_chi2_x.reserve(rows.size());
      trackmatch_external_pm_only_chi2_y.reserve(rows.size());
      trackmatch_external_pm_only_ndof_x.reserve(rows.size());
      trackmatch_external_pm_only_ndof_y.reserve(rows.size());
      trackmatch_external_dwg_only_expected_x.reserve(rows.size());
      trackmatch_external_dwg_only_expected_y.reserve(rows.size());
      trackmatch_external_dwg_only_tangent_x.reserve(rows.size());
      trackmatch_external_dwg_only_tangent_y.reserve(rows.size());
      trackmatch_external_dwg_only_chi2_x.reserve(rows.size());
      trackmatch_external_dwg_only_chi2_y.reserve(rows.size());
      trackmatch_external_dwg_only_ndof_x.reserve(rows.size());
      trackmatch_external_dwg_only_ndof_y.reserve(rows.size());
      trackmatch_true_frost_nearest_particle_id.reserve(rows.size());
      trackmatch_true_frost_nearest_position_x.reserve(rows.size());
      trackmatch_true_frost_nearest_position_y.reserve(rows.size());
      trackmatch_true_frost_nearest_tangent_x.reserve(rows.size());
      trackmatch_true_frost_nearest_tangent_y.reserve(rows.size());

      for (const auto &row : rows) {
        trackmatch_has_match.push_back(row.has_match);
        trackmatch_bm_track_id.push_back(row.bm_track_id);
        trackmatch_frost_match_bunch.push_back(row.frost_match_bunch);
        trackmatch_ninja_track_type.push_back(row.ninja_track_type);
        trackmatch_baby_mind_tangent_x.push_back(row.baby_mind_tangent_x);
        trackmatch_baby_mind_tangent_y.push_back(row.baby_mind_tangent_y);
        trackmatch_expected_x.push_back(row.expected_x);
        trackmatch_expected_y.push_back(row.expected_y);
        trackmatch_frost_nearest_x.push_back(row.frost_x);
        trackmatch_frost_nearest_y.push_back(row.frost_y);
        trackmatch_dx.push_back(row.dx);
        trackmatch_dy.push_back(row.dy);
        trackmatch_tangent_x.push_back(row.tangent_x);
        trackmatch_tangent_y.push_back(row.tangent_y);
        trackmatch_dtanx.push_back(row.dtanx);
        trackmatch_dtany.push_back(row.dtany);
        trackmatch_frost_is_hit.push_back(row.frost_is_hit);
        trackmatch_external_expected_x.push_back(row.external_expected_x);
        trackmatch_external_expected_y.push_back(row.external_expected_y);
        trackmatch_external_tangent_x.push_back(row.external_tangent_x);
        trackmatch_external_tangent_y.push_back(row.external_tangent_y);
        trackmatch_external_chi2_x.push_back(row.external_chi2_x);
        trackmatch_external_chi2_y.push_back(row.external_chi2_y);
        trackmatch_external_ndof_x.push_back(row.external_ndof_x);
        trackmatch_external_ndof_y.push_back(row.external_ndof_y);
        trackmatch_external_num_planes_proton_module_x.push_back(
          row.external_num_planes_proton_module_x);
        trackmatch_external_num_planes_proton_module_y.push_back(
          row.external_num_planes_proton_module_y);
        trackmatch_external_num_planes_downstream_wagasci_x.push_back(
          row.external_num_planes_downstream_wagasci_x);
        trackmatch_external_num_planes_downstream_wagasci_y.push_back(
          row.external_num_planes_downstream_wagasci_y);
        trackmatch_external_pm_only_expected_x.push_back(
          row.external_pm_only_expected_x);
        trackmatch_external_pm_only_expected_y.push_back(
          row.external_pm_only_expected_y);
        trackmatch_external_pm_only_tangent_x.push_back(
          row.external_pm_only_tangent_x);
        trackmatch_external_pm_only_tangent_y.push_back(
          row.external_pm_only_tangent_y);
        trackmatch_external_pm_only_chi2_x.push_back(
          row.external_pm_only_chi2_x);
        trackmatch_external_pm_only_chi2_y.push_back(
          row.external_pm_only_chi2_y);
        trackmatch_external_pm_only_ndof_x.push_back(
          row.external_pm_only_ndof_x);
        trackmatch_external_pm_only_ndof_y.push_back(
          row.external_pm_only_ndof_y);
        trackmatch_external_dwg_only_expected_x.push_back(
          row.external_dwg_only_expected_x);
        trackmatch_external_dwg_only_expected_y.push_back(
          row.external_dwg_only_expected_y);
        trackmatch_external_dwg_only_tangent_x.push_back(
          row.external_dwg_only_tangent_x);
        trackmatch_external_dwg_only_tangent_y.push_back(
          row.external_dwg_only_tangent_y);
        trackmatch_external_dwg_only_chi2_x.push_back(
          row.external_dwg_only_chi2_x);
        trackmatch_external_dwg_only_chi2_y.push_back(
          row.external_dwg_only_chi2_y);
        trackmatch_external_dwg_only_ndof_x.push_back(
          row.external_dwg_only_ndof_x);
        trackmatch_external_dwg_only_ndof_y.push_back(
          row.external_dwg_only_ndof_y);
        trackmatch_true_frost_nearest_particle_id.push_back(
          row.true_frost_nearest_particle_id);
        trackmatch_true_frost_nearest_position_x.push_back(
          row.true_frost_nearest_position_x);
        trackmatch_true_frost_nearest_position_y.push_back(
          row.true_frost_nearest_position_y);
        trackmatch_true_frost_nearest_tangent_x.push_back(
          row.true_frost_nearest_tangent_x);
        trackmatch_true_frost_nearest_tangent_y.push_back(
          row.true_frost_nearest_tangent_y);
      }

      match_info_out->Fill();
    }

    output_b2_file->cd();
    frost_match_out->Write("", TObject::kOverwrite);
    match_info_out->Write("", TObject::kOverwrite);

    if (datatype == B2DataType::kMonteCarlo) {
      CloneTreeIfExists(input_root.get(), output_b2_file, "frostmc");
      CloneTreeIfExists(input_root.get(), output_b2_file, "nRooTracker");
    }

    output_b2_file->Close();

    //efficiency
  if (total_tracks > 0) {
    double connection_efficiency = static_cast<double>(matched_tracks) / total_tracks;
    BOOST_LOG_TRIVIAL(info) << "Connection Efficiency: " << connection_efficiency * 100 << "%";
  } else {
    BOOST_LOG_TRIVIAL(info) << "No allowed tracks found, cannot calculate efficiency.";
  }

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  } catch (const std::invalid_argument &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Invalid argument error : " << error.what();
    std::exit(1);
  } catch (const std::out_of_range &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Out of range error : " << error.what();
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "==========FROST-Baby MIND Track Matching Finish==========";
  std::exit(0);

}
