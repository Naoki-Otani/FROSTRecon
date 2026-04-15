// system includes
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <memory>
#include <map>
#include <string>
#include <limits>
#include <cctype>

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

#include "TrackMatch.hpp"

namespace logging = boost::log;

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
};

struct NearestFrostPositionResult {
  bool found = false;
  double frost_position = B2_NON_INITIALIZED_VALUE;
  double diff = B2_NON_INITIALIZED_VALUE;
};

FrostMatchTrees OpenFrostMatchTrees(TFile *file) {
  FrostMatchTrees trees;
  trees.frost_match = dynamic_cast<TTree*>(file->Get("frost_match"));
  trees.match_info = dynamic_cast<TTree*>(file->Get("match_info"));
  if (!trees.frost_match || !trees.match_info) {
    throw std::runtime_error("Input file must contain frost_match and match_info trees");
  }
  return trees;
}

NearestFrostPositionResult FindNearestFrostPosition(const FrostEntryData &frost,
                                                    int bm_bunch,
                                                    int view,
                                                    double expected_position) {
  NearestFrostPositionResult result;

  const int frost_bunch_index = bm_bunch - 1;  // BM bunch: 1..8, FROST index: 0..7
  if (frost_bunch_index < 0 || frost_bunch_index >= 8) {
    return result;
  }

  const std::vector<std::vector<double>> *source = nullptr;
  if (view == B2View::kTopView) {
    source = frost.x_rec;
  } else if (view == B2View::kSideView) {
    source = frost.y_rec;
  } else {
    return result;
  }

  if (!source) return result;
  if (frost_bunch_index >= static_cast<int>(source->size())) return result;

  const std::vector<double> &positions = source->at(frost_bunch_index);
  if (positions.empty()) return result;

  double best_abs_diff = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < positions.size(); ++i) {
    const double diff = expected_position - positions.at(i);
    const double abs_diff = std::fabs(diff);
    if (!result.found || abs_diff < best_abs_diff) {
      result.found = true;
      result.frost_position = positions.at(i);
      result.diff = diff;
      best_abs_diff = abs_diff;
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
                                                 const FrostEntryData &frost, double z_shift) {
  FrostTrackCandidates candidates;

  const std::vector<double> hit_expected_position =
    CalculateExpectedPosition(ntbm, itrack, z_shift);
  const int bm_bunch = ntbm->GetBunch(itrack);
  const int frost_bunch_index = bm_bunch - 1;

  if (frost_bunch_index < 0 || frost_bunch_index >= 8) {
    BOOST_LOG_TRIVIAL(debug)
      << "Skip BM track " << itrack
      << " because BM bunch is out of FROST range: " << bm_bunch;
    return candidates;
  }

  const double dz_y =
    BABYMIND_POS_Z + BM_SECOND_LAYER_POS - NINJA_POS_Z_FROST - NINJA_FROST_POS_Z
    - (2 * B2View::kSideView - 1) * 10. + z_shift;
  const double dz_x =
    BABYMIND_POS_Z + BM_SECOND_LAYER_POS - NINJA_POS_Z_FROST - NINJA_FROST_POS_Z
    - (2 * B2View::kTopView - 1) * 10. + z_shift;

  if (frost.y_rec && frost_bunch_index < static_cast<int>(frost.y_rec->size())) {
    const std::vector<double> &ybunch = frost.y_rec->at(frost_bunch_index);
    for (std::size_t j = 0; j < ybunch.size(); ++j) {
      const double frost_y = ybunch.at(j);
      const double expected_y = hit_expected_position.at(B2View::kSideView);
      const double dy = expected_y - frost_y;
      if (std::fabs(dy) < TEMPORAL_ALLOWANCE[B2View::kSideView]) {
        candidates.frost_y_candidates.push_back(frost_y);
        candidates.expected_y_candidates.push_back(expected_y);
        candidates.difference_y_candidates.push_back(dy);
        candidates.tangent_y_candidates.push_back(dy / dz_y);
      }
    }
  }

  if (frost.x_rec && frost_bunch_index < static_cast<int>(frost.x_rec->size())) {
    const std::vector<double> &xbunch = frost.x_rec->at(frost_bunch_index);
    for (std::size_t j = 0; j < xbunch.size(); ++j) {
      const double frost_x = xbunch.at(j);
      const double expected_x = hit_expected_position.at(B2View::kTopView);
      const double dx = expected_x - frost_x;
      if (std::fabs(dx) < TEMPORAL_ALLOWANCE[B2View::kTopView]) {
        candidates.frost_x_candidates.push_back(frost_x);
        candidates.expected_x_candidates.push_back(expected_x);
        candidates.difference_x_candidates.push_back(dx);
        candidates.tangent_x_candidates.push_back(dx / dz_x);
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

void SetTruePositionAngle(const B2SpillSummary& spill_summary, NTBMSummary* ntbm_summary) {

  auto it_event = spill_summary.BeginTrueEvent();
  const auto *event = it_event.Next();
  auto &primary_vertex_summary = event->GetPrimaryVertex();

  auto it_emulsion = spill_summary.BeginEmulsion();
  TVector3 true_position;
  TVector3 true_direction;
  bool found_true_muon_in_tss = false;
  while (const auto *emulsion = it_emulsion.Next()) {
    if ( emulsion->GetParentTrackId() == 0 ) continue;
    if ( emulsion->GetParentTrackId() >= primary_vertex_summary.GetNumOutgoingTracks() )
      continue;
    // Get position of TSS downstream film position
    if (emulsion->GetFilmType() == B2EmulsionType::kShifter && emulsion->GetPlate() == 17) {
      int particle_id = emulsion->GetParentTrack().GetParticlePdg();
      if (!B2Pdg::IsMuonPlusOrMinus(particle_id)) continue;
      true_position = emulsion->GetAbsolutePosition().GetValue();
      true_direction = emulsion->GetTangent().GetValue();
      // Convert true position to FROST local coordinate consistently with
      // CalculateExpectedPosition() / FROST matching.
      true_position.SetX(true_position.X() + true_direction.X() * (NINJA_TSS_ATTACH_AC_THICK + 30.4)
			 - NINJA_POS_X - NINJA_FROST_POS_X);
      true_position.SetY(true_position.Y() + true_direction.Y() * (NINJA_TSS_ATTACH_AC_THICK + 10.4)
			 - NINJA_POS_Y - NINJA_FROST_POS_Y);
      true_position.SetZ(true_position.Z() - NINJA_POS_Z_FROST - NINJA_FROST_POS_Z);
      found_true_muon_in_tss = true;
      break;
    }
  }

  if ( !found_true_muon_in_tss ) return;

  std::vector<double> true_ninja_position;
  true_ninja_position.resize(2);
  true_ninja_position.at(B2View::kSideView) = true_position.Y();
  true_ninja_position.at(B2View::kTopView) = true_position.X();
  std::vector<double> true_ninja_tangent;
  true_ninja_tangent.resize(2);
  true_ninja_tangent.at(B2View::kSideView) = true_direction.Y();
  true_ninja_tangent.at(B2View::kTopView) = true_direction.X();

  for ( int icluster = 0; icluster < ntbm_summary->GetNumberOfFrostMatchEntries(); icluster++ ) {
    if ( ntbm_summary->GetBabyMindTrackId(icluster) >= 0 ) {
      ntbm_summary->SetNumberOfTrueParticles(icluster, 1);
      ntbm_summary->SetTrueParticleId(icluster, 0, (int)PDG_t::kMuonMinus);
      ntbm_summary->SetTruePosition(icluster, 0, true_ninja_position);
      ntbm_summary->SetTrueTangent(icluster, 0, true_ninja_tangent);
    } else {
      ntbm_summary->SetNumberOfTrueParticles(icluster, 0);
    }
  }

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

void TransferBabyMindTrackInfo(const B2SpillSummary &spill_summary, NTBMSummary *ntbm_summary, int datatype,
			       B2Dimension &dimension) {

  int itrack = 0;

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
      if ( nll_minus - nll_plus >= 4 )
	ntbm_summary->SetCharge(itrack, 1);
      else
	ntbm_summary->SetCharge(itrack, -1);
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
//一旦デバッグ用にオリジナルのB2トラック長を使用。
      // Use the original B2 track length here.
      // Some tracks contain invalid Baby MIND hit summaries that can crash
      // MyFuncCalculateTrackLength() when accessing hit->GetPlane().
      // ntbm_summary->SetTrackLengthTotal(itrack, track->GetTrackLengthTotal());


      itrack++;

    } // while track
  } // while vertex

}

void TransferMCInfo(const B2SpillSummary &spill_summary, NTBMSummary *ntbm_summary) {
  auto it_event = spill_summary.BeginTrueEvent();
  const auto *event = it_event.Next();
  ntbm_summary->SetNormalization(event->GetNormalization());
  ntbm_summary->SetTotalCrossSection(event->GetTotalCrossSection());
}

// main

int main(int argc, char *argv[]) {

  gErrorIgnoreLevel = kError;
  gRandom->SetSeed(1);

  B2Dimension dimension_((std::string)"/opt/wagasci_mc/WagasciMC/etc/wagasci/b2/geometry");

  logging::trivial::severity_level log_level = logging::trivial::info;
  if (argc == 6) {
    log_level = ParseLogLevel(argv[5]);
  }

  logging::core::get()->set_filter(
      logging::trivial::severity >= log_level);

  BOOST_LOG_TRIVIAL(info) << "==========FROST-Baby MIND Track Matching Start==========";

  if ( argc != 5 && argc != 6 ) {
    std::cerr << "Usage : " << argv[0]
              << " <input B2 file path>"
              << " <output ROOT file path>"
              << " <z shift>"
              << " <MC(0)/data(1)>"
              << " [trace|debug|info|warning|error|fatal]"
              << std::endl;
    std::exit(1);
  }

  try {
    B2Reader reader(argv[1]);
    std::unique_ptr<TFile> input_root(new TFile(argv[1], "READ"));
    if (!input_root || input_root->IsZombie()) {
      throw std::runtime_error("Cannot open input ROOT file");
    }
    FrostMatchTrees frost_trees = OpenFrostMatchTrees(input_root.get());

    std::vector<std::vector<double>> *x_rec = nullptr;
    std::vector<std::vector<double>> *y_rec = nullptr;
    std::vector<int> *is_hit = nullptr;
    frost_trees.frost_match->SetBranchAddress("x_rec", &x_rec);
    frost_trees.frost_match->SetBranchAddress("y_rec", &y_rec);
    frost_trees.frost_match->SetBranchAddress("is_hit", &is_hit);

    TTree *ntbm_tree = new TTree("ntbm", "NINJA BabyMIND Original Summary");
    ntbm_tree->SetDirectory(nullptr);
    NTBMSummary* my_ntbm = new NTBMSummary();
    ntbm_tree->Branch("NTBMSummary", &my_ntbm);

    double z_shift = std::stof(argv[3]);
    int datatype = std::stoi(argv[4]);

    int nspill = 0;
    std::vector<std::vector<TrackMatchRow>> spill_match_rows;

    //20241022 count matched tracks
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
      if ( datatype == B2DataType::kMonteCarlo)
	      TransferMCInfo(input_spill_summary, my_ntbm);

      FrostEntryData frost_entry;
      if (nspill < frost_trees.frost_match->GetEntries()) {
        frost_trees.frost_match->GetEntry(nspill);
        frost_entry.x_rec = x_rec;
        frost_entry.y_rec = y_rec;
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
	      TransferBabyMindTrackInfo(input_spill_summary, my_ntbm, datatype, dimension_);

        std::vector<TrackMatchRow> rows;

        for ( int ibmtrack = 0; ibmtrack < my_ntbm->GetNumberOfTracks(); ibmtrack++ ) {
          TrackMatchRow row;
          row.bm_track_id = ibmtrack;
          row.frost_match_bunch = my_ntbm->GetBunch(ibmtrack);
          row.ninja_track_type = my_ntbm->GetNinjaTrackType(ibmtrack);
          row.baby_mind_tangent_y = my_ntbm->GetBabyMindTangent(ibmtrack, B2View::kSideView);
          row.baby_mind_tangent_x = my_ntbm->GetBabyMindTangent(ibmtrack, B2View::kTopView);


          const std::vector<double> expected_position =
            CalculateExpectedPosition(my_ntbm, ibmtrack, z_shift);
          my_ntbm->SetExtrapolatedPosition(ibmtrack, expected_position);
          row.expected_y = expected_position.at(B2View::kSideView);
          row.expected_x = expected_position.at(B2View::kTopView);

          const int bm_bunch = my_ntbm->GetBunch(ibmtrack);
          const NearestFrostPositionResult nearest_y =
            FindNearestFrostPosition(frost_entry, bm_bunch, B2View::kSideView, row.expected_y);
          const NearestFrostPositionResult nearest_x =
            FindNearestFrostPosition(frost_entry, bm_bunch, B2View::kTopView, row.expected_x);

          const double dz_y =
            BABYMIND_POS_Z + BM_SECOND_LAYER_POS - NINJA_POS_Z_FROST - NINJA_FROST_POS_Z
            - (2 * B2View::kSideView - 1) * 10. + z_shift;
          const double dz_x =
            BABYMIND_POS_Z + BM_SECOND_LAYER_POS - NINJA_POS_Z_FROST - NINJA_FROST_POS_Z
            - (2 * B2View::kTopView - 1) * 10. + z_shift;

          if (nearest_y.found) {
            row.frost_y = nearest_y.frost_position;
            row.dy = nearest_y.diff;
            row.tangent_y = nearest_y.diff / dz_y;
          }
          if (nearest_x.found) {
            row.frost_x = nearest_x.frost_position;
            row.dx = nearest_x.diff;
            row.tangent_x = nearest_x.diff / dz_x;
          }

          const bool hit_expected = NinjaHitExpected(my_ntbm, ibmtrack, z_shift);
          if (hit_expected) {
            total_tracks++;
          }

          const FrostTrackCandidates candidates =
            CollectFrostTrackCandidates(my_ntbm, ibmtrack, frost_entry, z_shift);
          const FrostMatchCount match_count =
            CountAcceptedFrostMatches(candidates);

          AppendFrostTrackCandidates(my_ntbm, ibmtrack, bm_bunch, candidates);

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

        if ( datatype == B2DataType::kMonteCarlo &&
            my_ntbm->GetNumberOfFrostMatchEntries() > 0 )
          SetTruePositionAngle(input_spill_summary, my_ntbm);

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

    TTree *frost_match_out = frost_trees.frost_match->CloneTree(-1);
    frost_match_out->SetName("frost_match");
    frost_match_out->SetDirectory(output_b2_file);

    TTree *match_info_out = frost_trees.match_info->CloneTree(0);
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
    std::vector<Int_t> trackmatch_frost_is_hit;
    Double_t matchinfo_spill_pot = B2_NON_INITIALIZED_VALUE;
    Double_t matchinfo_bunch_pot[NUMBER_OF_BUNCHES];
    Int_t matchinfo_bsd_spill_number = B2_NON_INITIALIZED_VALUE;
    Double_t matchinfo_timestamp = B2_NON_INITIALIZED_VALUE;
    Int_t matchinfo_bsd_good_spill_flag = B2_NON_INITIALIZED_VALUE;
    Int_t matchinfo_wagasci_good_spill_flag = B2_NON_INITIALIZED_VALUE;
    Int_t matchinfo_detector_flags[8];

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
    match_info_out->Branch("trackmatch_frost_is_hit", &trackmatch_frost_is_hit);
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

    const Long64_t ninfo = frost_trees.match_info->GetEntries();
    for (Long64_t i = 0; i < ninfo; ++i) {
      frost_trees.match_info->GetEntry(i);
      if (i < frost_trees.frost_match->GetEntries()) {
        frost_trees.frost_match->GetEntry(i);
      } else {
        is_hit = nullptr;
      }
      const std::vector<TrackMatchRow> &rows =
        (i < static_cast<Long64_t>(spill_match_rows.size()))
          ? spill_match_rows.at(i)
          : std::vector<TrackMatchRow>{};

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
      trackmatch_frost_is_hit.clear();

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
      trackmatch_frost_is_hit.reserve(rows.size());

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

        Int_t frost_is_hit_value = -1;
        const int bunch = row.frost_match_bunch; // 1..8
        if (is_hit && bunch >= 1 &&
            bunch <= static_cast<int>(is_hit->size())) {
          frost_is_hit_value = is_hit->at(bunch - 1); // vector index is 0..7
        }
        trackmatch_frost_is_hit.push_back(frost_is_hit_value);
      }

      match_info_out->Fill();
    }

    output_b2_file->cd();
    frost_match_out->Write("", TObject::kOverwrite);
    match_info_out->Write("", TObject::kOverwrite);
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
