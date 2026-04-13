#ifndef NINJARECON_TRACKMATCH_HPP
#define NINJARECON_TRACKMATCH_HPP

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>

#include "B2SpillSummary.hh"
#include "B2HitSummary.hh"
#include "B2BeamSummary.hh"
#include "B2TrackSummary.hh"
#include "B2Dimension.hh"
#include "NTBMSummary.hh"

/**
 * Comparator for B2HitSummary vector sort in one Baby MIND B2TrackSummary
 * @param lhs left hand side object
 * @param rhs right hand side object
 * @return true if the objects should not be swapped
 */
bool CompareBabyMindHitInOneTrack(const B2HitSummary* lhs, const B2HitSummary *rhs);

/**
 * FROST entry information read from frost_match tree
 */
struct FrostEntryData {
  std::vector<std::vector<double>> *x_rec = nullptr;
  std::vector<std::vector<double>> *y_rec = nullptr;
};

/**
 * Open copied trees from HitConverter output
 */
struct FrostMatchTrees {
  TTree *frost_match = nullptr;
  TTree *match_info = nullptr;
};

/**
 * Count matched FROST candidates for one Baby MIND track
 * x_match_count : number of matched candidates in top view (X)
 * y_match_count : number of matched candidates in side view (Y)
 */
struct FrostMatchCount {
  int x_match_count = 0;
  int y_match_count = 0;
};

/**
 * Track-wise FROST candidate information for one Baby MIND track.
 * x/y candidates are stored independently because FROST does not provide
 * a unique 2D pairing between reconstructed x and y candidates.
 */
struct FrostTrackCandidates {
  std::vector<double> frost_y_candidates;
  std::vector<double> frost_x_candidates;
  std::vector<double> expected_y_candidates;
  std::vector<double> expected_x_candidates;
  std::vector<double> difference_y_candidates;
  std::vector<double> difference_x_candidates;
  std::vector<double> tangent_y_candidates;
  std::vector<double> tangent_x_candidates;
};


/**
 * Open HitConverter output trees
 * @param file input ROOT file
 * @return tree bundle
 */
FrostMatchTrees OpenFrostMatchTrees(TFile *file);

/**
 * Get position and error for one Baby MIND plane
 * @param position position list of Baby MIND hits
 * @param view view
 * @return position and error for one Baby MIND plane
 */
std::vector<std::vector<double> > CalcMergedOnePlanePositionAndError(std::vector<std::vector<double> > position, int view);

/**
 * Get position and error for Baby MIND planes
 * @param hits vector of Baby MIND B2HitSummary objects
 * @param datatype MC or real data
 * @return Baby MIND position and errors
 */
std::vector<std::vector<std::vector<std::vector<double> > > > GenerateMergedPositionAndErrors(std::vector<const B2HitSummary* > hits, int datatype);

/**
 * Fit Baby MIND hits with straight lines in each view
 * @param track reconstructed B2TrackSummary object
 * @param datatype MC or real data
 * @return fit parameters in each view
 * param.at(view).at(0) = slope, param.at(view).at(1) = intercept
 */
std::vector<std::vector<double> > FitBabyMind(const B2TrackSummary *track, int datatype);

/**
 * Get Baby MIND initial direction and position
 * @param track reconstructed B2TrackSummary object
 * @param datatype MC or real data
 * @return at(0) means y and at(1) does x directions
 * at(2) means y and at(3) does x positions
 */
std::vector<double> GetBabyMindInitialDirectionAndPosition(const B2TrackSummary *track, int datatype);

/**
 * Calculate hit expected position on the NINJA tracker position
 * @param ntbm NTBMSummary object of the spill in interest
 * @param itrack Baby MIND track id (incremented from 0 NINJA internally)
 * @param z_shift Difference of z distance from nominal
 * @return at(0) means y and at(1) means x
 */
std::vector<double> CalculateExpectedPosition(NTBMSummary *ntbm, int itrack, double z_shift);

/**
 * Check if the Baby MIND reconstructed track expected to have hits
 * in the NINJA tracker
 * @param ntbm NTBMSummary object of the spill in interest
 * @param itrack Baby MIND track id (incremented from 0 NINJA internally)
 * @param z_shift Difference of z distance from nominal
 * @return true if the track expected to have hits else false
 */
bool NinjaHitExpected(NTBMSummary *ntbm, int itrack, double z_shift);

/**
 * Collect accepted FROST candidates for one Baby MIND track.
 * FROST version directly uses x_rec/y_rec in frost_match tree.
 * Matching is performed only with the FROST bunch corresponding to the
 * Baby MIND bunch of the track.
 *
 * @param ntbm NTBMSummary object for this spill
 * @param itrack Baby MIND track id
 * @param frost one matched FROST entry
 * @param z_shift Difference of z distance from nominal
 * @return accepted x/y candidate information for this track
 */
FrostTrackCandidates CollectFrostTrackCandidates(NTBMSummary *ntbm, int itrack,
                                                 const FrostEntryData &frost, double z_shift);

/**
 * Count accepted x/y candidates stored in one FrostTrackCandidates object
 * @param candidates track-wise FROST candidate information
 * @return numbers of accepted X-view and Y-view candidates
 */
FrostMatchCount CountAcceptedFrostMatches(const FrostTrackCandidates &candidates);


/**
 * Append one BM-track-wise FROST candidate block to NTBMSummary
 * @param ntbm NTBMSummary object for this spill
 * @param bm_track_id Baby MIND track id
 * @param bunch Baby MIND bunch number corresponding to this track
 * @param candidates accepted FROST x/y candidates for this track
 */
void AppendFrostTrackCandidates(NTBMSummary *ntbm, int bm_track_id, int bunch,
                                const FrostTrackCandidates &candidates);

/**
 * Set TSS info as true position/angle information to evaluate
 * FROST performance with MC
 * @param spill_summary B2SpillSummary object
 * @param ntbm_summary NTBMSummary object
 */
void SetTruePositionAngle(const B2SpillSummary& spill_summary, NTBMSummary* ntbm_summary);

/**
 * Transfer Beam information from B2BeamSummary to NTBMSummary
 * @param spill_summary B2SpillSummary object
 * @param ntbm_summary NTBMSummary object
 */
void TransferBeamInfo(const B2SpillSummary& spill_summary, NTBMSummary* ntbm_summary);

/**
 * Transfer Baby MIND track info from B2TrackSummary to NTBMSummary
 * @param spill_summary B2SpillSummary object
 * @param ntbm_summary NTBMSummary objec
 * @param datatype MC or real data
 * @param dimension Baby MIND geometry informationt
 */
void TransferBabyMindTrackInfo(const B2SpillSummary& spill_summary, NTBMSummary *ntbm_summary,
                               int datatype, B2Dimension &dimension);

/**
 * Transfer MC normalization info from B2EventSummary to NTBMSummary
 * @param spill_summary B2SpillSummary object
 * @param ntbm_summary NTBMSummary object
 */
void TransferMCInfo(const B2SpillSummary& spill_summary, NTBMSummary *ntbm_summary);

/**
 * Calculate Baby MIND track length using detector material information
 * @param track reconstructed B2TrackSummary object
 * @param ax tangent in x-z
 * @param ay tangent in y-z
 * @param dimension Baby MIND geometry information
 * @return material-weighted track length
 */
double MyFuncCalculateTrackLength(const B2TrackSummary *track, double ax, double ay,
                                  B2Dimension &dimension);
#endif
