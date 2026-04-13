#ifndef NTBMSUMMARY_HH
#define NTBMSUMMARY_HH

#include <TObject.h>
#include <TVector2.h>

#include <vector>

#include "NTBMConst.hh"

#ifdef __ROOTCLING__
#pragma link off globals;
#pragma link off classes;
#pragma link off functions;
#pragma link C++ nestedclasses;
#pragma link C++ class NTBMSummary+;
#endif

/**
 * Class containing info for Baby MIND / FROST muon track matching
 */

class NTBMSummary: public TObject {

public :

  NTBMSummary();

  /**
   * Set the all members to zero
   * no memory is released.
   * @param option same as TObject option
   */
  void Clear(Option_t *option) override;

  /**
   * Print all the NTBM members value to the output stream
   * @param os output stream
   * @param obj object
   * @return output stream
   */
  friend std::ostream &operator<<(std::ostream &os, const NTBMSummary &obj);

  // Setter/Getter

  void SetEntryInDailyFile(int entry_in_daily_file);

  int GetEntryInDailyFile() const;

  void SetSpillPot(double spill_pot);

  double GetSpillPot() const;

  void SetBunchPot(int bunch, double bunch_pot);

  double GetBunchPot(int bunch) const;

  void SetBsdSpillNumber(int bsd_spill_number);

  int GetBsdSpillNumber() const;

  void SetTimestamp(double timestamp);

  double GetTimestamp() const;

  void SetBsdGoodSpillFlag(int bsd_good_spill_flag);

  int GetBsdGoodSpillFlag() const;

  void SetWagasciGoodSpillFlag(int wagasci_good_spill_flag);

  int GetWagasciGoodSpillFlag() const;

  void SetDetectorFlags(int detector, int detector_flag);

  int GetDetectorFlags(int detector) const;

  void SetNumberOfTracks(int number_of_tracks);

  int GetNumberOfTracks() const;

  void SetNinjaTrackType(int track, int ninja_track_type);

  int GetNinjaTrackType(int track) const;

  void SetMomentumType(int track, int momentum_type);

  int GetMomentumType(int track) const;

  void SetMomentum(int track, double momentum);

  double GetMomentum(int track) const;

  void SetMomentumError(int track, double momentum_error);

  double GetMomentumError(int track) const;

  void SetBabyMindPosition(int track, int view, double baby_mind_position);

  void SetBabyMindPosition(int track, std::vector<double> baby_mind_position);

  std::vector<double> GetBabyMindPosition(int track) const;

  double GetBabyMindPosition(int track, int view) const;

  void SetBabyMindPositionError(int track, int view, double baby_mind_position_error);

  void SetBabyMindPositionError(int track, std::vector<double> baby_mind_position_error);

  std::vector<double> GetBabyMindPositionError(int track) const;

  double GetBabyMindPositionError(int track, int view) const;

  void SetBabyMindTangent(int track, int view, double baby_mind_tangent);

  void SetBabyMindTangent(int track, std::vector<double> baby_mind_tangent);

  std::vector<double> GetBabyMindTangent(int track) const;

  double GetBabyMindTangent(int track, int view) const;

  void SetBabyMindTangentError(int track, int view, double baby_mind_tangent_error);

  void SetBabyMindTangentError(int track, std::vector<double> baby_mind_tangent_error);

  std::vector<double> GetBabyMindTangentError(int track) const;

  double GetBabyMindTangentError(int track, int view) const;

  void SetExtrapolatedPosition(int track, int view, double extrapolated_position);

  void SetExtrapolatedPosition(int track, std::vector<double> extrapolated_position);

  std::vector<double> GetExtrapolatedPosition(int track) const;

  double GetExtrapolatedPosition(int track, int view) const;

  void SetBabyMindMaximumPlane(int track, int baby_mind_maximum_plane);

  int GetBabyMindMaximumPlane(int track) const;

  void SetTrackLengthTotal(int track, double track_length_total);

  double GetTrackLengthTotal(int track) const;

  void SetCharge(int track, int charge);

  int GetCharge(int track) const;

  void SetDirection(int track, int direction);

  int GetDirection(int track) const;

  void SetBunch(int track, int bunch);

  int GetBunch(int track) const;

  void SetNumberOfFrostMatchEntries(int number_of_frost_match_entries);

  int GetNumberOfFrostMatchEntries() const;

  void SetBabyMindTrackId(int entry, int baby_mind_track_id);

  int GetBabyMindTrackId(int entry) const;

  void SetFrostMatchBunch(int entry, int frost_match_bunch);

  int GetFrostMatchBunch(int entry) const;

  void SetNumberOfFrostYCandidates(int entry, int number_of_frost_y_candidates);

  int GetNumberOfFrostYCandidates(int entry) const;

  void SetNumberOfFrostXCandidates(int entry, int number_of_frost_x_candidates);

  int GetNumberOfFrostXCandidates(int entry) const;

  void SetFrostYCandidates(int entry, std::vector<double> frost_y_candidates);

  std::vector<double> GetFrostYCandidates(int entry) const;

  void SetFrostXCandidates(int entry, std::vector<double> frost_x_candidates);

  std::vector<double> GetFrostXCandidates(int entry) const;

  void SetExpectedYCandidates(int entry, std::vector<double> expected_y_candidates);

  std::vector<double> GetExpectedYCandidates(int entry) const;

  void SetExpectedXCandidates(int entry, std::vector<double> expected_x_candidates);

  std::vector<double> GetExpectedXCandidates(int entry) const;

  void SetDifferenceYCandidates(int entry, std::vector<double> difference_y_candidates);

  std::vector<double> GetDifferenceYCandidates(int entry) const;

  void SetDifferenceXCandidates(int entry, std::vector<double> difference_x_candidates);

  std::vector<double> GetDifferenceXCandidates(int entry) const;

  void SetTangentYCandidates(int entry, std::vector<double> tangent_y_candidates);

  std::vector<double> GetTangentYCandidates(int entry) const;

  void SetTangentXCandidates(int entry, std::vector<double> tangent_x_candidates);

  std::vector<double> GetTangentXCandidates(int entry) const;

  void SetNumberOfTrueParticles(int cluster, int number_of_true_particles);

  int GetNumberOfTrueParticles(int cluster) const;

  void SetTrueParticleId(int cluster, int particle, int true_particle_id);

  void SetTrueParticleId(int cluster, std::vector<int> true_particle_id);

  std::vector<int> GetTrueParticleId(int cluster) const;

  int GetTrueParticleId(int cluster, int particle) const;

  void SetNormalization(double normalization);

  double GetNormalization() const;

  void SetTotalCrossSection(double total_cross_section);

  double GetTotalCrossSection() const;

  void SetTruePosition(int cluster, int particle, int view, double true_position);

  void SetTruePosition(int cluster, int particle, std::vector<double> true_position);

  void SetTruePosition(int cluster, std::vector<std::vector<double>> true_position);

  std::vector<std::vector<double>> GetTruePosition(int cluster) const;

  std::vector<double> GetTruePosition(int cluster, int particle) const;

  double GetTruePosition(int cluster, int particle, int view) const;

  void SetTrueTangent(int clcuster, int particle, int view, double true_tangent);

  void SetTrueTangent(int cluster, int particle, std::vector<double> true_tangent);

  void SetTrueTangent(int cluster, std::vector<std::vector<double>> true_tangent);

  std::vector<std::vector<double>> GetTrueTangent(int cluster) const;

  std::vector<double> GetTrueTangent(int cluster, int particle) const;

  double GetTrueTangent(int cluster, int particle, int view) const;

private :

  /////> File information/////
  int entry_in_daily_file_;

  /////> Beam information extracted from B2BeamSummary/////
  ///> Total POT of the spill
  double spill_pot_;
  ///> POT for each bunch
  double bunch_pot_[NUMBER_OF_BUNCHES];
  ///> BSD spill number
  int bsd_spill_number_;
  ///> Spill time stamp
  double timestamp_;
  ///> BSD spill flag
  int bsd_good_spill_flag_;
  ///> WAGASCI spill flag
  int wagasci_good_spill_flag_;
  ///> detector flags
  int detector_flags_[8];

  /////> Baby MIND information extracted from B2TrackSummary and B2ClusterSummary/////
  ///> track -> view(2)
  ///> Number of Baby MIND reconstructed tracks;
  int number_of_tracks_;
  ///> Track type (0:ECC interaction cand, 1:sand muon, 2:Upstream WAGASCI/PM -1:Downstream WAGASCI)
  std::vector<int> ninja_track_type_;
  ///> Momentum measurement type (0:range, 1:curvature)
  std::vector<int> momentum_type_;
  ///> Baby MIND reconstructed momentum
  std::vector<double> momentum_;
  ///> Baby MIND reconstructed momentum error
  std::vector<double> momentum_error_;
  ///> Baby MIND reconstructed position
  std::vector<std::vector<double>> baby_mind_position_;
  ///> Baby MIND reconstructed position error
  std::vector<std::vector<double>> baby_mind_position_error_;
  ///> Baby MIND reconstructed tangent
  std::vector<std::vector<double>> baby_mind_tangent_;
  ///> Baby MIND reconstructed tangent error
  std::vector<std::vector<double>> baby_mind_tangent_error_;
  ///> Baby MIND track-wise extrapolated position at FROST plane
  ///> track -> view(2)
  std::vector<std::vector<double>> extrapolated_position_;
  ///> Baby MIND maximum plane
  std::vector<int> baby_mind_maximum_plane_;
  ///> WAGASCI total material length
  std::vector<double> track_length_total_;
  ///> Baby MIND reconstructed charge (assuming muon +/-)
  std::vector<int> charge_;
  ///> Baby MIND reconstructed track direction (+/-)
  std::vector<int> direction_;
  ///> Bunch number where the track detected
  std::vector<int> bunch_;

  /////> BM-FROST matched entry information/////
  ///> one entry corresponds to one Baby MIND track having stored FROST candidates
  int number_of_frost_match_entries_;
  ///> Corresponding Baby MIND track id for each stored entry
  std::vector<int> baby_mind_track_id_;
  ///> Baby MIND bunch for each stored entry
  std::vector<int> frost_match_bunch_;
  ///> Number of stored FROST y candidates for each entry
  std::vector<int> number_of_frost_y_candidates_;
  ///> Number of stored FROST x candidates for each entry
  std::vector<int> number_of_frost_x_candidates_;
  ///> FROST y candidates for each stored entry
  std::vector<std::vector<double>> frost_y_candidates_;
  ///> FROST x candidates for each stored entry
  std::vector<std::vector<double>> frost_x_candidates_;
  ///> BM extrapolated y positions corresponding to frost_y_candidates_
  std::vector<std::vector<double>> expected_y_candidates_;
  ///> BM extrapolated x positions corresponding to frost_x_candidates_
  std::vector<std::vector<double>> expected_x_candidates_;
  ///> expected_y_candidates_ - frost_y_candidates_
  std::vector<std::vector<double>> difference_y_candidates_;
  ///> expected_x_candidates_ - frost_x_candidates_
  std::vector<std::vector<double>> difference_x_candidates_;
  ///> Tangent candidates in y-z from BM and FROST positions
  std::vector<std::vector<double>> tangent_y_candidates_;
  ///> Tangent candidates in x-z from BM and FROST positions
  std::vector<std::vector<double>> tangent_x_candidates_;

  /////> True particle information for MC/////
  ///> frost match entry -> true particle -> view(2)
  ///> Normalization factor from beam MC
  double normalization_;
  ///> Total cross section from NEUT
  double total_cross_section_;
  ///> Number of true particles making one cluster
  std::vector<int> number_of_true_particles_;
  ///> PDG particle id of true particles
  std::vector<std::vector<int>> true_particle_id_;
  ///> True position
  std::vector<std::vector<std::vector<double>>> true_position_;
  ///> True tangent
  std::vector<std::vector<std::vector<double>>> true_tangent_;

  ClassDefOverride(NTBMSummary, 15) // NT BM Summary
};

#endif
