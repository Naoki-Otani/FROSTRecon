#include "NTBMSummary.hh"
#include "NTBMConst.hh"

#include <ostream>

NTBMSummary::NTBMSummary() {
  NTBMSummary::Clear("C");
}

void NTBMSummary::Clear(Option_t *option) {
  entry_in_daily_file_ = -1;
  spill_pot_ = -1.;
  for ( double &i : bunch_pot_ )
    i = -1.;
  bsd_spill_number_ = 0;
  timestamp_ = 0.;
  bsd_good_spill_flag_ = -1;
  wagasci_good_spill_flag_ = -1;
  for (int i = 0; i < 8; i++)
    detector_flags_[i] = -1;
  number_of_tracks_ = 0;
  ninja_track_type_.clear();
  momentum_type_.clear();
  momentum_.clear();
  momentum_error_.clear();
  baby_mind_position_.clear();
  baby_mind_position_error_.clear();
  baby_mind_tangent_.clear();
  baby_mind_tangent_error_.clear();
  extrapolated_position_.clear();
  baby_mind_maximum_plane_.clear();
  track_length_total_.clear();
  charge_.clear();
  direction_.clear();
  bunch_.clear();
  number_of_frost_match_entries_ = 0;
  baby_mind_track_id_.clear();
  frost_match_bunch_.clear();
  number_of_frost_y_candidates_.clear();
  number_of_frost_x_candidates_.clear();
  frost_y_candidates_.clear();
  frost_x_candidates_.clear();
  expected_y_candidates_.clear();
  expected_x_candidates_.clear();
  difference_y_candidates_.clear();
  difference_x_candidates_.clear();
  tangent_y_candidates_.clear();
  tangent_x_candidates_.clear();
  normalization_ = 1.;
  total_cross_section_ = 1.;
  number_of_true_particles_.clear();
  true_particle_id_.clear();
  true_position_.clear();
  true_tangent_.clear();
  TObject::Clear(option);
}

std::ostream &operator<<(std::ostream &os, const NTBMSummary &obj) {
  os << "Entry in daily file = " << obj.entry_in_daily_file_ << "\n"
     << "Total POT of this spill = " << obj.spill_pot_ << "\n"
     << "POT for each bunch = ";
  for (int i = 0; i < NUMBER_OF_BUNCHES; i++) {
    os << i + 1 << " : " << obj.bunch_pot_[i];
    if (i != NUMBER_OF_BUNCHES - 1) os << ", ";
  }
  os << "\n"
     << "Timestamp = " << obj.timestamp_ << "\n"
     << "BSD good spill flag (good : 1, bad : 0) = " << obj.bsd_good_spill_flag_ << "\n"
     << "WAGASCI good spill flag (good : 1, bad : 0) = " << obj.wagasci_good_spill_flag_ << "\n"
     << "Detector flags (good : 1, bad : 0) = ";
  for (int i = 0; i < 8; i++) {
    os << i+1 << " : " << obj.detector_flags_[i];
    if (i != 7) os << ", ";
  }
  os << "\n"
     << "Number of Baby MIND tracks = " << obj.number_of_tracks_ << "\n"
     << "Baby MIND track type (ECC cand. : 0, Sand cand. : 1) = ";
  for (int i = 0; i < obj.number_of_tracks_; i++) {
    os << i + 1 << " : " << obj.ninja_track_type_.at(i);
    if(i != obj.number_of_tracks_ - 1) os << ", ";
  }
  os << "\n"
     << "Momentum measurement type (Baby MIND range : 0, curvature : 1) = ";
  for (int i = 0; i < obj.number_of_tracks_; i++) {
    os << i + 1 << " : " << obj.momentum_type_.at(i);
    if(i != obj.number_of_tracks_ - 1) os << ", ";
  }
  os << "\n"
     << "Absolute momentum = ";
  for (int i = 0; i < obj.number_of_tracks_; i++) {
    os << i + 1 << " : " << obj.momentum_.at(i) << " +/- "
       << obj.momentum_error_.at(i);
    if(i != obj.number_of_tracks_ - 1) os << ", ";
  }
  os << "\n"
     << "Baby MIND initial position = ";
  for (int i = 0; i < obj.number_of_tracks_; i++) {
    os << i + 1 << " : ( "
       << obj.baby_mind_position_.at(i).at(0) << " +/- "
       << obj.baby_mind_position_error_.at(i).at(0) << ", "
       << obj.baby_mind_position_.at(i).at(1) << "+/-"
       << obj.baby_mind_position_error_.at(i).at(1) << " )\n";
  }
  os << "\n"
     << "Baby MIND initial tangent = ";
  for (int i = 0; i < obj.number_of_tracks_; i++) {
    os << i + 1 << " : ( "
       << obj.baby_mind_tangent_.at(i).at(0) << " +/- "
       << obj.baby_mind_tangent_error_.at(i).at(0) << ", "
       << obj.baby_mind_tangent_.at(i).at(1) << " +/- "
       << obj.baby_mind_tangent_error_.at(i).at(1) <<" )\n";
  }
  os << "\n"
     << "Baby MIND extrapolated position at FROST plane = ";
  for (int i = 0; i < obj.number_of_tracks_; i++) {
    os << i + 1 << " : ( "
       << obj.extrapolated_position_.at(i).at(0) << ", "
       << obj.extrapolated_position_.at(i).at(1) << " )\n";
  }
  os << "\n"
     << "Baby MIND maximum plane = ";
  for (int i = 0; i < obj.number_of_tracks_; i++) {
    os << i + 1 << " : " << obj.baby_mind_maximum_plane_.at(i);
    if(i != obj.number_of_tracks_ - 1) os << ", ";
  }
  os << "\n"
     << "Track length total = ";
  for (int i = 0; i < obj.number_of_tracks_; i++) {
    os << i + 1 << " : " << obj.track_length_total_.at(i);
    if(i != obj.number_of_tracks_ - 1) os << ", ";
  }
  os << "\n"
     << "Charge = ";
  for (int i = 0; i < obj.number_of_tracks_; i++) {
    os << i + 1 << " : " << obj.charge_.at(i);
    if(i != obj.number_of_tracks_ - 1) os << ", ";
  }
  os << "\n"
     << "Direction = ";
  for (int i = 0; i < obj.number_of_tracks_; i++) {
    os << i + 1 << " : " << obj.direction_.at(i);
    if(i != obj.number_of_tracks_ - 1) os << ", ";
  }
  os << "\n"
     << "Bunch = ";
  for (int i = 0; i < obj.number_of_tracks_; i++) {
    os << i + 1 << " : " << obj.bunch_.at(i);
    if(i != obj.number_of_tracks_ - 1) os << ", ";
  }
  os << "\n"
     << "Number of BM-FROST matched entries = " << obj.number_of_frost_match_entries_ << "\n"
     << "Corresponding Baby MIND track ID = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : " << obj.baby_mind_track_id_.at(i);
    if(i != obj.number_of_frost_match_entries_ - 1) os << ", ";
  }
  os << "\n"
     << "Matched entry bunch = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : " << obj.frost_match_bunch_.at(i);
    if (i != obj.number_of_frost_match_entries_ - 1) os << ", ";
  }
  os << "\n"
     << "Number of FROST y candidates = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : " << obj.number_of_frost_y_candidates_.at(i);
    if (i != obj.number_of_frost_match_entries_ - 1) os << ", ";
  }
  os << "\n"
     << "Number of FROST x candidates = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : " << obj.number_of_frost_x_candidates_.at(i);
    if (i != obj.number_of_frost_match_entries_ - 1) os << ", ";
  }
  os << "\n"
     << "FROST y candidates = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : ( ";
    for (std::size_t j = 0; j < obj.frost_y_candidates_.at(i).size(); j++) {
      os << obj.frost_y_candidates_.at(i).at(j);
      if (j + 1 != obj.frost_y_candidates_.at(i).size()) os << ", ";
    }
    os << " )\n";
  }
  os << "\n"
     << "FROST x candidates = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : ( ";
    for (std::size_t j = 0; j < obj.frost_x_candidates_.at(i).size(); j++) {
      os << obj.frost_x_candidates_.at(i).at(j);
      if (j + 1 != obj.frost_x_candidates_.at(i).size()) os << ", ";
    }
    os << " )\n";
  }
  os << "\n"
     << "Expected y candidates = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : ( ";
    for (std::size_t j = 0; j < obj.expected_y_candidates_.at(i).size(); j++) {
      os << obj.expected_y_candidates_.at(i).at(j);
      if (j + 1 != obj.expected_y_candidates_.at(i).size()) os << ", ";
    }
    os << " )\n";
  }
  os << "\n"
     << "Expected x candidates = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : ( ";
    for (std::size_t j = 0; j < obj.expected_x_candidates_.at(i).size(); j++) {
      os << obj.expected_x_candidates_.at(i).at(j);
      if (j + 1 != obj.expected_x_candidates_.at(i).size()) os << ", ";
    }
    os << " )\n";
  }
  os << "\n"
     << "Difference y candidates = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : ( ";
    for (std::size_t j = 0; j < obj.difference_y_candidates_.at(i).size(); j++) {
      os << obj.difference_y_candidates_.at(i).at(j);
      if (j + 1 != obj.difference_y_candidates_.at(i).size()) os << ", ";
    }
    os << " )\n";
  }
  os << "\n"
     << "Difference x candidates = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : ( ";
    for (std::size_t j = 0; j < obj.difference_x_candidates_.at(i).size(); j++) {
      os << obj.difference_x_candidates_.at(i).at(j);
      if (j + 1 != obj.difference_x_candidates_.at(i).size()) os << ", ";
    }
    os << " )\n";
  }
  os << "\n"
     << "Tangent y candidates = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : ( ";
    for (std::size_t j = 0; j < obj.tangent_y_candidates_.at(i).size(); j++) {
      os << obj.tangent_y_candidates_.at(i).at(j);
      if (j + 1 != obj.tangent_y_candidates_.at(i).size()) os << ", ";
    }
    os << " )\n";
  }
  os << "\n"
     << "Tangent x candidates = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : ( ";
    for (std::size_t j = 0; j < obj.tangent_x_candidates_.at(i).size(); j++) {
      os << obj.tangent_x_candidates_.at(i).at(j);
      if (j + 1 != obj.tangent_x_candidates_.at(i).size()) os << ", ";
    }
    os << " )\n";
  }
  os << "\n"
     << "Normalization factor = " << obj.normalization_ << "\n"
     << "Total cross section = " << obj.total_cross_section_ << "\n";
  os << "\n"
     << "Number of true particles = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : "
       << obj.number_of_true_particles_.at(i);
    if(i != obj.number_of_frost_match_entries_ - 1) os << ", ";
  }
  os << "Particld PDG code = ";
  for (int i = 0; i < obj.number_of_frost_match_entries_;i++) {
    os << i + 1 << " : (";
    for (int j = 0; j < obj.number_of_true_particles_.at(i); j++) {
      os << obj.true_particle_id_.at(i).at(j);
      if (j != obj.number_of_true_particles_.at(i) - 1) os << ", ";
    }
    os << ")";
  }
  os << "True position = ";
  for (int i =0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : ";
      for (int j = 0; j < obj.number_of_true_particles_.at(i); j++) {
	os << j + 1 << " : ( "
	   << obj.true_position_.at(i).at(j).at(0) << ", "
	   << obj.true_position_.at(i).at(j).at(1) << ")\n";
      }
  }
  os << "True tangent = ";
  for (int i =0; i < obj.number_of_frost_match_entries_; i++) {
    os << i + 1 << " : ";
      for (int j = 0; j < obj.number_of_true_particles_.at(i); j++) {
	os << j + 1 << " : ( "
	   << obj.true_tangent_.at(i).at(j).at(0) << ", "
	   << obj.true_tangent_.at(i).at(j).at(1) << ")\n";
      }
  }

  return os;
}

void NTBMSummary::SetEntryInDailyFile(int entry_in_daily_file) {
  entry_in_daily_file_ = entry_in_daily_file;
}

int NTBMSummary::GetEntryInDailyFile() const {
  return entry_in_daily_file_;
}

void NTBMSummary::SetSpillPot(double spill_pot) {
  spill_pot_ = spill_pot;
}

double NTBMSummary::GetSpillPot() const {
  return spill_pot_;
}

void NTBMSummary::SetBunchPot(int bunch, double bunch_pot) {
  bunch_pot_[bunch] = bunch_pot;
}

double NTBMSummary::GetBunchPot(int bunch) const {
  if(bunch >= 0 && bunch < NUMBER_OF_BUNCHES)
    return bunch_pot_[bunch];
  else
    throw std::out_of_range("Bunch number out of range");
}

void NTBMSummary::SetBsdSpillNumber(int bsd_spill_number) {
  bsd_spill_number_ = bsd_spill_number;
}

int NTBMSummary::GetBsdSpillNumber() const {
  return bsd_spill_number_;
}

void NTBMSummary::SetTimestamp(double timestamp) {
  timestamp_ = timestamp;
}

double NTBMSummary::GetTimestamp() const {
  return timestamp_;
}

void NTBMSummary::SetBsdGoodSpillFlag(int bsd_good_spill_flag) {
  bsd_good_spill_flag_ = bsd_good_spill_flag;
}

int NTBMSummary::GetBsdGoodSpillFlag() const {
  return bsd_good_spill_flag_;
}

void NTBMSummary::SetWagasciGoodSpillFlag(int wagasci_good_spill_flag) {
  wagasci_good_spill_flag_ = wagasci_good_spill_flag;
}

int NTBMSummary::GetWagasciGoodSpillFlag() const {
  return wagasci_good_spill_flag_;
}

void NTBMSummary::SetDetectorFlags(int detector, int detector_flag) {
  if (detector >= 0 && detector < 8)
    detector_flags_[detector] = detector_flag;
  else
    throw std::out_of_range("Detector id out of range");
}

int NTBMSummary::GetDetectorFlags(int detector) const {
  if (detector >= 0 && detector < 8)
    return detector_flags_[detector];
  else
    throw std::out_of_range("Detector id out of range");
}

void NTBMSummary::SetNumberOfTracks(int number_of_tracks) {
  number_of_tracks_ = number_of_tracks;
  // Always set number of tracks before set other elements
  // related to Baby MIND
  ninja_track_type_.resize(number_of_tracks_);
  momentum_type_.resize(number_of_tracks_);
  momentum_.resize(number_of_tracks_);
  momentum_error_.resize(number_of_tracks_);
  baby_mind_position_.resize(number_of_tracks_);
  baby_mind_position_error_.resize(number_of_tracks_);
  baby_mind_tangent_.resize(number_of_tracks_);
  baby_mind_tangent_error_.resize(number_of_tracks_);
  extrapolated_position_.resize(number_of_tracks_);
  for(int i = 0; i < number_of_tracks_; i++) {
    baby_mind_position_.at(i).resize(2);
    baby_mind_position_error_.at(i).resize(2);
    baby_mind_tangent_.at(i).resize(2);
    baby_mind_tangent_error_.at(i).resize(2);
    extrapolated_position_.at(i).resize(2);
  }
  baby_mind_maximum_plane_.resize(number_of_tracks_);
  track_length_total_.resize(number_of_tracks_);
  charge_.resize(number_of_tracks_);
  direction_.resize(number_of_tracks_);
  bunch_.resize(number_of_tracks_);
}

int NTBMSummary::GetNumberOfTracks() const {
  return number_of_tracks_;
}

void NTBMSummary::SetNinjaTrackType(int track, int ninja_track_type) {
  ninja_track_type_.at(track) = ninja_track_type;
}

int NTBMSummary::GetNinjaTrackType(int track) const {
  if (track > number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return ninja_track_type_.at(track);
}

void NTBMSummary::SetMomentumType(int track, int momentum_type) {
  momentum_type_.at(track) = momentum_type;
}

int NTBMSummary::GetMomentumType(int track) const {
  if (track >= number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return momentum_type_.at(track);
}

void NTBMSummary::SetMomentum(int track, double momentum) {
  momentum_.at(track) = momentum;
}

double NTBMSummary::GetMomentum(int track) const {
  if (track >= number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return momentum_.at(track);
}

void NTBMSummary::SetMomentumError(int track, double momentum_error) {
  momentum_error_.at(track) = momentum_error;
}

double NTBMSummary::GetMomentumError(int track) const {
  if (track >= number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return momentum_error_.at(track);
}

void NTBMSummary::SetBabyMindPosition(int track, int view, double baby_mind_position) {
  if (view >= NUMBER_OF_VIEWS)
    throw std::out_of_range("View out of range");
  baby_mind_position_.at(track).at(view) = baby_mind_position;
}

void NTBMSummary::SetBabyMindPosition(int track, std::vector<double> baby_mind_position) {
  for (std::size_t view = 0; view < NUMBER_OF_VIEWS; view++)
    SetBabyMindPosition(track, view, baby_mind_position.at(view));
}

std::vector<double> NTBMSummary::GetBabyMindPosition(int track) const {
  if (track >= number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return baby_mind_position_.at(track);
}

double NTBMSummary::GetBabyMindPosition(int track, int view) const {
  if (view >= NUMBER_OF_VIEWS)
    throw std::out_of_range("View out of range");
  return GetBabyMindPosition(track).at(view);
}

void NTBMSummary::SetBabyMindPositionError(int track, int view, double baby_mind_position_error) {
  if (view >= NUMBER_OF_VIEWS)
    throw std::out_of_range("View out of range");
  baby_mind_position_error_.at(track).at(view) = baby_mind_position_error;
}

void NTBMSummary::SetBabyMindPositionError(int track, std::vector<double> baby_mind_position_error) {
  for (std::size_t view = 0; view < NUMBER_OF_VIEWS; view++)
    SetBabyMindPositionError(track, view, baby_mind_position_error.at(view));
}

std::vector<double> NTBMSummary::GetBabyMindPositionError(int track) const {
  if (track >= number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return baby_mind_position_error_.at(track);
}

double NTBMSummary::GetBabyMindPositionError(int track, int view) const {
  if (view >= NUMBER_OF_VIEWS)
    throw std::out_of_range("View out of range");
  return GetBabyMindPositionError(track).at(view);
}

void NTBMSummary::SetBabyMindTangent(int track, int view, double baby_mind_tangent) {
  if (view >= NUMBER_OF_VIEWS)
    throw std::out_of_range("View out of range");
  baby_mind_tangent_.at(track).at(view) = baby_mind_tangent;
}

void NTBMSummary::SetBabyMindTangent(int track, std::vector<double> baby_mind_tangent) {
  for(std::size_t view = 0; view < NUMBER_OF_VIEWS; view++)
    SetBabyMindTangent(track, view, baby_mind_tangent.at(view));
}

std::vector<double> NTBMSummary::GetBabyMindTangent(int track) const {
  if (track >= number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return baby_mind_tangent_.at(track);
}

double NTBMSummary::GetBabyMindTangent(int track, int view) const {
  if (view > NUMBER_OF_VIEWS)
    throw std::out_of_range("View our of range");
  return GetBabyMindTangent(track).at(view);
}

void NTBMSummary::SetBabyMindTangentError(int track, int view, double baby_mind_tangent_error) {
  if (view >= NUMBER_OF_VIEWS)
    throw std::out_of_range("View out of range");
  baby_mind_tangent_error_.at(track).at(view) = baby_mind_tangent_error;
}

void NTBMSummary::SetBabyMindTangentError(int track, std::vector<double> baby_mind_tangent_error) {
  for(std::size_t view = 0; view < NUMBER_OF_VIEWS; view++)
    SetBabyMindTangentError(track, view, baby_mind_tangent_error.at(view));
}

std::vector<double> NTBMSummary::GetBabyMindTangentError(int track) const {
  if (track >= number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return baby_mind_tangent_error_.at(track);
}

double NTBMSummary::GetBabyMindTangentError(int track, int view) const {
  if (view > NUMBER_OF_VIEWS)
    throw std::out_of_range("View our of range");
  return GetBabyMindTangentError(track).at(view);
}

void NTBMSummary::SetExtrapolatedPosition(int track, int view, double extrapolated_position) {
  if (view >= NUMBER_OF_VIEWS)
    throw std::out_of_range("View out of range");
  extrapolated_position_.at(track).at(view) = extrapolated_position;
}

void NTBMSummary::SetExtrapolatedPosition(int track, std::vector<double> extrapolated_position) {
  for (std::size_t view = 0; view < NUMBER_OF_VIEWS; view++)
    SetExtrapolatedPosition(track, view, extrapolated_position.at(view));
}

std::vector<double> NTBMSummary::GetExtrapolatedPosition(int track) const {
  if (track >= number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return extrapolated_position_.at(track);
}

double NTBMSummary::GetExtrapolatedPosition(int track, int view) const {
  if (view >= NUMBER_OF_VIEWS)
    throw std::out_of_range("View out of range");
  return GetExtrapolatedPosition(track).at(view);
}

void NTBMSummary::SetBabyMindMaximumPlane(int track, int baby_mind_maximum_plane) {
  baby_mind_maximum_plane_.at(track) = baby_mind_maximum_plane;
}

int NTBMSummary::GetBabyMindMaximumPlane(int track) const {
  if (track >= number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return baby_mind_maximum_plane_.at(track);
}

void NTBMSummary::SetTrackLengthTotal(int track, double track_length_total) {
  track_length_total_.at(track) = track_length_total;
}

double NTBMSummary::GetTrackLengthTotal(int track) const {
  if (track >= number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return track_length_total_.at(track);
}

void NTBMSummary::SetCharge(int track, int charge) {
  charge_.at(track) = charge;
}

int NTBMSummary::GetCharge(int track) const {
  if (track >= number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return charge_.at(track);
}

void NTBMSummary::SetDirection(int track, int direction) {
  direction_.at(track) = direction;
}

int NTBMSummary::GetDirection(int track) const {
  if (track >= number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return direction_.at(track);
}

void NTBMSummary::SetBunch(int track, int bunch) {
  bunch_.at(track) = bunch;
}

int NTBMSummary::GetBunch(int track) const {
  if (track >= number_of_tracks_)
    throw std::out_of_range("Number of track out of range");
  return bunch_.at(track);
}

void NTBMSummary::SetNumberOfFrostMatchEntries(int number_of_frost_match_entries) {
  number_of_frost_match_entries_ = number_of_frost_match_entries;
  baby_mind_track_id_.resize(number_of_frost_match_entries_);
  frost_match_bunch_.resize(number_of_frost_match_entries_);
  number_of_frost_y_candidates_.resize(number_of_frost_match_entries_);
  number_of_frost_x_candidates_.resize(number_of_frost_match_entries_);
  frost_y_candidates_.resize(number_of_frost_match_entries_);
  frost_x_candidates_.resize(number_of_frost_match_entries_);
  expected_y_candidates_.resize(number_of_frost_match_entries_);
  expected_x_candidates_.resize(number_of_frost_match_entries_);
  difference_y_candidates_.resize(number_of_frost_match_entries_);
  difference_x_candidates_.resize(number_of_frost_match_entries_);
  tangent_y_candidates_.resize(number_of_frost_match_entries_);
  tangent_x_candidates_.resize(number_of_frost_match_entries_);
  number_of_true_particles_.resize(number_of_frost_match_entries_);
  true_particle_id_.resize(number_of_frost_match_entries_);
  true_position_.resize(number_of_frost_match_entries_);
  true_tangent_.resize(number_of_frost_match_entries_);
}

int NTBMSummary::GetNumberOfFrostMatchEntries() const {
  return number_of_frost_match_entries_;
}

void NTBMSummary::SetBabyMindTrackId(int entry, int baby_mind_track_id) {
  baby_mind_track_id_.at(entry) = baby_mind_track_id;
}

int NTBMSummary::GetBabyMindTrackId(int entry) const {
  if (entry >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return baby_mind_track_id_.at(entry);
}

void NTBMSummary::SetFrostMatchBunch(int entry, int frost_match_bunch) {
  frost_match_bunch_.at(entry) = frost_match_bunch;
}

int NTBMSummary::GetFrostMatchBunch(int entry) const {
  if (entry >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return frost_match_bunch_.at(entry);
}

void NTBMSummary::SetNumberOfFrostYCandidates(int entry, int number_of_frost_y_candidates) {
  number_of_frost_y_candidates_.at(entry) = number_of_frost_y_candidates;
}

int NTBMSummary::GetNumberOfFrostYCandidates(int entry) const {
  if (entry >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return number_of_frost_y_candidates_.at(entry);
}

void NTBMSummary::SetNumberOfFrostXCandidates(int entry, int number_of_frost_x_candidates) {
  number_of_frost_x_candidates_.at(entry) = number_of_frost_x_candidates;
}

int NTBMSummary::GetNumberOfFrostXCandidates(int entry) const {
  if (entry >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return number_of_frost_x_candidates_.at(entry);
}

void NTBMSummary::SetFrostYCandidates(int entry, std::vector<double> frost_y_candidates) {
  frost_y_candidates_.at(entry) = frost_y_candidates;
}

std::vector<double> NTBMSummary::GetFrostYCandidates(int entry) const {
  if (entry >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return frost_y_candidates_.at(entry);
}

void NTBMSummary::SetFrostXCandidates(int entry, std::vector<double> frost_x_candidates) {
  frost_x_candidates_.at(entry) = frost_x_candidates;
}

std::vector<double> NTBMSummary::GetFrostXCandidates(int entry) const {
  if (entry >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return frost_x_candidates_.at(entry);
}

void NTBMSummary::SetExpectedYCandidates(int entry, std::vector<double> expected_y_candidates) {
  expected_y_candidates_.at(entry) = expected_y_candidates;
}

std::vector<double> NTBMSummary::GetExpectedYCandidates(int entry) const {
  if (entry >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return expected_y_candidates_.at(entry);
}

void NTBMSummary::SetExpectedXCandidates(int entry, std::vector<double> expected_x_candidates) {
  expected_x_candidates_.at(entry) = expected_x_candidates;
}

std::vector<double> NTBMSummary::GetExpectedXCandidates(int entry) const {
  if (entry >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return expected_x_candidates_.at(entry);
}

void NTBMSummary::SetDifferenceYCandidates(int entry, std::vector<double> difference_y_candidates) {
  difference_y_candidates_.at(entry) = difference_y_candidates;
}

std::vector<double> NTBMSummary::GetDifferenceYCandidates(int entry) const {
  if (entry >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return difference_y_candidates_.at(entry);
}

void NTBMSummary::SetDifferenceXCandidates(int entry, std::vector<double> difference_x_candidates) {
  difference_x_candidates_.at(entry) = difference_x_candidates;
}

std::vector<double> NTBMSummary::GetDifferenceXCandidates(int entry) const {
  if (entry >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return difference_x_candidates_.at(entry);
}

void NTBMSummary::SetTangentYCandidates(int entry, std::vector<double> tangent_y_candidates) {
  tangent_y_candidates_.at(entry) = tangent_y_candidates;
}

std::vector<double> NTBMSummary::GetTangentYCandidates(int entry) const {
  if (entry >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return tangent_y_candidates_.at(entry);
}

void NTBMSummary::SetTangentXCandidates(int entry, std::vector<double> tangent_x_candidates) {
  tangent_x_candidates_.at(entry) = tangent_x_candidates;
}

std::vector<double> NTBMSummary::GetTangentXCandidates(int entry) const {
  if (entry >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return tangent_x_candidates_.at(entry);
}

void NTBMSummary::SetNumberOfTrueParticles(int cluster, int number_of_true_particles) {
  number_of_true_particles_.at(cluster) = number_of_true_particles;
  true_particle_id_.at(cluster).resize(number_of_true_particles_.at(cluster));
  true_position_.at(cluster).resize(number_of_true_particles_.at(cluster));
  true_tangent_.at(cluster).resize(number_of_true_particles_.at(cluster));
  for ( int iparticle = 0; iparticle < number_of_true_particles_.at(cluster); iparticle++ ) {
    true_position_.at(cluster).at(iparticle).resize(2);
    true_tangent_.at(cluster).at(iparticle).resize(2);
  }
}

int NTBMSummary::GetNumberOfTrueParticles(int cluster) const {
  if (cluster >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return number_of_true_particles_.at(cluster);
}

void NTBMSummary::SetTrueParticleId(int cluster, int particle, int true_particle_id) {
  true_particle_id_.at(cluster).at(particle) = true_particle_id;
}

void NTBMSummary::SetTrueParticleId(int cluster, std::vector<int> true_particle_id) {
  for (int particle = 0; particle < number_of_true_particles_.at(cluster); particle++)
    SetTrueParticleId(cluster, particle, true_particle_id.at(particle));
}

std::vector<int> NTBMSummary::GetTrueParticleId(int cluster) const {
  return true_particle_id_.at(cluster);
}

int NTBMSummary::GetTrueParticleId(int cluster, int particle) const {
  return GetTrueParticleId(cluster).at(particle);
}

void NTBMSummary::SetNormalization(double normalization) {
  normalization_ = normalization;
}

double NTBMSummary::GetNormalization() const {
  return normalization_;
}

void NTBMSummary::SetTotalCrossSection(double total_cross_section) {
  total_cross_section_ = total_cross_section;
}

double NTBMSummary::GetTotalCrossSection() const {
  return total_cross_section_;
}

void NTBMSummary::SetTruePosition(int cluster, int particle, int view, double true_position) {
  if (view >= NUMBER_OF_VIEWS)
    throw std::out_of_range("View out of range");
  true_position_.at(cluster).at(particle).at(view) = true_position;
}

void NTBMSummary::SetTruePosition(int cluster, int particle, std::vector<double> true_position) {
  for(std::size_t view = 0; view < NUMBER_OF_VIEWS; view++)
    SetTruePosition(cluster, particle, view, true_position.at(view));
}

void NTBMSummary::SetTruePosition(int cluster, std::vector<std::vector<double>> true_position) {
  for(int particle = 0; particle < number_of_true_particles_.at(cluster); particle++)
    SetTruePosition(cluster, particle, true_position.at(particle));
}

std::vector<std::vector<double>> NTBMSummary::GetTruePosition(int cluster) const {
  if (cluster >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return true_position_.at(cluster);
}

std::vector<double> NTBMSummary::GetTruePosition(int cluster, int particle) const {
  if (particle >= number_of_true_particles_.at(cluster))
    throw std::out_of_range("Number of true particles our of range");
  return GetTruePosition(cluster).at(particle);
}

double NTBMSummary::GetTruePosition(int cluster, int particle, int view) const {
  if (view >= NUMBER_OF_VIEWS)
    throw std::out_of_range("View out of range");
  return GetTruePosition(cluster, particle).at(view);
}

void NTBMSummary::SetTrueTangent(int cluster, int particle, int view, double true_tangent) {
  if (view >= NUMBER_OF_VIEWS)
    throw std::out_of_range("View ourt of range");
  true_tangent_.at(cluster).at(particle).at(view) = true_tangent;
}

void NTBMSummary::SetTrueTangent(int cluster, int particle, std::vector<double> true_tangent) {
  for (std::size_t view = 0; view < NUMBER_OF_VIEWS; view++)
    SetTrueTangent(cluster, particle, view, true_tangent.at(view));
}

void NTBMSummary::SetTrueTangent(int cluster, std::vector<std::vector<double>> true_tangent) {
  for (int particle = 0; particle < number_of_true_particles_.at(cluster); particle++)
    SetTrueTangent(cluster, particle, true_tangent.at(particle));
}

std::vector<std::vector<double>> NTBMSummary::GetTrueTangent(int cluster) const {
  if (cluster >= number_of_frost_match_entries_)
    throw std::out_of_range("Number of frost match entry out of range");
  return true_tangent_.at(cluster);
}

std::vector<double> NTBMSummary::GetTrueTangent(int cluster, int particle) const {
  if (particle >= number_of_true_particles_.at(cluster))
    throw std::out_of_range("Number of true particles out of range");
  return GetTrueTangent(cluster).at(particle);
}

double NTBMSummary::GetTrueTangent(int cluster, int particle, int view) const {
  if (view >= NUMBER_OF_VIEWS)
    throw std::out_of_range("View out of range");
  return GetTrueTangent(cluster, particle).at(view);
}

ClassImp(NTBMSummary)
