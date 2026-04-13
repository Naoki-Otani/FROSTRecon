#ifndef NTBMCONST_HH
#define NTBMCONST_HH

//HitConverter for FROST
static const Int_t SPILL_MOD = 32768;      // 2^15
static const Int_t MAX_TIME_DIFF = 3000;   // sec allowance between BM and FROST unixtime

///> Default non initialized value
static const int NTBM_NON_INITIALIZED_VALUE = -2;

///> Number of J-PARC neutrino beam bunches
static const int NUMBER_OF_BUNCHES = 8;

///> Number of views
static const int NUMBER_OF_VIEWS = 2;

///> Baby MIND vertical scintillator overlap
static const double BM_VERTICAL_SCINTI_OVERLAP = 35.5; // mm
// Nominal position of the 2nd layer in the BM coordinate
// Offset + Scin_Mod_position.txt(2) + 5 cm (?) + 1.5 scintillators
static const double BM_SECOND_LAYER_POS = -2000. + 183. + 50. + 15.; // mm
// Scintillator position used in the WAGASCI is not well corrected for our usage
// plane more than 2nd should be much corrected.
static const double BM_SCI_CORRECTION = -31.5; // mm
static const double TEMPORAL_ALLOWANCE[2] = {200., 300.}; // mm
//static const double TEMPORAL_ALLOWANCE[2] = {500., 500.}; // mm 20241022

#endif
