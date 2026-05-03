#ifndef _INC_BMBASICRECON
#define _INC_BMBASICRECON

#include <vector>

#include <BMConst.hpp>

struct BMBasicRecon {
  // Baby MIND DATA STRUCTURE
  std::vector<double> mod;
  std::vector<double> view;
  std::vector<double> pln;
  std::vector<double> channel;
  std::vector<double> HG;
  std::vector<double> LHG;
  std::vector<double> RHG;
  std::vector<double> THG;
  std::vector<double> BHG;
  std::vector<double> Lgain;
  std::vector<double> Rgain;
  std::vector<double> Tgain;
  std::vector<double> Bgain;
  std::vector<double> Lpe;
  std::vector<double> Rpe;
  std::vector<double> Tpe;
  std::vector<double> Bpe;
  std::vector<double> LG;
  std::vector<double> Ltime;
  std::vector<double> Ftime;
  std::vector<double> Htime;
  std::vector<double> timedif;
  std::vector<double> bunch;

  int bm_event = 0;
  int year = 0;
  int date = 0;
  int mon = 0;
  int run = 0;
  int hour = 0;

  void Clear();
};
#endif
