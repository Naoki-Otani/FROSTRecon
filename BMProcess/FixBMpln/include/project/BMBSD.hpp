#ifndef _INC_BMBSD
#define _INC_BMBSD

#include <iostream>
#include <vector>
#include <cstdint>

#include <BMConst.hpp>

class BMBSD {

private:
public:
  // Baby MIND DATA STRUCTURE
  double pot;
  int unixtime;
  int goodspillflag;
  std::uint32_t bsd_spill_number;

  BMBSD();
  ~BMBSD();
};

#endif
