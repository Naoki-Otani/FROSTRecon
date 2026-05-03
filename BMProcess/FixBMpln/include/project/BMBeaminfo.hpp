#ifndef _INC_BMBEAMINFO
#define _INC_BMBEAMINFO

#include <iostream>
#include <vector>

#include <BMConst.hpp>

class BMBeaminfo {

private:
public:
  // Baby MIND DATA STRUCTURE
  int spillnum;
  int spillid;
  int extended_spillnum;

  BMBeaminfo();
  ~BMBeaminfo();
};

#endif
