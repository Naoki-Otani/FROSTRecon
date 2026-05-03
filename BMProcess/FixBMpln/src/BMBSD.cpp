#include <iostream>
#include <stdio.h>

#include "BMBSD.hpp"

using namespace std;

BMBSD::BMBSD() {
  pot = -1;
  unixtime = -1;
  goodspillflag = -999;
  bsd_spill_number = 0;
}

BMBSD::~BMBSD() {}
