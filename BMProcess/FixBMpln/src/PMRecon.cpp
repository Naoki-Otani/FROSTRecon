#include <iostream>

#include "PMRecon.hpp"

PMRecon::PMRecon() {
  mod = 0;
  view = 0;
  pln = 0;
  channel = 0;
  lope = 0;
  pe = 0;
  charge = 0;
  time = 0;
  bunch = 0;

  spill = -1;
  unixtime = -1;
}

PMRecon::~PMRecon() {}

void PMRecon::Clear() {
  mod = 0;
  view = 0;
  pln = 0;
  channel = 0;
  lope = 0;
  pe = 0;
  charge = 0;
  time = 0;
  bunch = 0;

  spill = -1;
  unixtime = -1;
}
