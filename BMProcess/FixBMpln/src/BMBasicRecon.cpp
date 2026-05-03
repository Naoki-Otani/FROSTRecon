#include <iostream>
#include <stdio.h>

#include "BMBasicRecon.hpp"

using namespace std;

void BMBasicRecon::Clear() {
  mod.clear();
  view.clear();
  pln.clear();
  channel.clear();
  HG.clear();
  LHG.clear();
  RHG.clear();
  THG.clear();
  BHG.clear();
  Lgain.clear();
  Rgain.clear();
  Tgain.clear();
  Bgain.clear();
  Lpe.clear();
  Rpe.clear();
  Tpe.clear();
  Bpe.clear();
  LG.clear();
  Ltime.clear();
  Ftime.clear();
  Htime.clear();
  timedif.clear();
  bunch.clear();

  bm_event = 0;
  year = 0;
  date = 0;
  mon = 0;
  run = 0;
  hour = 0;
}
