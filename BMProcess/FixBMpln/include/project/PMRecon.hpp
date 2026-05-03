#ifndef _INC_PMRecon
#define _INC_PMRecon

#include <iostream>
#include <vector>

#include <BMConst.hpp>

class PMRecon {
public:
  vector<double> *mod;
  vector<double> *view;
  vector<double> *pln;
  vector<double> *channel;
  vector<double> *lope;
  vector<double> *pe;
  vector<double> *charge;
  vector<double> *time;
  vector<double> *bunch;

  int spill;
  int unixtime;

  PMRecon();
  ~PMRecon();

  void Clear();
};

#endif
