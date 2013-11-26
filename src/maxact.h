#ifndef MAXACT_HPP
#define MAXACT_hpp

#include "util.h"
#include <algorithm>
#include <iostream>
#include <cmath>


#include "range.h"
typedef range_base<double> range;

double catt(const ct2x3 &ct, double theta);


enum Alternative {
  TWO_SIDE,
  GREATER,
  LESS
};



inline double catt_max(const ct2x3 &ct, bool useMax3 = true, Alternative alter=TWO_SIDE)
{
  double res;
  switch(alter){
  case TWO_SIDE:
    res = std::max(std::abs(catt(ct,0)), std::abs(catt(ct,1)));
    if(useMax3) res = std::max(std::abs(catt(ct,0.5)), res);
    break;
  case GREATER:
    res = std::max(catt(ct,0), catt(ct,1));
    if(useMax3) res = std::max(catt(ct,0.5), res);
    break;
  case LESS:
    res = std::min(catt(ct,0), catt(ct,1));
    if(useMax3) res = std::min(catt(ct,0.5), res);
    break;
  }
  return res;
}

double catt_p(const ct2x3 &ct, double theta, Alternative alter=TWO_SIDE);

class MaXact
{
 private:
  ct2x3 _ct;
  bool _useMax3;
  Alternative _alter;
  double _p, _d; // d is the test statistics
  double _dAdjust;
  range _rns[3];
  range _rbys[3]; //constrain by sum_col, sum_row
  double _logPathAll, _logConst;
  double proportionAcceptPoint();
  void calculateSmallP();
  void calculateP();
public:
  MaXact(const ct2x3 &ct, bool useMax3 = true, Alternative alter=TWO_SIDE);
  inline double p(){
    return _p;
  }
  inline double stat(){
    return _d;
  }
};


#endif
