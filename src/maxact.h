#ifndef MAXACT_HPP
#define MAXACT_hpp

#include "util.h"
#include <algorithm>
#include <iostream>
#include <cmath>



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

class maxacter
{
public:
  double p(const ct2x3 &ct, bool useMax3 = true, Alternative alter=TWO_SIDE);
};


#endif
