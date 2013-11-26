#include "maxact.h"
#include <cmath>
#include <vector>
#include <algorithm>

// use lgammafn from Rmath instead of lgamma, as it is not a function
// in C++98 although it is available in C++11
#include <Rmath.h>
#define lgamma lgammafn

//#include <boost/math/special_functions/gamma.hpp>

using std::ceil;
using std::floor;
using std::exp;




//A simple class range is used instead of boost::interval
//#include <boost/numeric/interval.hpp>
// using boost::numeric::interval_lib::policies;
// typedef boost::numeric::interval<double,
// 				 boost::numeric::interval_lib::policies<
// 				   boost::numeric::interval_lib::rounded_math<double>, 
// 				   boost::numeric::interval_lib::checking_no_nan<double> > > range;



using std::sqrt;

inline range getRange(const ct2x3 &ct, double theta, double d, enum Alternative alter=TWO_SIDE)
{
  double mid, wid;
  double tmp = ct.sum_col(2) + theta*ct.sum_col(1);
  mid = tmp*ct.sum_row(1)/ct.sum_all();
  wid = d*sqrt(ct.sum_row(0)*ct.sum_row(1)*
	       (ct.sum_col(2)+theta*theta*ct.sum_col(1) - tmp*tmp/ct.sum_all()))
    /ct.sum_all();
  range result;
  switch(alter){
  case TWO_SIDE:
    result=range(mid-wid, mid+wid);
    break;
  case GREATER:
    result=range(-1, mid+wid);
    break;
  case LESS:
    //wid and d have the same sign, so mid + wid NOT mid - wid
    result=range(mid+wid, ct.sum_row(1) + 1);   
  }
  return result;
}

double catt(const ct2x3 &ct, double theta)
{
  double t,var;
  double theta2 = theta * theta;
  //  double weights[3]={0,theta,1};
  double tmp = (ct.sum_col(2) + theta*ct.sum_col(1)) / ct.sum_all();
  t =
    (ct(1,2) + theta*ct(1,1)) / ct.sum_row(1)
    - tmp ;
  var = 
    ((ct.sum_col(2)+theta2*ct.sum_col(1))/ct.sum_all() -
     tmp*tmp)
    *ct.sum_row(0)/(ct.sum_row(1)*ct.sum_all());
  // cout<<((ct.sum_col(2)+theta2*ct.sum_col(1))/ct.sum_all() -
  // 	 pow((ct.sum_col(2)+theta*ct.sum_col(1))/ct.sum_all(),2))<<endl;
  // cout<<(((ct.sum_col(2)+theta2*ct.sum_col(1))/ct.sum_all() -
  // 	      pow((ct.sum_col(2)+theta*ct.sum_col(1))/ct.sum_all(),2))
  // 	 *ct.sum_row(0)/(ct.sum_row(1)*ct.sum_all()))<<endl;
  // cout<<"var "<<var<<endl;
  return t/sqrt(var);
}


//square of catt
inline double catt2(const ct2x3 &ct, double theta)
{
  double t,var;
  double theta2 = theta * theta;
  //  double weights[3]={0,theta,1};
  double tmp = (ct.sum_col(2) + theta*ct.sum_col(1)) / ct.sum_all();
  t =
    (ct(1,2) + theta*ct(1,1)) / ct.sum_row(1)
    - tmp ;
  var = 
    ((ct.sum_col(2)+theta2*ct.sum_col(1))/ct.sum_all() -
     tmp*tmp)
    *ct.sum_row(0)/(ct.sum_row(1)*ct.sum_all());
  return t*t/var;
}

//using boost::math::lgamma;

//log factorial
//the values are used repeately, so they are cached
class Lfactorial{
public:
  Lfactorial(unsigned cacheSize = 8192):
    _cacheSize(cacheSize),
    _cache(cacheSize)
  {
    for(unsigned i=0; i<_cacheSize; ++i){
      _cache[i] = lgamma(double(i+1));
    }
  }
  inline double operator()(unsigned z)
  {
    if(z<_cacheSize) {
      return _cache[z];
    }
    return lgamma(double(z+1));
  }
private:
  unsigned _cacheSize;
  std::vector<double> _cache;
};

Lfactorial lfactCached;

inline double lfactorial(unsigned z)
{
  return lgamma(double(z+1));
}

//double (*lf)(unsigned)=lfactorial;
//double (*lf)(unsigned)=lfactCached;
inline double lf(unsigned z)
{
  return lfactCached(z);
}

//the float number to release the tie, which might be caused by float error
//I'm not sure whether it is nessesary, but using it should be harmless
const double float_tie = 1e-10;

MaXact::MaXact(const ct2x3 &ct, bool useMax3, Alternative alter){
  _ct = ct;
  _useMax3 = useMax3;
  _alter = alter;
  _d = catt_max(_ct, _useMax3, _alter);
  _dAdjust = _d;
  if(_alter==LESS){
    _dAdjust += std::abs(_dAdjust)*float_tie;
  }else{
    _dAdjust -= std::abs(_dAdjust)*float_tie;
  }
  _rns[2] = getRange(_ct, 0, _dAdjust, _alter);   //s2
  _rns[1] = getRange(_ct, 0.5, _dAdjust, _alter); //s2 + 0.5*s1
  _rns[0] = getRange(_ct, 1, _dAdjust, _alter);   //s2 + s1
  using std::min;
  using std::max;
  _rbys[2] = range(max(0, _ct.sum_col(2) - _ct.sum_row(0)),
                  min(_ct.sum_col(2), _ct.sum_row(1)));
  _rbys[1] = range(max(0, _ct.sum_col(1) - _ct.sum_row(0)), 
                  min(_ct.sum_col(1), _ct.sum_row(1)));
  _rbys[0] = range(max(0, _ct.sum_col(0) - _ct.sum_row(0)),  
                  min(_ct.sum_col(0), _ct.sum_row(1)));
  
   _logPathAll = 
    lf(_ct.sum_all())
    - lf(_ct.sum_row(0))
    - lf(_ct.sum_row(1));
   _logConst =       
    lf(_ct.sum_col(2))
    +lf(_ct.sum_col(1))
    +lf(_ct.sum_col(0))
    - _logPathAll;
 
   if (proportionAcceptPoint() < 0.5){
     calculateP();
   } else{
     calculateSmallP();
   }
}

double MaXact::proportionAcceptPoint(){
  double totalPoint = 0;
  double acceptPoint = 0;
  range rn2 = _rns[2];
  for(int i = ceil(lower(_rbys[2])); i <= upper(_rbys[2]); ++i){
      range rn1;
      rn1 = _rns[0] - (double)i;
      if(_useMax3) rn1 = intersect(rn1, 2.0*(_rns[1]-(double)i));

      range rangeJ = intersect(_rbys[1], 
                               double(_ct.sum_row(1)-i) - _rbys[0]);
      rn1 = intersect(rn1, rangeJ);
      if (!empty(rangeJ)){
        // Number of points in the close interval rangeJ
        totalPoint += floor(upper(rangeJ)) - ceil(lower(rangeJ)) + 1;
      }
      if (!empty(rn1) && (i > lower(rn2)) && (i < upper(rn2))){
        // Number of points in the open interval rn1
        acceptPoint += ceil(upper(rn1)) - floor(lower(rn1)) - 1;
      }
  }
  if (totalPoint < 1){
    return 0.5;
  }else{
    return(acceptPoint / totalPoint);
  }
}

// When the test statistic is large and p is small, it is faster to
// sum the probabilty of events with large statistic
void MaXact::calculateSmallP()
{
  double sum_col01 = _ct.sum_col(0) + _ct.sum_col(1);
  double res = 0;
  range rn2 = _rns[2];
  for(int i = ceil(lower(_rbys[2])); i <= upper(_rbys[2]); ++i){
    if( i <= lower(rn2) || i >= upper(rn2) ){
      double tmp = _ct.sum_row(1) - i;
      res += exp(lf(_ct.sum_col(2)) - lf(i) - lf(_ct.sum_col(2) - i) +
                 lf(sum_col01) - lf(tmp) - lf(sum_col01 - tmp)
                 - _logPathAll);
    }else{
      double logPath2 = 
        //   lf(_ct.sum_col(2)) //moved into logConst
        - lf(i)
        - lf(_ct.sum_col(2)-i);
      range rn1;
      rn1 = _rns[0] - (double)i;
      if(_useMax3) rn1 = intersect(rn1, 2.0*(_rns[1]-(double)i));

      range rangeJ = intersect(_rbys[1], 
                               double(_ct.sum_row(1)-i) - _rbys[0]);
      for(int j = ceil(lower(rangeJ)); j <= upper(rangeJ); ++j){
        if( j > lower(rn1) && j < upper(rn1) ){
          j = ceil(upper(rn1));
          if (j > upper(rangeJ)) 
            break;
        }
        int k=_ct.sum_row(1)-i-j;
        double logPath1 =
          //lf(_ct.sum_col(1)) //moved into logConst
          - lf(j)
          - lf(_ct.sum_col(1)-j)
          //+ lf(_ct.sum_col(0)) //moved into logConst
          - lf(k)
          - lf(_ct.sum_col(0)-k);
        res += exp(logPath1 + logPath2 + _logConst);
      }
    }
  }

  //the result might be out of range slightly, because of the float error
  if(res>1){
    res=1;
  }else if(res<0){
    res=0;
  }
  _p = res;
}

void MaXact::calculateP()
{
  double res = 1;
  range rn2 = intersect(_rns[2], _rbys[2]);
  for (int i = ceil(lower(_rbys[2]));  i <= upper(_rbys[2]); ++i){
    if ( i <= lower(_rns[2])){
      // The largest number that should skip the loop
      i = floor(lower(_rns[2]));
      continue;
    }else if (i >= upper(_rns[2])){
      break;
    }
    double logPath2 = 
      //   lf(_ct.sum_col(2)) //moved into logConst
      - lf(i)
      - lf(_ct.sum_col(2)-i);
    range rn1 = intersect(_rbys[1], double(_ct.sum_row(1)-i)-_rbys[0]);
    // range from statistic
    range rns1 =  _rns[0] - (double)i;
    if(_useMax3) rns1 = intersect(rns1, 2.0*(_rns[1]-(double)i));
    // rn1 is a close range, while rns1 is an open range. If the
    // intersect is empty, it's conserve to skip the current loop
    if (empty(intersect(rn1, rns1))) continue;
    for(int j = ceil(lower(rn1));
        j <= upper(rn1);
        ++j){
      if ( j <= lower(rns1)){
        // The largest number that should skip the loop
        j = floor(lower(rns1));
        continue;
      }else if (j >= upper(rns1) ){
        break;
      }
      int k=_ct.sum_row(1)-i-j;
      double logPath1 =
	//lf(_ct.sum_col(1)) //moved into logConst
	- lf(j)
	- lf(_ct.sum_col(1)-j)
	//+ lf(_ct.sum_col(0)) //moved into logConst
	- lf(k)
	- lf(_ct.sum_col(0)-k);
      res -= exp(logPath1 + logPath2 + _logConst);
    }
  }

  //the result might be out of range slightly, because of the float error
  if(res>1){
    res=1;
  }else if(res<0){
    res=0;
  }
  _p = res;
}

double catt_p(const ct2x3 &ct, double theta, enum Alternative alter)
{
  range rnd; //constrain by catt
  double d = catt(ct, theta);
  if(alter==TWO_SIDE) d=std::abs(d);
  if(alter==LESS){
    d += std::abs(d)*float_tie;
  }else{
    d -= std::abs(d)*float_tie;
  }
  rnd = getRange(ct, theta, d, alter);   //s2 + theta*s1
  range rbys[3]; //constrain by sum_col, sum_row
  using std::min;
  rbys[2] = range(0, min(ct.sum_col(2), ct.sum_row(1)));
  rbys[1] = range(0, min(ct.sum_col(1), ct.sum_row(1)));
  rbys[0] = range(0, min(ct.sum_col(0), ct.sum_row(1)));
  
  double logPathAll = 
    lf(ct.sum_all())
    - lf(ct.sum_row(0))
    - lf(ct.sum_row(1));
  double logConst =       
    lf(ct.sum_col(2))
    +lf(ct.sum_col(1))
    +lf(ct.sum_col(0))
    -logPathAll;
  double res = 1;
  range rn2 =  rbys[2];
  if(empty(rn2)) return 1;
  for(int low=ceil(lower(rn2)), up=floor(upper(rn2)), i=low;
      i<=up;
      ++i){
    double logPath2 = 
      //   lf(ct.sum_col(2)) //moved into logConst
      - lf(i)
      - lf(ct.sum_col(2)-i);
    range rn1;
    rn1 = intersect(rbys[1], double(ct.sum_row(1)-i)-rbys[0]);
    if(empty(rn1)) continue;
    for(int low1=ceil(lower(rn1)), up1=floor(upper(rn1)), j=low1;
	j<=up1;
	++j){
      double tmp = theta*j + i;
      if( tmp < rnd.upper() && tmp > rnd.lower() ){
	int k=ct.sum_row(1)-i-j;
	double logPath1 =
	  //lf(ct.sum_col(1)) //moved into logConst
	  - lf(j)
	  - lf(ct.sum_col(1)-j)
	  //+ lf(ct.sum_col(0)) //moved into logConst
	  - lf(k)
	  - lf(ct.sum_col(0)-k);
	res -= exp(logPath1 + logPath2 + logConst);
      }
    }
  }

  //the result might be out of range slightly, because of the float error
  if(res>1){
    res=1;
  }else if(res<0){
    res=0;
  }
  return res;
}
