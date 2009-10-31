#include "maxact.h"
#include <cmath>
#include <vector>
#include <algorithm>
//#include <boost/math/special_functions/gamma.hpp>




#include "range.h"
typedef range_base<double> range;

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
    result=range(0, mid+wid);
    break;
  case LESS:
    //wid and d have the same sign, so mid + wid NOT mid - wid
    result=range(mid+wid, ct.sum_row(1));   
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
  Lfactorial(unsigned cacheSize = 2048):
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

double maxacter::p(const ct2x3 &ct, bool useMax3, enum Alternative alter)
{
  range rns[3]; //constrain by catt
  double d = catt_max(ct, useMax3, alter);
  if(alter==LESS){
    d += d*float_tie;
  }else{
    d -= d*float_tie;
  }
  rns[2] = getRange(ct, 0, d, alter);   //s2
  rns[1] = getRange(ct, 0.5, d, alter); //s2 + 0.5*s1
  rns[0] = getRange(ct, 1, d, alter);   //s2 + s1
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
  range rn2 = intersect(rns[2], rbys[2]);
  for(int low=ceil(lower(rn2)), up=floor(upper(rn2)), i=low;
      i<=up;
      ++i){
    double logPath2 = 
      //   lf(ct.sum_col(2)) //moved into logConst
      - lf(i)
      - lf(ct.sum_col(2)-i);
    range rn1;
    rn1 = rns[0] - (double)i;
    rn1 = intersect(rn1, rbys[1]);
    rn1 = intersect(rn1, double(ct.sum_row(1)-i)-rbys[0]);
    if(useMax3) rn1 = intersect(rn1, 2.0*(rns[1]-(double)i));
    if(empty(rn1)) continue;
    for(int low1=ceil(lower(rn1)), up1=floor(upper(rn1)), j=low1;
	j<=up1;
	++j){
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

  //the result might be out of range slightly, because of the float error
  if(res>1){
    res=1;
  }else if(res<0){
    res=0;
  }
  return res;
}

double catt_p(const ct2x3 &ct, double theta, enum Alternative alter)
{
  range rnd; //constrain by catt
  double d = catt(ct, theta);
  if(alter==TWO_SIDE) d=std::abs(d);
  if(alter==LESS){
    d += d*float_tie;
  }else{
    d -= d*float_tie;
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
      if( tmp<=rnd.upper() && tmp>=rnd.lower() ){
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
