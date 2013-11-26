#ifndef UTIL_HPP
#define UTIL_HPP
#include <iostream>

template<class T>
class tct2x3
{
public:
  //  tct2x3():_sum_cal(false){};
  tct2x3<T>& operator = (const tct2x3<T> &rh){
    std::copy(rh.data, rh.data + 6, data);
    sum_cal();
    return *this;
  }
  inline T& operator()(unsigned row, unsigned col){
    return data[row][col];
  }

  inline const T& operator()(unsigned row, unsigned col) const{
    return data[row][col];
  }

  
  inline const T& sum_row(unsigned ind) const{
    return _sum_row[ind];
  }
  inline const T& sum_col(unsigned ind) const{
    return _sum_col[ind];
  }
  inline const T& sum_all() const{
    return _sum_all;
  }


  //this function should be called before call sum_row / sum_col /sum_all
  inline void sum_cal(){
    _sum_row[0] = data[0][0]+data[0][1]+data[0][2];
    _sum_row[1] = data[1][0]+data[1][1]+data[1][2];
    _sum_col[0] = data[0][0]+data[1][0];
    _sum_col[1] = data[0][1]+data[1][1];
    _sum_col[2] = data[0][2]+data[1][2];
    _sum_all = _sum_row[0] + _sum_row[1];
  }
private:
  T data[2][3];
  T _sum_row[2];
  T _sum_col[3];
  T _sum_all;
  //  bool _sum_cal;
};

typedef tct2x3<int> ct2x3; 

template<class T>
std::ostream & operator<<(std::ostream &os, tct2x3<T> tc){
  for(unsigned i=0; i<2; ++i){
    for(unsigned j=0; j<3; ++j){
      os<<tc(i,j)<<'\t';
    }
    os<<"\n";
  }
  return os;
}

#endif
