#include "util.h"
#include "maxact.h"
extern "C" {
  maxacter maxact;
  void c_maxact_test(int *table, int *useMax3, int *alternative,double *stat, double *p){
    ct2x3 ct;
    ct(0,0)=table[0];
    ct(0,1)=table[1];
    ct(0,2)=table[2];
    ct(1,0)=table[3];
    ct(1,1)=table[4];
    ct(1,2)=table[5];
    ct.sum_cal();
    *stat = catt_max(ct, *useMax3, static_cast<Alternative>(*alternative));
    *p = maxact.p(ct, *useMax3, static_cast<Alternative>(*alternative));
  }
  
  void c_maxact(int *table, int *useMax3, int *alternative,double *stat){
    ct2x3 ct;
    ct(0,0)=table[0];
    ct(0,1)=table[1];
    ct(0,2)=table[2];
    ct(1,0)=table[3];
    ct(1,1)=table[4];
    ct(1,2)=table[5];
    ct.sum_cal();
    *stat = catt_max(ct, *useMax3, static_cast<Alternative>(*alternative));
  }
  
  void c_catt_test(int *table, double *theta, int *alternative, double *stat, double *p){
    ct2x3 ct;
    ct(0,0)=table[0];
    ct(0,1)=table[1];
    ct(0,2)=table[2];
    ct(1,0)=table[3];
    ct(1,1)=table[4];
    ct(1,2)=table[5];
    ct.sum_cal();
    *stat = catt(ct, *theta);
    *p = catt_p(ct, *theta, static_cast<Alternative>(*alternative));
  }


  void c_catt(int *table, double *theta, double *res){
    ct2x3 ct;
    ct(0,0)=table[0];
    ct(0,1)=table[1];
    ct(0,2)=table[2];
    ct(1,0)=table[3];
    ct(1,1)=table[4];
    ct(1,2)=table[5];
    ct.sum_cal();
    *res = catt(ct, *theta);
  }
}
