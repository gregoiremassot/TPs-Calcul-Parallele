#ifndef _CONJUGATE_GRADIENT_H
#define _CONJUGATE_GRADIENT_H

#include "linear_system.h"

class conjugate_gradient 
{
public:
  conjugate_gradient();
  
  virtual ~conjugate_gradient() { ; }
  
  void set_linear_system(linear_system* lLS) { SL = lLS; }
  void set_precision(double lprecision) { precision = lprecision; }
  void set_nbIterMax(int lnbItermax) { nbItermax = lnbItermax; }
  
  double inner_product(const double* vec1, const double* vec2) const;
  
  int solve();
  int solve_PCD();
    
  
protected:
  linear_system* SL;
  double precision;
  int nbItermax;
};

#endif