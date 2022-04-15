#ifndef codespec_HPP
#define codespec_HPP
#include <cmath>
#include <cassert>

extern const double test_E;
extern const double test_forces[9];
extern const double test_charges[3];
extern const double test_calc_esp[2];


void make_config_files(bool pbc=false);


void assert_all(double E, const double *forces, const double *charges, const double *calc_esp)
{
  assert(fabs(E - test_E) < 1e-5);
  for(int i = 0; i < 9; ++i)
  {
    assert(fabs(forces[i] - test_forces[i]) < 1e-5);
  }

  for(int i = 0; i < 3; ++i)
  {
    assert(fabs(charges[i] - test_charges[i]) < 1e-2);
  }

  for(int i = 0; i < 2; ++i)
  {
    assert(fabs(calc_esp[i] - test_calc_esp[i]) < 1e-2);
  }  

}


#endif //codespec_HPP
