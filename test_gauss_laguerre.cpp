/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I../include -o test_gauss_laguerre test_gauss_laguerre.cpp
*/

#include <cmath>
#include <iostream>
#include <iomanip>

#include "gauss_laguerre_integrate.h"

template<typename _Tp>
  struct my_f
  {
    _Tp a;
    _Tp b;
    _Tp c;
    _Tp d;
    _Tp e;

    _Tp
    operator()(_Tp x) const
    { return (((a * x + b) * x + c) * x + d) * x + e; }

    _Tp
    integ(_Tp x) const
    { return x * ((((a / 5 * x + b / 4) * x + c / 3) * x + d / 2) * x + e / 1); }
  };

template<typename _Tp>
  struct my_fc
  {
    _Tp a;
    _Tp b;
    _Tp c;

    _Tp
    operator()(_Tp x) const
    { return a * std::cos(b * x + c); }

    _Tp
    integ(_Tp x) const
    { return a * std::sin(b * x + c) / b; }
  };

template<typename _Tp>
  void
  test_laguerre_quad(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout.flags(std::ios::showpoint);
    auto width = 8 + std::cout.precision();

    std::cout << '\n';

    my_f<_Tp> params{ 0.5, 0.0, 0.0, -1.0, 1.3 };
    my_fc<_Tp> cparams{ 1.0, 1.0, 0.0 };

    //pfoo2 = __gnu_cxx::gauss_laguerre_prob_integrate<_Tp>(params, 3);
    //pfooc = __gnu_cxx::gauss_laguerre_prob_integrate<_Tp>(cparams, 29);

    auto xfoo2 = __gnu_cxx::gauss_laguerre_integrate<_Tp>(params, 3, 0.0);
    std::cout << "integral = " << std::setw(width) << xfoo2 << '\n';
    std::cout << "delta    = " << std::setw(width) << xfoo2 - 0.0L << '\n';
    auto xfooc = __gnu_cxx::gauss_laguerre_integrate<_Tp>(cparams, 29, 0.0);
    std::cout << "integral = " << std::setw(width) << xfooc << '\n';
    std::cout << "delta    = " << std::setw(width) << xfooc - 0.0L << '\n';
  }

int
main()
{
  std::cout << "\n\ndouble\n";
  test_laguerre_quad<double>();

  std::cout << "\n\nlong double\n";
  test_laguerre_quad<long double>();
}
