// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2016-2017 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <string>

#include "integration.h"

using namespace __gnu_cxx;

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_chebyshev_w(int n1, int n2, _Tp x)
  {
    const auto _S_eps = __gnu_cxx::__epsilon(x);
    const auto _S_inf = __gnu_cxx::__infinity(x);
    const auto _S_pi = __gnu_cxx::__const_pi(x);
    if (std::abs(x + _Tp{1}) < _S_eps)
      return (n1 + n2) & 1 ? -_S_inf : _S_inf;
    else
      return __gnu_cxx::chebyshev_w(n1, x)
	   * __gnu_cxx::chebyshev_w(n2, x)
	   * std::sqrt((_Tp{1} - x) / (_Tp{1} + x))
	   / _S_pi;
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_chebyshev_w()
  {
    const _Tp eps = std::numeric_limits<_Tp>::epsilon();

    int n1 = 0;
    for (; n1 <= 128; ++n1)
      {
	for (int n2 = 0; n2 <= n1; ++n2)
	  {
	    auto func = [n1, n2](_Tp x)->_Tp{return normalized_chebyshev_w(n1, n2, x);};
	    const _Tp integ_precision = _Tp{1000} * eps;
	    const _Tp comp_precision = _Tp{100000} * integ_precision;

	    auto [result, error]
//		= integrate(func, _Tp{-1} + 10000 * eps, _Tp{1}, integ_precision, _Tp{0});
//		= integrate_singular(func, _Tp{-1}, _Tp{1}, integ_precision, _Tp{0});
		= integrate_singular_endpoints(func,
				 _Tp{-1}, _Tp{1},
				 _Tp{-0.5}, _Tp{0.5}, 0, 0,
				 integ_precision, _Tp{0});
	    if (std::abs(delta<_Tp>(n1, n2) - result) > comp_precision)
	      {
		std::stringstream ss;
		ss.precision(std::numeric_limits<_Tp>::digits10);
		ss << std::showpoint << std::scientific;
		ss << "Integration failed at n1=" << n1 << ", n2=" << n2
		   << ", returning result " << result
		   << ", with error " << error
		   << " instead of the expected " << delta<_Tp>(n1, n2) << '\n';
		throw std::logic_error(ss.str());
	      }
	  }
	std::cout << "Integration successful for chebyshev_w polynomials up to n = " << n1
		  << '\n' << std::flush;
      }

    int ibot = n1 - 1;
    int itop = 2 * ibot;
    int del = 2;
    while (itop != ibot)
      {
	RESTART:
	for (int n2 = 0; n2 <= itop; n2 += del)
	  {
	    auto func = [n1 = itop, n2](_Tp x)->_Tp{return normalized_chebyshev_w(n1, n2, x);};
	    const _Tp integ_precision = _Tp{1000} * eps;
	    const _Tp comp_precision = _Tp{100000} * integ_precision;

	    auto [result, error]
//		= integrate(func, _Tp{-1} + 10000 * eps, _Tp{1}, integ_precision, _Tp{0});
//		= integrate_singular(func, _Tp{-1}, _Tp{1}, integ_precision, _Tp{0});
		= integrate_singular_endpoints(func,
				 _Tp{-1}, _Tp{1},
				 _Tp{-0.5}, _Tp{0.5}, 0, 0,
				 integ_precision, _Tp{0});

	    if (std::abs(delta<_Tp>(itop, n2) - result) > comp_precision)
	      {
		itop = (ibot + itop) / 2;
		goto RESTART;
	      }
	  }
	std::cout << "Integration successful for chebyshev_w polynomials up to n = " << itop
		  << '\n' << std::flush;
	ibot = itop;
	if (itop > 1000000)
	  {
	    std::cout << "\nGood enough!\n" << std::flush;
	    break;
	  }
	else if (itop <= std::numeric_limits<int>::max() / 2)
	  itop *= 2;
	else
	  break;
        del *= 2;
      }
  }

int
main()
{
  std::cout << "\n\nOrthonormality tests for float\n";
  try
    {
      test_chebyshev_w<float>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for double\n";
  try
    {
      test_chebyshev_w<double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for long double\n";
  try
    {
      test_chebyshev_w<long double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}
