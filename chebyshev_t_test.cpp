// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2016 Free Software Foundation, Inc.
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
#include <functional>
#include <stdexcept>
#include <sstream>
#include <string>

//#include "simple_integrate.h"
#include "integration.h"

using namespace __gnu_test;

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_chebyshev_t(int n1, int n2, _Tp x)
  {
    const auto _S_eps = __gnu_cxx::__epsilon(x);
    const auto _S_inf = std::numeric_limits<_Tp>::infinity();
    if (std::abs(x - _Tp{1}) < _S_eps)
      return _S_inf;
    else if (std::abs(x + _Tp{1}) < _S_eps)
      return (n1 + n2) & 1 ? -_S_inf : _S_inf;
    else
      return __gnu_cxx::chebyshev_t(n2, x)
	   * __gnu_cxx::chebyshev_t(n1, x)
	   / std::sqrt(_Tp{1} - x * x);
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_chebyshev_t()
  {
    const _Tp eps = std::numeric_limits<_Tp>::epsilon();

    // Neverending loop: runs until integration fails
    for (int n1 = 0; ; ++n1)
      {
	for (int n2 = 0; n2 <= n1; ++n2)
	  {
	    std::function<_Tp(_Tp)>
	      func([n1, n2](_Tp x)->_Tp{return normalized_chebyshev_t(n1, n2, x);});
	    _Tp integ_precision = _Tp{1000} * eps;
	    _Tp comp_precision = _Tp{10} * integ_precision;
	    _Tp integration_result, integration_error;

	    typedef std::pair<_Tp&,_Tp&> ret_type;
	    ret_type{integration_result, integration_error}
        	= integrate(func, _Tp{-1}, _Tp{1}, integ_precision, _Tp{0});
//        	= integrate_smooth(func, _Tp{-1}, _Tp{1}, integ_precision, _Tp{0});

            if (std::abs(delta<_Tp>(n1, n2) - integration_result) > comp_precision)
              {
        	std::stringstream ss;
        	ss.precision(-int(log10(eps)));
        	ss << "Integration failed at n1=" << n1 << ", n2=" << n2
        	   << ", returning result " << integration_result
        	   << " instead of the expected " << delta<_Tp>(n1, n2) << '\n';
        	throw std::logic_error(ss.str());
              }
	  }
	std::cout << "Integration successful for chebyshev_t polynomials up to n = " << n1
             << '\n';
      }
  }

int
main()
{
  try
    {
      test_chebyshev_t<double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  try
    {
      test_chebyshev_t<long double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}
