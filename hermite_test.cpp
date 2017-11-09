// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2011-2017 Free Software Foundation, Inc.
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
  normalized_hermite(int n1, int n2, _Tp x)
  {
    const auto _S_pi = __gnu_cxx::__const_pi(x);
    auto lnorm = _Tp{0.5} * (log(_S_pi) + _Tp(n1 + n2) * log(_Tp{2})
			  + __gnu_cxx::lfactorial<_Tp>(n1)
			  + __gnu_cxx::lfactorial<_Tp>(n2));
    return std::hermite(n2, x)
	 * std::exp(-x * x - lnorm)
	 * std::hermite(n1, x);
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_hermite()
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto infty = std::numeric_limits<_Tp>::infinity();
    const auto rel_precision = _Tp{1000} * eps;
    const auto abs_precision = _Tp{10} * rel_precision;

    int n1 = 0;
    for (; n1 <= 128; ++n1)
      {
	for (int n2 = 0; n2 <= n1; ++n2)
	  {
	    auto func = [n1, n2](_Tp x)
			-> _Tp
			{ return normalized_hermite(n1, n2, x); };
	    using func_t = decltype(func);

	    auto [result, error]
		//= integrate_singular(func, -infty, infty, rel_precision, _Tp{0});
		//= integrate(func, -infty, infty, rel_precision, _Tp{0});
		//= integrate_infinite(func, rel_precision, _Tp{0});
		//= integrate_singular_infinite(func, rel_precision, _Tp{0});
		//= integrate_oscillatory(inf_transform<func_t, _Tp>(func),
                //                        _Tp{0}, _Tp{1}, rel_precision, _Tp{0});
		= integrate_clenshaw_curtis(inf_transform<func_t, _Tp>(func),
                                            _Tp{0}, _Tp{1}, rel_precision, _Tp{0});

	    if (std::abs(delta<_Tp>(n1, n2) - result) > abs_precision)
	      {
		std::stringstream ss;
		ss.precision(std::numeric_limits<_Tp>::digits10);
		ss << std::showpoint << std::scientific;
		ss << "Integration failed at n1=" << n1 << ", n2=" << n2
		   << ", returning result " << result
		   << ", with error " << error
		   << " instead of the expected " << delta<_Tp>(n1, n2)
		   << " with absolute precision " << abs_precision << '\n';
		throw std::logic_error(ss.str());
	      }
	  }
	std::cout << "Integration successful for hermite polynomials up to n = " << n1
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
	    auto func = [n1 = itop, n2](_Tp x)
			-> _Tp
			{ return normalized_hermite(n1, n2, x); };

	    auto [result, error]
		= integrate_singular(func, -infty, infty, rel_precision, _Tp{0});

	    if (std::abs(delta<_Tp>(itop, n2) - result) > abs_precision)
	      {
		itop = (ibot + itop) / 2;
		goto RESTART;
	      }
	  }
	std::cout << "Integration successful for hermite polynomials up to n = " << itop
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
      test_hermite<float>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for double\n";
  try
    {
      test_hermite<double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for long double\n";
  try
    {
      test_hermite<long double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}
