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
  normalized_laguerre(int n1, int n2, _Tp x)
  {
    return std::exp(-x)
	 * std::laguerre(n1, x)
	 * std::laguerre(n2, x);
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_assoc_laguerre()
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto integ_prec = _Tp{1000} * eps;
    const auto cmp_prec = _Tp{10} * integ_prec;

    int n1 = 0;
    for (; n1 <= 128; ++n1)
      {
	for (int n2 = 0; n2 <= n1; ++n2)
	  {
	    auto func = [n1, n2](_Tp x)
			-> _Tp
			{ return normalized_laguerre<_Tp>(n1, n2, x); };

	    auto [result, error]
		= integrate_to_infinity(func, _Tp{0}, integ_prec, _Tp{0});

	    if (std::abs(delta<_Tp>(n1, n2) - result) > cmp_prec)
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
	std::cout << "Integration successful for laguerre polynomials up to n = " << n1
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
			{ return normalized_laguerre<_Tp>(n1, n2, x); };

	    auto [result, error]
		= integrate_to_infinity(func, _Tp{0}, integ_prec, _Tp{0});

	    if (std::abs(delta<_Tp>(itop, n2) - result) > cmp_prec)
	      {
		itop = (ibot + itop) / 2;
		goto RESTART;
	      }
	  }
	std::cout << "Integration successful for laguerre polynomials up to n = " << itop
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
      test_assoc_laguerre<float>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for double\n";
  try
    {
      test_assoc_laguerre<double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for long double\n";
  try
    {
      test_assoc_laguerre<long double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}
