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
#include <functional>
#include <stdexcept>
#include <sstream>
#include <string>

#include "integration.h"

using namespace __gnu_cxx;

// Function which should integrate to 1 for l1 == l2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_legendre(int l1, int l2, _Tp x)
  {
    return (_Tp(l1 + l2 + 1) / _Tp{2})
	 * std::legendre(l1,x)
	 * std::legendre(l2,x);
  }

template<typename _Tp>
  _Tp
  delta(int l1, int l2)
  { return l1 == l2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_legendre()
  {
    const _Tp eps = std::numeric_limits<_Tp>::epsilon();

    int l1 = 0;
    for (; l1 <= 128; ++l1)
      {
	for (int l2 = 0; l2 <= l1; ++l2)
	  {
	    std::function<_Tp(_Tp)> func(std::bind(&normalized_legendre<_Tp>, l1, l2,
					 std::placeholders::_1));
	    const _Tp integ_precision = _Tp{1000} * eps;
	    const _Tp comp_precision = _Tp{10} * integ_precision;

	    auto [result, error]
		= integrate_singular(func, _Tp{-1}, _Tp{1}, integ_precision, _Tp{0});

	    if (std::abs(delta<_Tp>(l1, l2) - result) > comp_precision)
	      {
		std::stringstream ss;
		ss.precision(std::numeric_limits<_Tp>::digits10);
		ss << std::showpoint << std::scientific;
		ss << "Integration failed at l1=" << l1 << ", l2=" << l2
		   << ", returning result " << result
		   << ", with error " << error
		   << " instead of the expected " << delta<_Tp>(l1, l2) << '\n';
		throw std::logic_error(ss.str());
	      }
	  }
	std::cout << "Integration successful for legendre polynomials up to l = " << l1
		  << '\n' << std::flush;
      }

    int ibot = l1 - 1;
    int itop = 2 * ibot;
    int del = 2;
    while (itop != ibot)
      {
	RESTART:
	for (int l2 = 0; l2 <= itop; l2 += del)
	  {
	    std::function<_Tp(_Tp)> func(std::bind(&normalized_legendre<_Tp>, itop, l2,
					 std::placeholders::_1));
	    const _Tp integ_precision = _Tp{1000} * eps;
	    const _Tp comp_precision = _Tp{10} * integ_precision;

	    auto [result, error]
		= integrate_singular(func, _Tp{-1}, _Tp{1}, integ_precision, _Tp{0});

	    if (std::abs(delta<_Tp>(itop, l2) - result) > comp_precision)
	      {
		itop = (ibot + itop) / 2;
		goto RESTART;
	      }
	  }
	std::cout << "Integration successful for legendre polynomials up to l = " << itop
		  << '\n' << std::flush;
	ibot = itop;
	if (itop <= std::numeric_limits<int>::max() / 2)
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
      test_legendre<float>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for double\n";
  try
    {
      test_legendre<double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for long double\n";
  try
    {
      test_legendre<long double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}
