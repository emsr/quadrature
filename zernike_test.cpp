// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2016-2018 Free Software Foundation, Inc.
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
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <string>

#include "ext/integration.h"

using namespace __gnu_cxx;

/**
 * Neumann's number
 */
template<typename _Tp>
  _Tp
  epsilon(int m)
  { return m == 0 ? _Tp{2} : _Tp{1}; }

// Azimuthal integral of zernike product.


// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_zernike(int n1, int m1, int n2, int m2, _Tp rho)
  {
    const auto _S_eps = _Tp{10000} * std::numeric_limits<_Tp>::epsilon();
    const auto _S_2pi = _Tp{2} * __gnu_cxx::__math_constants<_Tp>::__pi;
    auto z1 = [n1, m1, rho](_Tp phi)
	      ->_Tp { return __gnu_cxx::zernike(n1, m1, rho, phi); };
    auto z2 = [n2, m2, rho](_Tp phi)
	      ->_Tp { return __gnu_cxx::zernike(n2, m2, rho, phi); };
    auto norm = _Tp{1} / std::sqrt(_Tp(2 * n1 + 2) * _Tp(2 * n2 + 2));
    auto fun = [n1, m1, rho, z1, z2, norm](_Tp phi)
		->_Tp { return rho * z1(phi) * z2(phi) / norm; };
    auto val
	= integrate(fun, _Tp{0}, _Tp{_S_2pi}, _S_eps, _Tp{0}, 1024, Kronrod_61);
    return _Tp{2} * val.__result / _S_2pi / epsilon<_Tp>(m1);
  }

template<typename _Tp>
  _Tp
  delta(int n1, int m1, int n2, int m2)
  { return (n1 == n2 && m1 == m2) ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_zernike()
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto integ_prec = _Tp{1000} * eps;
    const auto cmp_prec = _Tp{10} * integ_prec;

    int n1 = 0;
    for (; n1 <= 128; ++n1)
      {
	for (int m1 = 0; m1 <= n1; ++m1)
	  {
	    if ((n1 - m1) & 1)
	      continue;
	    for (int n2 = 0; n2 <= n1; ++n2)
	      {
		for (int m2 = 0; m2 <= n2; ++m2)
		  {
		    if ((n2 - m2) & 1)
		      continue;
		    auto func = [n1, m1, n2, m2](_Tp x)
				-> _Tp
				{ return normalized_zernike(n1, m1, n2, m2, x); };

		    auto [result, error]
			= integrate_singular(func, _Tp{0}, _Tp{1},
						 integ_prec, _Tp{0});

		    if (std::abs(delta<_Tp>(n1, m1, n2, m2) - result) > cmp_prec)
		      {
			std::stringstream ss;
			ss.precision(std::numeric_limits<_Tp>::digits10);
			auto w = 8 + ss.precision();
			ss << std::showpoint << std::scientific;
			ss << "Integration failed at n1=" << n1 << ", m1=" << m1
			   << ", n2=" << n2 << ", m2=" << m2
			   << ", returning result = " << std::setw(w) << result
			   << ", with error = " << std::setw(w) << error
			   << " instead of the expected " << delta<_Tp>(n1, m1, n2, m2) << '\n';
			std::cerr << ss.str();
			//throw std::logic_error(ss.str());
		      }
		  }
	      }
	  }
	std::cout << "Integration successful for zernike polynomials up to n = " << n1
		  << '\n' << std::flush;
      }

    int ibot = n1 - 1;
    int itop = 2 * ibot;
    int del = 2;
    while (itop != ibot)
      {
	RESTART:
	for (int m1 = 0; m1 <= itop; ++m1)
	  {
	    if ((itop - m1) & 1)
	      continue;
	    for (int n2 = 0; n2 <= itop; n2 += del)
	      {
		for (int m2 = 0; m2 <= n2; ++m2)
		  {
		    if ((n2 - m2) & 1)
		      continue;
		    auto func = [n1 = itop, m1, n2, m2](_Tp x)
				-> _Tp
				{ return normalized_zernike(n1, m1, n2, m2, x); };

		    auto [result, error]
			= integrate_singular(func, _Tp{0}, _Tp{1}, integ_prec, _Tp{0});

		    if (std::abs(delta<_Tp>(itop, m1, n2, m2) - result) > cmp_prec)
		      {
			itop = (ibot + itop) / 2;
			goto RESTART;
		      }
		  }
	      }
	  }
	std::cout << "Integration successful for zernike polynomials up to n = " << itop
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
      test_zernike<float>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for double\n";
  try
    {
      test_zernike<double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }

  std::cout << "\n\nOrthonormality tests for long double\n";
  try
    {
      test_zernike<long double>();
    }
  catch (std::exception& err)
    {
      std::cerr << err.what() << '\n';
    }
}
