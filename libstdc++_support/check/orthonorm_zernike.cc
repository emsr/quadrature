// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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

#include <limits>
#include <cassert>
#include <ext/cmath>
#include "integration.h"

/**
 * Neumann's number
 */
template<typename _Tp>
  _Tp
  epsilon(int m)
  { return m == 0 ? _Tp{2} : _Tp{1}; }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
// This actually does the angular integral.
template<typename _Tp>
  _Tp
  norm_zernike(int n1, int m1, int n2, int m2, _Tp rho)
  {
    const auto _S_eps = _Tp{10000} * std::numeric_limits<_Tp>::epsilon();
    const auto _S_2pi = _Tp{2} * __gnu_cxx::__math_constants<_Tp>::__pi(rho);
    auto z1 = [n1, m1, rho](_Tp phi)
	      ->_Tp { return __gnu_cxx::zernike(n1, m1, rho, phi); };
    auto z2 = [n2, m2, rho](_Tp phi)
	      ->_Tp { return __gnu_cxx::zernike(n2, m2, rho, phi); };
    auto norm = _Tp{1} / std::sqrt(_Tp(2 * n1 + 2) * _Tp(2 * n2 + 2));
    auto fun = [n1, m1, rho, z1, z2, norm](_Tp phi)
	       ->_Tp { return rho * z1(phi) * z2(phi) / norm; };
    _Tp val;
    std::tie(val, std::ignore) = __gnu_cxx::integrate(fun, _Tp{0}, _Tp{_S_2pi},
					   _S_eps, _Tp{0}, 1024, QK_61);
    return _Tp{2} * val / _S_2pi / epsilon<_Tp>(m1);
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

    for (int n1 : {0, 2, 5})
      {
	for (int m1 = 0; m1 <= n1; ++m1)
	  {
	    if ((n1 - m1) & 1)
	      continue;
	    for (int n2 : {0, 2, 5})
	      {
		for (int m2 = 0; m2 <= n2; ++m2)
		  {
		    if ((n2 - m2) & 1)
		      continue;
		    auto func = [n1, m1, n2, m2](_Tp x)
				-> _Tp
				{ return norm_zernike(n1, m1, n2, m2, x); };

		    auto [result, error]
			= __gnu_cxx::integrate_singular(func, _Tp{0}, _Tp{1},
							 integ_prec, _Tp{0});

		    assert(std::abs(delta<_Tp>(n1, m1, n2, m2) - result)
				 < cmp_prec);
	      }
	  }
      }
  }

int
main()
{
  test_zernike<float>();

  test_zernike<double>();

  test_zernike<long double>();
}
