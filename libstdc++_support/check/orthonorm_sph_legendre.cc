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
#include <cmath>
#include "integration.h"

// Function which should integrate to 1 for l1 == l2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_sph_legendre(int l1, int m1, int l2, int m2, _Tp theta)
  {
    const auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
    return _Tp{2} * _S_pi * std::sin(theta)
	 * std::sph_legendre(l1, m1, x)
	 * std::sph_legendre(l2, m2, x);
  }

template<typename _Tp>
  _Tp
  delta(int l1, int l2)
  { return l1 == l2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_sph_legendre()
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
    const auto integ_prec = _Tp{1000} * eps;
    const auto cmp_prec = _Tp{10} * integ_prec;

    for (int n1 : {0, 2, 5})
      {
	for (int m1 = 0; m1 <= n1; ++m1)
	  {
	    for (int n2 : {0, 2, 5})
	      {
		for (int m2 = 0; m2 <= n2; ++m2)
		  {
		    auto func
			 = [l1, m1, l2, m2](_Tp theta)
			   -> _Tp
			   { return norm_sph_legendre(l1, m1, l2, m2, theta); };

		    auto [result, error]
			= __gnu_cxx::integrate(func, _Tp{0}, _S_pi,
						integ_prec, _Tp{0});

		    assert(std::abs(delta<_Tp>(n1, n2) - result) < cmp_prec);
		  }
	      }
	  }
      }
  }

int
main()
{
  test_sph_legendre<float>();

  test_sph_legendre<double>();

  test_sph_legendre<long double>();
}
