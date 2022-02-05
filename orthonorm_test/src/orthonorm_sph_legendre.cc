
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
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
#include <numbers>

#include <emsr/integration.h>
#include <emsr/math_constants.h>

// Function which should integrate to 1 for l1 == l2, 0 otherwise.
template<typename Tp>
  Tp
  norm_sph_legendre(int l1, int m1, int l2, int m2, Tp theta)
  {
    constexpr auto pi = emsr::pi_v<Tp>;
    return Tp{2} * pi * std::sin(theta)
	 * std::sph_legendre(l1, m1, theta)
	 * std::sph_legendre(l2, m2, theta);
  }

template<typename Tp>
  Tp
  delta(int l1, int l2)
  { return l1 == l2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  void
  test_sph_legendre()
  {
    constexpr auto eps = std::numeric_limits<Tp>::epsilon();
    constexpr auto pi = emsr::pi_v<Tp>;
    constexpr auto integ_prec = Tp{1000} * eps;
    constexpr auto cmp_prec = Tp{10} * integ_prec;

    for (int l1 : {0, 2, 5})
      {
	for (int m1 = 0; m1 <= l1; ++m1)
	  {
	    for (int l2 : {0, 2, 5})
	      {
		for (int m2 = 0; m2 <= l2; ++m2)
		  {
		    auto func
			 = [l1, m1, l2, m2](Tp theta)
			   -> Tp
			   { return norm_sph_legendre(l1, m1, l2, m2, theta); };

		    auto [result, error]
			= emsr::integrate(func, Tp{0}, pi,
						integ_prec, Tp{0});

		    assert(std::abs(delta<Tp>(l1, l2) - result) < cmp_prec);
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
