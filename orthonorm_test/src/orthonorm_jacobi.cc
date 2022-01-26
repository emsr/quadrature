// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
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

#include <emsr/integration.h>

// Try to manage the four-gamma ratio.
// alpha > -1, beta > -1.
template<typename Tp>
  Tp
  gamma_ratio(int n, Tp alpha, Tp beta)
  {
    constexpr auto eps = std::numeric_limits<Tp>::epsilon();
    if (std::abs(Tp(1) + alpha + beta) < eps)
      return Tp(0);
    else
      {
	auto gaman1 = std::tgamma(Tp(1) + alpha);
	auto gambn1 = std::tgamma(Tp(1) + beta);
	auto gamabn1 = std::tgamma(Tp(1) + alpha + beta);
	auto fact = gaman1 * gambn1 / gamabn1;
	for (int k = 1; k <= n; ++k)
	  fact *= (Tp(k) + alpha) * (Tp(k) + beta)
		/ (Tp(k) + alpha + beta) / Tp(k);
	return fact;
      }
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  normalized_jacobi(int n1, int n2, Tp alpha, Tp beta, Tp x)
  {
    constexpr auto eps = std::numeric_limits<Tp>::epsilon();
    if (std::abs(x - Tp{1}) < eps)
      return Tp{0};
    else if (std::abs(x + Tp{1}) < eps)
      return Tp{0};
    else
      {
	auto gam = gamma_ratio(n1, alpha, beta);
	auto norm = std::pow(Tp{2}, Tp{1} + alpha + beta)
		  * gam / (Tp(2 * n1 + 1) + alpha + beta);
	return std::pow(Tp{1} - x, alpha) * std::pow(Tp{1} + x, beta)
	     * emsr::jacobi(n1, alpha, beta, x)
	     * emsr::jacobi(n2, alpha, beta, x) / norm;
      }
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  void
  test_jacobi(Tp alpha, Tp beta)
  {
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto integ_prec = Tp{1000} * eps;
    const auto cmp_prec = Tp{10} * integ_prec;

    const bool singular = (alpha < Tp{0} || beta < Tp{0});

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [n1, n2, alpha, beta](Tp x)
			-> Tp
		     { return normalized_jacobi<Tp>(n1, n2, alpha, beta, x); };

	    auto [result, error]
		= singular
		? emsr::integrate_singular_endpoints(func,
					       Tp{-1}, Tp{1},
					       alpha, beta, 0, 0,
					       integ_prec, Tp{0})
		: emsr::integrate(func, Tp{-1}, Tp{1}, integ_prec, Tp{0});

	    assert(std::abs(delta<Tp>(n1, n2) - result) < cmp_prec);
	  }
      }
  }

int
main()
{
  test_jacobi<float>();

  test_jacobi<double>();

  test_jacobi<long double>();
}
