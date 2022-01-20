//
// Copyright (C) 2019 Free Software Foundation, Inc.
// Copyright (C) 2021-2022 Edward M. Smith-Rowland
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
//
// Ported from GSL by Jason Dick and Ed Smith-Rowland
// Originally written by Brian Gaugh
//
// This file implements an extrapolation table for use in integration schemes
// Based on gsl/integration/qelg.c

#ifndef QUADRATURE_POINT_H
#define QUADRATURE_POINT_H 1

#include <vector>

namespace emsr
{

  /**
   * A structure to store quadrature rules.
   */
  template<typename Tp>
    struct quadrature_point_t
    {
      Tp point;
      Tp weight;

      constexpr quadrature_point_t() = default;

      constexpr quadrature_point_t(Tp pt, Tp wt)
      : point(pt),
	weight(wt)
      { }
    };

  template<typename Tp>
    std::vector<quadrature_point_t<Tp>>
    legendre_zeros(unsigned int l);

  template<typename Tp>
    std::vector<quadrature_point_t<Tp>>
    laguerre_zeros(unsigned int n, Tp alpha1);

  template<typename Tp>
    std::vector<quadrature_point_t<Tp>>
    hermite_zeros(unsigned int n);

  template<typename Tp>
    std::vector<quadrature_point_t<Tp>>
    gegenbauer_zeros(unsigned int n, Tp lambda);

  template<typename Tp>
    std::vector<quadrature_point_t<Tp>>
    jacobi_zeros(unsigned int n, Tp alpha1, Tp beta1);

} // namespace emsr

#include <emsr/legendre_zeros.tcc>
#include <emsr/laguerre_zeros.tcc>
#include <emsr/hermite_zeros.tcc>
#include <emsr/gegenbauer_zeros.tcc>
#include <emsr/jacobi_zeros.tcc>

#endif // QUADRATURE_POINT_H
