//
// Copyright (C) 2011-2020 Free Software Foundation, Inc.
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
// This file implements a Gauss-Legendre table for use in integration schemes
// Based on gsl/integration/qelg.c

#ifndef GAUSS_LEGENDRE_TABLE_TCC
#define GAUSS_LEGENDRE_TABLE_TCC 1

#include <emsr/quadrature_point.h>

namespace emsr
{

  template<typename Tp>
      gauss_legendre_table<Tp>::gauss_legendre_table(std::size_t n)
      : order(n),
	point(nullptr),
	weight(nullptr),
	precomputed(false),
	i_precomp(-1),
	rule()
      {
	using prec_t = decltype(gauss_legendre_precomp[0]);
	auto prec_beg = gauss_legendre_precomp;
	auto prec_end = gauss_legendre_precomp + num_gauss_legendre_precomp;
	auto prec = std::find_if(prec_beg, prec_end,
				   [this, n](prec_t tab)
				   { return tab.order == this->order; });
	if (prec == prec_end)
	  this->rule = legendre_zeros<Tp>(this->order);
	else
	  {
	    this->precomputed = true;
	    this->i_precomp = prec - prec_beg;
	  }
      }

  /**
   * Routine to retrieve the i-th Gauss-Legendre point and weight.
   * Useful when the caller wishes to access the information stored in
   * the high-precision gauss_legendre_table struct.
   * Points are indexed and presented in increasing order to the caller.
   */
  template<typename Tp>
    std::tuple<Tp, Tp>
    gauss_legendre_table<Tp>::get_point(Tp lower, Tp upper, size_t i) const
    {
      const auto hwidth = (upper - lower) / Tp{2};
      const auto midpt = (lower + upper) / Tp{2};

      if (i >= this->order)
	throw std::domain_error("gauss_legendre_table: i must be less than n");

      // See comments above gauss_legendre_table for struct's x, w layout.
      // Simply unpack that layout into a sorted set of points, weights.
      Tp xi, wi;
      if (this->order & 1) // n is odd
	{
	  const auto k = int(i) - int(this->order) / 2;
	  const auto sign = (k < 0 ? -1 : +1);

	  xi = midpt + sign * hwidth * this->pt(sign * k);
	  wi =		  hwidth * this->wt(sign * k);
	}
      else if (/* n is even && */ i < this->order / 2)
	{
	  i = int(this->order) / 2 - 1 - int(i);
	  xi = midpt - hwidth * this->pt(i);
	  wi =	   hwidth * this->wt(i);
	}
      else // n is even and i >= n / 2
	{
	  i  -= this->order / 2;
	  xi = midpt + hwidth * this->pt(i);
	  wi =	   hwidth * this->wt(i);
	}

      return std::make_tuple(xi, wi);
    }

  template<typename Tp>
    Tp
    gauss_legendre_table<Tp>::pt(size_t i) const
    {
      if (this->precomputed)
	if (this->i_precomp == std::size_t(-1))
	  return this->point[i];
	else
	  return Tp(gauss_legendre_precomp[this->i_precomp].point[i]);
      else
	return this->rule[i + this->order / 2].point;
    }

  template<typename Tp>
    Tp
    gauss_legendre_table<Tp>::wt(size_t i) const
    {
      if (this->precomputed)
	if (this->i_precomp == std::size_t(-1))
	  return this->weight[i];
	else
	  return Tp(gauss_legendre_precomp[this->i_precomp].weight[i]);
      else
	return this->rule[i + this->order / 2].weight;
    }

} // namespace emsr

#endif // GAUSS_LEGENDRE_TABLE_TCC
