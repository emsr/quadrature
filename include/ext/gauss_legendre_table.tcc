// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2011-2020 Free Software Foundation, Inc.
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
//
// Ported from GSL by Jason Dick and Ed Smith-Rowland
// Originally written by Brian Gaugh
//
// This file implements a Gauss-Legendre table for use in integration schemes
// Based on gsl/integration/qelg.c

#ifndef GAUSS_LEGENDRE_TABLE_TCC
#define GAUSS_LEGENDRE_TABLE_TCC 1


#include <ext/quadrature_point.h>

namespace __gnu_cxx
{

  template<typename _Tp>
      gauss_legendre_table<_Tp>::gauss_legendre_table(std::size_t __n)
      : order(__n),
	point(nullptr),
	weight(nullptr),
	precomputed(false),
	i_precomp(-1),
	rule()
      {
	using __prec_t = decltype(gauss_legendre_precomp[0]);
	auto __prec_beg = gauss_legendre_precomp;
	auto __prec_end = gauss_legendre_precomp + num_gauss_legendre_precomp;
	auto __prec = std::find_if(__prec_beg, __prec_end,
				   [this, __n](__prec_t __tab)
				   { return __tab.order == this->order; });
	if (__prec == __prec_end)
	  this->rule = __legendre_zeros<_Tp>(this->order);
	else
	  {
	    this->precomputed = true;
	    this->i_precomp = __prec - __prec_beg;
	  }
      }

  /**
   * Routine to retrieve the i-th Gauss-Legendre point and weight.
   * Useful when the caller wishes to access the information stored in
   * the high-precision gauss_legendre_table struct.
   * Points are indexed and presented in increasing order to the caller.
   */
  template<typename _Tp>
    std::tuple<_Tp, _Tp>
    gauss_legendre_table<_Tp>::get_point(_Tp __lower, _Tp __upper, size_t __i) const
    {
      const auto __hwidth = (__upper - __lower) / _Tp{2};
      const auto __midpt = (__lower + __upper) / _Tp{2};

      if (__i >= this->order)
	std::__throw_domain_error("gauss_legendre_table: i must be less than n");

      // See comments above gauss_legendre_table for struct's x, w layout.
      // Simply unpack that layout into a sorted set of points, weights.
      _Tp __xi, __wi;
      if (this->order & 1) // n is odd
	{
	  const auto __k = int(__i) - int(this->order) / 2;
	  const auto sign = (__k < 0 ? -1 : +1);

	  __xi = __midpt + sign * __hwidth * this->pt(sign * __k);
	  __wi =		  __hwidth * this->wt(sign * __k);
	}
      else if (/* n is even && */ __i < this->order / 2)
	{
	  __i = int(this->order) / 2 - 1 - int(__i);
	  __xi = __midpt - __hwidth * this->pt(__i);
	  __wi =	   __hwidth * this->wt(__i);
	}
      else // n is even and __i >= n / 2
	{
	  __i  -= this->order / 2;
	  __xi = __midpt + __hwidth * this->pt(__i);
	  __wi =	   __hwidth * this->wt(__i);
	}

      return std::make_tuple(__xi, __wi);
    }

  template<typename _Tp>
    _Tp
    gauss_legendre_table<_Tp>::pt(size_t __i) const
    {
      if (this->precomputed)
	if (this->i_precomp == std::size_t(-1))
	  return this->point[__i];
	else
	  return _Tp(gauss_legendre_precomp[this->i_precomp].point[__i]);
      else
	return this->rule[__i + this->order / 2].__point;
    }

  template<typename _Tp>
    _Tp
    gauss_legendre_table<_Tp>::wt(size_t __i) const
    {
      if (this->precomputed)
	if (this->i_precomp == std::size_t(-1))
	  return this->weight[__i];
	else
	  return _Tp(gauss_legendre_precomp[this->i_precomp].weight[__i]);
      else
	return this->rule[__i + this->order / 2].__weight;
    }

} // namespace __gnu_cxx

#endif // GAUSS_LEGENDRE_TABLE_TCC
