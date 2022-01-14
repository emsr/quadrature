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
// Ported from GSL by Edward Smith-Rowland
// Originally written by Brian Gaugh
//
// Implements Gauss-Kronrod integration
// Based on gsl/integration/qk.c

#ifndef GAUSS_KRONROD_INTERGAL_H
#define GAUSS_KRONROD_INTERGAL_H 1

#include <type_traits>
#include <vector>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  enum Kronrod_Rule
  {
    Kronrod_15 = 15,
    Kronrod_21 = 21,
    Kronrod_31 = 31,
    Kronrod_41 = 41,
    Kronrod_51 = 51,
    Kronrod_61 = 61
  };

  /**
   * The return type for a Gauss-Kronrod rule.
   */
  template<typename _Tp, typename _RetTp>
    struct gauss_kronrod_integral_t
    {
      using _AreaTp = decltype(_RetTp{} * _Tp{});
      using _AbsAreaTp = decltype(std::abs(_AreaTp{}));

      /// Result of the integral.
      _AreaTp __result = _AreaTp{};
      /// Estimated error as difference between Gauss and Kronrod integrals.
      _AbsAreaTp __abserr = _AbsAreaTp{};
      /// Integral of absolute value of function.
      _AbsAreaTp __resabs = _AbsAreaTp{};
      /// Integral of absolute value of difference between function
      /// and weighted mean function value.
      _AbsAreaTp __resasc = _AbsAreaTp{};
    };

  template<typename _Tp>
    class gauss_kronrod_integral
    {
    public:

      explicit gauss_kronrod_integral(unsigned __gk_rule);

      template<typename _FuncTp>
	auto
	integrate(_FuncTp __func, _Tp __lower, _Tp __upper) const
	-> gauss_kronrod_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

      template<typename _FuncTp>
	auto
	operator()(_FuncTp __func, _Tp __lower, _Tp __upper) const
	-> gauss_kronrod_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
	{ return this->integrate(__func, __lower, __upper); }


      template<typename _FuncTp, typename _KronrodIter, typename _GaussIter>
	static auto
	_S_integrate(const _KronrodIter& __xgk,
		     const _GaussIter& __wg,
		     const _KronrodIter& __wgk,
		     _FuncTp __func, _Tp __lower, _Tp __upper)
	-> gauss_kronrod_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

    private:

      unsigned _M_rule = Kronrod_15;

      std::vector<_Tp> _M_x_kronrod;
      std::vector<_Tp> _M_w_gauss;
      std::vector<_Tp> _M_w_kronrod;
    };

  template<typename _Tp, typename _FuncTp>
    auto
    qk_integrate(_FuncTp __func, _Tp __lower, _Tp __upper,
		 Kronrod_Rule __qkintrule)
    -> gauss_kronrod_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // GAUSS_KRONROD_INTERGAL_H
