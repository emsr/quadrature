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
// Ported from GSL by Edward Smith-Rowland
// Originally written by Brian Gaugh
//
// Implements Gauss-Kronrod integration
// Based on gsl/integration/qk.c

#ifndef GAUSS_KRONROD_INTERGAL_H
#define GAUSS_KRONROD_INTERGAL_H 1

#include <type_traits>
#include <vector>

namespace emsr
{

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
  template<typename Tp, typename RetTp>
    struct gauss_kronrod_integral_t
    {
      using AreaTp = decltype(RetTp{} * Tp{});
      using AbsAreaTp = decltype(std::abs(AreaTp{}));

      /// Result of the integral.
      AreaTp result = AreaTp{};
      /// Estimated error as difference between Gauss and Kronrod integrals.
      AbsAreaTp abserr = AbsAreaTp{};
      /// Integral of absolute value of function.
      AbsAreaTp resabs = AbsAreaTp{};
      /// Integral of absolute value of difference between function
      /// and weighted mean function value.
      AbsAreaTp resasc = AbsAreaTp{};
    };

  template<typename Tp>
    class gauss_kronrod_integral
    {
    public:

      explicit gauss_kronrod_integral(unsigned gk_rule);

      template<typename FuncTp>
	auto
	integrate(FuncTp func, Tp lower, Tp upper) const
	-> gauss_kronrod_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

      template<typename FuncTp>
	auto
	operator()(FuncTp func, Tp lower, Tp upper) const
	-> gauss_kronrod_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
	{ return this->integrate(func, lower, upper); }


      template<typename FuncTp, typename KronrodIter, typename GaussIter>
	static auto
	s_integrate(const KronrodIter& xgk,
		     const GaussIter& wg,
		     const KronrodIter& wgk,
		     FuncTp func, Tp lower, Tp upper)
	-> gauss_kronrod_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

    private:

      unsigned m_rule = Kronrod_15;

      std::vector<Tp> m_x_kronrod;
      std::vector<Tp> m_w_gauss;
      std::vector<Tp> m_w_kronrod;
    };

  template<typename Tp, typename FuncTp>
    auto
    qk_integrate(FuncTp func, Tp lower, Tp upper,
		 Kronrod_Rule qkintrule)
    -> gauss_kronrod_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

} // namespace emsr

#endif // GAUSS_KRONROD_INTERGAL_H
