/* quadrature/simpson_integral.h
 *
 * Copyright (C) 2017-2020 Free Software Foundation, Inc.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef SIMPSON_INTEGRAL_H
#define SIMPSON_INTEGRAL_H 1

#include <type_traits>
#include <cstddef>
#include <limits>

namespace emsr
{

  /**
   * 
   */
  template<typename Tp, typename FuncTp>
    class composite_simpson_integral
    {
    public:

      using RetTp = std::invoke_result_t<FuncTp, Tp>;
      using AreaTp = decltype(RetTp{} * Tp{});
      using AbsAreaTp = decltype(std::abs(AreaTp{}));

      composite_simpson_integral(FuncTp fun, Tp a, Tp b,
				   std::size_t num_segs)
      : m_fun(fun), m_lower_lim(a), m_upper_lim(b),
	m_num_segs(num_segs), m_result()
      { }

      AreaTp operator()();

      template<typename FuncTp2>
	fixed_integral_t<Tp, std::invoke_result_t<FuncTp2, Tp>>
	integrate(FuncTp2 fun, Tp a, Tp b)
	{
	  composite_simpson_integral<FuncTp2, Tp>
	    trapi(fun, a, b, this->m_num_segments);
	  return {trapi()};
	}

    private:

      FuncTp m_fun;
      Tp m_lower_lim;
      Tp m_upper_lim;
      std::size_t m_num_segs;
      AreaTp m_result;
      AreaTp m_asymp_error;
    };

  /**
   * A globally adaptive recursive Simpson integrator.
   */
  template<typename Tp, typename FuncTp>
    class simpson_integral
    {
    public:

      using RetTp = std::invoke_result_t<FuncTp, Tp>;
      using AreaTp = decltype(RetTp{} * Tp{});
      using AbsAreaTp = decltype(std::abs(AreaTp{}));

      simpson_integral(FuncTp fun, Tp a, Tp b,
		       Tp abs_tol, Tp rel_tol)
      : m_fun(fun), m_lower_lim(a), m_upper_lim(b),
	m_abs_tol(std::abs(abs_tol)), m_rel_tol(std::abs(rel_tol)),
	m_result(), m_abs_error()
      { }

      AreaTp operator()();

      AbsAreaTp abs_error() const
      { return this->m_abs_error; }

      template<typename FuncTp2>
	adaptive_integral_t<Tp, std::invoke_result_t<FuncTp2, Tp>>
	integrate(FuncTp2 fun, Tp a, Tp b)
	{
	  simpson_integral<FuncTp2, Tp>
	    simpi(fun, a, b, this->m_abs_tol, this->m_rel_tol);
	  return {simpi(), simpi.abs_error() };
	}

      template<typename FuncTp2>
	adaptive_integral_t<Tp, std::invoke_result_t<FuncTp2, Tp>>
	operator()(FuncTp2 fun, Tp a, Tp b)
	{ return this->integrate(fun, a, b); }

    private:

      static constexpr auto s_max_iter = std::numeric_limits<Tp>::digits / 2;
      static constexpr auto s_min_delta
			   = std::sqrt(std::numeric_limits<Tp>::epsilon());

      AreaTp m_step();

      FuncTp m_fun;
      Tp m_lower_lim;
      Tp m_upper_lim;
      AbsAreaTp m_abs_tol;
      AbsAreaTp m_rel_tol;
      AreaTp m_result;
      AbsAreaTp m_abs_error;
      std::size_t m_iter = 0;
      std::size_t m_pow3 = 0;

    };

  template<typename Tp, typename FuncTp>
    inline adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_simpson(FuncTp __func, Tp __a, Tp __b,
		      Tp __max_abs_err, Tp __max_rel_err, int __max_iter)
    {
      simpson_integral<Tp, FuncTp>
	__simpi(__func, __a, __b, __max_abs_err, __max_rel_err, __max_iter);
      return {__simpi(), __simpi.abs_error()};
    }

} // namespace emsr

#endif // SIMPSON_INTEGRAL_H
