/* quadrature/midpoint_integral.h
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

#ifndef MIDPOINT_INTEGRAL_H
#define MIDPOINT_INTEGRAL_H 1

#include <type_traits>
#include <cstddef>
#include <limits>

namespace __gnu_cxx
{

  /**
   * 
   */
  template<typename Tp, typename FuncTp>
    class composite_midpoint_integral
    {
    public:

      using RetTp = std::invoke_result_t<FuncTp, Tp>;
      using AreaTp = decltype(RetTp{} * Tp{});
      using AbsAreaTp = decltype(std::abs(AreaTp{}));

      composite_midpoint_integral(FuncTp fun, Tp lower, Tp upper,
				  std::size_t num_segs)
      : m_fun(fun), m_lower_lim(lower), m_upper_lim(upper),
	m_num_segs(num_segs), m_result()
      { }

      AreaTp operator()();

      template<typename FuncTp2>
	fixed_integral_t<Tp, std::invoke_result_t<FuncTp2, Tp>>
	integrate(FuncTp2 fun, Tp lower, Tp upper)
	{
	  composite_midpoint_integral<FuncTp2, Tp>
	    trapi(fun, lower, upper, this->m_num_segments);
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
   * A globally adaptive recursive midpoint integrator.
   */
  template<typename Tp, typename FuncTp>
    class midpoint_integral
    {
    public:

      using RetTp = std::invoke_result_t<FuncTp, Tp>;
      using AreaTp = decltype(RetTp{} * Tp{});
      using AbsAreaTp = decltype(std::abs(AreaTp{}));

      midpoint_integral(FuncTp fun, Tp lower, Tp upper,
			Tp abs_tol, Tp rel_tol)
      : m_fun(fun), m_lower_lim(lower), m_upper_lim(upper),
	m_abs_tol(std::abs(abs_tol)), m_rel_tol(std::abs(rel_tol)),
	m_result(), m_abs_error()
      { }

      AreaTp operator()();

      AbsAreaTp abs_error() const
      { return this->m_abs_error; }

      template<typename FuncTp2>
	adaptive_integral_t<Tp, std::invoke_result_t<FuncTp2, Tp>>
	integrate(FuncTp2 fun, Tp lower, Tp upper)
	{
	  midpoint_integral<FuncTp2, Tp>
	    mpi(fun, lower, upper,
		  this->m_abs_tol, this->m_rel_tol);
	  return {mpi(), mpi.abs_error() };
	}

      template<typename FuncTp2>
	adaptive_integral_t<Tp, std::invoke_result_t<FuncTp2, Tp>>
	operator()(FuncTp2 fun, Tp lower, Tp upper)
	{ return this->integrate(fun, lower, upper); }

    private:

      static constexpr auto s_max_iter = std::numeric_limits<Tp>::digits / 2;
      static constexpr auto s_min_delta = std::sqrt(std::numeric_limits<Tp>::epsilon());

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
    integrate_midpoint(FuncTp func, Tp lower, Tp upper,
			Tp max_abs_err, Tp max_rel_err, int max_iter)
    {
      using integ_t = adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(lower) || std::isnan(upper)
          || std::isnan(max_abs_err) || std::isnan(max_rel_err))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}, absarea_t{}};
      else
	{
          midpoint_integral<Tp, FuncTp>
	    mpi(func, lower, upper, max_abs_err, max_rel_err, max_iter);
          return {mpi(), mpi.abs_error()};
	}
    }

} // namespace __gnu_cxx

#endif // MIDPOINT_INTEGRAL_H
