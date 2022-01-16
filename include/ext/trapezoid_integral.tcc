/* quadrature/trapezoid_integral.tcc
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

#ifndef TRAPEZOID_INTEGRAL_TCC
#define TRAPEZOID_INTEGRAL_TCC 1

#include <type_traits>
#include <cmath>

namespace __gnu_cxx
{

  /**
   * Integrate the function by naive subdivision.
   */
  template<typename Tp, typename FuncTp>
    typename composite_trapezoid_integral< Tp, FuncTp>::AreaTp
    composite_trapezoid_integral< Tp, FuncTp>::operator()()
    {
      const auto delta = (this->m_upper_lim - this->m_lower_lim)
			   / this->m_num_segs;

      auto sum = this->m_fun(this->m_lower_lim) / Tp{2};
      for (std::size_t j = 1; j < this->m_num_segs; ++j)
	sum += this->m_fun(this->m_lower_lim + j * delta);
      sum += this->m_fun(this->m_upper_lim) / Tp{2};

      this->m_result = sum * delta;

      return this->m_result;
    }

  /**
   * Integrate the function by globally adaptive binary subdivision.
   */
  template<typename Tp, typename FuncTp>
    typename trapezoid_integral< Tp, FuncTp>::AreaTp
    trapezoid_integral< Tp, FuncTp>::operator()()
    {
      auto sum_prev = this->m_step();
      for (std::size_t j = 1; j < s_max_iter; ++j)
	{
	  const auto sum = this->m_step();
	  this->m_abs_error = std::abs(sum - sum_prev);
	  if (this->m_abs_error < this->m_rel_tol * std::abs(sum))
	    return sum;
	  if (j > 6
	      && std::abs(sum) < this->m_rel_tol
	      && std::abs(sum_prev) < this->m_rel_tol )
	    return sum;
	  sum_prev = sum;
	}
      return sum_prev;
    }

  /**
   * Chances are, if the function returns Nan or inf, we stepped on a pole.
   */
  template<typename Tp, typename FuncTp>
    std::invoke_result_t<FuncTp, Tp>
    wrap_func(FuncTp func, Tp x)
    {
      auto y = func(x);
      if (std::isnan(y) || std::isinf(y))
	return Tp{0};
      else
	return y;
    }

  /**
   * 
   */
  template<typename Tp, typename FuncTp>
    typename trapezoid_integral< Tp, FuncTp>::AreaTp
    trapezoid_integral< Tp, FuncTp>::m_step()
    {
      if (this->m_iter == 0)
	{
	  this->m_iter = 1;
	  this->m_result = (this->m_upper_lim - this->m_lower_lim)
			  * (wrap_func(this->m_fun, this->m_lower_lim)
			   + wrap_func(this->m_fun, this->m_upper_lim))
		       / Tp{2};
	  this->m_pow2 = 1;
	}
      else
	{
	  ++this->m_iter;
	  const auto del = (this->m_upper_lim - this->m_lower_lim)
			   / this->m_pow2;
	  if (std::abs(del) < s_min_delta)
	    return this->m_result;
	  auto x = this->m_lower_lim + del / Tp{2};
	  auto sum = AreaTp{};
	  for (std::size_t j = 0; j < this->m_pow2; ++j, x += del)
	    sum += wrap_func(this->m_fun, x);
	  this->m_result = (this->m_result + del * sum) / Tp{2};
	  this->m_pow2 *= 2;
	}
      return this->m_result;
    }

} // namespace __gnu_cxx

#endif // TRAPEZOID_INTEGRAL_TCC
