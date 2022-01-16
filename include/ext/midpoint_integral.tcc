/* quadrature/midpoint_integral.tcc
 *
 * Copyright (C) 2018-2020 Free Software Foundation, Inc.
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

#ifndef MIDPOINT_INTEGRAL_TCC
#define MIDPOINT_INTEGRAL_TCC 1

#include <cmath>

namespace __gnu_cxx
{

  /**
   * Integrate the function by naive subdivision.
   */
  template<typename Tp, typename FuncTp>
    typename composite_midpoint_integral< Tp, FuncTp>::AreaTp
    composite_midpoint_integral< Tp, FuncTp>::operator()()
    {
      const auto delta = (this->m_upper_lim - this->m_lower_lim)
			   / this->m_num_segs;

      auto sum = this->m_fun(this->m_lower_lim + delta / Tp{2}) / Tp{2};
      for (std::size_t j = 1; j < this->m_num_segs; ++j)
	sum += this->m_fun(this->m_lower_lim + Tp(j + 0.5) * delta);

      this->m_result = sum * delta;

      return this->m_result;
    }

  /**
   * Integrate the function by naive subdivision.
   */
  template<typename Tp, typename FuncTp>
    typename midpoint_integral< Tp, FuncTp>::AreaTp
    midpoint_integral< Tp, FuncTp>::operator()()
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
   * 
   */
  template<typename Tp, typename FuncTp>
    typename midpoint_integral< Tp, FuncTp>::AreaTp
    midpoint_integral< Tp, FuncTp>::m_step()
    {
      const auto a = this->m_lower_lim;
      const auto b = this->m_upper_lim;
      if (this->m_iter == 0)
	{
	  this->m_iter = 1;
	  this->m_pow3 = 1;
	  const auto m = (a + b) / Tp{2};
	  this->m_result = (b - a) * this->m_fun(m);
	}
      else
	{
	  ++this->m_iter;
	  const auto del = (b - a) / Tp(3 * this->m_pow3);
	  if (std::abs(del) < s_min_delta)
	    return this->m_result;
	  const auto ddel = Tp{2} * del;
	  auto m = a + del / Tp{2};
	  auto sum = AreaTp{};
	  for (auto j = 1u; j <= this->m_pow3; ++j)
	    {
	      sum += this->m_fun(m);
	      m += ddel;
	      sum += this->m_fun(m);
	      m += del;
	    }
	  this->m_result += (b - a) * sum / Tp(this->m_pow3);
	  this->m_result /= Tp{3};
	  this->m_pow3 *= 3;
	}
      return this->m_result;
    }

} // namespace __gnu_cxx

#endif // MIDPOINT_INTEGRAL_TCC
