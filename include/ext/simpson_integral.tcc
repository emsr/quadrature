//
// Copyright (C) 2017-2020 Free Software Foundation, Inc.
// Copyright (C) 2021-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//

#ifndef SIMPSON_INTEGRAL_TCC
#define SIMPSON_INTEGRAL_TCC 1

#include <cmath>

namespace emsr
{

  /**
   * Integrate the function by naive subdivision.
   */
  template<typename Tp, typename FuncTp>
    typename composite_simpson_integral<Tp, FuncTp>::AreaTp
    composite_simpson_integral<Tp, FuncTp>::operator()()
    {
      const auto delta = (this->m_upper_lim - this->m_lower_lim)
			   / this->m_num_segs / 2;

      auto x = this->m_lower_lim;
      auto sum = this->m_fun(x);
      x += delta;
      sum += Tp{2} * this->m_fun(x);
      for (std::size_t j = 1; j < this->m_num_segs; ++j)
	{
	  x += delta;
	  sum += Tp{4} * this->m_fun(x);
	  x += delta;
	  sum += Tp{2} * this->m_fun(x);
	}
      sum += this->m_fun(this->m_upper_lim);

      this->m_result = sum * delta / Tp{6};

      return this->m_result;
    }

  /**
   * Integrate the function by globally adaptive subdivision.
   */
  template<typename Tp, typename FuncTp>
    typename simpson_integral<Tp, FuncTp>::AreaTp
    simpson_integral<Tp, FuncTp>::operator()()
    {
      auto simp_prev = this->m_step();
      auto sum_prev = simp_prev;
      for (std::size_t j = 1; j <= s_max_iter; ++j)
	{
	  const auto simp = this->m_step();
	  const auto sum = (Tp{4} * simp - simp_prev) / Tp{3};
	  this->m_abs_error = std::abs(sum - sum_prev);
	  if (this->m_abs_error < this->m_rel_tol * std::abs(sum))
	    return sum;
	  if (j > 6
		&& std::abs(sum) < this->m_rel_tol
		&& std::abs(sum_prev) < this->m_rel_tol  )
	    return sum;
	  sum_prev = sum;
	  simp_prev = simp;
	}

      return sum_prev;
    }

  /**
   * 
   */
  template<typename Tp, typename FuncTp>
    typename simpson_integral< Tp, FuncTp>::AreaTp
    simpson_integral< Tp, FuncTp>::m_step()
    {
      const auto a = this->m_lower_lim;
      const auto b = this->m_upper_lim;
      if (this->m_iter == 0)
	{
	  this->m_iter = 1;
	  this->m_pow3 = 1;
	  const auto m = (a + b) / Tp{2};
	  auto midp = this->m_fun(m);
	  auto trap = (this->m_fun(a) + this->m_fun(b)) / Tp{2};
	  this->m_result = (b - a) * (trap + Tp{2} * midp) / Tp{3};
	}
      else
	{
	  ++this->m_iter;
	  const auto del = (b - a) / Tp(3 * this->m_pow3);
	  if (std::abs(del) < s_min_delta)
	    return this->m_result;
	  const auto ddel = Tp{2} * del;
	  auto t = a + del;
	  auto trap = AreaTp{};
	  auto m = a + del / Tp{2};
	  auto midp = AreaTp{};
	  for (auto j = 1u; j <= this->m_pow3; ++j)
	    {
	      trap += this->m_fun(t);
	      t += del;
	      trap += this->m_fun(t);
	      t += ddel;
	      midp += this->m_fun(m);
	      m += ddel;
	      midp += this->m_fun(m);
	      m += del;
	    }
	  this->m_result += (b - a) * (trap + Tp{2} * midp)
			   / Tp(3 * this->m_pow3);
	  this->m_result /= Tp{3};
	  this->m_pow3 *= 3;
	}
      return this->m_result;
    }

} // namespace emsr

#endif // SIMPSON_INTEGRAL_TCC
