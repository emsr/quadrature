/* quadrature/simpson_integral.tcc
 *
 * Copyright (C) 2017-2018 Free Software Foundation, Inc.
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

#ifndef SIMPSON_INTEGRAL_TCC
#define SIMPSON_INTEGRAL_TCC 1

#include <cmath>

namespace __gnu_cxx
{

  /**
   * Integrate the function by naive subdivision.
   */
  template<typename _Tp, typename _FuncTp>
    typename composite_simpson_integral<_Tp, _FuncTp>::_AreaTp
    composite_simpson_integral<_Tp, _FuncTp>::operator()()
    {
      const auto __delta = (this->_M_upper_lim - this->_M_lower_lim)
			   / this->_M_num_segs / 2;

      auto __x = this->_M_lower_lim;
      auto __sum = this->_M_fun(__x);
      __x += __delta;
      __sum += _Tp{2} * this->_M_fun(__x);
      for (std::size_t __j = 1; __j < this->_M_num_segs; ++__j)
	{
	  __x += __delta;
	  __sum += _Tp{4} * this->_M_fun(__x);
	  __x += __delta;
	  __sum += _Tp{2} * this->_M_fun(__x);
	}
      __sum += this->_M_fun(this->_M_upper_lim);

      this->_M_result = __sum * __delta / _Tp{6};

      return this->_M_result;
    }

  /**
   * Integrate the function by globally adaptive subdivision.
   */
  template<typename _Tp, typename _FuncTp>
    typename simpson_integral<_Tp, _FuncTp>::_AreaTp
    simpson_integral<_Tp, _FuncTp>::operator()()
    {
      auto __simp_prev = this->_M_step();
      auto __sum_prev = __simp_prev;
      for (std::size_t __j = 1; __j <= _S_max_iter; ++__j)
	{
	  const auto __simp = this->_M_step();
	  const auto __sum = (_Tp{4} * __simp - __simp_prev) / _Tp{3};
	  this->_M_abs_error = std::abs(__sum - __sum_prev);
	  if (this->_M_abs_error < this->_M_rel_tol * std::abs(__sum))
	    return __sum;
	  if (__j > 6
		&& std::abs(__sum) < this->_M_rel_tol
		&& std::abs(__sum_prev) < this->_M_rel_tol  )
	    return __sum;
	  __sum_prev = __sum;
	  __simp_prev = __simp;
	}

      return __sum_prev;
    }

  /**
   * 
   */
  template<typename _Tp, typename _FuncTp>
    typename simpson_integral< _Tp, _FuncTp>::_AreaTp
    simpson_integral< _Tp, _FuncTp>::_M_step()
    {
      const auto __a = this->_M_lower_lim;
      const auto __b = this->_M_upper_lim;
      if (this->_M_iter == 0)
	{
	  this->_M_iter = 1;
	  this->_M_pow3 = 1;
	  const auto __m = (__a + __b) / _Tp{2};
	  auto __midp = this->_M_fun(__m);
	  auto __trap = (this->_M_fun(__a) + this->_M_fun(__b)) / _Tp{2};
	  this->_M_result = (__b - __a) * (__trap + _Tp{2} * __midp) / _Tp{3};
	}
      else
	{
	  ++this->_M_iter;
	  const auto __del = (__b - __a) / _Tp(3 * this->_M_pow3);
	  if (std::abs(__del) < _S_min_delta)
	    return this->_M_result;
	  const auto __ddel = _Tp{2} * __del;
	  auto __t = __a + __del;
	  auto __trap = _AreaTp{};
	  auto __m = __a + __del / _Tp{2};
	  auto __midp = _AreaTp{};
	  for (auto __j = 1u; __j <= this->_M_pow3; ++__j)
	    {
	      __trap += this->_M_fun(__t);
	      __t += __del;
	      __trap += this->_M_fun(__t);
	      __t += __ddel;
	      __midp += this->_M_fun(__m);
	      __m += __ddel;
	      __midp += this->_M_fun(__m);
	      __m += __del;
	    }
	  this->_M_result += (__b - __a) * (__trap + _Tp{2} * __midp)
			   / _Tp(3 * this->_M_pow3);
	  this->_M_result /= _Tp{3};
	  this->_M_pow3 *= 3;
	}
      return this->_M_result;
    }

} // namespace __gnu_cxx

#endif // SIMPSON_INTEGRAL_TCC
