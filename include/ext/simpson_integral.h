
/* quadrature/simpson_integral.h
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

#ifndef SIMPSON_INTEGRAL_H
#define SIMPSON_INTEGRAL_H 1

#include "trapezoid_integral.h"

namespace __gnu_cxx
{

template<typename _Func, typename _Tp>
  class simpson_integral : trapezoid_integral<_Func, _Tp>
  {
  public:

    simpson_integral(_Func __fun, _Tp __a, _Tp __b, _Tp __tol)
    : trapezoid_integral(__fun, __a, __b, __tol)
    { }

    _Tp
    operator()()
    {
      auto __trap_prev = this->_M_step();
      auto __sum_prev = __trap_prev;
      for (std::size_t __j = 1; __j <= _S_max_iter; ++__j)
	{
          const auto __trap = this->_M_step();
	  const auto __sum = (_Tp{4} * __trap - __trap_prev) / _Tp{3};
	  this->_M_abs_error = std::abs(__sum - __sum_prev);
	  if (this->_M_abs_error < this->_M_rel_tol * std::abs(__sum))
	    return __sum;
	  if (__j > 6
		&&  std::abs(__sum) < this->_M_rel_tol
		&& std::abs(__sum_prev) < this->_M_rel_tol  )
	    return __sum;
	  __sum_prev = __sum;
	  __trap_prev = __trap;
	}

      return __sum_prev;
    }

  };

} // namespace __gnu_cxx

#endif // SIMPSON_INTEGRAL_H
