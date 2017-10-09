/* quadrature/midpoint_integral.h
 *
 * Copyright (C) 2017 Free Software Foundation, Inc.
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

template<typename _Func, typename _Tp>
  class midpoint_integral
  {
  public:

    using _RetTp = std::invoke_result_t<_Func, _Tp>;
    using _AreaTp = decltype(_RetTp{} * _Tp{});

    midpoint_integral(_Func __fun, _Tp __a, _Tp __b, _Tp __tol)
    : _M_fun(__fun), _M_lower_lim(__a), _M_upper_lim(__b),
      _M_rel_tol(std::abs(__tol)), _M_result(), _M_abs_error()
    { }

    _AreaTp operator()();

    _AreaTp abs_error() const
    { return this->_M_abs_error; }

  private:

    static constexpr auto _S_max_iter = std::numeric_limits<_Tp>::digits / 2;
    static constexpr auto _S_min_delta
			 = std::sqrt(std::numeric_limits<_Tp>::epsilon());

    _AreaTp _M_step();

    _Func _M_fun;
    _Tp _M_lower_lim;
    _Tp _M_upper_lim;
    _AreaTp _M_rel_tol;
    _AreaTp _M_result;
    _AreaTp _M_abs_error;
    std::size_t _M_iter = 0;
    std::size_t _M_pow3 = 0;
  };

} // namespace __gnu_cxx

#include "midpoint_integral.tcc"

#endif // MIDPOINT_INTEGRAL_H
