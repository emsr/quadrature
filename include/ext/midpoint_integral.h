/* quadrature/midpoint_integral.h
 *
 * Copyright (C) 2017-2019 Free Software Foundation, Inc.
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
  template<typename _Tp, typename _FuncTp>
    class composite_midpoint_integral
    {
    public:

      using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
      using _AreaTp = decltype(_RetTp{} * _Tp{});
      using _AbsAreaTp = decltype(std::abs(_AreaTp{}));

      composite_midpoint_integral(_FuncTp __fun, _Tp __a, _Tp __b,
				   std::size_t __num_segs)
      : _M_fun(__fun), _M_lower_lim(__a), _M_upper_lim(__b),
	_M_num_segs(__num_segs), _M_result()
      { }

      _AreaTp operator()();

      template<typename _FuncTp2>
	fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp2, _Tp>>
	integrate(_FuncTp2 __fun, _Tp __a, _Tp __b)
	{
	  composite_midpoint_integral<_FuncTp2, _Tp>
	    __trapi(__fun, __a, __b, this->_M_num_segments);
	  return {__trapi()};
	}

    private:

      _FuncTp _M_fun;
      _Tp _M_lower_lim;
      _Tp _M_upper_lim;
      std::size_t _M_num_segs;
      _AreaTp _M_result;
      _AreaTp _M_asymp_error;
    };

  /**
   * A globally adaptive recursive midpoint integrator.
   */
  template<typename _Tp, typename _FuncTp>
    class midpoint_integral
    {
    public:

      using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
      using _AreaTp = decltype(_RetTp{} * _Tp{});
      using _AbsAreaTp = decltype(std::abs(_AreaTp{}));

      midpoint_integral(_FuncTp __fun, _Tp __a, _Tp __b,
			_Tp __abs_tol, _Tp __rel_tol)
      : _M_fun(__fun), _M_lower_lim(__a), _M_upper_lim(__b),
	_M_abs_tol(std::abs(__abs_tol)), _M_rel_tol(std::abs(__rel_tol)),
	_M_result(), _M_abs_error()
      { }

      _AreaTp operator()();

      _AbsAreaTp abs_error() const
      { return this->_M_abs_error; }

      template<typename _FuncTp2>
	adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp2, _Tp>>
	integrate(_FuncTp2 __fun, _Tp __a, _Tp __b)
	{
	  midpoint_integral<_FuncTp2, _Tp>
	    __mpi(__fun, __a, __b,
		  this->_M_abs_tol, this->_M_rel_tol);
	  return {__mpi(), __mpi.abs_error() };
	}

      template<typename _FuncTp2>
	adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp2, _Tp>>
	operator()(_FuncTp2 __fun, _Tp __a, _Tp __b)
	{ return this->integrate(__fun, __a, __b); }

    private:

      static constexpr auto _S_max_iter = std::numeric_limits<_Tp>::digits / 2;
      static constexpr auto _S_min_delta
			   = std::sqrt(std::numeric_limits<_Tp>::epsilon());

      _AreaTp _M_step();

      _FuncTp _M_fun;
      _Tp _M_lower_lim;
      _Tp _M_upper_lim;
      _AbsAreaTp _M_abs_tol;
      _AbsAreaTp _M_rel_tol;
      _AreaTp _M_result;
      _AbsAreaTp _M_abs_error;
      std::size_t _M_iter = 0;
      std::size_t _M_pow3 = 0;
    };

  template<typename _Tp, typename _FuncTp>
    inline adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_midpoint(_FuncTp __func, _Tp __a, _Tp __b,
			_Tp __max_abs_err, _Tp __max_rel_err, int __max_iter)
    {
      midpoint_integral<_Tp, _FuncTp>
	__mpi(__func, __a, __b, __max_abs_err, __max_rel_err, __max_iter);
      return {__mpi(), __mpi.abs_error()};
    }

} // namespace __gnu_cxx

#endif // MIDPOINT_INTEGRAL_H
