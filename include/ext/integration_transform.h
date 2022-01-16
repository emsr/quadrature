// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2016-2020 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
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

#ifndef INTEGRATION_TRANSFORM_H
#define INTEGRATION_TRANSFORM_H 1

#include <limits>

namespace __gnu_cxx
{
  //  Is this useful?
  template<typename _Tp, typename _FuncTp>
    struct mapper
    {
      _FuncTp _M_func;

      mapper(_FuncTp __func)
      : _M_func(__func)
      { }
    };

  /**
   * Map a function defined on @f$ x \in (-\infty, +\infty) @f$
   * to one defined on @f$ t \in (0, 1] @f$:
   * $f[
   *    \int_{-\infty}^{+\infty}dx f(x) = \int_{0}^{1}dt \frac{f(t)}{t^2}
   * $f]
   * $f[
   *    x(t) = \frac{1-t}{t}
   * $f]
   * Note: @f$ g(t) @f$ actually returns @f$ [f(-x) + f(+x)]dx/dt @f$.
   * Note: @f$ x:0->\infty @f$ as @f$ t:1->0 @f$.
   */
  template<typename _Tp, typename _FuncTp>
    struct map_minf_pinf_symm
    {
      _FuncTp _M_func;

      map_minf_pinf_symm(_FuncTp __func)
      : _M_func(__func)
      { }

      std::invoke_result_t<_FuncTp, _Tp>
      operator()(_Tp __t) const
      {
	if (__t == _Tp{0})
	  return _M_func(-std::numeric_limits<_Tp>::infinity());
	else if (__t == _Tp{1})
	  return _M_func(+std::numeric_limits<_Tp>::infinity());
	else
	  {
	    const auto __x = (_Tp{1} - __t) / __t;
	    const auto __y = _M_func(+__x) + _M_func(-__x);
	    return __y / __t / __t;
	  }
      }
    };

  /**
   * Map a function defined on @f$ x \in (-\infty, +\infty) @f$
   * to one defined on @f$ t \in (0, 1) @f$:
   * $f[
   *    \int_{-\infty}^{+\infty}dx f(x) = \int_{0}^{1}dt \frac{f(t)}{t^2}
   * $f]
   * $f[
   *    x(t) = -\frac{1}{t} + \frac{1}{1-t}
   * $f]
   */
  template<typename _Tp, typename _FuncTp>
    struct map_minf_pinf
    {
      _FuncTp _M_func;

      map_minf_pinf(_FuncTp __func)
      : _M_func(__func)
      { }

      std::invoke_result_t<_FuncTp, _Tp>
      operator()(_Tp __t) const
      {
	if (__t == _Tp{0})
	  return _M_func(-std::numeric_limits<_Tp>::infinity());
	else if (__t == _Tp{1})
	  return _M_func(+std::numeric_limits<_Tp>::infinity());
	else
	  {
	    const auto __inv_t = _Tp{1} / __t;
	    const auto __inv_1mt = _Tp{1} / (_Tp{1} - __t);
	    const auto __x = -__inv_t + __inv_1mt;
	    const auto __y = _M_func(__x);
	    return __y * (__inv_t * __inv_t + __inv_1mt * __inv_1mt);
	  }
      }
    };

  /**
   * Map a function defined on @f$ x \in (-\infty, b] @f$
   * to one defined on @f$ t \in (0, 1] @f$:
   * $f[
   *    \int_{-\infty}^{b}dx f(x) = \int_{0}^{1}dt \frac{f(t)}{t^2}
   * $f]
   * $f[
   *    x(t) = b - \frac{1-t}{t}
   * $f]
   */
  template<typename _Tp, typename _FuncTp>
    struct map_minf_b
    {
      _FuncTp _M_func;
      const _Tp _M_b;

      map_minf_b(_FuncTp __func, _Tp __b)
      : _M_func(__func),
	_M_b(__b)
      { }

      std::invoke_result_t<_FuncTp, _Tp>
      operator()(_Tp __t) const
      {
	if (__t == _Tp{0})
	  return _M_func(-std::numeric_limits<_Tp>::infinity());
	else
	  {
	    const auto __inv_t = _Tp{1} / __t;
	    const auto __x = _M_b - (_Tp{1} - __t) * __inv_t;
	    const auto __y = _M_func(__x);
	    return __y * __inv_t * __inv_t;
	  }
      }
    };

  /**
   * Map a function defined on @f$ x \in [a, +\infty) @f$
   * to one defined on @f$ t \in (0, 1] @f$:
   * $f[
   *    \int_{a}^{+\infty}dx f(x) = \int_{0}^{1}dt \frac{f(t)}{(1-t)^2}
   * $f]
   * $f[
   *    x(t) = a + \frac{t}{1-t}
   * $f]
   */
  template<typename _Tp, typename _FuncTp>
    struct map_a_pinf
    {
      _FuncTp _M_func;
      const _Tp _M_a;

      map_a_pinf(_FuncTp __func, _Tp __a)
      : _M_func(__func),
	_M_a(__a)
      { }

      std::invoke_result_t<_FuncTp, _Tp>
      operator()(_Tp __t) const
      {
	if (__t == _Tp{1})
	  return _M_func(+std::numeric_limits<_Tp>::infinity());
	else
	  {
	    const auto __inv_1mt = _Tp{1} / (_Tp{1} - __t);
	    const auto __x = _M_a + __t * __inv_1mt;
	    const auto __y = _M_func(__x);
	    return __y * __inv_1mt * __inv_1mt;
	  }
      }
    };

} // namespace __gnu_cxx

#endif // INTEGRATION_TRANSFORM_H
