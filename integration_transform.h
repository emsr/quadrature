// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2016-2018 Free Software Foundation, Inc.
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

namespace __gnu_cxx
{
  //  Is this useful?
  template<typename _FuncTp, typename _Tp>
    struct mapper
    {
      _FuncTp _M_func;

      mapper(_FuncTp __func)
      : _M_func(__func)
      { }
    };

  /**
   * Map a function defined on @f$ x \in (-\infty, +\infty) @f$
   * to one defined on @f$ t \in (0, 1] @f$.
   * Note: g(t) actually returns [f(-x) + f(+x)]dx/dt.
   */
  template<typename _FuncTp, typename _Tp>
    struct map_minf_pinf
    {
      _FuncTp _M_func;

      map_minf_pinf(_FuncTp __func)
      : _M_func(__func)
      { }

      _Tp
      operator()(_Tp __t) const
      {
	auto __x = (_Tp{1} - __t) / __t;
	auto __y = _M_func(__x) + _M_func(-__x);
	return (__y / __t) / __t;
      }
    };

  /**
   * Map a function defined on @f$ x \in (-\infty, b] @f$
   * to one defined on @f$ t \in (0, 1] @f$.
   */
  template<typename _FuncTp, typename _Tp>
    struct map_minf_b
    {
      _FuncTp _M_func;
      const _Tp _M_b;

      map_minf_b(_FuncTp __func, _Tp __b)
      : _M_func(__func),
	_M_b(__b)
      { }

      _Tp
      operator()(_Tp __t) const
      {
	auto __x = _M_b - (_Tp{1} - __t) / __t;
	auto __y = _M_func(__x);
	return (__y / __t) / __t;
      }
    };

  /**
   * Map a function defined on @f$ x \in [a, +\infty) @f$
   * to one defined on @f$ t \in (0, 1] @f$.
   */
  template<typename _FuncTp, typename _Tp>
    struct map_a_pinf
    {
      _FuncTp _M_func;
      const _Tp _M_a;

      map_a_pinf(_FuncTp __func, _Tp __a)
      : _M_func(__func),
	_M_a(__a)
      { }

      _Tp
      operator()(_Tp __t) const
      {
	_Tp __x = _M_a + (_Tp{1} - __t) / __t;
	_Tp __y = _M_func(__x);
	return (__y / __t) / __t;
      }
    };

} // namespace __gnu_cxx

#endif // INTEGRATION_TRANSFORM_H
