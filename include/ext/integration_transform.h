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

namespace emsr
{
  //  Is this useful?
  template<typename Tp, typename FuncTp>
    struct mapper
    {
      FuncTp m_func;

      mapper(FuncTp func)
      : m_func(func)
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
  template<typename Tp, typename FuncTp>
    struct map_minf_pinf_symm
    {
      FuncTp m_func;

      map_minf_pinf_symm(FuncTp func)
      : m_func(func)
      { }

      std::invoke_result_t<FuncTp, Tp>
      operator()(Tp t) const
      {
	if (t == Tp{0})
	  return m_func(-std::numeric_limits<Tp>::infinity());
	else if (t == Tp{1})
	  return m_func(+std::numeric_limits<Tp>::infinity());
	else
	  {
	    const auto x = (Tp{1} - t) / t;
	    const auto y = m_func(+x) + m_func(-x);
	    return y / t / t;
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
  template<typename Tp, typename FuncTp>
    struct map_minf_pinf
    {
      FuncTp m_func;

      map_minf_pinf(FuncTp func)
      : m_func(func)
      { }

      std::invoke_result_t<FuncTp, Tp>
      operator()(Tp t) const
      {
	if (t == Tp{0})
	  return m_func(-std::numeric_limits<Tp>::infinity());
	else if (t == Tp{1})
	  return m_func(+std::numeric_limits<Tp>::infinity());
	else
	  {
	    const auto inv_t = Tp{1} / t;
	    const auto inv_1mt = Tp{1} / (Tp{1} - t);
	    const auto x = -inv_t + inv_1mt;
	    const auto y = m_func(x);
	    return y * (inv_t * inv_t + inv_1mt * inv_1mt);
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
  template<typename Tp, typename FuncTp>
    struct map_minf_b
    {
      FuncTp m_func;
      const Tp m_b;

      map_minf_b(FuncTp func, Tp b)
      : m_func(func),
	m_b(b)
      { }

      std::invoke_result_t<FuncTp, Tp>
      operator()(Tp t) const
      {
	if (t == Tp{0})
	  return m_func(-std::numeric_limits<Tp>::infinity());
	else
	  {
	    const auto inv_t = Tp{1} / t;
	    const auto x = m_b - (Tp{1} - t) * inv_t;
	    const auto y = m_func(x);
	    return y * inv_t * inv_t;
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
  template<typename Tp, typename FuncTp>
    struct map_a_pinf
    {
      FuncTp m_func;
      const Tp m_a;

      map_a_pinf(FuncTp func, Tp a)
      : m_func(func),
	m_a(a)
      { }

      std::invoke_result_t<FuncTp, Tp>
      operator()(Tp t) const
      {
	if (t == Tp{1})
	  return m_func(+std::numeric_limits<Tp>::infinity());
	else
	  {
	    const auto inv_1mt = Tp{1} / (Tp{1} - t);
	    const auto x = m_a + t * inv_1mt;
	    const auto y = m_func(x);
	    return y * inv_1mt * inv_1mt;
	  }
      }
    };

} // namespace emsr

#endif // INTEGRATION_TRANSFORM_H
