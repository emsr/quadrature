/* integration/glfixed_integrate.tcc
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * Copyright (C) 2016-2020 Free Software Foundation, Inc.
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

#ifndef GLFIXED_INTEGRATE_TCC
#define GLFIXED_INTEGRATE_TCC 1

#include <type_traits>

#include <ext/gauss_legendre_table.h>

namespace emsr
{

  template<typename Tp, typename FuncTp>
    decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
    glfixed_integrate(const gauss_legendre_table<Tp>& t,
		      FuncTp func,
		      Tp lower, Tp upper)
    {
      using RetTp = std::invoke_result_t<FuncTp, Tp>;
      using AreaTp = decltype(RetTp{} * Tp{});

      const int n = t.order;

      auto m = (n + 1) >> 1;
      auto A = (upper - lower) / Tp{2};
      auto B = (upper + lower) / Tp{2};

      if (n & 1) // n is odd.
	{
	  auto sum = t.wt(0) * func(B);
	  for (int i = 1; i < m; ++i)
	    {
	      auto Ax = A * t.pt(i);
	      sum += t.wt(i) * (func(B + Ax) + func(B - Ax));
	    }
	  return A * sum;
	}
      else // n is even.
	{
	  auto sum = AreaTp{0};
	  for (int i = 0; i < m; ++i)
	    {
	      auto Ax = A * t.pt(i);
	      sum += t.wt(i) * (func(B + Ax) + func(B - Ax));
	    }
	  return A * sum;
	}
    }

} // namespace emsr

#endif // GLFIXED_INTEGRATE_TCC
