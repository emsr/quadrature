//
// Copyright (C) 2011-2020 Free Software Foundation, Inc.
// Copyright (C) 2021-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.

#ifndef GAUSS_HERMITE_INTEGRATE_H
#define GAUSS_HERMITE_INTEGRATE_H 1

#include <type_traits>
#include <vector>

#include <ext/quadrature_point.h>

namespace emsr
{

template<typename Tp, typename Func>
  decltype(std::invoke_result_t<Func, Tp>{} * Tp{})
  gauss_hermite_integrate(Func func, unsigned int n)
  {
    using RetTp = std::invoke_result_t<Func, Tp>;
    using AreaTp = decltype(RetTp{} * Tp{});

    if(n == 0)
      throw std::domain_error("gauss_hermite_integrate: "
    				"Hermite order must be greater than 0");
    else
     {
	const auto rule = hermite_zeros<Tp>(n);
	auto sum = AreaTp{};
	for (const auto& pt : rule)
	  {
	    const auto x = pt.point;
	    const auto w = pt.weight;
	    sum += w * func(x);
	  }
	return sum;
      }
  }

} // namespace emsr

#endif // GAUSS_HERMITE_INTEGRATE_H
