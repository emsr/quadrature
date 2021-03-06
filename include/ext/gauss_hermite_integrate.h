// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2011-2020 Free Software Foundation, Inc.
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

#ifndef GAUSS_HERMITE_INTEGRATE_H
#define GAUSS_HERMITE_INTEGRATE_H 1

#include <vector>
#include <ext/quadrature_point.h>

namespace __gnu_cxx
{

template<typename _Tp, typename _Func>
  decltype(std::invoke_result_t<_Func, _Tp>{} * _Tp{})
  gauss_hermite_integrate(_Func __func, unsigned int __n)
  {
    using _RetTp = std::invoke_result_t<_Func, _Tp>;
    using _AreaTp = decltype(_RetTp{} * _Tp{});

    if(__n == 0)
      std::__throw_domain_error("gauss_hermite_integrate: "
    				"Hermite order must be greater than 0");
    else
     {
	const auto __rule = __hermite_zeros<_Tp>(__n);
	auto __sum = _AreaTp{};
	for (const auto& __pt : __rule)
	  {
	    const auto __x = __pt.__point;
	    const auto __w = __pt.__weight;
	    __sum += __w * __func(__x);
	  }
	return __sum;
      }
  }

} // namespace __gnu_cxx

#endif // GAUSS_HERMITE_INTEGRATE_H
