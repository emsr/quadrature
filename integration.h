// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2011-2018 Free Software Foundation, Inc.
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


#ifndef INTEGRATION_H
#define INTEGRATION_H 1

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Return true if the supplied absolute and relative tolerances
   * can be achieved.
   */
  template<typename _Tp>
    constexpr bool
    valid_tolerances(_Tp __max_abs_err, _Tp __max_rel_err)
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      return __max_abs_err <= _Tp{0}
	     && (__max_rel_err < _Tp{50} * _S_eps
		  || __max_rel_err < 0.5e-28);
      // I don't understand the etymology of this last number.
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#include "qk_integrate.tcc"
#include "qag_integrate.tcc"
#include "qags_integrate.tcc"
#include "qng_integrate.tcc"
#include "qagp_integrate.tcc"
#include "qcheb_integrate.tcc"
#include "qawc_integrate.tcc"
#include "qaws_integrate.tcc"
#include "qawo_integrate.tcc"
#include "qawf_integrate.tcc"
#include "glfixed_integrate.tcc"
#include "cquad_integrate.tcc"

#include "gauss_quadrature.h"

#include "integrate.h"

#endif // INTEGRATION_H
