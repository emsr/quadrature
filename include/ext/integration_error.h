// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2018-2020 Free Software Foundation, Inc.
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
// Used to interpret the error estimate of the integration routines
// Based on gsl/integration/err.c

#ifndef INTEGRATION_ERROR_H
#define INTEGRATION_ERROR_H 1

#include <cmath>
#include <limits>
#include <sstream>
#include <string_view>

namespace __gnu_cxx
{

  enum
  {
    NO_ERROR,
    MAX_ITER_ERROR,
    ROUNDOFF_ERROR,
    SINGULAR_ERROR,
    EXTRAP_ROUNDOFF_ERROR,
    DIVERGENCE_ERROR,
    MAX_SUBDIV_ERROR,
    TOLERANCE_ERROR,
    UNKNOWN_ERROR
  };

  template<typename AreaTp, typename AbsAreaTp>
    class integration_error : public std::runtime_error
    {
      AreaTp _M_result;
      AbsAreaTp _M_abserr;
      int _M_errcode;

    public:

      integration_error(const char* __what, int __errcode,
			AreaTp __result, AbsAreaTp __abserr)
      : std::runtime_error(__what),
	_M_result(__result),
	_M_abserr(__abserr),
	_M_errcode(__errcode)
      { }

      int
      error_code() const
      { return _M_errcode; }

      AreaTp
      result() const
      { return this->_M_result; }

      AbsAreaTp
      abserr() const
      { return this->_M_abserr; }
    };

  /**
   * Throws appropriate error if errcode nonzero
   */
  template<typename AreaTp, typename AbsAreaTp>
    void
    check_error(std::string_view __func, int __errcode,
		AreaTp __result, AbsAreaTp __abserr)
    {
      std::ostringstream msg;
      msg << __func << ": ";

      if (__errcode > 2)
	--__errcode;
      switch(__errcode)
	{
	case NO_ERROR:
	  return;
	case MAX_ITER_ERROR:
	  msg << "Number of iterations was insufficient";
	  break;
	case ROUNDOFF_ERROR:
	  msg << "Cannot reach tolerance because of roundoff error";
	  break;
	case SINGULAR_ERROR:
	  msg << "Bad integrand behavior found in the integration interval";
	  break;
	case EXTRAP_ROUNDOFF_ERROR:
	  msg << "Roundoff error detected in the extrapolation";
	  break;
	case DIVERGENCE_ERROR:
	  msg << "Integral is divergent, or slowly convergent";
	  break;
	case MAX_SUBDIV_ERROR:
	  msg << "Maximum number of subdivisions reached";
	  break;
        case TOLERANCE_ERROR:
	  msg << "Cannot reach tolerance with maximum order rule";
	  break;
	default:
	  msg << "Could not integrate function";
	}

      throw integration_error(msg.str().c_str(), __errcode, __result, __abserr);
    }

  /**
   * 
   */
  template<typename AreaTp, typename AbsAreaTp>
    AbsAreaTp
    rescale_error(AreaTp __err,
		  const AbsAreaTp __result_abs, const AbsAreaTp __result_asc)
    {
      const auto _S_eps = std::numeric_limits<AbsAreaTp>::epsilon();
      const auto _S_min = std::numeric_limits<AbsAreaTp>::min();

      AbsAreaTp __abserr = std::abs(__err);
      if (__result_asc != AbsAreaTp{} && __abserr != AbsAreaTp{})
	{
	  auto __scale = std::pow((AbsAreaTp{200} * __abserr / __result_asc), AbsAreaTp{1.5});

	  if (__scale < AbsAreaTp{1})
	    __abserr = __result_asc * __scale;
	  else
	    __abserr = __result_asc;
	}
      if (__result_abs > _S_min / (AbsAreaTp{50} * _S_eps))
	{
	  auto __min_err = AbsAreaTp{50} * _S_eps * __result_abs;

	  if (__min_err > __abserr)
	    __abserr = __min_err;
	}

      return __abserr;
    }

} // namespace __gnu_cxx

#endif // INTEGRATION_ERROR_H
