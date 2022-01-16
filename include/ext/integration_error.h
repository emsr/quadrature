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

namespace emsr
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
      AreaTp m_result;
      AbsAreaTp m_abserr;
      int m_errcode;

    public:

      integration_error(const char* what, int errcode,
			AreaTp result, AbsAreaTp abserr)
      : std::runtime_error(what),
	m_result(result),
	m_abserr(abserr),
	m_errcode(errcode)
      { }

      int
      error_code() const
      { return m_errcode; }

      AreaTp
      result() const
      { return this->m_result; }

      AbsAreaTp
      abserr() const
      { return this->m_abserr; }
    };

  /**
   * Throws appropriate error if errcode nonzero
   */
  template<typename AreaTp, typename AbsAreaTp>
    void
    check_error(std::string_view func, int errcode,
		AreaTp result, AbsAreaTp abserr)
    {
      std::ostringstream msg;
      msg << func << ": ";

      if (errcode > 2)
	--errcode;
      switch(errcode)
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

      throw integration_error(msg.str().c_str(), errcode, result, abserr);
    }

  /**
   * 
   */
  template<typename AreaTp, typename AbsAreaTp>
    AbsAreaTp
    rescale_error(AreaTp err,
		  const AbsAreaTp result_abs, const AbsAreaTp result_asc)
    {
      const auto s_eps = std::numeric_limits<AbsAreaTp>::epsilon();
      const auto s_min = std::numeric_limits<AbsAreaTp>::min();

      AbsAreaTp abserr = std::abs(err);
      if (result_asc != AbsAreaTp{} && abserr != AbsAreaTp{})
	{
	  auto scale = std::pow((AbsAreaTp{200} * abserr / result_asc), AbsAreaTp{1.5});

	  if (scale < AbsAreaTp{1})
	    abserr = result_asc * scale;
	  else
	    abserr = result_asc;
	}
      if (result_abs > s_min / (AbsAreaTp{50} * s_eps))
	{
	  auto min_err = AbsAreaTp{50} * s_eps * result_abs;

	  if (min_err > abserr)
	    abserr = min_err;
	}

      return abserr;
    }

} // namespace emsr

#endif // INTEGRATION_ERROR_H
