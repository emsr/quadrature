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
//
// Ported from GSL by Edward Smith-Rowland
// Originally written by Brian Gaugh
//
// Implements integration using a recursive locally-adaptive algorithm
// Based on gsl/integration/qag.c

#ifndef QAG_INTEGRATE_TCC
#define QAG_INTEGRATE_TCC 1

#include <type_traits>
#include <utility>
#include <limits>
#include <string>
#include <stdexcept>

#include <ext/integration_workspace.h>

namespace __gnu_cxx
{

  /**
   * Integrates a function from finite a to finite b using an adaptive
   * quadrature rule.  The integration domain is recursively subdivided
   * prioritized by the segment with the greatest absolute error
   * estimate.  Subdivision continues until a prescribed error target is reached
   * or until a maximum number of divisions is performed.
   *
   * Once either the absolute or relative error limit is reached,
   * qag_integrate() returns
   *
   * @tparam FuncTp     A function type that takes a single real scalar
   *                     argument and returns a real scalar.
   * @tparam Tp         A real type for the limits of integration and the step.
   * @tparam Integrator A non-adaptive integrator that is able to return
   *                     an error estimate in addition to the result.
   *
   * @param[in] workspace The workspace that manages adaptive quadrature
   * @param[in] func The single-variable function to be integrated
   * @param[in] lower The lower limit of integration
   * @param[in] upper The upper limit of integration
   * @param[in] max_abs_err The limit on absolute error
   * @param[in] max_rel_err The limit on relative error
   * @param[in] quad The quadrature stepper taking a function object
   *                   and two integration limits
   *
   * @return A tuple with the first value being the integration result,
   *	     and the second value being the estimated error.
   */
  template<typename Tp, typename FuncTp,
	   typename Integrator = gauss_kronrod_integral<Tp>>
    auto
    qag_integrate(integration_workspace<Tp,
		  std::invoke_result_t<FuncTp, Tp>>& workspace,
		  FuncTp func,
		  Tp lower, Tp upper,
		  Tp max_abs_err, Tp max_rel_err,
		  Integrator quad = gauss_kronrod_integral<Tp>(Kronrod_21))
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      const auto max_iter = workspace.capacity();
      // Try to adjust tests for varing precision.
      const auto s_rel_err = std::pow(Tp{10},
				 -std::numeric_limits<Tp>::digits / Tp{10});

      if (!valid_tolerances(max_abs_err, max_rel_err))
	{
	  std::ostringstream msg;
	  msg << "qag_integrate: Tolerance cannot be achieved with given "
		   "absolute (" << max_abs_err << ") and relative ("
		<< max_rel_err << ") error limits.";
	  throw std::runtime_error(msg.str().c_str());
	}

      auto [result0, abserr0, resabs0, resasc0]
	= quad(func, lower, upper);

      auto tolerance = std::max(max_abs_err, max_rel_err * std::abs(result0));

      // Compute roundoff tolerance.
      const auto round_off = Tp{10} * tolerance * resabs0;

      if (abserr0 <= round_off && abserr0 > tolerance)
	throw integration_error("qag_integrate: "
				"Cannot reach tolerance because "
				"of roundoff error on first attempt",
				ROUNDOFF_ERROR, result0, abserr0);
      else if ((abserr0 <= tolerance && abserr0 != resasc0)
		|| abserr0 == Tp{0})
	return {result0, abserr0};
      else if (max_iter == 1)
	throw integration_error("qag_integrate: "
				"A maximum of one iteration was insufficient",
				MAX_ITER_ERROR, result0, abserr0);

      workspace.clear();
      workspace.append(lower, upper, result0, abserr0);

      auto area = result0;
      auto errsum = abserr0;
      int error_type = NO_ERROR;
      std::size_t iteration = 1;

      int roundoff_type1 = 0, roundoff_type2 = 0;
      do
	{
	  // Bisect the subinterval with the largest error estimate
	  const auto& curr = workspace.retrieve();

	  const auto a1 = curr.lower_lim;
	  const auto mid = (curr.lower_lim + curr.upper_lim) / Tp{2};
	  const auto a2 = mid;
	  const auto b2 = curr.upper_lim;

	  auto [area1, error1, resabs1, resasc1]
	    = quad(func, a1, mid);

	  auto [area2, error2, resabs2, resasc2]
	    = quad(func, a2, b2);

	  const auto area12 = area1 + area2;
	  const auto error12 = error1 + error2;
	  const auto delta = area12 - curr.result;

	  area += delta;
	  errsum += error12 - curr.abs_error;

	  tolerance = std::max(max_abs_err,
				 max_rel_err * std::abs(area));

	  if (resasc1 != error1 && resasc2 != error2)
	    {
	      if (std::abs(delta) <= s_rel_err * std::abs(area12)
		  && error12 >= Tp{0.99} * curr.abs_error)
		++roundoff_type1;
	      if (iteration >= 10 && error12 > curr.abs_error)
		++roundoff_type2;
	    }

	  if (errsum > tolerance)
	    {
	      if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
		error_type = ROUNDOFF_ERROR;

	      // Set error flag in the case of bad integrand behaviour at
	      // a point of the integration range.
	      if (workspace.subinterval_too_small(a1, a2, b2))
		error_type = SINGULAR_ERROR;
	    }

	  workspace.split(mid, area1, error1, area2, error2);

	  ++iteration;
	}
      while (iteration < max_iter
	     && !error_type
	     && errsum > tolerance);

      auto result = workspace.total_integral();
      auto abserr = errsum;

      if (errsum <= tolerance)
	return {result, abserr};

      if (error_type == NO_ERROR && iteration >= max_iter)
	error_type = MAX_ITER_ERROR;

      check_error(__func__, error_type, result, abserr);
      throw integration_error("qag_integrate: Unknown error.",
			      UNKNOWN_ERROR, result, abserr);
    }

} // namespace __gnu_cxx

#endif // QAG_INTEGRATE_TCC
