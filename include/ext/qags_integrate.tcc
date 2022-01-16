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
// Implements integration using a recursive Gauss-Kronrod integration
// with singularities
// Based on gsl/integration/qags.c

#ifndef QAGS_INTEGRATE_TCC
#define QAGS_INTEGRATE_TCC 1

#include <type_traits>
#include <utility>
#include <limits>
#include <tuple>
#include <stdexcept>

#include <ext/integration_error.h>
#include <ext/integration_transform.h>
#include <ext/integration_workspace.h>
#include <ext/extrapolation_table.h>

namespace __gnu_cxx
{

  /**
   * Adaptively integrate a potentially singular function from a to b
   * using a recursive Gauss-Kronrod algorithm.
   *
   * @tparam FuncTp     A function type that takes a single real scalar
   *                     argument and returns a real scalar.
   * @tparam Tp         A real type for the limits of integration and the step.
   * @tparam Integrator A non-adaptive integrator that is able to return
   *                     an error estimate in addition to the result.
   *
   * @param[in] workspace The workspace that manages adaptive quadrature
   * @param[in] func The single-variable function to be integrated
   * @param[in] pts The sorted array of points including the integration
   *                  limits and intermediate discontinuities/singularities
   * @param[in] max_abs_err The limit on absolute error
   * @param[in] max_rel_err The limit on relative error
   * @param[in] quad The quadrature stepper taking a function object
   *                   and two integration limits
   */
  template<typename Tp, typename FuncTp,
	   typename Integrator = gauss_kronrod_integral<Tp>>
    auto
    qags_integrate(integration_workspace<Tp,
			std::invoke_result_t<FuncTp, Tp>>& workspace,
		   FuncTp func,
		   Tp lower, Tp upper,
		   Tp max_abs_err, Tp max_rel_err,
		   Integrator quad = gauss_kronrod_integral<Tp>(Kronrod_15))
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using AreaTp = std::invoke_result_t<FuncTp, Tp>;
      using AbsAreaTp = decltype(std::abs(AreaTp{}));

      const auto s_max = std::numeric_limits<Tp>::max();
      const auto max_iter = workspace.capacity();
      // Try to adjust tests for varing precision.
      const auto s_rel_err = std::pow(Tp{10},
				 -std::numeric_limits<Tp>::digits / Tp{10});
      bool extrapolate = false;
      bool allow_extrapolation = true;

      if (!valid_tolerances(max_abs_err, max_rel_err))
	{
	  std::ostringstream msg;
	  msg << "qags_integrate: Tolerance cannot be achieved with given "
		   "absolute (" << max_abs_err << ") and relative ("
		<< max_rel_err << ") error limits.";
	  throw std::runtime_error(msg.str().c_str());
	}

      workspace.clear();

      // Perform the first integration.

      auto [result0, abserr0, resabs0, resasc0]
	= quad(func, lower, upper);

      workspace.append(lower, upper, result0, abserr0);

      auto tolerance = std::max(max_abs_err,
				  max_rel_err * std::abs(result0));

      // Compute roundoff tolerance.
      const auto round_off = Tp{10} * tolerance * resabs0;

      if (abserr0 <= round_off && abserr0 > tolerance)
	throw integration_error("qags_integrate: "
				"Cannot reach tolerance because "
				"of roundoff error on first attempt",
				ROUNDOFF_ERROR, result0, abserr0);
      else if ((abserr0 <= tolerance && abserr0 != resasc0)
		|| abserr0 == Tp{0})
	return {result0, abserr0};
      else if (max_iter == 1)
	throw integration_error("qags_integrate: "
				"A maximum of one iteration was insufficient",
				MAX_ITER_ERROR, result0, abserr0);

      extrapolation_table<AreaTp, AbsAreaTp> table;
      table.append(result0);

      auto res_ext = result0;
      auto err_ext = s_max;

      auto area = result0;
      auto errsum = abserr0;
      auto iteration = 1u;
      auto ktmin = 0u;
      auto ertest = Tp{0};
      auto error_over_large_intervals = Tp{0};
      auto reseps = AreaTp{0}, abseps = AbsAreaTp{0}, correc = AbsAreaTp{0};
      int error_type = NO_ERROR, error_type2 = NO_ERROR;
      int roundoff_type1 = 0, roundoff_type2 = 0, roundoff_type3 = 0;
      do
	{
	  // Bisect the subinterval with the largest error estimate.
	  const auto& curr = workspace.retrieve();
	  const auto current_depth = workspace.curr_depth() + 1;
	  const auto a1 = curr.lower_lim;
	  const auto mid = (curr.lower_lim + curr.upper_lim) / Tp{2};
	  const auto a2 = mid;
	  const auto b2 = curr.upper_lim;

	  ++iteration;

	  auto [area1, error1, resabs1, resasc1]
	    = quad(func, a1, mid);

	  auto [area2, error2, resabs2, resasc2]
	    = quad(func, a2, b2);

	  const auto area12 = area1 + area2;
	  const auto error12 = error1 + error2;
	  const auto last_e_i = curr.abs_error;
	  const auto delta = area12 - curr.result;

	  // Improve previous approximations to the integral and test for
	  // accuracy.

	  area += delta;
	  errsum += error12 - curr.abs_error;

	  tolerance = std::max(max_abs_err,
				 max_rel_err * std::abs(area));

	  if (resasc1 != error1 && resasc2 != error2)
	    {

	      if (std::abs(delta) <= s_rel_err * std::abs(area12)
		  && error12 >= Tp{0.99} * curr.abs_error)
		{
		  if (!extrapolate)
		    ++roundoff_type1;
		  else
		    ++roundoff_type2;
		}
	      if (iteration > 10 && error12 > curr.abs_error)
		++roundoff_type3;
	    }

	  // Test for roundoff and eventually set error flag.
	  if (roundoff_type1 + roundoff_type2 >= 10
	      || roundoff_type3 >= 20)
	    error_type = ROUNDOFF_ERROR;

	  if (roundoff_type2 >= 5)
	    error_type2 = MAX_ITER_ERROR;

	  // Set error flag in the case of bad integrand behaviour at
	  // a point of the integration range.

	  if (workspace.subinterval_too_small(a1, a2, b2))
	    error_type = EXTRAP_ROUNDOFF_ERROR;

	  // Split the current interval in two.
	  workspace.split(mid, area1, error1, area2, error2);

	  if (errsum <= tolerance)
	    {
	      const auto result = workspace.total_integral();
	      check_error(__func__, error_type, result, errsum);
	      return {result, errsum};
	    }

	  if (error_type != NO_ERROR)
	    break;

	  if (iteration >= max_iter - 1)
	    {
	      error_type = MAX_ITER_ERROR;
	      break;
	    }

	  if (iteration == 2) // Set up variables on first iteration.
	    {
	      error_over_large_intervals = errsum;
	      ertest = tolerance;
	      table.append(area);
	      continue;
	    }

	  if (!allow_extrapolation)
	    continue;

	  error_over_large_intervals -= last_e_i;

	  if (current_depth < workspace.max_depth())
	    error_over_large_intervals += error12;

	  if (!extrapolate)
	    {
	      // Test whether the interval to be bisected next is the
	      // smallest interval.
	      if (workspace.large_interval())
		continue;
	      else
		{
		  extrapolate = true;
		  workspace.increment_curr_index();
		}
	    }

	  // The smallest interval has the largest error.  Before
	  // bisecting decrease the sum of the errors over the larger
	  // intervals (error_over_large_intervals) and perform
	  // extrapolation.

	  if (error_type2 == NO_ERROR
	      && error_over_large_intervals > ertest)
	    if (workspace.increment_curr_index())
	      continue;

	  // Perform extrapolation.
	  table.append(area);
	  std::tie(reseps, abseps) = table.qelg();

	  ++ktmin;
	  if (ktmin > 5 && err_ext < 0.001 * errsum)
	    error_type = DIVERGENCE_ERROR;

	  if (abseps < err_ext)
	    {
	      ktmin = 0;
	      err_ext = abseps;
	      res_ext = reseps;
	      correc = error_over_large_intervals;
	      ertest = std::max(max_abs_err,
				  max_rel_err * std::abs(reseps));
	      if (err_ext <= ertest)
		break;
	    }

	  // Prepare bisection of the smallest interval.
	  if (table.get_nn() == 1)
	    allow_extrapolation = false;

	  if (error_type == DIVERGENCE_ERROR)
	    break;

	  // Work on interval with largest error.
	  workspace.reset_curr_index();
	  extrapolate = false;
	  error_over_large_intervals = errsum;
	}
      while (iteration < max_iter);

      auto result = res_ext;
      auto abserr = err_ext;

      if (err_ext == s_max)
	{
	  const auto result = workspace.total_integral();
	  check_error(__func__, error_type, result, errsum);
	  return {result, errsum};
	}

      if (error_type != NO_ERROR || error_type2 != NO_ERROR)
	{
	  if (error_type2 != NO_ERROR)
	    err_ext += correc;

	  if (error_type == NO_ERROR)
	    error_type = SINGULAR_ERROR;

	  if (res_ext != Tp{0} && area != Tp{0})
	    {
	      if (err_ext / std::abs(res_ext) > errsum / std::abs(area))
		{
		  const auto result = workspace.total_integral();
		  check_error(__func__, error_type, result, errsum);
		  return {result, errsum};
		}
	    }
	  else if (err_ext > errsum)
	    {
	      const auto result = workspace.total_integral();
	      check_error(__func__, error_type, result, errsum);
	      return {result, errsum};
	    }
	  else if (area == Tp{0})
	    {
	      check_error(__func__, error_type, result, errsum);
	      throw std::runtime_error("qags_integrate: Unknown error.");
	    }
	}

      // Test on divergence.
      auto positive_integrand = test_positivity(result0, resabs0);
      auto max_area = std::max(std::abs(res_ext), std::abs(area));
      if (!positive_integrand && max_area < Tp{0.01} * resabs0)
	{
	  check_error(__func__, error_type, area, errsum);
	  throw std::runtime_error("qags_integrate: Unknown error.");
	}

      auto ratio = std::abs(res_ext / area); // FIXME: Added abs for complex type issues.
      if (ratio < Tp{0.01} || ratio > Tp{100}
	  || errsum > std::abs(area))
	error_type = UNKNOWN_ERROR;

      if (error_type == NO_ERROR)
	return {result, abserr};

      check_error(__func__, error_type, decltype(result)(), decltype(abserr)());
      throw integration_error("qags_integrate: Unknown error.",
			      UNKNOWN_ERROR, result, abserr);
    }

  /**
   * Integrate a potentially singular function defined over (-\infty, +\infty).
   */
  template<typename Tp, typename FuncTp>
    auto
    qagi_integrate(integration_workspace<Tp,
			std::invoke_result_t<FuncTp, Tp>>& workspace,
		   FuncTp func,
		   Tp max_abs_err, Tp max_rel_err)
    -> adaptive_integral_t<Tp, decltype(map_minf_pinf<Tp, FuncTp>(func)(Tp{}))>
    {
      using FuncTp2 = decltype(map_minf_pinf<Tp, FuncTp>(func));
      return qags_integrate<Tp, FuncTp2>(workspace,
				map_minf_pinf<Tp, FuncTp>(func),
				Tp{0}, Tp{1}, max_abs_err, max_rel_err);
    }

  /**
   * Integrate a potentially singular symmetric function
   * defined over (-\infty, +\infty).
   */
  template<typename Tp, typename FuncTp>
    auto
    qagis_integrate(integration_workspace<Tp,
			std::invoke_result_t<FuncTp, Tp>>& workspace,
		    FuncTp func,
		    Tp max_abs_err, Tp max_rel_err)
    -> adaptive_integral_t<Tp, decltype(map_minf_pinf_symm<Tp, FuncTp>(func)(Tp{}))>
    {
      using FuncTp2 = decltype(map_minf_pinf_symm<Tp, FuncTp>(func));
      return qags_integrate<Tp, FuncTp2>(workspace,
				map_minf_pinf_symm<Tp, FuncTp>(func),
				Tp{0}, Tp{1}, max_abs_err, max_rel_err);
    }

  /**
   * Integrate a potentially singular function defined over (-\infty, b].
   */
  template<typename Tp, typename FuncTp>
    auto
    qagil_integrate(integration_workspace<Tp,
			std::invoke_result_t<FuncTp, Tp>>& workspace,
		    FuncTp func, Tp upper,
		    Tp max_abs_err, Tp max_rel_err)
    -> adaptive_integral_t<Tp, decltype(map_minf_b<Tp, FuncTp>(func, upper)(Tp{}))>
    {
      using FuncTp2 = decltype(map_minf_b<Tp, FuncTp>(func, Tp{}));
      return qags_integrate<Tp, FuncTp2>(workspace,
				map_minf_b<Tp, FuncTp>(func, upper),
				Tp{0}, Tp{1}, max_abs_err, max_rel_err);
    }

  /**
   * Integrate a potentially singular function defined over [a, +\infty).
   */
  template<typename Tp, typename FuncTp>
    auto
    qagiu_integrate(integration_workspace<Tp,
			std::invoke_result_t<FuncTp, Tp>>& workspace,
		    FuncTp func, Tp lower,
		    Tp max_abs_err, Tp max_rel_err)
    -> adaptive_integral_t<Tp, decltype(map_a_pinf<Tp, FuncTp>(func, lower)(Tp{}))>
    {
      using FuncTp2 = decltype(map_a_pinf<Tp, FuncTp>(func, Tp{}));
      return qags_integrate<Tp, FuncTp2>(workspace,
				map_a_pinf<Tp, FuncTp>(func, lower),
				Tp{0}, Tp{1}, max_abs_err, max_rel_err);
    }

} // namespace __gnu_cxx

#endif // QAGS_INTEGRATE_TCC
