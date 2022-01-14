// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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
// Ported from GSL by Edward Smith-Rowland
// Originally written by Brian Gaugh
//
// Implements integration using a recursive Gauss-Kronrod algorithm
// Based on gsl/integration/qagp.c

#ifndef QAGP_INTEGRATE_TCC
#define QAGP_INTEGRATE_TCC 1

#include <type_traits>
#include <tuple>

#include <ext/integration_workspace.h>
#include <ext/extrapolation_table.h>

namespace __gnu_cxx
{

  /**
   * Adaptively integrate a function with known singular/discontinuous points.
   *
   * @tparam _FuncTp     A function type that takes a single real scalar
   *                     argument and returns a real scalar.
   * @tparam _Tp         A real type for the limits of integration and the step.
   * @tparam _Integrator A non-adaptive integrator that is able to return
   *                     an error estimate in addition to the result.
   *
   * @param[in] __workspace The workspace that manages adaptive quadrature
   * @param[in] __func The single-variable function to be integrated
   * @param[in] __pts The sorted array of points including the integration
   *                  limits and intermediate discontinuities/singularities
   * @param[in] __max_abs_err The limit on absolute error
   * @param[in] __max_rel_err The limit on relative error
   * @param[in] __quad The quadrature stepper taking a function object
   *                   and two integration limits
   */
  template<typename _Tp, typename _FuncTp,
	   typename _Integrator = gauss_kronrod_integral<_Tp>>
    auto
    qagp_integrate(integration_workspace<_Tp,
			std::invoke_result_t<_FuncTp, _Tp>>& __workspace,
		   _FuncTp __func,
		   std::vector<_Tp> __pts,
		   _Tp __max_abs_err, _Tp __max_rel_err,
		   _Integrator __quad = gauss_kronrod_integral<_Tp>(Kronrod_21))
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    {
      const auto _S_max = std::numeric_limits<_Tp>::max();
      const auto __max_iter = __workspace.capacity();
      const auto __n_ivals = __pts.size() - 1;
      // Try to adjust tests for varing precision.
      const auto _S_rel_err = std::pow(_Tp{10},
				 -std::numeric_limits<_Tp>::digits / _Tp{10});

      bool __extrapolate = false;
      bool __allow_extrapolation = true;

      if (!valid_tolerances(__max_abs_err, __max_rel_err))
	{
	  std::ostringstream __msg;
	  __msg << "qagp_integrate: Tolerance cannot be achieved with given "
		   "absolute (" << __max_abs_err << ") and relative ("
		<< __max_rel_err << ") error limits.";
	  std::__throw_runtime_error(__msg.str().c_str());
	}

      if (__pts.size() > __workspace.capacity())
	std::__throw_runtime_error("qagp_integrate: "
				 "number of points exceeds size of workspace");

      // Check that the integration range and break points
      // are in ascending order.

      if (!std::is_sorted(std::begin(__pts), std::end(__pts)))
	std::__throw_runtime_error("qagp_integrate: "
				   "points are not in an ascending sequence");

      __workspace.clear();

      // Perform the first integration.
      auto __result0 = _Tp{0};
      auto __abserr0 = _Tp{0};
      auto __resabs0 = _Tp{0};
      for (std::size_t __i = 0; __i < __n_ivals; ++__i)
	{
	  const auto __lower = __pts[__i];
	  const auto __upper = __pts[__i + 1];

	  auto [__area0, __error0, __resabs0, __resasc0]
	    = __quad(__func, __lower, __upper);

	  __result0 += __area0;
	  __abserr0 += __error0;
	  __resabs0 += __resabs0;
	  std::size_t __level = (__error0 == __resasc0 && __error0 != _Tp{0})
				? 1 : 0;
	  __workspace.append(__lower, __upper, __area0, __error0, __level);
	}

      // Compute the initial error estimate.
      auto __errsum = _Tp{0};
      for (std::size_t __i = 0; __i < __n_ivals; ++__i)
	{
	  if (__workspace.depth(__i) == 1)
	    {
	      __workspace.set_abs_error(__i, __abserr0);
	      __workspace.set_depth(__i, 0);
	    }
	  __errsum += __workspace.abs_error(__i);
	}	

      // We must re-sort because the errors were reassigned.
      __workspace.sort_error();

      auto __tolerance = std::max(__max_abs_err,
				  __max_rel_err * std::abs(__result0));

      // Compute roundoff tolerance.
      const auto __round_off = _Tp{10} * __tolerance * __resabs0;

      if (__abserr0 <= __round_off && __abserr0 > __tolerance)
	throw integration_error("qagp_integrate: "
				"Cannot reach tolerance because "
				"of roundoff error on first attempt",
				ROUNDOFF_ERROR, __result0, __abserr0);
      else if (__abserr0 <= __tolerance)
	return {__result0, __abserr0};
      else if (__max_iter == 1)
	throw integration_error("qagp_integrate: "
				"A maximum of one iteration was insufficient",
				MAX_ITER_ERROR, __result0, __abserr0);

      extrapolation_table<_Tp> __table;
      __table.append(__result0);

      auto __res_ext = __result0;
      auto __err_ext = _S_max;

      auto __area = __result0;
      auto __iteration = __n_ivals - 1;
      auto __ktmin = 0u;
      auto __ertest = __tolerance;
      auto __error_over_large_intervals = __errsum;
      auto __reseps = _Tp{0}, __abseps = _Tp{0}, __correc = _Tp{0};
      int __error_type = NO_ERROR, __error_type2 = NO_ERROR;
      int __roundoff_type1 = 0, __roundoff_type2 = 0, __roundoff_type3 = 0;
      do
	{
	  // Bisect the subinterval with the largest error estimate.
	  const auto& __curr = __workspace.retrieve();

	  const auto __current_depth = __workspace.curr_depth() + 1;

	  const auto __a1 = __curr.__lower_lim;
	  const auto __mid = (__curr.__lower_lim + __curr.__upper_lim) / _Tp{2};
	  const auto __a2 = __mid;
	  const auto __b2 = __curr.__upper_lim;

	  ++__iteration;

	  auto [__area1, __error1, __resabs1, __resasc1]
	    = __quad(__func, __a1, __mid);

	  auto [__area2, __error2, __resabs2, __resasc2]
	    = __quad(__func, __a2, __b2);

	  const auto __area12 = __area1 + __area2;
	  const auto __error12 = __error1 + __error2;
	  const auto __last_e_i = __curr.__abs_error;
	  const auto __delta = __area12 - __curr.__result;

	  // Improve previous approximations to the integral and test for
	  // accuracy.

	  __area += __delta;
	  __errsum += __error12 - __curr.__abs_error;

	  __tolerance = std::max(__max_abs_err,
				 __max_rel_err * std::abs(__area));

	  if (__resasc1 != __error1 && __resasc2 != __error2)
	    {
	      if (std::abs(__delta) <= _S_rel_err * std::abs(__area12)
		  && __error12 >= _Tp{0.99} * __curr.__abs_error)
		{
		  if (!__extrapolate)
		    ++__roundoff_type1;
		  else
		    ++__roundoff_type2;
		}
	      if (__iteration > 10 && __error12 > __curr.__abs_error)
		++__roundoff_type3;
	    }

	  // Test for roundoff and eventually set error flag.

	  if (__roundoff_type1 + __roundoff_type2 >= 10
	   || __roundoff_type3 >= 20)
	    __error_type = ROUNDOFF_ERROR;

	  if (__roundoff_type2 >= 5)
	    __error_type2 = MAX_ITER_ERROR;

	  // Set error flag in the case of bad integrand behaviour at
	  // a point of the integration range.

	  if (__workspace.subinterval_too_small(__a1, __a2, __b2))
	    __error_type = EXTRAP_ROUNDOFF_ERROR;

	  // Split the current interval in two.
	  __workspace.split(__mid, __area1, __error1, __area2, __error2);

	  if (__errsum <= __tolerance)
	    {
	      const auto __result = __workspace.total_integral();
	      check_error(__func__, __error_type, __result, __errsum);
	      return {__result, __errsum};
	    }

	  if (__error_type != NO_ERROR)
	    break;

	  if (__iteration >= __max_iter - 1)
	    {
	      __error_type = MAX_ITER_ERROR;
	      break;
	    }

	  if (!__allow_extrapolation)
	    continue;

	  __error_over_large_intervals -= __last_e_i;

	  if (__current_depth < __workspace.max_depth())
	    __error_over_large_intervals += __error12;

	  if (!__extrapolate)
	    {
	      // Test whether the interval to be bisected next is the
	      // smallest interval.
	      if (__workspace.large_interval())
		continue;
	      else
		{
		  __extrapolate = true;
		  __workspace.increment_curr_index();
		}
	    }

	  // The smallest interval has the largest error.  Before
	  // bisecting decrease the sum of the errors over the larger
	  // intervals (error_over_large_intervals) and perform
	  // extrapolation.

	  if (__error_type2 == NO_ERROR
	   && __error_over_large_intervals > __ertest)
	    if (__workspace.increment_curr_index())
	      continue;

	  // Perform extrapolation.
	  __table.append(__area);

	  if (__table.get_nn() < 3)
	    {
	      __workspace.reset_curr_index();
	      __extrapolate = false;
	      __error_over_large_intervals = __errsum;
	      continue;
	    }

	  std::tie(__reseps, __abseps) = __table.qelg();

	  ++__ktmin;
	  if (__ktmin > 5 && __err_ext < 0.001 * __errsum)
	    __error_type = DIVERGENCE_ERROR;

	  if (__abseps < __err_ext)
	    {
	      __ktmin = 0;
	      __err_ext = __abseps;
	      __res_ext = __reseps;
	      __correc = __error_over_large_intervals;
	      __ertest = std::max(__max_abs_err,
				  __max_rel_err * std::abs(__reseps));
	      if (__err_ext <= __ertest)
		break;
	    }

	  // Prepare bisection of the smallest interval.
	  if (__table.get_nn() == 1)
	    __allow_extrapolation = false;

	  if (__error_type == DIVERGENCE_ERROR)
	    break;

	  // Work on interval with largest error.
	  __workspace.reset_curr_index();
	  __extrapolate = false;
	  __error_over_large_intervals = __errsum;
	}
      while (__iteration < __max_iter);

      auto __result = __res_ext;
      auto __abserr = __err_ext;

      if (__err_ext == _S_max)
	{
	  const auto __result = __workspace.total_integral();
	  check_error(__func__, __error_type, __result, __errsum);
	  return {__result, __errsum};
	}

      if (__error_type != NO_ERROR || __error_type2 != NO_ERROR)
	{
	  if (__error_type2 != NO_ERROR)
	    __err_ext += __correc;

	  if (__error_type == NO_ERROR)
	    __error_type = SINGULAR_ERROR;

	  if (__result != _Tp{0} && __area != _Tp{0})
	    {
	      if (__err_ext / std::abs(__res_ext) > __errsum / std::abs(__area))
		{
		  const auto __result = __workspace.total_integral();
		  check_error(__func__, __error_type, __result, __errsum);
		  return {__result, __errsum};
		}
	    }
	  else if (__err_ext > __errsum)
	    {
	      const auto __result = __workspace.total_integral();
	      check_error(__func__, __error_type, __result, __errsum);
	      return {__result, __errsum};
	    }
	  else if (__area == _Tp{0})
	    {
	      check_error(__func__, __error_type, __result, __abserr);
	      std::__throw_runtime_error("qagp_integrate: Unknown error.");
	    }
	}

      // Test on divergence.
      auto __positive_integrand = __test_positivity(__result0, __resabs0);
      auto __max_area = std::max(std::abs(__res_ext), std::abs(__area));
      if (!__positive_integrand && __max_area < _Tp{0.01} * __resabs0)
	{
	  check_error(__func__, __error_type, __result, __abserr);
	  std::__throw_runtime_error("qagp_integrate: Unknown error.");
	}

      auto __ratio = __res_ext / __area;
      if (__ratio < _Tp{0.01} || __ratio > _Tp{100}
	  || __errsum > std::abs(__area))
	__error_type = UNKNOWN_ERROR;

      if (__error_type == NO_ERROR)
	return {__result, __abserr};

      check_error(__func__, __error_type, __result, __abserr);
      throw integration_error("qagp_integrate: Unknown error.",
			      UNKNOWN_ERROR, __result, __abserr);
    }

} // namespace __gnu_cxx

#endif // QAGP_INTEGRATE_H
