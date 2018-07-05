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
//
// Ported from GSL by Edward Smith-Rowland
// Originally written by Brian Gaugh
//
// This file implements Gauss-Kronrod integration with singularities
// Based on gsl/integration/qags.c

#ifndef QAGS_INTEGRATE_TCC
#define QAGS_INTEGRATE_TCC 1

#include <utility>
#include <limits>
#include <tuple>
#include <stdexcept>

#include "integration_error.h"
#include "integration_transform.h"
#include "qk_integrate.tcc"
#include "integration_workspace.h"
#include "extrapolation_table.h"

namespace __gnu_cxx
{

  /**
   * Integrate potentially singular function from a to b using recursive
   * Gauss-Kronrod algorithm.
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
   * @param[in] __max_iter The maximum number of integration steps allowed
   * @param[in] __max_abs_err The limit on absolute error
   * @param[in] __max_rel_err The limit on relative error
   * @param[in] __quad The quadrature stepper taking a function object
   *                   and two integration limits
   */
  template<typename _Tp, typename _FuncTp, typename _Integrator>
    std::tuple<_Tp, _Tp>
    qags_integrate(integration_workspace<_Tp>& __workspace,
		   _FuncTp __func,
		   _Tp __lower, _Tp __upper,
		   _Tp __max_abs_err, _Tp __max_rel_err,
		   _Integrator __quad)
    {
      const auto _S_max = std::numeric_limits<_Tp>::max();
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto __max_iter = __workspace.capacity();
      // Try to adjust tests for varing precision.
      const auto _S_rel_err = std::pow(_Tp{10},
				 -std::numeric_limits<_Tp>::digits / _Tp{10});
#ifdef VERBOSE_DEBUG
      std::cerr.precision(8);
      std::cerr << "\nC++\n";
#endif
      bool __extrapolate = false;
      bool __allow_extrapolation = true;

      if (valid_tolerances(__max_abs_err, __max_rel_err))
	std::__throw_runtime_error("qags_integrate: "
				   "Tolerance cannot be achieved "
				   "with given absolute "
				   "and relative error limits.");

      __workspace.clear();

      // Perform the first integration.

      _Tp __result0, __abserr0, __resabs0, __resasc0;
      std::tie(__result0, __abserr0, __resabs0, __resasc0)
	  = __quad(__func, __lower, __upper);
#ifdef VERBOSE_DEBUG
      std::cerr << "  result0 = " << __result0
		<< "  abserr0 = " << __abserr0
		<< "  resabs0 = " << __resabs0
		<< "  resasc0 = " << __resasc0
		<< "  abserr0 == resasc0 : "
		<< (__abserr0 == __resasc0) << '\n';
#endif

      __workspace.append(__lower, __upper, __result0, __abserr0);
#ifdef VERBOSE_DEBUG
      dump_ws(__workspace, "qags", "first quad");
#endif

      auto __tolerance = std::max(__max_abs_err,
				  __max_rel_err * std::abs(__result0));

      if (__abserr0 <= 100 * _S_eps * __resabs0
	  && __abserr0 > __tolerance)
	__throw_integration_error("qags_integrate: "
				  "cannot reach tolerance because "
				  "of roundoff error on first attempt",
				  ROUNDOFF_ERROR, __result0, __abserr0);
      else if ((__abserr0 <= __tolerance && __abserr0 != __resasc0)
		|| __abserr0 == _Tp{0})
	return std::make_tuple(__result0, __abserr0);
      else if (__max_iter == 1)
	__throw_integration_error("qags_integrate: "
				  "a maximum of one iteration was insufficient",
				  MAX_ITER_ERROR, __result0, __abserr0);

      extrapolation_table<_Tp> __table;
      __table.append(__result0);

      auto __res_ext = __result0;
      auto __err_ext = _S_max;

      auto __area = __result0;
      auto __errsum = __abserr0;
      auto __iteration = 1u;
      auto __ktmin = 0u;
      auto __ertest = _Tp{0};
      auto __error_over_large_intervals = _Tp{0};
      auto __reseps = _Tp{0}, __abseps = _Tp{0}, __correc = _Tp{0};
      int __error_type = NO_ERROR, __error_type2 = NO_ERROR;
      int __roundoff_type1 = 0, __roundoff_type2 = 0, __roundoff_type3 = 0;
      do
	{
	  // Bisect the subinterval with the largest error estimate.
	  const auto& __curr = __workspace.retrieve();

	  const auto __current_depth = __workspace.curr_depth() + 1;
#ifdef VERBOSE_DEBUG
	  std::cerr << "current_depth = " << __current_depth << '\n';
#endif

	  const auto __a1 = __curr.__lower_lim;
	  const auto __mid = (__curr.__lower_lim + __curr.__upper_lim) / _Tp{2};
	  const auto __a2 = __mid;
	  const auto __b2 = __curr.__upper_lim;

	  ++__iteration;

	  _Tp __area1, __error1, __resabs1, __resasc1;
	  std::tie(__area1, __error1, __resabs1, __resasc1)
	    = __quad(__func, __a1, __mid);

	  _Tp __area2, __error2, __resabs2, __resasc2;
	  std::tie(__area2, __error2, __resabs2, __resasc2)
	    = __quad(__func, __a2, __b2);

	  const auto __area12 = __area1 + __area2;
	  const auto __error12 = __error1 + __error2;
	  const auto __last_e_i = __curr.__abs_error;

	  // Improve previous approximations to the integral and test for
	  // accuracy.

	  __errsum += __error12 - __curr.__abs_error;
	  __area += __area12 - __curr.__result;

	  __tolerance = std::max(__max_abs_err,
				 __max_rel_err * std::abs(__area));

	  if (__resasc1 != __error1 && __resasc2 != __error2)
	    {
	      const auto __delta = __curr.__result - __area12;

	      if (std::abs(__delta) <= _S_rel_err * std::abs(__area12)
		  && __error12 >= 0.99 * __curr.__abs_error)
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
#ifdef VERBOSE_DEBUG
	  dump_ws(__workspace, "qags", "split");
	  std::cerr << "errsum    = " << __errsum << '\n';
	  std::cerr << "tolerance = " << __tolerance << '\n';
#endif

	  if (__errsum <= __tolerance)
	    {
	      const auto __result = __workspace.total_integral();
	      __check_error(__func__, __error_type, __result, __errsum);
	      return std::make_tuple(__result, __errsum);
	    }

	  if (__error_type != NO_ERROR)
	    break;

	  if (__iteration >= __max_iter - 1)
	    {
	      __error_type = MAX_ITER_ERROR;
	      break;
	    }

	  if (__iteration == 2) // Set up variables on first iteration.
	    {
	      __error_over_large_intervals = __errsum;
	      __ertest = __tolerance;
	      __table.append(__area);
	      continue;
	    }

	  if (!__allow_extrapolation)
	    continue;

	  __error_over_large_intervals -= __last_e_i;

	  if (__current_depth < __workspace.max_depth())
	    __error_over_large_intervals += __error12;
#ifdef VERBOSE_DEBUG
	  std::cerr << "last_e_i    = " << __last_e_i << "\n";
	  std::cerr << "error12     = " << __error12 << "\n";
	  std::cerr << "eoli        = " << __error_over_large_intervals << "\n";
	  std::cerr << "extrapolate = " << __extrapolate << "\n";
	  std::cerr << "wkli        = " << __workspace.large_interval()<< "\n";
	  std::cerr << "error_type2 : " << __error_type2
	   << "  error_over_large_intervals = " << __error_over_large_intervals
	   << "  ertest = " << __ertest << "\n";
	  std::cerr << "extrapolate : " << __extrapolate << "\n";
	  std::cerr << "large_interval : " << __workspace.large_interval()
		    << "\n";
#endif

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
#ifdef VERBOSE_DEBUG
	  dump_ws(__workspace, "qagp", "extrap");
	  std::cerr << "nn          = " << __table.get_nn() << "\n";
#endif

	  std::tie(__reseps, __abseps) = __table.qelg();
#ifdef VERBOSE_DEBUG
	  std::cerr << "reseps      = " << __reseps
		    << "  abseps      = " << __abseps << "\n";
#endif

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
#ifdef VERBOSE_DEBUG
	  dump_ws(__workspace, "qagp", "stop extrap");
#endif
	}
      while (__iteration < __max_iter);

      auto __result = __res_ext;
      auto __abserr = __err_ext;

      if (__err_ext == _S_max)
	{
	  const auto __result = __workspace.total_integral();
	  __check_error(__func__, __error_type, __result, __errsum);
	  return std::make_tuple(__result, __errsum);
	}

      if (__error_type != NO_ERROR || __error_type2 != NO_ERROR)
	{
	  if (__error_type2 != NO_ERROR)
	    __err_ext += __correc;

	  if (__error_type == NO_ERROR)
	    __error_type = SINGULAR_ERROR;

	  if (__res_ext != _Tp{0} && __area != _Tp{0})
	    {
	      if (__err_ext / std::abs(__res_ext) > __errsum / std::abs(__area))
		{
		  const auto __result = __workspace.total_integral();
		  __check_error(__func__, __error_type, __result, __errsum);
		  return std::make_tuple(__result, __errsum);
		}
	    }
	  else if (__err_ext > __errsum)
	    {
	      const auto __result = __workspace.total_integral();
	      __check_error(__func__, __error_type, __result, __errsum);
	      return std::make_tuple(__result, __errsum);
	    }
	  else if (__area == _Tp{0})
	    {
	      __check_error<_Tp>(__func__, __error_type);
	      std::__throw_runtime_error("qags_integrate: Unknown error.");
	    }
	}

      // Test on divergence.
      auto __positive_integrand = __test_positivity(__result0, __resabs0);
      auto __max_area = std::max(std::abs(__res_ext), std::abs(__area));
      if (!__positive_integrand && __max_area < 0.01 * __resabs0)
	{
	  __check_error<_Tp>(__func__, __error_type);
	  std::__throw_runtime_error("qags_integrate: Unknown error.");
	}

      auto __ratio = __res_ext / __area;
      if (__ratio < 0.01 || __ratio > 100.0 || __errsum > std::abs(__area))
	__error_type = UNKNOWN_ERROR;

      if (__error_type == NO_ERROR)
	return std::make_tuple(__result, __abserr);

      __check_error<_Tp>(__func__, __error_type);
      __throw_integration_error("qags_integrate: Unknown error.",
				UNKNOWN_ERROR, __result, __abserr);
    }

  /**
   * Specialize to Gauss-Kronrod integration step with default 21-point rule.
   *
   * @param[in] __workspace The workspace that manages adaptive quadrature
   * @param[in] __func The single-variable function to be integrated
   * @param[in] __lower The lower limit of integration
   * @param[in] __upper The upper limit of integration
   * @param[in] __max_iter The maximum number of integration steps allowed
   * @param[in] __max_abs_err The limit on absolute error
   * @param[in] __max_rel_err The limit on relative error
   * @param[in] __max_iter The maximum number of iterations
   * @param[in] __qkintrule The size of the Gauss-Kronrod integration scheme
   */
  template<typename _Tp, typename _FuncTp>
    std::tuple<_Tp, _Tp>
    qags_integrate(integration_workspace<_Tp>& __workspace,
		  _FuncTp __func,
		  _Tp __lower, _Tp __upper,
		  _Tp __max_abs_err, _Tp __max_rel_err,
                  Kronrod_Rule __qk_rule = QK_21)
    {
      auto __quad
	= [__qk_rule]
	  (_FuncTp __func, _Tp __lower, _Tp __upper)
	  -> std::tuple<_Tp, _Tp, _Tp, _Tp>
	  { return qk_integrate(__func, __lower, __upper, __qk_rule); };

      return qags_integrate(__workspace, __func, __lower, __upper,
			    __max_abs_err, __max_rel_err, __quad);
    }

  /**
   * Integrate a potentially singular function defined over (-\infty, +\infty).
   */
  template<typename _Tp, typename _FuncTp>
    std::tuple<_Tp, _Tp>
    qagi_integrate(integration_workspace<_Tp>& __workspace,
		   _FuncTp __func,
		   _Tp __max_abs_err, _Tp __max_rel_err)
    {
      return qags_integrate(__workspace, map_minf_pinf<_Tp, _FuncTp>(__func),
			    _Tp{0}, _Tp{1}, __max_abs_err, __max_rel_err);
    }

  /**
   * Integrate a potentially singular function defined over (-\infty, b].
   */
  template<typename _Tp, typename _FuncTp>
    std::tuple<_Tp, _Tp>
    qagil_integrate(integration_workspace<_Tp>& __workspace,
		    _FuncTp __func, _Tp __upper,
		    _Tp __max_abs_err, _Tp __max_rel_err)
    {
      return qags_integrate(__workspace, map_minf_b(__func, __upper),
			    _Tp{0}, _Tp{1}, __max_abs_err, __max_rel_err);
    }

  /**
   * Integrate a potentially singular function defined over [a, +\infty).
   */
  template<typename _Tp, typename _FuncTp>
    std::tuple<_Tp, _Tp>
    qagiu_integrate(integration_workspace<_Tp>& __workspace,
		    _FuncTp __func, _Tp __lower,
		    _Tp __max_abs_err, _Tp __max_rel_err)
    {
      return qags_integrate(__workspace, map_a_pinf(__func, __lower),
			    _Tp{0}, _Tp{1}, __max_abs_err, __max_rel_err);
    }

} // namespace __gnu_cxx

#endif // QAGS_INTEGRATE_TCC
