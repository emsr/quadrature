// -*- C++ -*-
// Integration utilities for C++.
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
// Implements integration using a recursive locally-adaptive algorithm
// Based on gsl/integration/qag.c

#ifndef QAG_INTEGRATE_TCC
#define QAG_INTEGRATE_TCC 1

#include <utility>
#include <limits>
#include <string>
#include <stdexcept>

#include "integration_workspace.h"

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
   * @tparam _FuncTp     A function type that takes a single real scalar
   *                     argument and returns a real scalar.
   * @tparam _Tp         A real type for the limits of integration and the step.
   * @tparam _Integrator A non-adaptive integrator that is able to return
   *                     an error estimate in addition to the result.
   *
   * @param[in] __workspace The workspace that manages adaptive quadrature
   * @param[in] __func The single-variable function to be integrated
   * @param[in] __lower The lower limit of integration
   * @param[in] __upper The upper limit of integration
   * @param[in] __max_iter The maximum number of integration steps allowed
   * @param[in] __max_abs_err The limit on absolute error
   * @param[in] __max_rel_err The limit on relative error
   * @param[in] __max_iter The maximum number of iterations
   * @param[in] __quad The quadrature stepper taking a function object
   *                   and two integration limits
   *
   * @return A tuple with the first value being the integration result,
   *	     and the second value being the estimated error.
   */
  template<typename _Tp, typename _FuncTp,
	   typename _Integrator = gauss_kronrod_integral<_Tp>>
    auto
    qag_integrate(integration_workspace<_Tp,
			std::invoke_result_t<_FuncTp, _Tp>>& __workspace,
		  _FuncTp __func,
		  _Tp __lower, _Tp __upper,
		  _Tp __max_abs_err, _Tp __max_rel_err,
		  std::size_t __max_iter,
		  _Integrator __quad = gauss_kronrod_integral<_Tp>(Kronrod_21))
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      // Try to adjust tests for varing precision.
      const auto _S_rel_err = std::pow(_Tp{10},
				 -std::numeric_limits<_Tp>::digits / _Tp{10});

      if (!valid_tolerances(__max_abs_err, __max_rel_err))
	{
	  std::ostringstream __msg;
	  __msg << "qag_integrate: Tolerance cannot be achieved with given "
		   "absolute (" << __max_abs_err << ") and relative ("
		<< __max_rel_err << ") error limits.";
	  std::__throw_runtime_error(__msg.str().c_str());
	}

      auto [__result0, __abserr0, __resabs0, __resasc0]
	= __quad(__func, __lower, __upper);

      // Test on accuracy
      auto __tolerance = std::max(__max_abs_err,
				  __max_rel_err * std::abs(__result0));

      // Compute roundoff
      const auto __round_off = _Tp{100} * _S_eps * __resabs0;

      if (__abserr0 <= __round_off && __abserr0 > __tolerance)
	__throw_integration_error("qag_integrate: "
				  "cannot reach tolerance because "
				  "of roundoff error on first attempt",
				  ROUNDOFF_ERROR, __result0, __abserr0);
      else if ((__abserr0 <= __tolerance && __abserr0 != __resasc0)
		|| __abserr0 == _Tp{0})
	return {__result0, __abserr0};
      else if (__max_iter == 1)
	__throw_integration_error("qag_integrate: "
				  "a maximum of one iteration was insufficient",
				  MAX_ITER_ERROR, __result0, __abserr0);

      __workspace.clear();
      __workspace.append(__lower, __upper, __result0, __abserr0);

      auto __area = __result0;
      auto __errsum = __abserr0;
      int __error_type = NO_ERROR;
      std::size_t __iteration = 1;

      int __roundoff_type1 = 0, __roundoff_type2 = 0;
      do
	{
	  // Bisect the subinterval with the largest error estimate
	  const auto& __curr = __workspace.retrieve();

	  const auto __a1 = __curr.__lower_lim;
	  const auto __mid = (__curr.__lower_lim + __curr.__upper_lim) / _Tp{2};
	  const auto __a2 = __mid;
	  const auto __b2 = __curr.__upper_lim;

	  auto [__area1, __error1, __resabs1, __resasc1]
	    = __quad(__func, __a1, __mid);

	  auto [__area2, __error2, __resabs2, __resasc2]
	    = __quad(__func, __a2, __b2);

	  const auto __area12 = __area1 + __area2;
	  const auto __error12 = __error1 + __error2;

	  __area += __area12 - __curr.__result;
	  __errsum += __error12 - __curr.__abs_error;

	  if (__resasc1 != __error1 && __resasc2 != __error2)
	    {
	      const auto __delta = __curr.__result - __area12;

	      if (std::abs(__delta) <= _S_rel_err * std::abs(__area12)
		 && __error12 >= 0.99 * __curr.__abs_error)
		++__roundoff_type1;
	      if (__iteration >= 10 && __error12 > __curr.__abs_error)
		++__roundoff_type2;
	    }

	  __tolerance = std::max(__max_abs_err,
				 __max_rel_err * std::abs(__area));

	  if (__errsum > __tolerance)
	    {
	      if (__roundoff_type1 >= 6 || __roundoff_type2 >= 20)
		__error_type = ROUNDOFF_ERROR;

	      // Set error flag in the case of bad integrand behaviour at
	      // a point of the integration range.
	      if (__workspace.subinterval_too_small(__a1, __a2, __b2))
		__error_type = SINGULAR_ERROR;
	    }

	  __workspace.split(__mid, __area1, __error1, __area2, __error2);

	  ++__iteration;
	}
      while (__iteration < __max_iter
	     && !__error_type
	     && __errsum > __tolerance);

      auto __result = __workspace.total_integral();
      auto __abserr = __errsum;

      if (__errsum <= __tolerance)
	return {__result, __abserr};

      if (__error_type == NO_ERROR && __iteration >= __max_iter)
	__error_type = MAX_ITER_ERROR;

      __check_error<_Tp>(__func__, __error_type, __result, __abserr);
      __throw_integration_error("qag_integrate: Unknown error.",
				UNKNOWN_ERROR, __result, __abserr);
    }

} // namespace __gnu_cxx

#endif // QAG_INTEGRATE_TCC
