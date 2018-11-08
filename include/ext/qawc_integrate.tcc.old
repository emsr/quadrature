// qawc_integrate.tcc
//
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
// Copyright (C) 2016-2018 Free Software Foundation, Inc.
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
// Ported from GSL by Ed Smith-Rowland
// Originally written by Brian Gaugh
//
// Implements qawc integration algorithm
// Based on gsl/integration/qawc.c

#ifndef QAWC_INTEGRATE_TCC
#define QAWC_INTEGRATE_TCC 1

#include <array>

#include "integration_workspace.h"

namespace __gnu_cxx
{

  template<typename _Tp, typename _FuncTp>
    auto
    qc25c(_FuncTp __func, _Tp __lower, _Tp __upper, _Tp __center)
    -> std::tuple<decltype(_Tp{} * __func(_Tp{})), _Tp, bool>;

  template<typename _Tp>
    std::vector<_Tp>
    compute_moments(std::size_t __N, _Tp __cc);

  /**
   * Adaptive integration for Cauchy principal values:
   * @f[
   *   I = \int_a^b dx f(x) / (x - c)
   * @f]
   * The adaptive bisection algorithm of QAG is used, with modifications
   * to ensure that subdivisions do not occur at the singular point x = c.
   * When a subinterval contains the point x = c or is close to it then
   * a special 25-point modified Clenshaw-Curtis rule is used to control
   * the singularity. Further away from the singularity the algorithm uses
   * a user-supplied integration rule. 
   */
  template<typename _Tp, typename _FuncTp>
    auto
    qawc_integrate(integration_workspace<_Tp,
			std::invoke_result_t<_FuncTp, _Tp>>& __workspace,
		   _FuncTp __func,
		   _Tp __lower, _Tp __upper, _Tp __center,
		   _Tp __max_abs_err, _Tp __max_rel_err)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    {
      auto __result = _Tp{};
      auto __abserr = _Tp{};

      const auto __limit = __workspace.capacity();
      // Try to adjust tests for varing precision.
      const auto _M_rel_err = std::pow(_Tp{10},
				 -std::numeric_limits<_Tp>::digits / _Tp{10});

      int __sign = 1;
      if (__upper < __lower)
	{
	  std::swap(__lower, __upper);
	  __sign = -1;
	}

      if (!valid_tolerances(__max_abs_err, __max_rel_err))
	{
	  std::ostringstream __msg;
	  __msg << "qawc_integrate: Tolerance cannot be achieved with given "
		   "absolute (" << __max_abs_err << ") and relative ("
		<< __max_rel_err << ") error limits.";
	  std::__throw_runtime_error(__msg.str().c_str());
	}

      if (__center == __lower || __center == __upper)
	std::__throw_runtime_error ("qawc_integrate: "
				    "Cannot integrate with singularity "
				    "on endpoint.");

      __workspace.clear();

      // Perform the first integration.
      auto [__result0, __abserr0, __err_reliable]
	= qc25c(__func, __lower, __upper, __center);

      __workspace.append(__lower, __upper, __result0, __abserr0);

      // Test on accuracy; Use 0.01 relative error as an extra safety
      // margin on the first iteration (ignored for subsequent iterations).
      auto __tolerance = std::max(__max_abs_err,
				  __max_rel_err * std::abs( __result0));
      if (__abserr0 < __tolerance && __abserr0 < 0.01 * std::abs(__result0))
	return {__sign * __result0, __abserr0};
      else if (__limit == 1)
	__throw_integration_error("qawc_integrate: "
				  "a maximum of one iteration was insufficient",
				 MAX_ITER_ERROR, __sign * __result0, __abserr0);

      auto __area = __result0;
      auto __errsum = __abserr0;
      auto __iteration = 1u;
      int __error_type = NO_ERROR;
      int __roundoff_type1 = 0, __roundoff_type2 = 0;
      do
	{
	  // Bisect the subinterval with the largest error estimate.
	  const auto& __curr = __workspace.retrieve();

	  const auto __a1 = __curr.__lower_lim;
	  const auto __b2 = __curr.__upper_lim;
	  auto __mid = (__curr.__lower_lim + __curr.__upper_lim) / _Tp{2};
	  if (__center > __a1 && __center <= __mid)
	    __mid = (__center + __b2) / _Tp{2};
	  else if (__center > __mid && __center < __b2)
	    __mid = (__a1 + __center) / _Tp{2};
	  const auto __a2 = __mid;

	  auto [__area1, __error1, __err_reliable1]
	    = qc25c(__func, __a1, __mid, __center);

	  auto [__area2, __error2, __err_reliable2]
	    = qc25c(__func, __a2, __b2, __center);

	  const auto __area12 = __area1 + __area2;
	  const auto __error12 = __error1 + __error2;

	  __errsum += __error12 - __curr.__abs_error;
	  __area += __area12 - __curr.__result;

	  if (__err_reliable1 && __err_reliable2)
	    {
	      const auto __delta = __curr.__result - __area12;

	      if (std::abs (__delta) <= _M_rel_err * std::abs (__area12)
	    	   && __error12 >= 0.99 * __curr.__abs_error)
		++__roundoff_type1;
	      if (__iteration >= 10 && __error12 > __curr.__abs_error)
		++__roundoff_type2;
	    }

	  __tolerance = std::max(__max_abs_err, __max_rel_err * std::abs(__area));
	  if (__errsum > __tolerance)
	    {
	      if (__roundoff_type1 >= 6 || __roundoff_type2 >= 20)
		__error_type = ROUNDOFF_ERROR;

	      // set error flag in the case of bad integrand behaviour at
	      // a point of the integration range
	      if (__workspace.subinterval_too_small(__a1, __a2, __b2))
		__error_type = SINGULAR_ERROR;
	    }

	  __workspace.split(__mid, __area1, __error1, __area2, __error2);

	  ++__iteration;
	}
      while (__iteration < __limit && !__error_type && __errsum > __tolerance);

      __result = __sign * __workspace.total_integral();
      __abserr = __errsum;

      if (__iteration == __limit)
	__error_type = MAX_SUBDIV_ERROR;

      if (__errsum <= __tolerance)
	return {__result, __abserr};

      if (__error_type == NO_ERROR)
	return {__result, __abserr};

      __check_error<_Tp>(__func__, __error_type, __result, __abserr);
      __throw_integration_error("qawc_integrate: Unknown error.",
				UNKNOWN_ERROR, __result, __abserr);
    }

  /**
   * Adaptive integration for Cauchy principal values:
   * @f[
   *   I = \int_a^b dx f(x) / (x - c)
   * @f]
   * The adaptive bisection algorithm of QAG is used, with modifications
   * to ensure that subdivisions do not occur at the singular point x = c.
   * When a subinterval contains the point x = c or is close to it then
   * a special 25-point modified Clenshaw-Curtis rule is used to control
   * the singularity. Further away from the singularity the algorithm uses
   * the 15-point Gauss-Kronrod rule. 
   */
  template<typename _Tp, typename _FuncTp>
    auto
    qc25c(_FuncTp __func, _Tp __lower, _Tp __upper, _Tp __center)
    -> std::tuple<decltype(_Tp{} * __func(_Tp{})), _Tp, bool>
    {
      const auto __cc = (_Tp{2} * __center - __upper - __lower)
		      / (__upper - __lower);

      _Tp __result, __abserr;
      bool __err_reliable;

      if (std::abs(__cc) > 1.1)
	{
	  auto __func_cauchy = [__func, __center](_Tp __x)
				-> _Tp
				{ return __func(__x) / (__x - __center); };

	  auto [__result, __abserr, __resabs, __resasc]
	    = qk_integrate(__func_cauchy, __lower, __upper, Kronrod_15);

	  if (__abserr == __resasc)
	    __err_reliable = false;
	  else
	    __err_reliable = true;

	  return std::make_tuple(__result, __abserr, __err_reliable);
	}
      else
	{
	  auto [__cheb12, __cheb24] = qcheb_integrate(__func, __lower, __upper);
	  const auto __moment = compute_moments(__cheb24.size(), __cc);

	  auto __res12 = _Tp{0};
	  for (size_t __i = 0u; __i < __cheb12.size(); ++__i)
	    __res12 += __cheb12[__i] * __moment[__i];

	  auto __res24 = _Tp{0};
	  for (size_t __i = 0u; __i < __cheb24.size(); ++__i)
	    __res24 += __cheb24[__i] * __moment[__i];

	  __result = __res24;
	  __abserr = std::abs(__res24 - __res12);
	  __err_reliable = false;

	  return std::make_tuple(__result, __abserr, __err_reliable);
	}
    }

  /**
   * Compute Clenshaw-Curtis moments.
   * An iterator range version would be nicer I think.
   */
  template<typename _Tp>
    std::vector<_Tp>
    compute_moments(std::size_t __N, _Tp __cc)
    {
      std::vector<_Tp> __moment(__N);

      auto __a0 = std::log(std::abs((_Tp{1} - __cc) / (_Tp{1} + __cc)));
      auto __a1 = _Tp{2} + __a0 * __cc;

      __moment[0] = __a0;
      __moment[1] = __a1;

      for (size_t __k = 2; __k < __N; ++__k)
	{
	  _Tp __a2;

	  if ((__k % 2) == 0)
	    __a2 = _Tp{2} * __cc * __a1 - __a0;
	  else
	    {
	      const auto __km1 = _Tp(__k - 1);
	      __a2 = _Tp{2} * __cc * __a1
		   - __a0
		   - _Tp{4} / (__km1 * __km1 - _Tp{1});
	    }

	  __moment[__k] = __a2;

	  __a0 = __a1;
	  __a1 = __a2;
	}

      return __moment;
    }

} // namespace

#endif // QAWC_INTEGRATE_TCC
