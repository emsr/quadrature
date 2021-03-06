// integration/qawo_integrate.tcc
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
// Ported from GSL by Ed Smith-Rowland
// Originally written by Brian Gaugh
//
// Implements qawo integration algorithm
// Based on gsl/integration/qawo.c

#ifndef QAWO_INTEGRATE_TCC
#define QAWO_INTEGRATE_TCC 1

#include <cmath>

#include <ext/oscillatory_integration_table.h>

namespace __gnu_cxx
{

  template<typename _Tp, typename _FuncTp>
    auto
    qc25f(oscillatory_integration_table<_Tp>& __wf,
	  _FuncTp __func, _Tp __lower, _Tp __upper,
	  std::size_t __depth)
    -> gauss_kronrod_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

  template<typename _Tp, typename _FuncTp>
    auto
    qawo_integrate(integration_workspace<_Tp,
			std::invoke_result_t<_FuncTp, _Tp>>& __workspace,
		   oscillatory_integration_table<_Tp>& __wf,
		   _FuncTp __func,
		   const _Tp __lower,
		   const _Tp __max_abs_err, const _Tp __max_rel_err)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    {
      const auto _S_max = std::numeric_limits<_Tp>::max();
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto __limit = __workspace.capacity();
      // Try to adjust tests for varing precision. 10^{-5.3} for double.
      const auto _S_rel_err = std::pow(_Tp{10},
				 -std::numeric_limits<_Tp>::digits / _Tp{10});

      bool __extrapolate = false;
      bool __extall = false;
      bool __allow_extrapolation = true;

      __workspace.clear();
      //__wf.clear();
      extrapolation_table<_Tp> __table;

      const auto __upper = __lower + __wf.get_length();
      const auto __abs_omega = std::abs(__wf.omega);

      auto __result = _Tp{0};
      auto __abserr = _Tp{0};

      if (!valid_tolerances(__max_abs_err, __max_rel_err))
	{
	  std::ostringstream __msg;
	  __msg << "qawo_integrate: Tolerance cannot be achieved with given "
		   "absolute (" << __max_abs_err << ") and relative ("
		<< __max_rel_err << ") error limits.";
	  std::__throw_runtime_error(__msg.str().c_str());
	}

      // Perform the first integration.
      auto [__result0, __abserr0, __resabs0, __resasc0]
	= qc25f(__wf, __func, __lower, __upper, 0);

      __workspace.append(__lower, __upper, __result0, __abserr0);

      auto __tolerance = std::max(__max_abs_err,
				  __max_rel_err * std::abs(__result0));

      if (__abserr0 <= _Tp{100} * _S_eps * __resabs0
	  && __abserr0 > __tolerance)
	__throw_integration_error("qawo_integrate: "
				  "Cannot reach tolerance because "
				  "of roundoff error on first attempt",
				  ROUNDOFF_ERROR, __result0, __abserr0);
      else if ((__abserr0 <= __tolerance && __abserr0 != __resasc0)
		|| __abserr0 == _Tp{0})
	return {__result0, __abserr0};
      else if (__limit == 1)
	__throw_integration_error("qawo_integrate: "
				  "A maximum of one iteration was insufficient",
				  MAX_ITER_ERROR, __result0, __abserr0);

      if (0.5 * __abs_omega * std::abs(__upper - __lower) <= _Tp{2})
	{
	  __table.append(__result0);
	  __extall = true;
	}

      auto __res_ext = __result0;
      auto __err_ext = _S_max;

      auto __area = __result0;
      auto __errsum = __abserr0;
      auto __iteration = 1u;
      auto __ktmin = 0u;
      auto __ertest = _Tp{0};
      auto __error_over_large_intervals = _Tp{0};
      auto __correc = _Tp{0};
      int __roundoff_type1 = 0, __roundoff_type2 = 0, __roundoff_type3 = 0;
      int __error_type = NO_ERROR, __error_type2 = NO_ERROR;
      do
	{
	  // Bisect the subinterval with the largest error estimate.
	  const auto& __curr = __workspace.retrieve();

	  const auto __current_depth = __workspace.curr_depth() + 1;

	  if (__current_depth >= __wf.n)
	    {
	      __error_type = -1; // exceeded limit of table.
	      break;
	    }

	  const auto __a1 = __curr.__lower_lim;
	  const auto __mid = (__curr.__lower_lim + __curr.__upper_lim) / _Tp{2};
	  const auto __a2 = __mid;
	  const auto __b2 = __curr.__upper_lim;

	  ++__iteration;

	  auto [__area1, __error1, __resabs1, __resasc1]
	    = qc25f(__wf, __func, __a1, __mid, __current_depth);

	  auto [__area2, __error2, __resabs2, __resasc2]
	    = qc25f(__wf, __func, __a2, __b2, __current_depth);

	  const auto __area12 = __area1 + __area2;
	  const auto __error12 = __error1 + __error2;
	  const auto __last_e_i = __curr.__abs_error;

	  // Improve previous approximations to the integral and test
	  // for accuracy.
	  __area += __area12 - __curr.__result;
	  __errsum += __error12 - __curr.__abs_error;

	  __tolerance = std::max(__max_abs_err,
				 __max_rel_err * std::abs(__area));

	  if (__resasc1 != __error1 && __resasc2 != __error2)
	    {
	      const auto __delta = __curr.__result - __area12;

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
	      __result = __workspace.total_integral();
	      __abserr = __errsum;
	      __check_error<_Tp>(__func__, __error_type, __result, __abserr);
	      return {__result, __abserr};
	    }

	  if (__error_type != NO_ERROR)
	    break;

	  if (__iteration >= __limit - 1)
	    {
	      __error_type = MAX_ITER_ERROR;
	      break;
	    }

	  // Set up variables on first iteration.
	  if (__iteration == 2 && __extall)
	    {
	      __error_over_large_intervals = __errsum;
	      __ertest = __tolerance;
	      __table.append(__area);
	      continue;
	    }

	  if (!__allow_extrapolation)
	    continue;

	  if (__extall)
	    {
	      __error_over_large_intervals -= __last_e_i;

	      if (__current_depth < __workspace.max_depth())
		__error_over_large_intervals += __error12;

	      if (__extrapolate)
		if (__error_type2 != NO_ERROR
		 && __error_over_large_intervals > __ertest)
		  if (__workspace.increment_curr_index())
		    continue;
	    }

	  if (__workspace.large_interval())
	    continue;

	  if (__extall)
	    {
	      __extrapolate = true;
	      __workspace.increment_curr_index();
	    }
	  else
	    {
	      // Test whether the interval to be bisected next is the
	      // smallest interval.
	      auto __width = __workspace.upper_lim() - __workspace.lower_lim();

	      if (0.25 * std::abs(__width) * __abs_omega > _Tp{2})
		continue;

	      __extall = true;
	      __error_over_large_intervals = __errsum;
	      __ertest = __tolerance;
	      continue;
	    }

	  if (__error_type2 != NO_ERROR
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

	  _Tp __reseps, __abseps;
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
      while (__iteration < __limit);

      __result = __res_ext;
      __abserr = __err_ext;

      if (__err_ext == _S_max)
	{
	  __result = __workspace.total_integral();
	  __abserr = __errsum;
	  __check_error<_Tp>(__func__, __error_type, __result, __abserr);
	  return {__result, __abserr};
	}

      if (__error_type != NO_ERROR || __error_type2 != NO_ERROR)
	{
	  if (__error_type2 != NO_ERROR)
	    __err_ext += __correc;

	  if (__error_type == NO_ERROR)
	    __error_type = SINGULAR_ERROR;

	  if (__result != 0 && __area != 0)
	    {
	      if (__err_ext / std::abs(__res_ext) > __errsum / std::abs(__area))
		{
		  __result = __workspace.total_integral();
		  __abserr = __errsum;
		  __check_error<_Tp>(__func__, __error_type,
				     __result, __abserr);
		  return {__result, __abserr};
		}
	    }
	  else if (__err_ext > __errsum)
	    {
	      __result = __workspace.total_integral();
	      __abserr = __errsum;
	      __check_error<_Tp>(__func__, __error_type, __result, __abserr);
	      return {__result, __abserr};
	    }
	  else if (__area == _Tp{0})
	    {
	      __check_error<_Tp>(__func__, __error_type, __result, __abserr);
	      return {__result, __abserr};
	    }
	}

      // Test on divergence.
      bool __positive_integrand = __test_positivity(__result0, __resabs0);
      auto __max_area = std::max(std::abs(__res_ext), std::abs(__area));
      if (!__positive_integrand && __max_area < _Tp{0.01} * __resabs0)
	{
	  __check_error<_Tp>(__func__, __error_type, __result, __abserr);
	  return {__result, __abserr};
	}

      auto __ratio = __res_ext / __area;
      if (__ratio < _Tp{0.01} || __ratio > _Tp{100}
	  || __errsum > std::abs(__area))
	__error_type = UNKNOWN_ERROR;

      __check_error<_Tp>(__func__, __error_type, __result, __abserr);
      return {__result, __abserr};
    }


  template<typename _Tp, typename _FuncTp>
    auto
    qc25f(oscillatory_integration_table<_Tp>& __wf,
	  _FuncTp __func, _Tp __lower, _Tp __upper,
	  std::size_t __depth)
    -> gauss_kronrod_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    {
      const auto _S_max = std::numeric_limits<_Tp>::max();
      const auto __center = ( __lower + __upper) / _Tp{2};
      const auto __half_length = (__upper - __lower) / _Tp{2};
      const auto __omega = __wf.omega;

      const auto __par = __omega * __half_length;

      if (std::abs(__par) < _Tp{2})
	{
	  if (__wf.circfun == oscillatory_integration_table<_Tp>::INTEG_SINE)
	    {
	      auto __wfun = [__func, __omega](_Tp __x)
			    ->_Tp
			    { return std::sin(__omega * __x) * __func(__x); };
	      return qk_integrate(__wfun, __lower, __upper, Kronrod_15);
	    }
	  else
	    {
	      auto __wfun = [__func, __omega](_Tp __x)
			    ->_Tp
			    { return std::cos(__omega * __x) * __func(__x); };
	      return qk_integrate(__wfun, __lower, __upper, Kronrod_15);
	    }
	}
      else
	{
	  auto __chout = qcheb_integrate(__func, __lower, __upper);
	  const auto& __cheb12 = __chout.__cheb12;
	  const auto& __cheb24 = __chout.__cheb24;

	  if (__depth >= __wf.n)
	    {
	      // Table overflow should not happen, check before calling.
	      std::__throw_runtime_error("Table overflow in internal function");
	    }

	  const auto& __moment = __wf.get_moments(__depth);

	  auto __res12_cos = __cheb12[12] * __moment[12];
	  auto __res12_sin = _Tp{0};
	  for (int __i = 0; __i < 6; ++__i)
	    {
	      const std::size_t __k = 10 - 2 * __i;
	      __res12_cos += __cheb12[__k] * __moment[__k];
	      __res12_sin += __cheb12[__k + 1] * __moment[__k + 1];
	    }

	  auto __res24_cos = __cheb24[24] * __moment[24];
	  auto __res24_sin = _Tp{0};
	  auto __result_abs = std::abs(__cheb24[24]);
	  for (int __i = 0; __i < 12; ++__i)
	    {
	      const std::size_t __k = 22 - 2 * __i;
	      __res24_cos += __cheb24[__k] * __moment[__k];
	      __res24_sin += __cheb24[__k + 1] * __moment[__k + 1];
	      __result_abs += std::abs(__cheb24[__k])
			    + std::abs(__cheb24[__k + 1]);
	    }

	  const auto __est_cos = std::abs(__res24_cos - __res12_cos);
	  const auto __est_sin = std::abs(__res24_sin - __res12_sin);

	  const auto __c = __half_length * std::cos(__center * __omega);
	  const auto __s = __half_length * std::sin(__center * __omega);

	  _Tp __result, __abserr, __resabs, __resasc;
	  if (__wf.circfun == oscillatory_integration_table<_Tp>::INTEG_SINE)
	    {
	      __result = __c * __res24_sin + __s * __res24_cos;
	      __abserr = std::abs(__c * __est_sin) + std::abs(__s * __est_cos);
	    }
	  else
	    {
	      __result = __c * __res24_cos - __s * __res24_sin;
	      __abserr = std::abs(__c * __est_cos) + std::abs(__s * __est_sin);
	    }

	  __resabs = __result_abs * __half_length;
	  __resasc = _S_max;

	  return {__result, __abserr, __resabs, __resasc};
	}
    }

} // namespace __gnu_cxx

#endif // QAWO_INTEGRATE_TCC

