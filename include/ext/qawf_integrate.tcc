// integration/qawf.c
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
// Implements qawf integration algorithm
// Based on gsl/integration/qawf.c

#ifndef QAWF_INTEGRATE_TCC
#define QAWF_INTEGRATE_TCC 1

#include <cmath>

#include <ext/integration_workspace.h>
#include <ext/oscillatory_integration_table.h>

namespace __gnu_cxx
{

  /**
   * This function attempts to compute a Fourier integral of the function f
   * over the semi-infinite interval [a,+\infty)
   */
  template<typename _Tp, typename _FuncTp>
    auto
    qawf_integrate(integration_workspace<_Tp,
			std::invoke_result_t<_FuncTp, _Tp>>& __workspace,
		   integration_workspace<_Tp,
			std::invoke_result_t<_FuncTp, _Tp>>& __cycle_workspace,
		   oscillatory_integration_table<_Tp>& __wf,
		   _FuncTp __func,
		   _Tp __lower, _Tp __max_abs_err)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    {
      std::size_t __ktmin = 0;
      std::size_t __iteration = 0;

      extrapolation_table<_Tp> __table;

      _Tp __cycle;
      auto __omega = __wf.omega;

      const _Tp __p = 0.9;
      _Tp __factor = 1;
      _Tp __initial_eps, __eps;
      int __error_type = NO_ERROR;

      int __status = 0;
      auto __result = _Tp{0};
      auto __abserr = _Tp{0};

      const auto _S_max = std::numeric_limits<_Tp>::max();
      const auto __limit = __workspace.capacity();

      __workspace.clear();
      __cycle_workspace.clear();
      //__wf.clear();

      // Test on accuracy.
      if (__max_abs_err <= _Tp{0})
	std::__throw_domain_error("absolute tolerance epsabs must be positive") ;

      if (__omega == _Tp{0})
	{
	  if (__wf.circfun == oscillatory_integration_table<_Tp>::INTEG_SINE)
	    // The function sin(w x) f(x) is always zero for w = 0.
	    return {_Tp{0}, _Tp{0}};
	  else
	    // The function cos(w x) f(x) is always f(x) for w = 0.
	    return qagiu_integrate(__cycle_workspace, __func, __lower, __max_abs_err,
				   _Tp{0});
	}

      if (__max_abs_err * (_Tp{1} - __p) > std::numeric_limits<_Tp>::min())
	__eps = __max_abs_err * (_Tp{1} - __p);
      else
	__eps = __max_abs_err;

      __initial_eps = __eps;

      auto __area = _Tp{0};
      auto __errsum = _Tp{0};
      auto __res_ext = _Tp{0};
      auto __err_ext = _S_max;
      auto __correc = _Tp{0};
      auto __total_error = _Tp{0};
      auto __truncation_error = _Tp{0};

      __cycle = (2 * std::floor(std::abs(__omega)) + 1)
	      * M_PI / std::abs(__omega);

      __wf.set_length(__cycle);

      for (__iteration = 0; __iteration < __limit; ++__iteration)
	{
	  const auto __a1 = __lower + __iteration * __cycle;
	  const auto __b1 = __a1 + __cycle;

	  const auto __max_abs_err1 = __eps * __factor;

	  auto __out1 = qawo_integrate(__cycle_workspace, __wf, __func, __a1,
			     __max_abs_err1, _Tp{0});
	  auto __area1 = __out1.__result;
	  auto __error1 = __out1.__abserr;

	  __workspace.append(__a1, __b1, __area1, __error1);

	  __factor *= __p;

	  __area += __area1;
	  __errsum += __error1;

	  // Estimate the truncation error as 50 times the final term.
	  __truncation_error = 50 * std::abs(__area1);
	  __total_error = __errsum + __truncation_error;
	  if (__total_error < __max_abs_err && __iteration > 4)
	    goto compute_result;

	  if (__error1 > __correc)
	    __correc = __error1;

	  if (__status)
	    __eps = std::max(__initial_eps, __correc * (_Tp{1} - __p));

	  if (__status && __total_error < 10 * __correc && __iteration > 3)
	    goto compute_result;

	  __table.append(__area);

	  if (__table.get_nn() < 2)
	    continue;

	  _Tp __reseps, __erreps;
	  std::tie(__reseps, __erreps) = __table.qelg();

	  ++__ktmin;
	  if (__ktmin >= 15 && __err_ext < _Tp{0.001L} * __total_error)
	    __error_type = EXTRAP_ROUNDOFF_ERROR;

	  if (__erreps < __err_ext)
	    {
	      __ktmin = 0;
	      __err_ext = __erreps;
	      __res_ext = __reseps;

	      if (__err_ext + 10 * __correc <= __max_abs_err)
		break;
	      if (__err_ext <= __max_abs_err && 10 * __correc >= __max_abs_err)
		break;
	    }
	}

      if (__iteration == __limit)
	__error_type = MAX_ITER_ERROR;

      if (__err_ext == _S_max)
	goto compute_result;

      __err_ext += 10 * __correc;

      __result = __res_ext;
      __abserr = __err_ext;

      if (__error_type == NO_ERROR)
	return {__result, __abserr};

      if (__res_ext != _Tp{0} && __area != _Tp{0})
	{
	  if (__err_ext / std::abs(__res_ext) > __errsum / std::abs(__area))
	    goto compute_result;
	}
      else if (__err_ext > __errsum)
	goto compute_result;
      else if (__area == _Tp{0})
	goto return_error;

      if (__error_type == EXTRAP_ROUNDOFF_ERROR)
	__err_ext += __truncation_error;

      goto return_error;

    compute_result:

      __result = __area;
      __abserr = __total_error;

      if (__error_type == NO_ERROR)
	return {__result, __abserr};

    return_error:

      __check_error(__func__, __error_type, __result, __abserr);
      __throw_integration_error("qawf_integrate: Unknown error.", UNKNOWN_ERROR,
				__result, __abserr);
    }

} // namespace __gn_test

#endif // QAWF_INTEGRATE_TCC
