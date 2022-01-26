//
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
// Copyright (C) 2016-2020 Free Software Foundation, Inc.
// Copyright (C) 2021-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
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

#include <stdexcept>
#include <type_traits>
#include <cmath>

#include <emsr/integration_workspace.h>
#include <emsr/oscillatory_integration_table.h>

namespace emsr
{

  /**
   * This function attempts to compute a Fourier integral of the function f
   * over the semi-infinite interval [a,+\infty)
   */
  template<typename Tp, typename FuncTp>
    auto
    qawf_integrate(integration_workspace<Tp,
			std::invoke_result_t<FuncTp, Tp>>& workspace,
		   integration_workspace<Tp,
			std::invoke_result_t<FuncTp, Tp>>& cycle_workspace,
		   oscillatory_integration_table<Tp>& wf,
		   FuncTp func,
		   Tp lower, Tp max_abs_err)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using AreaTp = std::invoke_result_t<FuncTp, Tp>;
      using AbsAreaTp = decltype(std::abs(AreaTp{}));

      std::size_t ktmin = 0;
      std::size_t iteration = 0;

      extrapolation_table<AreaTp, AbsAreaTp> table;

      Tp cycle;
      auto omega = wf.omega;

      const Tp p = 0.9;
      Tp factor = 1;
      Tp initial_eps, eps;
      int error_type = NO_ERROR;

      int status = 0;
      auto result = Tp{0};
      auto abserr = Tp{0};

      const auto s_max = std::numeric_limits<Tp>::max();
      const auto limit = workspace.capacity();

      workspace.clear();
      cycle_workspace.clear();
      //wf.clear();

      // Test on accuracy.
      if (max_abs_err <= Tp{0})
	throw std::domain_error("absolute tolerance epsabs must be positive") ;

      if (omega == Tp{0})
	{
	  if (wf.circfun == oscillatory_integration_table<Tp>::INTEG_SINE)
	    // The function sin(w x) f(x) is always zero for w = 0.
	    return {Tp{0}, Tp{0}};
	  else
	    // The function cos(w x) f(x) is always f(x) for w = 0.
	    return qagiu_integrate(cycle_workspace, func, lower, max_abs_err,
				   Tp{0});
	}

      if (max_abs_err * (Tp{1} - p) > std::numeric_limits<Tp>::min())
	eps = max_abs_err * (Tp{1} - p);
      else
	eps = max_abs_err;

      initial_eps = eps;

      auto area = Tp{0};
      auto errsum = Tp{0};
      auto res_ext = Tp{0};
      auto err_ext = s_max;
      auto correc = Tp{0};
      auto total_error = Tp{0};
      auto truncation_error = Tp{0};

      cycle = (2 * std::floor(std::abs(omega)) + 1)
	      * M_PI / std::abs(omega);

      wf.set_length(cycle);

      for (iteration = 0; iteration < limit; ++iteration)
	{
	  const auto a1 = lower + iteration * cycle;
	  const auto b1 = a1 + cycle;

	  const auto max_abs_err1 = eps * factor;

	  auto out1 = qawo_integrate(cycle_workspace, wf, func, a1,
			     max_abs_err1, Tp{0});
	  auto area1 = out1.result;
	  auto error1 = out1.abserr;

	  workspace.append(a1, b1, area1, error1);

	  factor *= p;

	  area += area1;
	  errsum += error1;

	  // Estimate the truncation error as 50 times the final term.
	  truncation_error = 50 * std::abs(area1);
	  total_error = errsum + truncation_error;
	  if (total_error < max_abs_err && iteration > 4)
	    goto compute_result;

	  if (error1 > correc)
	    correc = error1;

	  if (status)
	    eps = std::max(initial_eps, correc * (Tp{1} - p));

	  if (status && total_error < 10 * correc && iteration > 3)
	    goto compute_result;

	  table.append(area);

	  if (table.get_nn() < 2)
	    continue;

	  Tp reseps, erreps;
	  std::tie(reseps, erreps) = table.qelg();

	  ++ktmin;
	  if (ktmin >= 15 && err_ext < Tp{0.001L} * total_error)
	    error_type = EXTRAP_ROUNDOFF_ERROR;

	  if (erreps < err_ext)
	    {
	      ktmin = 0;
	      err_ext = erreps;
	      res_ext = reseps;

	      if (err_ext + 10 * correc <= max_abs_err)
		break;
	      if (err_ext <= max_abs_err && 10 * correc >= max_abs_err)
		break;
	    }
	}

      if (iteration == limit)
	error_type = MAX_ITER_ERROR;

      if (err_ext == s_max)
	goto compute_result;

      err_ext += 10 * correc;

      result = res_ext;
      abserr = err_ext;

      if (error_type == NO_ERROR)
	return {result, abserr};

      if (res_ext != Tp{0} && area != Tp{0})
	{
	  if (err_ext / std::abs(res_ext) > errsum / std::abs(area))
	    goto compute_result;
	}
      else if (err_ext > errsum)
	goto compute_result;
      else if (area == Tp{0})
	goto return_error;

      if (error_type == EXTRAP_ROUNDOFF_ERROR)
	err_ext += truncation_error;

      goto return_error;

    compute_result:

      result = area;
      abserr = total_error;

      if (error_type == NO_ERROR)
	return {result, abserr};

    return_error:

      check_error(__func__, error_type, result, abserr);
      throw integration_error("qawf_integrate: Unknown error.", UNKNOWN_ERROR,
			      result, abserr);
    }

} // namespace emsr

#endif // QAWF_INTEGRATE_TCC
