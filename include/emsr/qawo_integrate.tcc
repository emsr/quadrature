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
// Implements qawo integration algorithm
// Based on gsl/integration/qawo.c

#ifndef QAWO_INTEGRATE_TCC
#define QAWO_INTEGRATE_TCC 1

#include <stdexcept>
#include <type_traits>
#include <cmath>

#include <emsr/oscillatory_integration_table.h>

namespace emsr
{

  template<typename Tp, typename FuncTp>
    auto
    qc25f(oscillatory_integration_table<Tp>& wf,
	  FuncTp func, Tp lower, Tp upper,
	  std::size_t depth)
    -> gauss_kronrod_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

  template<typename Tp, typename FuncTp>
    auto
    qawo_integrate(integration_workspace<Tp,
			std::invoke_result_t<FuncTp, Tp>>& workspace,
		   oscillatory_integration_table<Tp>& wf,
		   FuncTp func,
		   const Tp lower,
		   const Tp max_abs_err, const Tp max_rel_err)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using AreaTp = std::invoke_result_t<FuncTp, Tp>;
      using AbsAreaTp = decltype(std::abs(AreaTp{}));

      const auto s_max = std::numeric_limits<Tp>::max();
      const auto s_eps = std::numeric_limits<Tp>::epsilon();
      const auto limit = workspace.capacity();
      // Try to adjust tests for varing precision. 10^{-5.3} for double.
      const auto s_rel_err = std::pow(Tp{10},
				 -std::numeric_limits<Tp>::digits / Tp{10});

      bool extrapolate = false;
      bool extall = false;
      bool allow_extrapolation = true;

      workspace.clear();
      //wf.clear();
      extrapolation_table<AreaTp, AbsAreaTp> table;

      const auto upper = lower + wf.get_length();
      const auto abs_omega = std::abs(wf.omega);

      auto result = Tp{0};
      auto abserr = Tp{0};

      if (!valid_tolerances(max_abs_err, max_rel_err))
	{
	  std::ostringstream msg;
	  msg << "qawo_integrate: Tolerance cannot be achieved with given "
		   "absolute (" << max_abs_err << ") and relative ("
		<< max_rel_err << ") error limits.";
	  throw std::runtime_error(msg.str().c_str());
	}

      // Perform the first integration.
      auto [result0, abserr0, resabs0, resasc0]
	= qc25f(wf, func, lower, upper, 0);

      workspace.append(lower, upper, result0, abserr0);

      auto tolerance = std::max(max_abs_err,
				  max_rel_err * std::abs(result0));

      if (abserr0 <= Tp{100} * s_eps * resabs0
	  && abserr0 > tolerance)
	throw integration_error("qawo_integrate: "
				"Cannot reach tolerance because "
				"of roundoff error on first attempt",
				ROUNDOFF_ERROR, result0, abserr0);
      else if ((abserr0 <= tolerance && abserr0 != resasc0)
		|| abserr0 == Tp{0})
	return {result0, abserr0};
      else if (limit == 1)
	throw integration_error("qawo_integrate: "
				"A maximum of one iteration was insufficient",
				MAX_ITER_ERROR, result0, abserr0);

      if (0.5 * abs_omega * std::abs(upper - lower) <= Tp{2})
	{
	  table.append(result0);
	  extall = true;
	}

      auto res_ext = result0;
      auto err_ext = s_max;

      auto area = result0;
      auto errsum = abserr0;
      auto iteration = 1u;
      auto ktmin = 0u;
      auto ertest = Tp{0};
      auto error_over_large_intervals = Tp{0};
      auto correc = Tp{0};
      int roundoff_type1 = 0, roundoff_type2 = 0, roundoff_type3 = 0;
      int error_type = NO_ERROR, error_type2 = NO_ERROR;
      do
	{
	  // Bisect the subinterval with the largest error estimate.
	  const auto& curr = workspace.retrieve();

	  const auto current_depth = workspace.curr_depth() + 1;

	  if (current_depth >= wf.n)
	    {
	      error_type = -1; // exceeded limit of table.
	      break;
	    }

	  const auto a1 = curr.lower_lim;
	  const auto mid = (curr.lower_lim + curr.upper_lim) / Tp{2};
	  const auto a2 = mid;
	  const auto b2 = curr.upper_lim;

	  ++iteration;

	  auto [area1, error1, resabs1, resasc1]
	    = qc25f(wf, func, a1, mid, current_depth);

	  auto [area2, error2, resabs2, resasc2]
	    = qc25f(wf, func, a2, b2, current_depth);

	  const auto area12 = area1 + area2;
	  const auto error12 = error1 + error2;
	  const auto last_e_i = curr.abs_error;

	  // Improve previous approximations to the integral and test
	  // for accuracy.
	  area += area12 - curr.result;
	  errsum += error12 - curr.abs_error;

	  tolerance = std::max(max_abs_err,
				 max_rel_err * std::abs(area));

	  if (resasc1 != error1 && resasc2 != error2)
	    {
	      const auto delta = curr.result - area12;

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
	      result = workspace.total_integral();
	      abserr = errsum;
	      check_error(__func__, error_type, result, abserr);
	      return {result, abserr};
	    }

	  if (error_type != NO_ERROR)
	    break;

	  if (iteration >= limit - 1)
	    {
	      error_type = MAX_ITER_ERROR;
	      break;
	    }

	  // Set up variables on first iteration.
	  if (iteration == 2 && extall)
	    {
	      error_over_large_intervals = errsum;
	      ertest = tolerance;
	      table.append(area);
	      continue;
	    }

	  if (!allow_extrapolation)
	    continue;

	  if (extall)
	    {
	      error_over_large_intervals -= last_e_i;

	      if (current_depth < workspace.max_depth())
		error_over_large_intervals += error12;

	      if (extrapolate)
		if (error_type2 != NO_ERROR
		 && error_over_large_intervals > ertest)
		  if (workspace.increment_curr_index())
		    continue;
	    }

	  if (workspace.large_interval())
	    continue;

	  if (extall)
	    {
	      extrapolate = true;
	      workspace.increment_curr_index();
	    }
	  else
	    {
	      // Test whether the interval to be bisected next is the
	      // smallest interval.
	      auto width = workspace.upper_lim() - workspace.lower_lim();

	      if (0.25 * std::abs(width) * abs_omega > Tp{2})
		continue;

	      extall = true;
	      error_over_large_intervals = errsum;
	      ertest = tolerance;
	      continue;
	    }

	  if (error_type2 != NO_ERROR
	   && error_over_large_intervals > ertest)
	    if (workspace.increment_curr_index())
	      continue;

	  // Perform extrapolation.

	  table.append(area);

	  if (table.get_nn() < 3)
	    {
	      workspace.reset_curr_index();
	      extrapolate = false;
	      error_over_large_intervals = errsum;
	      continue;
	    }

	  Tp reseps, abseps;
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
	      ertest = std::max(max_abs_err, max_rel_err * std::abs(reseps));
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
      while (iteration < limit);

      result = res_ext;
      abserr = err_ext;

      if (err_ext == s_max)
	{
	  result = workspace.total_integral();
	  abserr = errsum;
	  check_error(__func__, error_type, result, abserr);
	  return {result, abserr};
	}

      if (error_type != NO_ERROR || error_type2 != NO_ERROR)
	{
	  if (error_type2 != NO_ERROR)
	    err_ext += correc;

	  if (error_type == NO_ERROR)
	    error_type = SINGULAR_ERROR;

	  if (result != 0 && area != 0)
	    {
	      if (err_ext / std::abs(res_ext) > errsum / std::abs(area))
		{
		  result = workspace.total_integral();
		  abserr = errsum;
		  check_error(__func__, error_type, result, abserr);
		  return {result, abserr};
		}
	    }
	  else if (err_ext > errsum)
	    {
	      result = workspace.total_integral();
	      abserr = errsum;
	      check_error(__func__, error_type, result, abserr);
	      return {result, abserr};
	    }
	  else if (area == Tp{0})
	    {
	      check_error(__func__, error_type, result, abserr);
	      return {result, abserr};
	    }
	}

      // Test on divergence.
      bool positive_integrand = test_positivity(result0, resabs0);
      auto max_area = std::max(std::abs(res_ext), std::abs(area));
      if (!positive_integrand && max_area < Tp{0.01} * resabs0)
	{
	  check_error(__func__, error_type, result, abserr);
	  return {result, abserr};
	}

      auto ratio = res_ext / area;
      if (ratio < Tp{0.01} || ratio > Tp{100}
	  || errsum > std::abs(area))
	error_type = UNKNOWN_ERROR;

      check_error(__func__, error_type, result, abserr);
      return {result, abserr};
    }


  template<typename Tp, typename FuncTp>
    auto
    qc25f(oscillatory_integration_table<Tp>& wf,
	  FuncTp func, Tp lower, Tp upper,
	  std::size_t depth)
    -> gauss_kronrod_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      const auto s_max = std::numeric_limits<Tp>::max();
      const auto center = ( lower + upper) / Tp{2};
      const auto half_length = (upper - lower) / Tp{2};
      const auto omega = wf.omega;

      const auto par = omega * half_length;

      if (std::abs(par) < Tp{2})
	{
	  if (wf.circfun == oscillatory_integration_table<Tp>::INTEG_SINE)
	    {
	      auto wfun = [func, omega](Tp x)
			    ->Tp
			    { return std::sin(omega * x) * func(x); };
	      return qk_integrate(wfun, lower, upper, Kronrod_15);
	    }
	  else
	    {
	      auto wfun = [func, omega](Tp x)
			    ->Tp
			    { return std::cos(omega * x) * func(x); };
	      return qk_integrate(wfun, lower, upper, Kronrod_15);
	    }
	}
      else
	{
	  auto chout = qcheb_integrate(func, lower, upper);
	  const auto& cheb12 = chout.cheb12;
	  const auto& cheb24 = chout.cheb24;

	  if (depth >= wf.n)
	    {
	      // Table overflow should not happen, check before calling.
	      throw std::runtime_error("Table overflow in internal function");
	    }

	  const auto& moment = wf.get_moments(depth);

	  auto res12_cos = cheb12[12] * moment[12];
	  auto res12_sin = Tp{0};
	  for (int i = 0; i < 6; ++i)
	    {
	      const std::size_t k = 10 - 2 * i;
	      res12_cos += cheb12[k] * moment[k];
	      res12_sin += cheb12[k + 1] * moment[k + 1];
	    }

	  auto res24_cos = cheb24[24] * moment[24];
	  auto res24_sin = Tp{0};
	  auto result_abs = std::abs(cheb24[24]);
	  for (int i = 0; i < 12; ++i)
	    {
	      const std::size_t k = 22 - 2 * i;
	      res24_cos += cheb24[k] * moment[k];
	      res24_sin += cheb24[k + 1] * moment[k + 1];
	      result_abs += std::abs(cheb24[k])
			  + std::abs(cheb24[k + 1]);
	    }

	  const auto est_cos = std::abs(res24_cos - res12_cos);
	  const auto est_sin = std::abs(res24_sin - res12_sin);

	  const auto c = half_length * std::cos(center * omega);
	  const auto s = half_length * std::sin(center * omega);

	  Tp result, abserr, resabs, resasc;
	  if (wf.circfun == oscillatory_integration_table<Tp>::INTEG_SINE)
	    {
	      result = c * res24_sin + s * res24_cos;
	      abserr = std::abs(c * est_sin) + std::abs(s * est_cos);
	    }
	  else
	    {
	      result = c * res24_cos - s * res24_sin;
	      abserr = std::abs(c * est_cos) + std::abs(s * est_sin);
	    }

	  resabs = result_abs * half_length;
	  resasc = s_max;

	  return {result, abserr, resabs, resasc};
	}
    }

} // namespace emsr

#endif // QAWO_INTEGRATE_TCC

