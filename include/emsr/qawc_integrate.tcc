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
// Implements qawc integration algorithm
// Based on gsl/integration/qawc.c

#ifndef QAWC_INTEGRATE_TCC
#define QAWC_INTEGRATE_TCC 1

#include <stdexcept>
#include <type_traits>
#include <array>

namespace emsr
{

  template<typename Tp, typename FuncTp,
	   typename Integrator = gauss_kronrod_integral<Tp>>
    auto
    qc25c(FuncTp func, Tp lower, Tp upper, Tp center,
	  Integrator quad = gauss_kronrod_integral<Tp>(Kronrod_15))
    -> std::tuple<decltype(Tp{} * func(Tp{})), Tp, bool>;

  template<typename Tp>
    std::vector<Tp>
    compute_moments(std::size_t N, Tp cc);

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
   * a user-supplied integration rule (default 15-point Gauss-Kronrod). 
   */
  template<typename Tp, typename FuncTp,
	   typename Integrator = gauss_kronrod_integral<Tp>>
    auto
    qawc_integrate(integration_workspace<Tp,
			std::invoke_result_t<FuncTp, Tp>>& workspace,
		   FuncTp func,
		   Tp lower, Tp upper, Tp center,
		   Tp max_abs_err, Tp max_rel_err,
		   Integrator quad = gauss_kronrod_integral<Tp>(Kronrod_15))
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using RetTp = std::invoke_result_t<FuncTp, Tp>;
      using AreaTp = decltype(RetTp{} * Tp{});
      using AbsAreaTp = decltype(std::abs(AreaTp{}));

      auto result = AreaTp{};
      auto abserr = AbsAreaTp{};

      const auto limit = workspace.capacity();
      // Try to adjust tests for varing precision.
      const auto m_rel_err = std::pow(Tp{10.0},
				 -std::numeric_limits<Tp>::digits / Tp{10.0});

      int sign = 1;
      if (upper < lower)
	{
	  std::swap(lower, upper);
	  sign = -1;
	}

      if (!valid_tolerances(max_abs_err, max_rel_err))
	{
	  std::ostringstream msg;
	  msg << "qawc_integrate: Tolerance cannot be achieved with given "
		   "absolute (" << max_abs_err << ") and relative ("
		<< max_rel_err << ") error limits.";
	  throw std::runtime_error(msg.str().c_str());
	}

      if (center == lower || center == upper)
	throw std::runtime_error ("qawc_integrate: "
				    "Cannot integrate with singularity "
				    "on endpoint.");

      workspace.clear();

      // Perform the first integration.
      auto [result0, abserr0, err_reliable]
	= qc25c(func, lower, upper, center, quad);

      workspace.append(lower, upper, result0, abserr0);

      // Test on accuracy; Use 0.01 relative error as an extra safety
      // margin on the first iteration (ignored for subsequent iterations).
      auto tolerance = std::max(max_abs_err,
				  max_rel_err * std::abs( result0));
      if (abserr0 < tolerance && abserr0
	  < Tp{0.01} * std::abs(result0))
	return {sign * result0, abserr0};
      else if (limit == 1)
	throw integration_error("qawc_integrate: "
				"A maximum of one iteration was insufficient",
				MAX_ITER_ERROR, sign * result0, abserr0);

      auto area = result0;
      auto errsum = abserr0;
      auto iteration = 1u;
      int error_type = NO_ERROR;
      int roundoff_type1 = 0, roundoff_type2 = 0;
      do
	{
	  // Bisect the subinterval with the largest error estimate.
	  const auto& curr = workspace.retrieve();

	  const auto a1 = curr.lower_lim;
	  const auto b2 = curr.upper_lim;
	  auto mid = (curr.lower_lim + curr.upper_lim) / Tp{2};
	  if (center > a1 && center <= mid)
	    mid = (center + b2) / Tp{2};
	  else if (center > mid && center < b2)
	    mid = (a1 + center) / Tp{2};
	  const auto a2 = mid;

	  auto [area1, error1, err_reliable1]
	    = qc25c(func, a1, mid, center, quad);

	  auto [area2, error2, err_reliable2]
	    = qc25c(func, a2, b2, center, quad);

	  const auto area12 = area1 + area2;
	  const auto error12 = error1 + error2;

	  errsum += error12 - curr.abs_error;
	  area += area12 - curr.result;

	  if (err_reliable1 && err_reliable2)
	    {
	      const auto delta = curr.result - area12;

	      if (std::abs (delta) <= m_rel_err * std::abs (area12)
	    	   && error12 >= 0.99 * curr.abs_error)
		++roundoff_type1;
	      if (iteration >= 10 && error12 > curr.abs_error)
		++roundoff_type2;
	    }

	  tolerance = std::max(max_abs_err, max_rel_err * std::abs(area));
	  if (errsum > tolerance)
	    {
	      if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
		error_type = ROUNDOFF_ERROR;

	      // set error flag in the case of bad integrand behaviour at
	      // a point of the integration range
	      if (workspace.subinterval_too_small(a1, a2, b2))
		error_type = SINGULAR_ERROR;
	    }

	  workspace.split(mid, area1, error1, area2, error2);

	  ++iteration;
	}
      while (iteration < limit && !error_type && errsum > tolerance);

      result = sign * workspace.total_integral();
      abserr = errsum;

      if (iteration == limit)
	error_type = MAX_SUBDIV_ERROR;

      if (errsum <= tolerance)
	return {result, abserr};

      if (error_type == NO_ERROR)
	return {result, abserr};

      check_error(__func__, error_type, result, abserr);
      throw integration_error("qawc_integrate: Unknown error.",
			      UNKNOWN_ERROR, result, abserr);
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
   * a user-supplied integration rule. 
   */
  template<typename Tp, typename FuncTp,
	   typename Integrator = gauss_kronrod_integral<Tp>>
    auto
    qc25c(FuncTp func, Tp lower, Tp upper, Tp center,
	  Integrator quad)
    -> std::tuple<decltype(Tp{} * func(Tp{})), Tp, bool>
    {
      using RetTp = std::invoke_result_t<FuncTp, Tp>;
      using AreaTp = decltype(RetTp{} * Tp{});
      using AbsAreaTp = decltype(std::abs(AreaTp{}));
      bool err_reliable;

      auto result = AreaTp{};
      auto abserr = AbsAreaTp{};

      const auto cc = (Tp{2} * center - upper - lower)
		      / (upper - lower);

      if (std::abs(cc) > Tp{1.1})
	{
	  auto func_cauchy = [func, center](Tp x)
				-> Tp
				{ return func(x) / (x - center); };

	  auto [result, abserr, resabs, resasc]
	    = quad(func_cauchy, lower, upper);

	  if (abserr == resasc)
	    err_reliable = false;
	  else
	    err_reliable = true;

	  return std::make_tuple(result, abserr, err_reliable);
	}
      else
	{
	  auto [cheb12, cheb24] = qcheb_integrate(func, lower, upper);
	  const auto moment = compute_moments(cheb24.size(), cc);

	  auto res12 = AreaTp{0};
	  for (size_t i = 0u; i < cheb12.size(); ++i)
	    res12 += cheb12[i] * moment[i];

	  auto res24 = AreaTp{0};
	  for (size_t i = 0u; i < cheb24.size(); ++i)
	    res24 += cheb24[i] * moment[i];

	  result = res24;
	  abserr = std::abs(res24 - res12);
	  err_reliable = false;

	  return std::make_tuple(result, abserr, err_reliable);
	}
    }

  /**
   * Compute Clenshaw-Curtis moments.
   * An iterator range version would be nicer I think.
   */
  template<typename Tp>
    std::vector<Tp>
    compute_moments(std::size_t N, Tp cc)
    {
      std::vector<Tp> moment(N);

      auto a0 = std::log(std::abs((Tp{1} - cc) / (Tp{1} + cc)));
      auto a1 = Tp{2} + a0 * cc;

      moment[0] = a0;
      moment[1] = a1;

      for (size_t k = 2; k < N; ++k)
	{
	  Tp a2;

	  if ((k % 2) == 0)
	    a2 = Tp{2} * cc * a1 - a0;
	  else
	    {
	      const auto km1 = Tp(k - 1);
	      a2 = Tp{2} * cc * a1
		   - a0
		   - Tp{4} / (km1 * km1 - Tp{1});
	    }

	  moment[k] = a2;

	  a0 = a1;
	  a1 = a2;
	}

      return moment;
    }

} // namespace emsr

#endif // QAWC_INTEGRATE_TCC
