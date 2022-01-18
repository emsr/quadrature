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
// Implements QAWS integration algorithm
// Based on gsl/integration/qaws.c

#ifndef QAWS_INTEGRATE_TCC
#define QAWS_INTEGRATE_TCC 1

#include <type_traits>
#include <array>

#include <ext/integration_workspace.h>
#include <ext/qaws_integration_table.h>

namespace emsr
{

  template<typename Tp, typename FuncTp>
    struct fn_qaws;

  template<typename AreaTp>
    struct compute_result_t
    {
      AreaTp res12;
      AreaTp res24;
    };

  template<typename Tp, typename RetTp>
    auto
    compute_result(const std::array<Tp, 25>& r,
		   const std::array<RetTp, 13>& cheb12,
		   const std::array<RetTp, 25>& cheb24)
    -> compute_result_t<decltype(Tp{} * RetTp{})>;

 template<typename Tp, typename FuncTp,
	  typename Integrator = gauss_kronrod_integral<Tp>>
    std::tuple<Tp, Tp, bool>
    qc25s(qaws_integration_table<Tp>& t,
	  FuncTp func, Tp lower, Tp upper, Tp a1, Tp mid,
	  Integrator quad = gauss_kronrod_integral<Tp>(Kronrod_15));

  /**
   * The singular weight function is defined by:
   * @f[
   *    W(x) = (x-a)^\alpha (b-x)^\beta log^\mu (x-a) log^\nu (b-x)
   * @f]
   * where @f$ \alpha > -1 @f$, @f$ \beta > -1 @f$,
   * and @f$ \mu = 0 @f$, 1, @f$ \nu = 0, 1 @f$.
   *
   * The weight function can take four different forms depending
   * on the values of \mu and \nu,
   * @f[
   *    W(x) = (x-a)^\alpha (b-x)^\beta                   (\mu = 0, \nu = 0)
   *    W(x) = (x-a)^\alpha (b-x)^\beta log(x-a)          (\mu = 1, \nu = 0)
   *    W(x) = (x-a)^\alpha (b-x)^\beta log(b-x)          (\mu = 0, \nu = 1)
   *    W(x) = (x-a)^\alpha (b-x)^\beta log(x-a) log(b-x) (\mu = 1, \nu = 1)
   * @f]
   *
   * The QAWS algorithm is designed for integrands with algebraic-logarithmic
   * singularities at the end-points of an integration region.
   *
   * In order to work efficiently the algorithm requires a precomputed table
   * of Chebyshev moments.
   */
  template<typename Tp, typename FuncTp,
	   typename Integrator = gauss_kronrod_integral<Tp>>
    auto
    qaws_integrate(integration_workspace<Tp,
			std::invoke_result_t<FuncTp, Tp>>& workspace,
		   qaws_integration_table<Tp>& table,
		   FuncTp func,
		   Tp lower, Tp upper,
		   Tp max_abs_err, Tp max_rel_err,
		   Integrator quad = gauss_kronrod_integral<Tp>(Kronrod_15))
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      // Try to adjust tests for varing precision.
      const auto m_rel_err = std::pow(Tp{10},
				 -std::numeric_limits<Tp>::digits / Tp{10});

      if (upper <= lower)
	throw std::runtime_error("qaws_integrate: "
				   "Limits must form an ascending sequence");
      if (!valid_tolerances(max_abs_err, max_rel_err))
	{
	  std::ostringstream msg;
	  msg << "qaws_integrate: Tolerance cannot be achieved with given "
		   "absolute (" << max_abs_err << ") and relative ("
		<< max_rel_err << ") error limits.";
	  throw std::runtime_error(msg.str().c_str());
	}

      const auto limit = workspace.capacity();
      workspace.clear();

      // Perform the first integration.
      Tp result0, abserr0;
      {
	const auto a1 = lower;
	const auto mid = (lower + upper) / Tp{2};
	const auto b2 = upper;

        auto [area1, error1, err_reliable1]
	  = qc25s(table, func, lower, upper, a1, mid, quad);
	workspace.append(a1, mid, area1, error1);

	auto [area2, error2, err_reliable2]
	  = qc25s(table, func, lower, upper, mid, b2, quad);
	workspace.append(mid, b2, area2, error2);

	result0 = area1 + area2;
	abserr0 = error1 + error2;
      }

      // Test on accuracy; Use 0.01 relative error as an extra safety
      // margin on the first iteration (ignored for subsequent iterations).
      auto tolerance = std::max(max_abs_err, max_rel_err * std::abs(result0));
      if (abserr0 < tolerance && abserr0
	  < Tp{0.01} * std::abs(result0))
	return {result0, abserr0};
      else if (limit == 1)
	throw integration_error("qaws_integrate: "
				"A maximum of one iteration was insufficient",
				MAX_ITER_ERROR, result0, abserr0);

      auto area = result0;
      auto errsum = abserr0;
      auto iteration = 2u;
      int error_type = NO_ERROR;
      int roundoff_type1 = 0, roundoff_type2 = 0;
      do
	{
	  // Bisect the subinterval with the largest error estimate.
	  const auto& curr = workspace.retrieve();

	  const auto a1 = curr.lower_lim;
	  const auto mid = (curr.lower_lim + curr.upper_lim) / Tp{2};
	  const auto b2 = curr.upper_lim;

	  auto [area1, error1, err_reliable1]
	    = qc25s(table, func, lower, upper, a1, mid, quad);

	  auto [area2, error2, err_reliable2]
	    = qc25s(table, func, lower, upper, mid, b2, quad);

	  const auto area12 = area1 + area2;
	  const auto error12 = error1 + error2;

	  errsum += error12 - curr.abs_error;
	  area += area12 - curr.result;

	  if (err_reliable1 && err_reliable2)
	    {
	      const auto delta = curr.result - area12;

	      if (std::abs (delta) <= m_rel_err * std::abs(area12)
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
	      if (workspace.subinterval_too_small(a1, mid, b2))
		error_type = SINGULAR_ERROR;
	    }

	  workspace.split(mid, area1, error1, area2, error2);

	  ++iteration;
	}
      while (iteration < limit && !error_type && errsum > tolerance);

      const auto result = workspace.total_integral();
      const auto abserr = errsum;

      if (iteration == limit)
	error_type = MAX_SUBDIV_ERROR;

      if (errsum <= tolerance)
	return {result, abserr};

      if (error_type == NO_ERROR)
	return {result, abserr};

      check_error(__func__, error_type, result, abserr);
      throw integration_error("qaws_integrate: Unknown error.",
			      UNKNOWN_ERROR, result, abserr);
    }

  /**
   *
   */
  template<typename Tp, typename FuncTp,
	   typename Integrator = gauss_kronrod_integral<Tp>>
    std::tuple<Tp, Tp, bool>
    qc25s(qaws_integration_table<Tp>& t,
	  FuncTp func, Tp lower, Tp upper, Tp a1, Tp b1,
	  Integrator quad)
    {
      fn_qaws<Tp, FuncTp> fqaws(&t, func, lower, upper);

      if (a1 == lower && (t.alpha != Tp{0} || t.mu != 0))
	{
	  const auto factor
	    = std::pow(0.5 * (b1 - a1), t.alpha + Tp{1});

	  auto f = [fqaws](Tp x)
		     -> Tp { return fqaws.eval_right(x); };
	  auto chout = qcheb_integrate(f, a1, b1);
	  const auto& cheb12 = chout.cheb12;
	  const auto& cheb24 = chout.cheb24;

	  if (t.mu == 0)
	    {
	      const auto u = factor;

	      auto [res12, res24]
		= compute_result(t.ri, cheb12, cheb24);

	      const auto result = u * res24;
	      const auto abserr = std::abs(u * (res24 - res12));
	      return std::make_tuple(result, abserr, false);
	    }
	  else
	    {
	      const auto u = factor * std::log(b1 - a1);
	      const auto v = factor;

	      auto [res12a, res24a]
		= compute_result(t.ri, cheb12, cheb24);
	      auto [res12b, res24b]
		= compute_result(t.rg, cheb12, cheb24);

	      const auto result = u * res24a + v * res24b;
	      const auto abserr = std::abs(u * (res24a - res12a))
				  + std::abs(v * (res24b - res12b));
	      return std::make_tuple(result, abserr, false);
	    }
	}
      else if (b1 == upper && (t.beta != Tp{0} || t.nu != 0))
	{
	  auto factor = std::pow(0.5 * (b1 - a1), t.beta + Tp{1});

	  auto f = [fqaws](Tp x)
		     -> Tp { return fqaws.eval_left(x); };
	  auto chout = qcheb_integrate(f, a1, b1);
	  const auto& cheb12 = chout.cheb12;
	  const auto& cheb24 = chout.cheb24;

	  if (t.nu == 0)
	    {
	      const auto u = factor;

	      auto [res12, res24]
		= compute_result(t.rj, cheb12, cheb24);

	      const auto result = u * res24;
	      const auto abserr = std::abs(u * (res24 - res12));
	      return std::make_tuple(result, abserr, false);
	    }
	  else
	    {
	      const auto u = factor * std::log(b1 - a1);
	      const auto v = factor;

	      auto [res12a, res24a]
		= compute_result(t.rj, cheb12, cheb24);
	      auto [res12b, res24b]
		= compute_result(t.rh, cheb12, cheb24);

	      const auto result = u * res24a + v * res24b;
	      const auto abserr = std::abs(u * (res24a - res12a))
				  + std::abs(v * (res24b - res12b));
	      return std::make_tuple(result, abserr, false);
	    }
	}
      else
	{
	  auto f = [fqaws](Tp x)
		     ->Tp
		     { return fqaws.eval_middle(x); };

	  auto [result, abserr, resabs, resasc]
	    = quad(f, a1, b1);

	  bool err_reliable;
	  if (abserr == resasc)
	    err_reliable = false;
	  else
	    err_reliable = true;

	  return std::make_tuple(result, abserr, err_reliable);
	}
    }

  /*
   *
   */
  template<typename Tp, typename FuncTp>
    struct fn_qaws
    {
      using RetTp = std::invoke_result_t<FuncTp, Tp>;
      using AreaTp = decltype(RetTp{} * Tp{});

      const qaws_integration_table<Tp>* table;
      FuncTp func;
      Tp a;
      Tp b;

      fn_qaws(const qaws_integration_table<Tp>* tab,
	      FuncTp func, Tp a_in, Tp b_in)
      : table(tab),
	func(func), a(a_in), b(b_in)
      { }

      RetTp eval_middle(Tp) const;
      RetTp eval_left(Tp) const;
      RetTp eval_right(Tp) const;
    };

  template<typename Tp, typename FuncTp>
    std::invoke_result_t<FuncTp, Tp>
    fn_qaws<Tp, FuncTp>::eval_middle(Tp x) const
    {
      auto factor = Tp{1};

      if (this->table->alpha != Tp{0})
	factor *= std::pow(x - this->a, this->table->alpha);

      if (table->mu == 1)
	factor *= std::log(x - this->a);

      if (this->table->beta != Tp{0})
	factor *= std::pow(this->b - x, this->table->beta);

      if (table->nu == 1)
	factor *= std::log(this->b - x);

      return factor * this->func(x);
    }

  template<typename Tp, typename FuncTp>
    std::invoke_result_t<FuncTp, Tp>
    fn_qaws<Tp, FuncTp>::eval_left(Tp x) const
    {
      auto factor = Tp{1};

      if (this->table->alpha != Tp{0})
	factor *= std::pow(x - this->a, this->table->alpha);

      if (this->table->mu == 1)
	factor *= std::log(x - this->a);

      return factor * this->func(x);
    }

  template<typename Tp, typename FuncTp>
    std::invoke_result_t<FuncTp, Tp>
    fn_qaws<Tp, FuncTp>::eval_right(Tp x) const
    {
      auto factor = Tp{1};

      if (this->table->beta != Tp{0})
	factor *= std::pow(this->b - x, this->table->beta);

      if (this->table->nu == 1)
	factor *= std::log(this->b - x);

      return factor * this->func(x);
    }

  template<typename Tp, typename RetTp>
    auto
    compute_result(const std::array<Tp, 25>& r,
		   const std::array<RetTp, 13>& cheb12,
		   const std::array<RetTp, 25>& cheb24)
    -> compute_result_t<decltype(Tp{} * RetTp{})>
    {
      using AreaTp = decltype(RetTp{} * Tp{});

      auto res12 = AreaTp{};
      for (size_t i = 0; i < cheb12.size(); ++i)
	res12 += r[i] * cheb12[i];

      auto res24 = AreaTp{};
      for (size_t i = 0; i < cheb24.size(); ++i)
	res24 += r[i] * cheb24[i];

      return compute_result_t<AreaTp>{res12, res24};
    }

} // namespace emsr

#endif // QAWS_INTEGRATE_TCC
