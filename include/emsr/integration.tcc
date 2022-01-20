//
// Copyright (C) 2011-2020 Free Software Foundation, Inc.
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

#ifndef INTEGRATION_TCC
#define INTEGRATION_TCC 1

#include <type_traits>

namespace emsr
{

  template<typename Tp>
    constexpr
    error_tolerance_t<Tp>::
    error_tolerance_t(Tp max_abs_err, Tp max_rel_err,
			 unsigned min_num_passes)
    : m_max_abs_err{std::abs(max_abs_err)},
      m_max_rel_err{std::abs(max_rel_err)},
      m_min_num_passes(min_num_passes == 0 ? 1 : min_num_passes),
      m_num_passes(0)
    {
      if (!s_valid_tolerances(this->m_max_abs_err, this->m_max_rel_err))
	{
	  std::ostringstream msg;
	  msg << "integral_error_t: Integration tolerance cannot be achieved"
		   " with given absolute (" << this->m_max_abs_err
		<< ") and relative (" << this->m_max_rel_err
		<< ") error limits.";
	  throw std::runtime_error(msg.str().c_str());
	}
    }

  /**
   * Integrate a smooth function from a to b.
   *
   * Higher-order rules converge more rapidly for most functions,
   * but may slow convergence for less well-behaved ones.
   *
   * @param func The function to be integrated.
   * @param lower The lower limit of integration.
   * @param upper The upper limit of integration.
   * @param max_abs_error The absolute error limit.
   * @param max_rel_error The relative error limit.
   * @param max_iter is the maximum number of iterations allowed
   * @param qkintrule is the Gauss-Kronrod integration rule.
   * @return A structure containing the integration result and the error.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate(FuncTp func,
	      Tp lower, Tp upper,
	      Tp max_abs_error,
	      Tp max_rel_error,
	      std::size_t max_iter,
	      Kronrod_Rule qkintrule)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(lower) || std::isnan(upper)
          || std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}, absarea_t{}};
      else
	{
          integration_workspace<Tp, std::invoke_result_t<FuncTp, Tp>>
	    workspace(max_iter);
          return qag_integrate(workspace, func,
			       lower, upper,
			       max_abs_error, max_rel_error,
			       gauss_kronrod_integral<Tp>(qkintrule));
	}
    }

  /**
   * Integrates a smooth function from -infinity to +infinity.
   *
   * @param func The function to be integrated.
   * @param max_abs_error The absolute error limit.
   * @param max_rel_error The relative error limit.
   * @param max_iter is the maximum number of iterations allowed
   * @param qkintrule is the Gauss-Kronrod integration rule.
   * @return A structure containing the integration result and the error.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_minf_pinf(FuncTp func,
			Tp max_abs_error,
			Tp max_rel_error,
			std::size_t max_iter,
			Kronrod_Rule qkintrule)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else
	{
          integration_workspace<Tp, std::invoke_result_t<FuncTp, Tp>>
	    workspace(max_iter);
          return qag_integrate(workspace,
			       map_minf_pinf<Tp, FuncTp>(func),
			       Tp{0}, Tp{1},
			       max_abs_error, max_rel_error,
			       gauss_kronrod_integral<Tp>(qkintrule));
	}
    }

  /**
   * Integrate a smooth function from -infinity to finite b.
   *
   * @param func The function to be integrated.
   * @param upper The upper limit of integration.
   * @param max_abs_error The absolute error limit.
   * @param max_rel_error The relative error limit.
   * @param max_iter is the maximum number of iterations allowed
   * @param qkintrule is the Gauss-Kronrod integration rule.
   * @return A structure containing the integration result and the error.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_minf_upper(FuncTp func,
			 Tp upper,
			 Tp max_abs_error,
			 Tp max_rel_error,
			 std::size_t max_iter,
			 Kronrod_Rule qkintrule)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(upper)
          || std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else
	{
          integration_workspace<Tp, std::invoke_result_t<FuncTp, Tp>>
	    workspace(max_iter);
          return qag_integrate(workspace,
			       map_minf_b(func, upper),
			       Tp{0}, Tp{1},
			       max_abs_error, max_rel_error,
			       gauss_kronrod_integral<Tp>(qkintrule));
	}
    }

  /**
   * Integrate a smooth function from finite lower limit to +infinity.
   *
   * @param func The function to be integrated.
   * @param upper The upper limit of integration.
   * @param max_abs_error The absolute error limit.
   * @param max_rel_error The relative error limit.
   * @param max_iter is the maximum number of iterations allowed
   * @param qkintrule is the Gauss-Kronrod integration rule.
   * @return A structure containing the integration result and the error.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_lower_pinf(FuncTp func,
			 Tp lower,
			 Tp max_abs_error,
			 Tp max_rel_error,
			 std::size_t max_iter,
			 Kronrod_Rule qkintrule)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(lower)
          || std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else
	{
          integration_workspace<Tp, std::invoke_result_t<FuncTp, Tp>>
	    workspace(max_iter);
          return qag_integrate(workspace,
			       map_a_pinf(func, lower),
			       Tp{0}, Tp{1},
			       max_abs_error, max_rel_error,
			       gauss_kronrod_integral<Tp>(qkintrule));
	}
    }

  /**
   *  Adaptive Gauss-Kronrod integration optimized for
   *  discontinuous or singular functions
   *
   * @param func The function to be integrated.
   * @param lower The lower limit of integration.
   * @param upper The upper limit of integration.
   * @param max_abs_error The absolute error limit.
   * @param max_rel_error The relative error limit.
   * @param max_iter is the maximum number of iterations allowed
   * @return A structure containing the integration result and the error.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_kronrod_singular(FuncTp func,
			       Tp lower, Tp upper,
			       Tp max_abs_error,
			       Tp max_rel_error,
			       std::size_t max_iter)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(lower) || std::isnan(upper)
          || std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}, absarea_t{}};
      else
	{
          integration_workspace<Tp, std::invoke_result_t<FuncTp, Tp>>
	    workspace(max_iter);
          return qags_integrate(workspace, func, lower, upper,
			        max_abs_error, max_rel_error);
	}
    }

  /**
   * Integrate a potentially singular function from -infinity to +infinity
   *
   * @param func The function to be integrated.
   * @param max_abs_error The absolute error limit.
   * @param max_rel_error The relative error limit.
   * @param max_iter is the maximum number of iterations allowed
   * @return A structure containing the integration result and the error.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_singular_minf_pinf(FuncTp func,
				 Tp max_abs_error,
				 Tp max_rel_error,
				 std::size_t max_iter)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else
	{
          integration_workspace<Tp, std::invoke_result_t<FuncTp, Tp>>
	    workspace(max_iter);
          return qags_integrate(workspace,
			        map_minf_pinf<Tp, FuncTp>(func),
			        Tp{0}, Tp{1},
			        max_abs_error, max_rel_error);
	}
    }

  /**
   * Integrate a potentially singular function from -infinity to b
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_singular_minf_upper(FuncTp func, Tp upper,
				  Tp max_abs_error,
				  Tp max_rel_error,
				  std::size_t max_iter)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(upper)
          || std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else
	{
          integration_workspace<Tp, std::invoke_result_t<FuncTp, Tp>>
	    workspace(max_iter);
          return qags_integrate(workspace,
			        map_minf_b(func, upper),
			        Tp{0}, Tp{1},
			        max_abs_error, max_rel_error);
	}
    }

  /**
   * Integrates a potentially singular function from a to +infinity
   *
   * @param func The function to be integrated.
   * @param lower The lower limit of integration.
   * @param max_abs_error The absolute error limit.
   * @param max_rel_error The relative error limit.
   * @param max_iter is the maximum number of iterations allowed
   * @return A structure containing the integration result and the error.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_singular_lower_pinf(FuncTp func, Tp lower,
				  Tp max_abs_error,
				  Tp max_rel_error,
				  std::size_t max_iter)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(lower)
          || std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else
	{
          integration_workspace<Tp, std::invoke_result_t<FuncTp, Tp>>
	    workspace(max_iter);
          return qags_integrate(workspace,
			        map_a_pinf(func, lower),
			        Tp{0}, Tp{1},
			        max_abs_error, max_rel_error);
	}
    }

  /**
   * Integrate a potentially singular function.
   *
   * @param func The function to be integrated.
   * @param lower The lower limit of integration.
   * @param upper The upper limit of integration.
   * @param max_abs_error The absolute error limit.
   * @param max_rel_error The relative error limit.
   * @param max_iter is the maximum number of iterations allowed
   * @return A structure containing the integration result and the error.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_singular(FuncTp func,
		       Tp lower, Tp upper,
		       Tp max_abs_error,
		       Tp max_rel_error,
		       std::size_t max_iter)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using std::isnan; // To avoid ambiguous overload
      using RetTp = std::invoke_result_t<FuncTp, Tp>;

      using integ_t = adaptive_integral_t<Tp, RetTp>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      constexpr auto infty = std::numeric_limits<Tp>::infinity();

      if (std::isnan(lower) || std::isnan(upper)
          || std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}, absarea_t{}};
      else if (lower == -infty)
	{
	  if (upper == infty) // Integration from -inf to +inf
	    return integrate_minf_pinf(func, max_abs_error,
				      max_rel_error, max_iter);
	  else if (upper == -infty)
	    throw std::runtime_error("integrate_singular: "
				       "Attempt to integrate from -infinity"
				       " to -infinity");
	  else // Integration from -inf to finite value
	    return integrate_minf_upper(func, upper,
					max_abs_error,
					max_rel_error,
					max_iter);
	}
      else if (lower == infty)
	{
	  if (upper == infty)
	    throw std::runtime_error("integrate_singular: "
				       "Attempted to integrate from +infinity"
				       " to +infinity");
	  else if (upper == -infty) // Integration from +inf to -inf,
	    {		     // Call integrate_infinite() and flip sign
	      adaptive_integral_t<Tp, RetTp> res
		= integrate_minf_pinf(func, max_abs_error,
				      max_rel_error, max_iter);
	      return {-res.result, res.abserr};
	    }
	  else // Integration from +inf to finite value,
	    { // Call integrate_singular_to_infinity and flip sign
	      adaptive_integral_t<Tp, RetTp>
		res = integrate_lower_pinf(func, upper,
					     max_abs_error,
					     max_rel_error,
					     max_iter);
	      return {-res.result, res.abserr};
	    }
	}
      else // a is finite
	{
	  if (upper == infty)
	    return integrate_lower_pinf(func, lower,
					max_abs_error,
					max_rel_error, max_iter);
	  else if (upper == -infty)
	    { // Call integrate_from_infinity and flip sign
	      adaptive_integral_t<Tp, RetTp>
		res = integrate_minf_upper(func, lower,
					     max_abs_error,
					     max_rel_error,
					     max_iter);
	      return {-res.result, res.abserr};
	    }
	  else // Both a and b finite, call integrate_singular
	    return integrate_kronrod_singular(func, lower, upper,
					      max_abs_error, max_rel_error,
					      max_iter);
	}
    }

  /**
   * Integrate an oscillatory function.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_oscillatory(FuncTp func,
			  Tp lower, Tp upper,
			  Tp max_abs_error,
			  Tp max_rel_error,
			  std::size_t max_iter)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(lower) || std::isnan(upper)
          || std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}, absarea_t{}};
      else
	{
          using iwksp_t
	    = integration_workspace<Tp, std::invoke_result_t<FuncTp, Tp>>;
          using oitab_t = oscillatory_integration_table<Tp>;

          iwksp_t w(max_iter);
          oitab_t wo(upper, Tp{1}, oitab_t::INTEG_SINE, max_iter);
          return qawo_integrate(w, wo, func, lower,
			        max_abs_error, max_rel_error);
	}
    }

  /**
   * Adaptively integrate a function with known singular/discontinuous points.
   *
   * @tparam FuncTp     A function type that takes a single real scalar
   *                     argument and returns a real scalar.
   * @tparam Tp         A real type for the limits of integration and the step.
   */
  template<typename FuncTp, typename FwdIter, typename Tp>//, typename Integrator>
    auto
    integrate_multisingular(FuncTp func,
			    FwdIter ptbeg, FwdIter ptend,
			    Tp max_abs_error,
			    Tp max_rel_error,
			    std::size_t max_iter)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      //using pt_t = decltype(*ptbeg);
      //using ret_t = decltype(func(pt_t()));
      //using area_t = decltype(pt_t() * func(pt_t()));
      //using err_t = decltype(std::abs(area_t()));

      if (std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else
	{
          integration_workspace<Tp, std::invoke_result_t<FuncTp, Tp>>
	    wk(max_iter);
          return qagp_integrate(wk, func, std::vector(ptbeg, ptend),
			        max_abs_error, max_rel_error);
	}
    }

  /**
   * Integrate a function using an adaptive Clenshaw-Curtis algorithm.
   *
   * @tparam FuncTp     A function type that takes a single real scalar
   *                     argument and returns a real scalar.
   * @tparam Tp         A real type for the limits of integration and the step.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_clenshaw_curtis(FuncTp func,
			      Tp lower, Tp upper,
			      Tp max_abs_error,
			      Tp max_rel_error,
			      std::size_t max_iter)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(lower) || std::isnan(upper)
          || std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}, absarea_t{}};
      else
	{
          using RetTp = std::invoke_result_t<FuncTp, Tp>;

          cquad_workspace<Tp, RetTp> ws(max_iter);

          return cquad_integrate(ws, func, lower, upper,
			         max_abs_error, max_rel_error);
	}
    }

  /**
   * Adaptively integrate a function using a recursive Gauss-Kronrod quadrature
   * called the Patterson algorithm.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_patterson(FuncTp func,
			Tp lower, Tp upper,
			Tp max_abs_error,
			Tp max_rel_error)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(lower) || std::isnan(upper)
          || std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}, absarea_t{}};
      else
	{
          return qng_integrate(func, lower, upper,
			       max_abs_error, max_rel_error);
	}
    }

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
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_singular_endpoints(FuncTp func,
				 Tp lower, Tp upper,
				 Tp alpha, Tp beta,
				 int mu, int nu,
				 Tp max_abs_error,
				 Tp max_rel_error,
				 std::size_t max_iter)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(lower) || std::isnan(upper)
          || std::isnan(alpha) || std::isnan(beta)
          || std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}, absarea_t{}};
      else
	{
          integration_workspace<Tp, std::invoke_result_t<FuncTp, Tp>>
	    wksp(max_iter);
          qaws_integration_table<Tp> tab(alpha, beta, mu, nu);

          return qaws_integrate(wksp, tab, func, lower, upper,
			        max_abs_error, max_rel_error);
	}
    }

  /**
   * Integrate the principal value of a function with a Cauchy singularity.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_cauchy_principal_value(FuncTp func,
				     Tp lower, Tp upper, Tp center,
				     Tp max_abs_error, Tp max_rel_error,
				     std::size_t max_iter)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      if (std::isnan(lower) || std::isnan(upper) || std::isnan(center)
          || std::isnan(max_abs_error) || std::isnan(max_rel_error))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}, absarea_t{}};
      else
	{
          integration_workspace<Tp, std::invoke_result_t<FuncTp, Tp>>
	    wksp(max_iter);

          return  qawc_integrate(wksp, func, lower, upper, center,
			         max_abs_error, max_rel_error);
	}
    }

} // namespace emsr

#include <emsr/gauss_kronrod_integral.tcc>
#include <emsr/qag_integrate.tcc>
#include <emsr/qags_integrate.tcc>
#include <emsr/qng_integrate.tcc>
#include <emsr/qagp_integrate.tcc>
#include <emsr/qcheb_integrate.tcc>
#include <emsr/qawc_integrate.tcc>
#include <emsr/qaws_integrate.tcc>
#include <emsr/qawo_integrate.tcc>
#include <emsr/qawf_integrate.tcc>
#include <emsr/glfixed_integrate.tcc>
#include <emsr/cquad_integrate.tcc>

#include <emsr/gauss_quadrature.h>

#endif // INTEGRATION_TCC
