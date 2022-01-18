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


#ifndef INTEGRATION_H
#define INTEGRATION_H 1

#include <limits>
#include <tuple>
#include <complex> // For complex abs
#include <type_traits>

#include <ext/quadrature_point.h>
#include <ext/gauss_kronrod_integral.h>

namespace emsr
{

  /**
   * The return type for a fixed integral rule.
   * Fixed integral types could return errors too!
   * Globally vs. locally adaptive integrals need to be distinguished.
   * abs(diff) vs. diff(abs) is also a thing.
   */
  template<typename Tp, typename RetTp>
    struct fixed_integral_t
    {
      using AreaTp = decltype(RetTp{} * Tp{});

      /// Result of the integral.
      AreaTp result = AreaTp{};
    };

  /**
   * The return type for an adaptive integral rule.
   */
  template<typename Tp, typename RetTp>
    struct adaptive_integral_t
    {
      using AreaTp = decltype(RetTp{} * Tp{});
      using AbsAreaTp = decltype(std::abs(AreaTp{}));

      /// Result of the integral.
      AreaTp    result = AreaTp{};
      /// Absolute value of estimated error.
      AbsAreaTp abserr = AbsAreaTp{};
    };

} // namespace emsr

#include <ext/trapezoid_integral.h>
#include <ext/midpoint_integral.h>
#include <ext/simpson_integral.h>
#include <ext/gauss_quadrature.h>

namespace emsr
{

  /**
   * An error tolerance model for integrals.
   */
  template<typename Tp>
    struct error_tolerance_t
    {
      constexpr
      error_tolerance_t(Tp max_abs_err, Tp max_rel_err,
			unsigned min_num_passes);

      /// Maximum absolute error tolerance.
      Tp m_max_abs_err;
      /// Maximum relative error tolerance.
      Tp m_max_rel_err;
      /// Current tolerance.
      Tp m_tolerance = this->tolerance(Tp{1});
      /// Minimum number of consecutive passes before result is OK.
      unsigned m_min_num_passes;
      /// Current number of passes.
      unsigned m_num_passes;

      /// Set and return the tolerance given an initial result.
      template<typename ResTp>
	constexpr Tp
	tolerance(ResTp result)
	{
	  this->m_tolerance
		  = std::max(this->m_max_abs_err,
			     Tp(this->m_max_rel_err * std::abs(result)));
	  return this->m_tolerance;
	}

      /// Return whether convergence is acceptable tolerance
      /// given two subsequent results.
      template<typename ResTp>
	bool
	test(ResTp curr_result, ResTp prev_result)
	{
	  const auto del = std::abs(curr_result - prev_result);
	  if (del < this->m_max_abs_err
	      || del < this->m_max_rel_err * std::abs(curr_result))
	    ++this->m_num_passes;
	  else
	    this->m_num_passes = 0;
	  return this->m_num_passes >= this->m_min_num_passes;
	}

      /// Return the current tolerance.
      constexpr Tp
      tolerance() const
      { return this->m_tolerance; }

      /// Test for valid error tolerances.
      static constexpr bool
      s_valid_tolerances(Tp max_abs_err, Tp max_rel_err)
      {
	constexpr auto s_eps = std::numeric_limits<Tp>::epsilon();
	return !(max_abs_err <= Tp{0}
	      && (max_rel_err < Tp{50} * s_eps
	       || max_rel_err < 0.5e-28));
	// I don't understand the etymology of this last number.
      }
    };

  // Hack for now...
  template<typename Tp>
    inline bool
    valid_tolerances(Tp max_abs_err, Tp max_rel_err)
    {
      return error_tolerance_t<Tp>::
	     s_valid_tolerances(max_abs_err, max_rel_err);
    }

} // namespace emsr

namespace emsr
{

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
	      std::size_t max_iter = 1024,
	      Kronrod_Rule qkintrule = Kronrod_21)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

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
			std::size_t max_iter = 1024,
			Kronrod_Rule qkintrule = Kronrod_21)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

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
			 std::size_t max_iter = 1024,
			 Kronrod_Rule qkintrule = Kronrod_21)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

  template<typename Tp, typename FuncTp>
    inline auto
    integrate_lower_minf(FuncTp func,
			 Tp lower,
			 Tp max_abs_error,
			 Tp max_rel_error,
			 std::size_t max_iter = 1024,
			 Kronrod_Rule qkintrule = Kronrod_21)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      return -integrate_minf_upper(func, lower,
				   max_abs_error, max_rel_error,
				   max_iter, qkintrule);
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
			 std::size_t max_iter = 1024,
			 Kronrod_Rule qkintrule = Kronrod_21)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

  template<typename Tp, typename FuncTp>
    inline auto
    integrate_pinf_upper(FuncTp func,
			 Tp upper,
			 Tp max_abs_error,
			 Tp max_rel_error,
			 std::size_t max_iter = 1024,
			 Kronrod_Rule qkintrule = Kronrod_21)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      return -integrate_lower_pinf(func, upper,
				   max_abs_error, max_rel_error,
				   max_iter, qkintrule);
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
			       std::size_t max_iter = 1024)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

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
				 std::size_t max_iter = 1024)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

  /**
   * Integrate a potentially singular function from -infinity to b
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_singular_minf_upper(FuncTp func, Tp upper,
				  Tp max_abs_error,
				  Tp max_rel_error,
				  std::size_t max_iter = 1024)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

  template<typename Tp, typename FuncTp>
    inline auto
    integrate_singular_lower_minf(FuncTp func, Tp lower,
				  Tp max_abs_error,
				  Tp max_rel_error,
				  std::size_t max_iter = 1024)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      return -integrate_singular_minf_upper(func, lower,
					    max_abs_error, max_rel_error,
					    max_iter);
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
				  std::size_t max_iter = 1024)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

  template<typename Tp, typename FuncTp>
    inline auto
    integrate_singular_pinf_upper(FuncTp func, Tp upper,
				  Tp max_abs_error,
				  Tp max_rel_error,
				  std::size_t max_iter = 1024)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      return -integrate_singular_lower_pinf(func, upper,
					    max_abs_error, max_rel_error,
					    max_iter);
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
		       std::size_t max_iter = 1024)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

  /**
   * Integrate an oscillatory function.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_oscillatory(FuncTp func,
			  Tp lower, Tp upper,
			  Tp max_abs_error,
			  Tp max_rel_error,
			  std::size_t max_iter = 1024)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

  /**
   * Adaptively integrate a function with known singular/discontinuous points.
   *
   * @tparam FuncTp     A function type that takes a single real scalar
   *                     argument and returns a real scalar.
   * @tparam Tp         A real type for the limits of integration and the step.
   */
  template<typename FuncTp, typename FwdIter, typename Tp>
    auto
    integrate_multisingular(FuncTp func,
			    FwdIter ptbeg, FwdIter ptend,
			    Tp max_abs_error,
			    Tp max_rel_error,
			    std::size_t max_iter = 1024)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

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
			      std::size_t max_iter = 1024)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

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
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

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
				 std::size_t max_iter = 1024)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

  /**
   * Integrate the principal value of a function with a Cauchy singularity.
   */
  template<typename Tp, typename FuncTp>
    auto
    integrate_cauchy_principal_value(FuncTp func,
				     Tp lower, Tp upper, Tp center,
				     Tp max_abs_err, Tp max_rel_err,
				     std::size_t max_iter = 1024)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>;

  /**
   * @f[
   *    \int_{-1}^{+1}f(x)dx
   * @f]
   * Making the change of variables:
   * @f[
   *    x = tanh\left[\frac{\pi}{2}sinh(u)\right],
   *   dx = \frac{\pi}{2}
   *        \frac{cosh(u)}{cosh^2\left[\frac{\pi}{2}sinh(u)\right]}du
   * @f]
   * gives the following integral:
   * @f[
   *    \int_{-\infty}^{+\infty}f(tanh\left[\frac{\pi}{2}sinh(u)\right])
   *     = \sum_{k=-n}^{+n} 
   * @f]
   */
  template<typename Tp, typename FuncTp>
    adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_tanh_sinh(FuncTp func, Tp a, Tp b,
			Tp max_abs_err, Tp max_rel_err,
			int max_iter = 4);

  /**
   * @f[
   *    \int_{-\infty}^{+\infty}f(x)dx
   * @f]
   * Making the change of variables:
   * @f[
   *    x = sinh\left[\frac{\pi}{2}sinh(u)\right],
   *   dx = \frac{\pi}{2} cosh(u)
   *        cosh\left[\frac{\pi}{2}sinh(u)\right]du
   * @f]
   * gives the following integral:
   * @f[
   *    \int_{-\infty}^{+\infty}f(sinh\left[\frac{\pi}{2}sinh(u)\right])
   *     = \sum_{k=-n}^{+n} 
   * @f]
   */
  template<typename Tp, typename FuncTp>
    adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_sinh_sinh(FuncTp func,
			Tp max_abs_err, Tp max_rel_err,
			int max_iter = 8);

  /**
   * @f[
   *    \int_{0}^{+\infty}f(x)dx
   * @f]
   * Making the change of variables:
   * @f[
   *    x = exp\left[\frac{\pi}{2}sinh(u)\right],
   *   dx = \frac{\pi}{2} cosh(u)
   *        exp\left[\frac{\pi}{2}sinh(u)\right]du
   * @f]
   * gives the following integral:
   * @f[
   *    \int_{0}^{+\infty}f(exp\left[\frac{\pi}{2}sinh(u)\right])
   *     = \sum_{k=-n}^{+n} 
   * @f]
   *
   * This function allows a non-zero lower limit @c a.
   *
   * @param  func  The function to be integrated.
   * @param  a  The lower limit of the semi-infinite integral.
   */
  template<typename Tp, typename FuncTp>
    adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_exp_sinh(FuncTp func, Tp a,
			Tp max_abs_err, Tp max_rel_err,
			int max_iter = 4);

  template<typename Tp, typename FuncTp>
    adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_trapezoid(FuncTp func, Tp a, Tp b,
			Tp max_abs_err, Tp max_rel_err,
			int max_iter);

  template<typename Tp, typename FuncTp>
    adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_midpoint(FuncTp func, Tp a, Tp b,
			Tp max_abs_err, Tp max_rel_err,
			int max_iter);

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_legendre(int n,
				   FuncTp func, Tp a, Tp b);

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_chebyshev_t(int n,
				      FuncTp func, Tp a, Tp b);

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_chebyshev_u(int n,
				      FuncTp func, Tp a, Tp b);

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_chebyshev_v(int n,
				      FuncTp func, Tp a, Tp b);

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_chebyshev_w(int n,
				      FuncTp func, Tp a, Tp b);

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_gegenbauer(int n, Tp lambda,
				     FuncTp func, Tp a, Tp b);

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_jacobi(int n, Tp alf, Tp bet,
				 FuncTp func, Tp a, Tp b);

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_laguerre(int n, Tp alf,
				   FuncTp func, Tp a, Tp b);

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_hermite(int n, Tp alf,
				  FuncTp func, Tp a, Tp b);

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_exponential(int n, Tp alf,
				      FuncTp func, Tp a, Tp b);

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_rational(int n, Tp alf, Tp bet,
				   FuncTp func, Tp a, Tp b);

} // namespace emsr

#include <ext/trapezoid_integral.tcc>
#include <ext/midpoint_integral.tcc>
#include <ext/simpson_integral.tcc>
#include <ext/gauss_kronrod_integral.tcc>
#include <ext/qag_integrate.tcc>
#include <ext/qags_integrate.tcc>
#include <ext/qng_integrate.tcc>
#include <ext/qagp_integrate.tcc>
#include <ext/qcheb_integrate.tcc>
#include <ext/qawc_integrate.tcc>
#include <ext/qaws_integrate.tcc>
#include <ext/qawo_integrate.tcc>
#include <ext/qawf_integrate.tcc>
#include <ext/glfixed_integrate.tcc>
#include <ext/cquad_integrate.tcc>
#include <ext/double_exp_integrate.tcc>
#include <ext/gauss_quadrature.tcc>

#include <ext/integration.tcc>

#endif // INTEGRATION_H
