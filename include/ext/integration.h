// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2011-2020 Free Software Foundation, Inc.
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


#ifndef INTEGRATION_H
#define INTEGRATION_H 1

#include <limits>
#include <tuple>
#include <complex> // For complex abs
#include <type_traits>

#include <ext/quadrature_point.h>
#include <ext/gauss_kronrod_integral.h>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * The return type for a fixed integral rule.
   * Fixed integral types could return errors too!
   * Globally vs. locally adaptive integrals need to be distinguished.
   * abs(diff) vs. diff(abs) is also a thing.
   */
  template<typename _Tp, typename _RetTp>
    struct fixed_integral_t
    {
      using _AreaTp = decltype(_RetTp{} * _Tp{});

      /// Result of the integral.
      _AreaTp __result = _AreaTp{};
    };

  /**
   * The return type for an adaptive integral rule.
   */
  template<typename _Tp, typename _RetTp>
    struct adaptive_integral_t
    {
      using _AreaTp = decltype(_RetTp{} * _Tp{});
      using _AbsAreaTp = decltype(std::abs(_AreaTp{}));

      /// Result of the integral.
      _AreaTp    __result = _AreaTp{};
      /// Absolute value of estimated error.
      _AbsAreaTp __abserr = _AbsAreaTp{};
    };

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#include <ext/trapezoid_integral.h>
#include <ext/midpoint_integral.h>
#include <ext/simpson_integral.h>
#include <ext/gauss_quadrature.h>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * An error tolerance model for integrals.
   */
  template<typename _Tp>
    struct error_tolerance_t
    {
      constexpr
      error_tolerance_t(_Tp __max_abs_err, _Tp __max_rel_err,
			unsigned __min_num_passes);

      /// Maximum absolute error tolerance.
      _Tp _M_max_abs_err;
      /// Maximum relative error tolerance.
      _Tp _M_max_rel_err;
      /// Current tolerance.
      _Tp _M_tolerance = this->tolerance(_Tp{1});
      /// Minimum number of consecutive passes before result is OK.
      unsigned _M_min_num_passes;
      /// Current number of passes.
      unsigned _M_num_passes;

      /// Set and return the tolerance given an initial result.
      template<typename _ResTp>
	constexpr _Tp
	tolerance(_ResTp __result)
	{
	  this->_M_tolerance
		  = std::max(this->_M_max_abs_err,
			     _Tp(this->_M_max_rel_err * std::abs(__result)));
	  return this->_M_tolerance;
	}

      /// Return whether convergence is acceptable tolerance
      /// given two subsequent results.
      template<typename _ResTp>
	bool
	test(_ResTp __curr_result, _ResTp __prev_result)
	{
	  const auto __del = std::abs(__curr_result - __prev_result);
	  if (__del < this->_M_max_abs_err
	      || __del < this->_M_max_rel_err * std::abs(__curr_result))
	    ++this->_M_num_passes;
	  else
	    this->_M_num_passes = 0;
	  return this->_M_num_passes >= this->_M_min_num_passes;
	}

      /// Return the current tolerance.
      constexpr _Tp
      tolerance() const
      { return this->_M_tolerance; }

      /// Test for valid error tolerances.
      static constexpr bool
      _S_valid_tolerances(_Tp __max_abs_err, _Tp __max_rel_err)
      {
	constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
	return !(__max_abs_err <= _Tp{0}
	      && (__max_rel_err < _Tp{50} * _S_eps
	       || __max_rel_err < 0.5e-28));
	// I don't understand the etymology of this last number.
      }
    };

  // Hack for now...
  template<typename _Tp>
    inline bool
    valid_tolerances(_Tp __max_abs_err, _Tp __max_rel_err)
    {
      return error_tolerance_t<_Tp>::
	     _S_valid_tolerances(__max_abs_err, __max_rel_err);
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

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
  template<typename _Tp, typename _FuncTp>
    auto
    integrate(_FuncTp __func,
	      _Tp __lower, _Tp __upper,
	      _Tp __max_abs_error,
	      _Tp __max_rel_error,
	      std::size_t __max_iter = 1024,
	      Kronrod_Rule __qkintrule = Kronrod_21)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

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
  template<typename _Tp, typename _FuncTp>
    auto
    integrate_minf_pinf(_FuncTp __func,
			_Tp __max_abs_error,
			_Tp __max_rel_error,
			std::size_t __max_iter = 1024,
			Kronrod_Rule __qkintrule = Kronrod_21)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

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
  template<typename _Tp, typename _FuncTp>
    auto
    integrate_minf_upper(_FuncTp __func,
			 _Tp __upper,
			 _Tp __max_abs_error,
			 _Tp __max_rel_error,
			 std::size_t __max_iter = 1024,
			 Kronrod_Rule __qkintrule = Kronrod_21)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

  template<typename _Tp, typename _FuncTp>
    inline auto
    integrate_lower_minf(_FuncTp __func,
			 _Tp __lower,
			 _Tp __max_abs_error,
			 _Tp __max_rel_error,
			 std::size_t __max_iter = 1024,
			 Kronrod_Rule __qkintrule = Kronrod_21)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    {
      return -integrate_minf_upper(__func, __lower,
				   __max_abs_error, __max_rel_error,
				   __max_iter, __qkintrule);
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
  template<typename _Tp, typename _FuncTp>
    auto
    integrate_lower_pinf(_FuncTp __func,
			 _Tp __lower,
			 _Tp __max_abs_error,
			 _Tp __max_rel_error,
			 std::size_t __max_iter = 1024,
			 Kronrod_Rule __qkintrule = Kronrod_21)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

  template<typename _Tp, typename _FuncTp>
    inline auto
    integrate_pinf_upper(_FuncTp __func,
			 _Tp __upper,
			 _Tp __max_abs_error,
			 _Tp __max_rel_error,
			 std::size_t __max_iter = 1024,
			 Kronrod_Rule __qkintrule = Kronrod_21)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    {
      return -integrate_lower_pinf(__func, __upper,
				   __max_abs_error, __max_rel_error,
				   __max_iter, __qkintrule);
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
  template<typename _Tp, typename _FuncTp>
    auto
    integrate_kronrod_singular(_FuncTp __func,
			       _Tp __lower, _Tp __upper,
			       _Tp __max_abs_error,
			       _Tp __max_rel_error,
			       std::size_t __max_iter = 1024)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

  /**
   * Integrate a potentially singular function from -infinity to +infinity
   *
   * @param func The function to be integrated.
   * @param max_abs_error The absolute error limit.
   * @param max_rel_error The relative error limit.
   * @param max_iter is the maximum number of iterations allowed
   * @return A structure containing the integration result and the error.
   */
  template<typename _Tp, typename _FuncTp>
    auto
    integrate_singular_minf_pinf(_FuncTp __func,
				 _Tp __max_abs_error,
				 _Tp __max_rel_error,
				 std::size_t __max_iter = 1024)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

  /**
   * Integrate a potentially singular function from -infinity to b
   */
  template<typename _Tp, typename _FuncTp>
    auto
    integrate_singular_minf_upper(_FuncTp __func, _Tp __upper,
				  _Tp __max_abs_error,
				  _Tp __max_rel_error,
				  std::size_t __max_iter = 1024)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

  template<typename _Tp, typename _FuncTp>
    inline auto
    integrate_singular_lower_minf(_FuncTp __func, _Tp __lower,
				  _Tp __max_abs_error,
				  _Tp __max_rel_error,
				  std::size_t __max_iter = 1024)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    {
      return -integrate_singular_minf_upper(__func, __lower,
					    __max_abs_error, __max_rel_error,
					    __max_iter);
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
  template<typename _Tp, typename _FuncTp>
    auto
    integrate_singular_lower_pinf(_FuncTp __func, _Tp __lower,
				  _Tp __max_abs_error,
				  _Tp __max_rel_error,
				  std::size_t __max_iter = 1024)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

  template<typename _Tp, typename _FuncTp>
    inline auto
    integrate_singular_pinf_upper(_FuncTp __func, _Tp __upper,
				  _Tp __max_abs_error,
				  _Tp __max_rel_error,
				  std::size_t __max_iter = 1024)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    {
      return -integrate_singular_lower_pinf(__func, __upper,
					    __max_abs_error, __max_rel_error,
					    __max_iter);
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
  template<typename _Tp, typename _FuncTp>
    auto
    integrate_singular(_FuncTp __func,
		       _Tp __lower, _Tp __upper,
		       _Tp __max_abs_error,
		       _Tp __max_rel_error,
		       std::size_t __max_iter = 1024)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

  /**
   * Integrate an oscillatory function.
   */
  template<typename _Tp, typename _FuncTp>
    auto
    integrate_oscillatory(_FuncTp __func,
			  _Tp __lower, _Tp __upper,
			  _Tp __max_abs_error,
			  _Tp __max_rel_error,
			  std::size_t __max_iter = 1024)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

  /**
   * Adaptively integrate a function with known singular/discontinuous points.
   *
   * @tparam _FuncTp     A function type that takes a single real scalar
   *                     argument and returns a real scalar.
   * @tparam _Tp         A real type for the limits of integration and the step.
   */
  template<typename _FuncTp, typename _FwdIter, typename _Tp>
    auto
    integrate_multisingular(_FuncTp __func,
			    _FwdIter __ptbeg, _FwdIter __ptend,
			    _Tp __max_abs_error,
			    _Tp __max_rel_error,
			    std::size_t __max_iter = 1024)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

  /**
   * Integrate a function using an adaptive Clenshaw-Curtis algorithm.
   *
   * @tparam _FuncTp     A function type that takes a single real scalar
   *                     argument and returns a real scalar.
   * @tparam _Tp         A real type for the limits of integration and the step.
   */
  template<typename _Tp, typename _FuncTp>
    auto
    integrate_clenshaw_curtis(_FuncTp __func,
			      _Tp __lower, _Tp __upper,
			      _Tp __max_abs_error,
			      _Tp __max_rel_error,
			      std::size_t __max_iter = 1024)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

  /**
   * Adaptively integrate a function using a recursive Gauss-Kronrod quadrature
   * called the Patterson algorithm.
   */
  template<typename _Tp, typename _FuncTp>
    auto
    integrate_patterson(_FuncTp __func,
			_Tp __lower, _Tp __upper,
			_Tp __max_abs_error,
			_Tp __max_rel_error)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

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
  template<typename _Tp, typename _FuncTp>
    auto
    integrate_singular_endpoints(_FuncTp __func,
				 _Tp __lower, _Tp __upper,
				 _Tp __alpha, _Tp __beta,
				 int __mu, int __nu,
				 _Tp __max_abs_error,
				 _Tp __max_rel_error,
				 std::size_t __max_iter = 1024)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

  /**
   * Integrate the principal value of a function with a Cauchy singularity.
   */
  template<typename _Tp, typename _FuncTp>
    auto
    integrate_cauchy_principal_value(_FuncTp __func,
				     _Tp __lower, _Tp __upper, _Tp __center,
				     _Tp __max_abs_err, _Tp __max_rel_err,
				     std::size_t __max_iter = 1024)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>;

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
  template<typename _Tp, typename _FuncTp>
    adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_tanh_sinh(_FuncTp __func, _Tp __a, _Tp __b,
			_Tp __max_abs_err, _Tp __max_rel_err,
			int __max_iter = 4);

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
  template<typename _Tp, typename _FuncTp>
    adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_sinh_sinh(_FuncTp __func,
			_Tp __max_abs_err, _Tp __max_rel_err,
			int __max_iter = 8);

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
  template<typename _Tp, typename _FuncTp>
    adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_exp_sinh(_FuncTp __func, _Tp __a,
			_Tp __max_abs_err, _Tp __max_rel_err,
			int __max_iter = 4);

  template<typename _Tp, typename _FuncTp>
    adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_trapezoid(_FuncTp __func, _Tp __a, _Tp __b,
			_Tp __max_abs_err, _Tp __max_rel_err,
			int __max_iter);

  template<typename _Tp, typename _FuncTp>
    adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_midpoint(_FuncTp __func, _Tp __a, _Tp __b,
			_Tp __max_abs_err, _Tp __max_rel_err,
			int __max_iter);

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_legendre(int __n,
				   _FuncTp __func, _Tp __a, _Tp __b);

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_chebyshev_t(int __n,
				      _FuncTp __func, _Tp __a, _Tp __b);

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_chebyshev_u(int __n,
				      _FuncTp __func, _Tp __a, _Tp __b);

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_chebyshev_v(int __n,
				      _FuncTp __func, _Tp __a, _Tp __b);

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_chebyshev_w(int __n,
				      _FuncTp __func, _Tp __a, _Tp __b);

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_gegenbauer(int __n, _Tp __lambda,
				     _FuncTp __func, _Tp __a, _Tp __b);

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_jacobi(int __n, _Tp __alf, _Tp __bet,
				 _FuncTp __func, _Tp __a, _Tp __b);

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_laguerre(int __n, _Tp __alf,
				   _FuncTp __func, _Tp __a, _Tp __b);

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_hermite(int __n, _Tp __alf,
				  _FuncTp __func, _Tp __a, _Tp __b);

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_exponential(int __n, _Tp __alf,
				      _FuncTp __func, _Tp __a, _Tp __b);

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_rational(int __n, _Tp __alf, _Tp __bet,
				   _FuncTp __func, _Tp __a, _Tp __b);

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

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
