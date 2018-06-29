// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2011-2018 Free Software Foundation, Inc.
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

#ifndef INTEGRATE_H
#define INTEGRATE_H 1

#include <limits>

namespace __gnu_cxx
{

  /**
   * Integrate a smooth function from a to b.
   *
   * Higher-order rules converge more rapidly for most functions,
   * but may slow convergence for less well-behaved ones.
   *
   * @param __func The function to be integrated.
   * @param __lower The lower limit of integration.
   * @param __upper The upper limit of integration.
   * @param __max_abs_error The absolute error limit.
   * @param __max_rel_error The relative error limit.
   * @param __max_iter is the maximum number of iterations allowed
   * @param __qkintrule is the Gauss-Kronrod integration rule.
   * @return A structure containing the integration result and the error.
   */
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate(_FuncTp __func,
	      _Tp __lower, _Tp __upper,
	      _Tp __max_abs_error,
	      _Tp __max_rel_error,
	      std::size_t __max_iter = 1024,
	      Kronrod_Rule __qkintrule = QK_61)
    {
      integration_workspace<_Tp> __workspace(__max_iter);
      return qag_integrate(__workspace, __func,
			   __lower, __upper,
			   __max_abs_error, __max_rel_error,
			   __max_iter,
			   __qkintrule);
    }

  /**
   * Integrates a smooth function from -infinity to +infinity.
   *
   * @param __func The function to be integrated.
   * @param __max_abs_error The absolute error limit.
   * @param __max_rel_error The relative error limit.
   * @param __max_iter is the maximum number of iterations allowed
   * @param __qkintrule is the Gauss-Kronrod integration rule.
   * @return A structure containing the integration result and the error.
   */
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_infinite(_FuncTp __func,
		       _Tp __max_abs_error,
		       _Tp __max_rel_error,
		       std::size_t __max_iter = 1024,
		       Kronrod_Rule __qkintrule = QK_61)
    {
      integration_workspace<_Tp> __workspace(__max_iter);
      return qag_integrate(__workspace,
			   map_minf_pinf<_FuncTp, _Tp>(__func),
			   _Tp{0}, _Tp{1},
			   __max_abs_error, __max_rel_error,
			   __max_iter,
			   __qkintrule);
    }

  /**
   * Integrate a smooth function from -infinity to finite b.
   *
   * @param __func The function to be integrated.
   * @param __upper The upper limit of integration.
   * @param __max_abs_error The absolute error limit.
   * @param __max_rel_error The relative error limit.
   * @param __max_iter is the maximum number of iterations allowed
   * @param __qkintrule is the Gauss-Kronrod integration rule.
   * @return A structure containing the integration result and the error.
   */
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_from_infinity(_FuncTp __func,
			    _Tp __upper,
			    _Tp __max_abs_error,
			    _Tp __max_rel_error,
			    std::size_t __max_iter = 1024,
			    Kronrod_Rule __qkintrule = QK_61)
    {
      integration_workspace<_Tp> __workspace(__max_iter);
      return qag_integrate(__workspace,
			   map_minf_b(__func, __upper),
			   _Tp{0}, _Tp{1},
			   __max_abs_error, __max_rel_error,
			   __max_iter,
			   __qkintrule);
    }

  /**
   * Integrate a smooth function from finite a to +infinity.
   *
   * @param __func The function to be integrated.
   * @param __upper The upper limit of integration.
   * @param __max_abs_error The absolute error limit.
   * @param __max_rel_error The relative error limit.
   * @param __max_iter is the maximum number of iterations allowed
   * @param __qkintrule is the Gauss-Kronrod integration rule.
   * @return A structure containing the integration result and the error.
   */
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_to_infinity(_FuncTp __func,
			  _Tp __lower,
			  _Tp __max_abs_error,
			  _Tp __max_rel_error,
			  std::size_t __max_iter = 1024,
			  Kronrod_Rule __qkintrule = QK_61)
    {
      integration_workspace<_Tp> __workspace(__max_iter);
      return qag_integrate(__workspace,
			   map_a_pinf(__func, __lower),
			   _Tp{0}, _Tp{1},
			   __max_abs_error, __max_rel_error,
			   __max_iter,
			   __qkintrule);
    }

  /**
   *  Adaptive Gauss-Kronrod integration optimized for
   *  discontinuous or singular functions
   *
   * @param __func The function to be integrated.
   * @param __lower The lower limit of integration.
   * @param __upper The upper limit of integration.
   * @param __max_abs_error The absolute error limit.
   * @param __max_rel_error The relative error limit.
   * @param __max_iter is the maximum number of iterations allowed
   * @return A structure containing the integration result and the error.
   */
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_kronrod_singular(_FuncTp __func,
			       _Tp __lower, _Tp __upper,
			       _Tp __max_abs_error,
			       _Tp __max_rel_error,
			       std::size_t __max_iter = 1024)
    {
      integration_workspace<_Tp> __workspace(__max_iter);
      return qags_integrate(__workspace, __func, __lower, __upper,
			    __max_abs_error, __max_rel_error);
    }

  /**
   * Integrate a potentially singular function from -infinity to +infinity
   *
   * @param __func The function to be integrated.
   * @param __max_abs_error The absolute error limit.
   * @param __max_rel_error The relative error limit.
   * @param __max_iter is the maximum number of iterations allowed
   * @return A structure containing the integration result and the error.
   */
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_singular_infinite(_FuncTp __func,
				_Tp __max_abs_error,
				_Tp __max_rel_error,
				std::size_t __max_iter = 1024)
    {
      integration_workspace<_Tp> __workspace(__max_iter);
      return qags_integrate(__workspace,
			    map_minf_pinf<_FuncTp, _Tp>(__func),
			    _Tp{0}, _Tp{1},
			    __max_abs_error, __max_rel_error);
    }

  /**
   * Integrate a potentially singular function from -infinity to b
   */
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_singular_from_infinity(_FuncTp __func, _Tp __upper,
				     _Tp __max_abs_error,
				     _Tp __max_rel_error,
				     std::size_t __max_iter = 1024)
    {
      integration_workspace<_Tp> __workspace(__max_iter);
      return qags_integrate(__workspace,
			    map_minf_b(__func, __upper),
			    _Tp{0}, _Tp{1},
			    __max_abs_error, __max_rel_error);
    }

  /**
   * Integrates a potentially singular function from a to +infinity
   *
   * @param __func The function to be integrated.
   * @param __lower The lower limit of integration.
   * @param __max_abs_error The absolute error limit.
   * @param __max_rel_error The relative error limit.
   * @param __max_iter is the maximum number of iterations allowed
   * @return A structure containing the integration result and the error.
   */
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_singular_to_infinity(_FuncTp __func, _Tp __lower,
				   _Tp __max_abs_error,
				   _Tp __max_rel_error,
				   std::size_t __max_iter = 1024)
    {
      integration_workspace<_Tp> __workspace(__max_iter);
      return qags_integrate(__workspace,
			    map_a_pinf(__func, __lower),
			    _Tp{0}, _Tp{1},
			    __max_abs_error, __max_rel_error);
    }

  /**
   * Integrate a potentially singular function.
   *
   * @param __func The function to be integrated.
   * @param __lower The lower limit of integration.
   * @param __upper The upper limit of integration.
   * @param __max_abs_error The absolute error limit.
   * @param __max_rel_error The relative error limit.
   * @param __max_iter is the maximum number of iterations allowed
   * @return A structure containing the integration result and the error.
   */
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_singular(_FuncTp __func,
		       _Tp __lower, _Tp __upper,
		       _Tp __max_abs_error,
		       _Tp __max_rel_error,
		       std::size_t __max_iter = 1024)
    {
      using std::isnan; // To avoid ambiguous overload
      const _Tp __infty = std::numeric_limits<_Tp>::infinity();
      const _Tp __NaN = std::numeric_limits<_Tp>::quiet_NaN();

      if (isnan(__lower) || isnan(__upper))
	return std::make_tuple(__NaN, __NaN);
      else if (__lower == -__infty)
	{
	  if (__upper == __infty) // Integration from -inf to +inf
	    return integrate_singular_infinite(__func, __max_abs_error,
				      __max_rel_error, __max_iter);
	  else if (__upper == -__infty)
	    std::__throw_runtime_error("integrate_singular: "
				       "Attempt to integrate from -infinity"
				       " to -infinity");
	  else // Integration from -inf to finite value
	    return integrate_singular_from_infinity(__func, __upper,
						    __max_abs_error,
						    __max_rel_error,
						    __max_iter);
	}
      else if (__lower == __infty)
	{
	  if (__upper == __infty)
	    std::__throw_runtime_error("integrate_singular: "
				       "Attempted to integrate from +infinity"
				       " to +infinity");
	  else if (__upper == -__infty) // Integration from +inf to -inf,
	    {		     // Call integrate_infinite() and flip sign
	      std::tuple<_Tp, _Tp> __res
		= integrate_singular_infinite(__func, __max_abs_error,
				     __max_rel_error, __max_iter);
	      return std::make_tuple(-std::get<0>(__res), std::get<1>(__res));
	    }
	  else // Integration from +inf to finite value,
	    { // Call integrate_singular_to_infinity and flip sign
	      std::tuple<_Tp, _Tp>
		__res = integrate_singular_to_infinity(__func, __upper,
						       __max_abs_error,
						       __max_rel_error,
						       __max_iter);
	      return std::make_tuple(-std::get<0>(__res), std::get<1>(__res));
	    }
	}
      else // a is finite
	{
	  if (__upper == __infty)
	    return integrate_singular_to_infinity(__func, __lower,
						  __max_abs_error,
						  __max_rel_error, __max_iter);
	  else if (__upper == -__infty)
	    { // Call integrate_from_infinity and flip sign
	      std::tuple<_Tp, _Tp>
		__res = integrate_singular_from_infinity(__func, __lower,
							 __max_abs_error,
							 __max_rel_error,
							 __max_iter);
	      return std::make_tuple(-std::get<0>(__res), std::get<1>(__res));
	    }
	  else // Both a and b finite, call integrate_singular
	    return integrate_kronrod_singular(__func, __lower, __upper,
					      __max_abs_error, __max_rel_error,
					      __max_iter);
	}
    }

  /**
   * Integrate an oscillatory function.
   */
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_oscillatory(_FuncTp __func,
			  _Tp __lower, _Tp __upper,
			  _Tp __max_abs_error,
			  _Tp __max_rel_error,
			  std::size_t __max_iter = 1024)
    {
      using __iwksp_t = integration_workspace<_Tp>;
      using __oitab_t = oscillatory_integration_table<_Tp>;

      __iwksp_t __w(__max_iter);
      __oitab_t __wo(__upper, _Tp{1}, __oitab_t::INTEG_SINE, __max_iter);
      return qawo_integrate(__w, __wo, __func, __lower,
			    __max_abs_error, __max_rel_error);
    }

  /**
   * Adaptively integrate a function with known singular/discontinuous points.
   *
   * @tparam _FuncTp     A function type that takes a single real scalar
   *                     argument and returns a real scalar.
   * @tparam _Tp         A real type for the limits of integration and the step.
   */
  template<typename _FuncTp, typename _FwdIter, typename _Tp>//, typename _Integrator>
    inline std::tuple<_Tp, _Tp>
    integrate_multisingular(_FuncTp __func,
			    _FwdIter __ptbeg, _FwdIter __ptend,
			    _Tp __max_abs_error,
			    _Tp __max_rel_error,
			    std::size_t __max_iter = 1024)
    {
      //using __pt_t = decltype(*__ptbeg);
      //using __ret_t = decltype(__func(__pt_t()));
      //using __area_t = decltype(__pt_t() * __func(__pt_t()));
      //using __err_t = decltype(std::abs(__area_t()));

      integration_workspace<_Tp> __wk(__max_iter);
      return qagp_integrate(__wk, __func, std::vector(__ptbeg, __ptend),
			    __max_abs_error, __max_rel_error);
    }

  /**
   * Integrate a function using an adaptive Clenshaw-Curtis algorithm.
   *
   * @tparam _FuncTp     A function type that takes a single real scalar
   *                     argument and returns a real scalar.
   * @tparam _Tp         A real type for the limits of integration and the step.
   */
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_clenshaw_curtis(_FuncTp __func,
			      _Tp __lower, _Tp __upper,
			      _Tp __max_abs_error,
			      _Tp __max_rel_error,
			      std::size_t __max_iter = 1024)
    {
      //using __pt_t = decltype(*__ptbeg);
      //using __ret_t = decltype(__func(__pt_t()));
      //using __area_t = decltype(__pt_t() * __func(__pt_t()));
      //using __err_t = decltype(std::abs(__area_t()));

      cquad_workspace<_Tp> __ws(__max_iter);

      return cquad_integrate(__ws, __func, __lower, __upper,
			     __max_abs_error, __max_rel_error);
    }

  /**
   * Adaptively integrate a function using a recursive Gauss-Kronrod quadrature
   * called the Patterson algorithm.
   */
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_patterson(_FuncTp __func,
			_Tp __lower, _Tp __upper,
			_Tp __max_abs_error,
			_Tp __max_rel_error)
    {
      return qng_integrate(__func, __lower, __upper,
			   __max_abs_error, __max_rel_error);
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
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_singular_endpoints(_FuncTp __func,
				 _Tp __lower, _Tp __upper,
				 _Tp __alpha, _Tp __beta,
				 int __mu, int __nu,
				 _Tp __max_abs_error,
				 _Tp __max_rel_error,
				 std::size_t __max_iter = 1024)
    {
      integration_workspace<_Tp> __wksp(__max_iter);
      qaws_integration_table<_Tp> __tab(__alpha, __beta, __mu, __nu);

      return qaws_integrate(__wksp, __tab, __func, __lower, __upper,
			    __max_abs_error, __max_rel_error);
    }

  /**
   * Integrate the principal value of a function with a Cauchy singularity.
   */
  template<typename _FuncTp, typename _Tp>
    inline std::tuple<_Tp, _Tp>
    integrate_cauchy_principal_value(_FuncTp __func,
				     _Tp __lower, _Tp __upper, _Tp __center,
				     _Tp __max_abs_err, _Tp __max_rel_err,
				     std::size_t __max_iter = 1024)
    {
      integration_workspace<_Tp> __wksp(__max_iter);

      return  qawc_integrate(__wksp, __func, __lower, __upper, __center,
			     __max_abs_err, __max_rel_err);
    }


} // namespace __gnu_cxx

#endif // INTEGRATE_H
