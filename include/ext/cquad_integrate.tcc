/* integration/cquad_integrate.tcc
 *
 * Copyright (C) 2010 Pedro Gonnet
 * Copyright (C) 2016-2019 Free Software Foundation, Inc.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
// Ported from GSL by Ed Smith-Rowland
// Originally written by Pedro Gonnet
//
// This file implements the cquad integration scheme.
// Based on structs in gsl/integration/cquad.c

#ifndef CQUAD_INTEGRATE_TCC
#define CQUAD_INTEGRATE_TCC 1

#include <ext/cquad_const.tcc>
#include <ext/cquad_workspace.h>

namespace __gnu_cxx
{

  /**
   * Compute the product of the fx with one of the inverse
   * Vandermonde-like matrices.
   */
  template<typename _RetTp>
    void
    _Vinvfx(const std::array<_RetTp, 33>& __fx,
	    _RetTp* __coeff, const int __depth)
    {
      switch (__depth)
	{
	case 0:
	  for (int __i = 0; __i <= 4; ++__i)
	    {
	      __coeff[__i] = _RetTp{0};
	      for (int __j = 0; __j <= 4; ++__j)
		__coeff[__i] += _RetTp(V1inv[__i * 5 + __j]) * __fx[__j * 8];
	    }
	  break;

	case 1:
	  for (int __i = 0; __i <= 8; ++__i)
	    {
	      __coeff[__i] = _RetTp{0};
	      for (int __j = 0; __j <= 8; ++__j)
		__coeff[__i] +=  _RetTp(V2inv[__i * 9 + __j]) * __fx[__j * 4];
	    }
	  break;

	case 2:
	  for (int __i = 0; __i <= 16; ++__i)
	    {
	      __coeff[__i] = _RetTp{0};
	      for (int __j = 0; __j <= 16; ++__j)
		__coeff[__i] +=  _RetTp(V3inv[__i * 17 + __j]) * __fx[__j * 2];
	    }
	  break;

	case 3:
	  for (int __i = 0; __i <= 32; ++__i)
	    {
	      __coeff[__i] = _RetTp{0};
	      for (int __j = 0; __j <= 32; ++__j)
		__coeff[__i] +=  _RetTp(V4inv[__i * 33 + __j]) * __fx[__j];
	    }
	  break;
	}
    }

  /**
   * Downdate the interpolation given by the n coefficients c
   * by removing the nodes with indices in NaNs.
   */
  template<typename _Tp>
    void
    downdate(_Tp* __coeff, std::ptrdiff_t __n, std::ptrdiff_t __depth,
	     std::ptrdiff_t* __NaN, std::ptrdiff_t __num_NaNs)
    {
      constexpr std::ptrdiff_t __bidx[4] = { 0, 6, 16, 34 };
      _Tp __b_new[34], __alpha;

      for (std::ptrdiff_t __i = 0; __i <= __n + 1; ++__i)
	__b_new[__i] = bee[__bidx[__depth] + __i];
      for (std::ptrdiff_t __i = 0; __i < __num_NaNs; ++__i)
	{
	  __b_new[__n + 1] = __b_new[__n + 1] / _Tp(Lalpha[__n]);
	  __b_new[__n] = (__b_new[__n] + _Tp(xi[__NaN[__i]]) * __b_new[__n + 1])
		       / _Tp(Lalpha[__n - 1]);
	  for (std::ptrdiff_t __j = __n - 1; __j > 0; --__j)
	    __b_new[__j] = (__b_new[__j] + _Tp(xi[__NaN[__i]]) * __b_new[__j + 1]
			   - _Tp(Lgamma[__j + 1]) * __b_new[__j + 2])
			 / _Tp(Lalpha[__j - 1]);
	  for (std::ptrdiff_t __j = 0; __j <= __n; ++__j)
	    __b_new[__j] = __b_new[__j + 1];
	  __alpha = __coeff[__n] / __b_new[__n];
	  for (std::ptrdiff_t __j = 0; __j < __n; ++__j)
	    __coeff[__j] -= __alpha * __b_new[__j];
	  __coeff[__n] = _Tp{0};
	  --__n;
	}
    }

  /**
   * CQUAD is a new doubly-adaptive general-purpose quadrature routine
   * which can handle most types of singularities, non-numerical function
   * values such as Inf or NaN, as well as some divergent integrals.
   * It generally requires more function evaluations than
   * the other integration routines in this library, yet fails
   * less often for difficult integrands.
   *
   * The underlying algorithm uses a doubly-adaptive scheme
   * in which Clenshaw-Curtis quadrature rules of increasing degree
   * are used to compute the integral in each interval.
   * The L_2-norm of the difference between the underlying interpolatory
   * polynomials of two successive rules is used as an error estimate.
   * The interval is subdivided if the difference between two successive
   * rules is too large or a rule of maximum degree has been reached.
   *
   * The CQUAD algorithm divides the integration region into subintervals,
   * and in each iteration, the subinterval with the largest estimated error
   * is processed. The algorithm uses Clenshaw-Curtis quadrature rules
   * of degree 4, 8, 16 and 32 over 5, 9, 17 and 33 nodes respectively.
   * Each interval is initialized with the lowest-degree rule.
   * When an interval is processed, the next-higher degree rule is evaluated
   * and an error estimate is computed based on the L_2-norm of the difference
   * between the underlying interpolating polynomials of both rules.
   * If the highest-degree rule has already been used, or the interpolatory
   * polynomials differ significantly, the interval is bisected. 
   */
  template<typename _Tp, typename _FuncTp>
    auto
    cquad_integrate(cquad_workspace<_Tp, std::invoke_result_t<_FuncTp, _Tp>>& __ws,
		    _FuncTp __func,
		    _Tp __a, _Tp __b,
		    _Tp __epsabs, _Tp __epsrel)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    {
      using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
      using _AreaTp = decltype(_RetTp{} * _Tp{});
      //using _AbsAreaTp = decltype(std::abs(_AreaTp{}));
      using _AbsRetTp = decltype(std::abs(_RetTp{}));

      // Some constants that we will need.
      constexpr std::ptrdiff_t __n[4] = { 4, 8, 16, 32 };
      constexpr std::ptrdiff_t __skip[4] = { 8, 4, 2, 1 };
      constexpr std::ptrdiff_t __idx[4] = { 0, 5, 14, 31 };
      constexpr std::ptrdiff_t __ndiv_max = 20;
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
      constexpr auto _S_inf = std::numeric_limits<_Tp>::infinity();
      constexpr _Tp _S_sqrt2 = M_SQRT2;
      constexpr auto __w = _S_sqrt2 / _Tp{2};

      // Actual variables (as opposed to constants above).
      bool __split;
      std::ptrdiff_t __num_NaNs, __NaN[32];
      _AbsRetTp __nc, __ncdiff;

      // Check for unreasonable accuracy demands.
      if (__epsabs < _Tp{0} || __epsrel < _Tp{0})
	std::__throw_domain_error("tolerances may not be negative");
      if (__epsabs <= _Tp{0} && __epsrel < _S_eps)
	std::__throw_domain_error("unreasonable accuracy requirement");

      // Create the first interval.
      __ws.clear();
      cquad_interval<_Tp, _RetTp> __iv;
      auto __m = (__a + __b) / _Tp{2};
      auto __h = (__b - __a) / _Tp{2};
      __num_NaNs = 0;
      for (std::ptrdiff_t __i = 0; __i <= __n[3]; ++__i)
	{
	  __iv.fx[__i] = __func(__m + _Tp(xi[__i]) * __h);
	  if (std::isinf(__iv.fx[__i]) || std::isnan(__iv.fx[__i]))
	    {
	      __NaN[__num_NaNs++] = __i;
	      __iv.fx[__i] = _RetTp{0};
	    }
	}
      _Vinvfx(__iv.fx, &(__iv._M_coeff[__idx[0]]), 0);
      _Vinvfx(__iv.fx, &(__iv._M_coeff[__idx[3]]), 3);
      _Vinvfx(__iv.fx, &(__iv._M_coeff[__idx[2]]), 2);
      for (std::ptrdiff_t __i = 0; __i < __num_NaNs; ++__i)
	__iv.fx[__NaN[__i]] = _S_NaN;
      __iv._M_lower_lim = __a;
      __iv._M_upper_lim = __b;
      __iv.depth = 3;
      __iv.rdepth = 1;
      __iv.ndiv = 0;
      __iv._M_result = _Tp{2} * __h * __iv._M_coeff[__idx[3]] * __w;
      __nc = _AbsRetTp{0};
      for (std::ptrdiff_t __i = __n[2] + 1; __i <= __n[3]; ++__i)
	{
	  const auto __temp = std::abs(__iv._M_coeff[__idx[3] + __i]);
	  __nc += __temp * __temp;
	}
      __ncdiff = __nc;
      for (std::ptrdiff_t __i = 0; __i <= __n[2]; ++__i)
	{
	  const auto __temp = std::abs(__iv._M_coeff[__idx[2] + __i]
				     - __iv._M_coeff[__idx[3] + __i]);
	  __ncdiff += __temp * __temp;
	  __nc += std::abs(__iv._M_coeff[__idx[3] + __i])
		* std::abs(__iv._M_coeff[__idx[3] + __i]);
	}
      __ncdiff = std::sqrt(__ncdiff);
      __nc = std::sqrt(__nc);
      __iv._M_abs_error = __ncdiff * _Tp{2} * __h;
      if (__ncdiff / __nc > _Tp{0.1} && __iv._M_abs_error < _Tp{2} * __h * __nc)
	__iv._M_abs_error = _Tp{2} * __h * __nc;
      __ws.push(__iv);

      // Main loop...
      auto __igral = __iv._M_result;
      auto __igral_final = _AreaTp{0};
      auto __err = __iv._M_abs_error;
      auto __err_final = _Tp{0};
      while (__ws.size() > 0 && __err > _Tp{0} &&
	     !(__err <= std::abs(__igral) * __epsrel || __err <= __epsabs)
	     && !(__err_final > std::abs(__igral) * __epsrel
		  && __err - __err_final < std::abs(__igral) * __epsrel)
	     && !(__err_final > __epsabs && __err - __err_final < __epsabs))
	{
	  // Put our finger on the interval with the largest error.
	  auto& __iv = __ws.top();
	  __m = (__iv._M_lower_lim + __iv._M_upper_lim) / _Tp{2};
	  __h = (__iv._M_upper_lim - __iv._M_lower_lim) / _Tp{2};

	  // Should we try to increase the degree?
	  if (__iv.depth < 3)
	    {
	      // Keep tabs on some variables.
	      auto __depth = ++__iv.depth;

	      // Get the new (missing) function values.
	      for (std::ptrdiff_t __i = __skip[__depth];
			__i <= 32; __i += 2 * __skip[__depth])
		__iv.fx[__i] = __func(__m + _Tp(xi[__i]) * __h);
	      __num_NaNs = 0;
	      for (std::ptrdiff_t __i = 0; __i <= 32; __i += __skip[__depth])
		if (std::isinf(__iv.fx[__i]) || std::isnan(__iv.fx[__i]))
		  {
		    __NaN[__num_NaNs++] = __i;
		    __iv.fx[__i] = _RetTp{0};
		  }

	      // Compute the new coefficients.
	      _Vinvfx(__iv.fx, &(__iv._M_coeff[__idx[__depth]]), __depth);

	      // Downdate any NaNs.
	      if (__num_NaNs > 0)
		{
		  downdate(&(__iv._M_coeff[__idx[__depth]]),
			   __n[__depth], __depth, __NaN, __num_NaNs);
		  for (std::ptrdiff_t __i = 0; __i < __num_NaNs; ++__i)
		    __iv.fx[__NaN[__i]] = _S_NaN;
		}

	      // Compute the error estimate.
	      __nc = _AbsRetTp{0};
	      for (std::ptrdiff_t __i = __n[__depth - 1] + 1;
		   __i <= __n[__depth]; ++__i)
		{
		  const auto __temp = std::abs(__iv._M_coeff[__idx[__depth] + __i]);
		  __nc += __temp * __temp;
		}
	      __ncdiff = __nc;
	      for (std::ptrdiff_t __i = 0; __i <= __n[__depth - 1]; ++__i)
		{
		  const auto __temp
		    = std::abs(__iv._M_coeff[__idx[__depth - 1] + __i]
			     - __iv._M_coeff[__idx[__depth] + __i]);
		  __ncdiff += __temp * __temp;
		  __nc += abs(__iv._M_coeff[__idx[__depth] + __i])
			* abs(__iv._M_coeff[__idx[__depth] + __i]);
		}
	      __ncdiff = std::sqrt(__ncdiff);
	      __nc = std::sqrt(__nc);
	      __iv._M_abs_error = __ncdiff * _Tp{2} * __h;

	      // Compute the local integral.
	      __iv._M_result = _Tp{2} * __h * __w
			     * __iv._M_coeff[__idx[__depth]];

	      // Split the interval prematurely?
	      __split = (__nc > _Tp{0} && __ncdiff / __nc > _Tp{0.1});
	    }
	  else // Maximum degree reached, just split.
	    __split = true;

	  // Should we drop this interval?
	  if ((__m + __h * _Tp(xi[0])) >= (__m + __h * _Tp(xi[1]))
	      || (__m + __h * _Tp(xi[31])) >= (__m + __h * _Tp(xi[32]))
	      || __iv._M_abs_error < std::abs(__iv._M_result) * _S_eps * 10)
	    {
	      // Keep this interval's contribution.
	      __err_final += __iv._M_abs_error;
	      __igral_final += __iv._M_result;

	      __ws.pop();
	    }
	  else if (__split) // Do we need to split this interval?
	    {
	      // Some values we will need often...
	      auto __depth = __iv.depth;

	      // Generate the interval on the left.
	      cquad_interval<_Tp, _RetTp> __ivl;
	      __ivl._M_lower_lim = __iv._M_lower_lim;
	      __ivl._M_upper_lim = __m;
	      __ivl.depth = 0;
	      __ivl.rdepth = __iv.rdepth + 1;
	      __ivl.fx[0] = __iv.fx[0];
	      __ivl.fx[32] = __iv.fx[16];
	      for (std::ptrdiff_t __i = __skip[0]; __i < 32; __i += __skip[0])
		__ivl.fx[__i] = __func((__ivl._M_lower_lim + __ivl._M_upper_lim)
			      / _Tp{2} + _Tp(xi[__i]) * __h / _Tp{2});
	      __num_NaNs = 0;
	      for (std::ptrdiff_t __i = 0; __i <= 32; __i += __skip[0])
		{
		  if (std::isinf(__ivl.fx[__i]) || std::isnan(__ivl.fx[__i]))
		    {
		      __NaN[__num_NaNs++] = __i;
		      __ivl.fx[__i] = _RetTp{0};
		    }
		}
	      _Vinvfx(__ivl.fx, __ivl._M_coeff, 0);
	      if (__num_NaNs > 0)
		{
		  downdate(__ivl._M_coeff, __n[0], 0, __NaN, __num_NaNs);
		  for (std::ptrdiff_t __i = 0; __i < __num_NaNs; ++__i)
		    __ivl.fx[__NaN[__i]] = _S_NaN;
		}
	      for (std::ptrdiff_t __i = 0; __i <= __n[__depth]; ++__i)
		{
		  __ivl._M_coeff[__idx[__depth] + __i] = _Tp{0};
		  for (std::ptrdiff_t __j = __i; __j <= __n[__depth]; ++__j)
		    __ivl._M_coeff[__idx[__depth] + __i]
			+= _Tp(Tleft[__i * 33 + __j])
			   * __iv._M_coeff[__idx[__depth] + __j];
		}
	      __ncdiff = _AbsRetTp{0};
	      for (std::ptrdiff_t __i = 0; __i <= __n[0]; ++__i)
		{
		  const auto __temp
		    = std::abs(__ivl._M_coeff[__i]
			     - __ivl._M_coeff[__idx[__depth] + __i]);
		  __ncdiff += __temp * __temp;
		}
	      for (std::ptrdiff_t __i = __n[0] + 1; __i <= __n[__depth]; ++__i)
		{
		  const auto __temp
		    = std::abs(__ivl._M_coeff[__idx[__depth] + __i]);
		  __ncdiff += __temp * __temp;
		}
	      __ncdiff = std::sqrt(__ncdiff);
	      __ivl._M_abs_error = __ncdiff * __h;

	      // Check for divergence.
	      __ivl.ndiv = __iv.ndiv
			 + (std::abs(__iv._M_coeff[0]) > 0
			 && std::abs(__ivl._M_coeff[0] / __iv._M_coeff[0]) > 2);
	      if (__ivl.ndiv > __ndiv_max && 2 * __ivl.ndiv > __ivl.rdepth)
		{
		  //const auto __result = std::copysign(_S_inf, __igral);
		  return {_AreaTp(_S_inf), _S_inf};
		}

	      // Compute the local integral.
	      __ivl._M_result = __h * __w * __ivl._M_coeff[0];

	      // Generate the interval on the right.
	      cquad_interval<_Tp, _RetTp> __ivr;
	      __ivr._M_lower_lim = __m;
	      __ivr._M_upper_lim = __iv._M_upper_lim;
	      __ivr.depth = 0;
	      __ivr.rdepth = __iv.rdepth + 1;
	      __ivr.fx[0] = __iv.fx[16];
	      __ivr.fx[32] = __iv.fx[32];
	      for (std::ptrdiff_t __i = __skip[0]; __i < 32; __i += __skip[0])
		__ivr.fx[__i] = __func((__ivr._M_lower_lim + __ivr._M_upper_lim)
			      / _Tp{2} + _Tp(xi[__i]) * __h / _Tp{2});
	      __num_NaNs = 0;
	      for (std::ptrdiff_t __i = 0; __i <= 32; __i += __skip[0])
		{
		  if (std::isinf(__ivr.fx[__i]) || std::isnan(__ivr.fx[__i]))
		    {
		      __NaN[__num_NaNs++] = __i;
		      __ivr.fx[__i] = _RetTp{0};
		    }
		}
	      _Vinvfx (__ivr.fx, __ivr._M_coeff, 0);
	      if (__num_NaNs > 0)
		{
		  downdate(__ivr._M_coeff, __n[0], 0, __NaN, __num_NaNs);
		  for (std::ptrdiff_t __i = 0; __i < __num_NaNs; ++__i)
		    __ivr.fx[__NaN[__i]] = _S_NaN;
		}
	      for (std::ptrdiff_t __i = 0; __i <= __n[__depth]; ++__i)
		{
		  __ivr._M_coeff[__idx[__depth] + __i] = _Tp{0};
		  for (std::ptrdiff_t __j = __i; __j <= __n[__depth]; ++__j)
		    __ivr._M_coeff[__idx[__depth] + __i]
			+= _Tp(Tright[__i * 33 + __j])
			 * __iv._M_coeff[__idx[__depth] + __j];
		}
	      __ncdiff = _AbsRetTp{0};
	      for (std::ptrdiff_t __i = 0; __i <= __n[0]; ++__i)
		{
		  const auto __temp
		    = std::abs(__ivr._M_coeff[__i]
			    - __ivr._M_coeff[__idx[__depth] + __i]);
		  __ncdiff += __temp * __temp;
		}
	      for (std::ptrdiff_t __i = __n[0] + 1; __i <= __n[__depth]; ++__i)
		{
		  const auto __temp
		    = std::abs(__ivr._M_coeff[__idx[__depth] + __i]);
		  __ncdiff += __temp * __temp;
		}
	      __ncdiff = std::sqrt(__ncdiff);
	      __ivr._M_abs_error = __ncdiff * __h;

	      // Check for divergence.
	      __ivr.ndiv = __iv.ndiv
			 + (std::abs(__iv._M_coeff[0]) > 0
			 && std::abs(__ivr._M_coeff[0] / __iv._M_coeff[0]) > 2);
	      if (__ivr.ndiv > __ndiv_max && 2 * __ivr.ndiv > __ivr.rdepth)
		{
		  //const auto __result = std::copysign(_S_inf, __igral);
		  return {_AreaTp(_S_inf), _S_inf};
		}

	      // Compute the local integral.
	      __ivr._M_result = __h * __w * __ivr._M_coeff[0];

	      __ws.pop();
	      __ws.push(__ivl);
	      __ws.push(__ivr);
	    }
	  else // Otherwise, just fix-up the heap.
	    __ws.update();

	  // Collect the value of the integral and error.
	  __igral = __igral_final + __ws.total_integral();
	  __err = __err_final + __ws.total_error();
	}

      return {__igral, __err};
    }

} // namespace __gnu_cxx

#endif // CQUAD_INTEGRATE_TCC
