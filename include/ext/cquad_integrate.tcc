//
// Copyright (C) 2010 Pedro Gonnet
// Copyright (C) 2016-2020 Free Software Foundation, Inc.
// Copyright (C) 2021-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
///
// Ported from GSL by Ed Smith-Rowland
// Originally written by Pedro Gonnet
//
// This file implements the cquad integration scheme.
// Based on structs in gsl/integration/cquad.c

#ifndef CQUAD_INTEGRATE_TCC
#define CQUAD_INTEGRATE_TCC 1

#include <type_traits>

#include <ext/cquad_const.tcc>
#include <ext/cquad_workspace.h>
#include <ext/complex_util.h> // isinf/isnan for complex

namespace emsr
{

  /**
   * Compute the product of the fx with one of the inverse
   * Vandermonde-like matrices.
   */
  template<typename RetTp>
    void
    Vinvfx(const std::array<RetTp, 33>& fx,
	   RetTp* coeff, const int depth)
    {
      switch (depth)
	{
	case 0:
	  for (int i = 0; i <= 4; ++i)
	    {
	      coeff[i] = RetTp{0};
	      for (int j = 0; j <= 4; ++j)
		coeff[i] += RetTp(V1inv[i * 5 + j]) * fx[j * 8];
	    }
	  break;

	case 1:
	  for (int i = 0; i <= 8; ++i)
	    {
	      coeff[i] = RetTp{0};
	      for (int j = 0; j <= 8; ++j)
		coeff[i] +=  RetTp(V2inv[i * 9 + j]) * fx[j * 4];
	    }
	  break;

	case 2:
	  for (int i = 0; i <= 16; ++i)
	    {
	      coeff[i] = RetTp{0};
	      for (int j = 0; j <= 16; ++j)
		coeff[i] +=  RetTp(V3inv[i * 17 + j]) * fx[j * 2];
	    }
	  break;

	case 3:
	  for (int i = 0; i <= 32; ++i)
	    {
	      coeff[i] = RetTp{0};
	      for (int j = 0; j <= 32; ++j)
		coeff[i] +=  RetTp(V4inv[i * 33 + j]) * fx[j];
	    }
	  break;
	}
    }

  /**
   * Downdate the interpolation given by the n coefficients c
   * by removing the nodes with indices in NaNs.
   */
  template<typename Tp>
    void
    downdate(Tp* coeff, std::ptrdiff_t n, std::ptrdiff_t depth,
	     std::ptrdiff_t* NaN, std::ptrdiff_t num_NaNs)
    {
      constexpr std::ptrdiff_t bidx[4] = { 0, 6, 16, 34 };
      Tp b_new[34], alpha;

      for (std::ptrdiff_t i = 0; i <= n + 1; ++i)
	b_new[i] = bee[bidx[depth] + i];
      for (std::ptrdiff_t i = 0; i < num_NaNs; ++i)
	{
	  b_new[n + 1] = b_new[n + 1] / Tp(Lalpha[n]);
	  b_new[n] = (b_new[n] + Tp(xi[NaN[i]]) * b_new[n + 1])
		       / Tp(Lalpha[n - 1]);
	  for (std::ptrdiff_t j = n - 1; j > 0; --j)
	    b_new[j] = (b_new[j] + Tp(xi[NaN[i]]) * b_new[j + 1]
			   - Tp(Lgamma[j + 1]) * b_new[j + 2])
			 / Tp(Lalpha[j - 1]);
	  for (std::ptrdiff_t j = 0; j <= n; ++j)
	    b_new[j] = b_new[j + 1];
	  alpha = coeff[n] / b_new[n];
	  for (std::ptrdiff_t j = 0; j < n; ++j)
	    coeff[j] -= alpha * b_new[j];
	  coeff[n] = Tp{0};
	  --n;
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
  template<typename Tp, typename FuncTp>
    auto
    cquad_integrate(cquad_workspace<Tp, std::invoke_result_t<FuncTp, Tp>>& ws,
		    FuncTp func,
		    Tp a, Tp b,
		    Tp epsabs, Tp epsrel)
    -> adaptive_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      using RetTp = std::invoke_result_t<FuncTp, Tp>;
      using AreaTp = decltype(RetTp{} * Tp{});
      //using AbsAreaTp = decltype(std::abs(AreaTp{}));
      using AbsRetTp = decltype(std::abs(RetTp{}));

      // Some constants that we will need.
      constexpr std::ptrdiff_t n[4] = { 4, 8, 16, 32 };
      constexpr std::ptrdiff_t skip[4] = { 8, 4, 2, 1 };
      constexpr std::ptrdiff_t idx[4] = { 0, 5, 14, 31 };
      constexpr std::ptrdiff_t ndiv_max = 20;
      constexpr auto s_eps = std::numeric_limits<Tp>::epsilon();
      constexpr auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
      constexpr auto s_inf = std::numeric_limits<Tp>::infinity();
      constexpr Tp s_sqrt2 = M_SQRT2;
      constexpr auto w = s_sqrt2 / Tp{2};

      // Actual variables (as opposed to constants above).
      bool split;
      std::ptrdiff_t num_NaNs, NaN[32];
      AbsRetTp nc, ncdiff;

      // Check for unreasonable accuracy demands.
      if (epsabs < Tp{0} || epsrel < Tp{0})
	throw std::domain_error("tolerances may not be negative");
      if (epsabs <= Tp{0} && epsrel < s_eps)
	throw std::domain_error("unreasonable accuracy requirement");

      // Create the first interval.
      ws.clear();
      cquad_interval<Tp, RetTp> iv;
      auto m = (a + b) / Tp{2};
      auto h = (b - a) / Tp{2};
      num_NaNs = 0;
      for (std::ptrdiff_t i = 0; i <= n[3]; ++i)
	{
	  iv.fx[i] = func(m + Tp(xi[i]) * h);
	  if (std::isinf(iv.fx[i]) || std::isnan(iv.fx[i]))
	    {
	      NaN[num_NaNs++] = i;
	      iv.fx[i] = RetTp{0};
	    }
	}
      Vinvfx(iv.fx, &(iv.m_coeff[idx[0]]), 0);
      Vinvfx(iv.fx, &(iv.m_coeff[idx[3]]), 3);
      Vinvfx(iv.fx, &(iv.m_coeff[idx[2]]), 2);
      for (std::ptrdiff_t i = 0; i < num_NaNs; ++i)
	iv.fx[NaN[i]] = s_NaN;
      iv.m_lower_lim = a;
      iv.m_upper_lim = b;
      iv.depth = 3;
      iv.rdepth = 1;
      iv.ndiv = 0;
      iv.m_result = Tp{2} * h * iv.m_coeff[idx[3]] * w;
      nc = AbsRetTp{0};
      for (std::ptrdiff_t i = n[2] + 1; i <= n[3]; ++i)
	{
	  const auto temp = std::abs(iv.m_coeff[idx[3] + i]);
	  nc += temp * temp;
	}
      ncdiff = nc;
      for (std::ptrdiff_t i = 0; i <= n[2]; ++i)
	{
	  const auto temp = std::abs(iv.m_coeff[idx[2] + i]
				     - iv.m_coeff[idx[3] + i]);
	  ncdiff += temp * temp;
	  nc += std::abs(iv.m_coeff[idx[3] + i])
		* std::abs(iv.m_coeff[idx[3] + i]);
	}
      ncdiff = std::sqrt(ncdiff);
      nc = std::sqrt(nc);
      iv.m_abs_error = ncdiff * Tp{2} * h;
      if (ncdiff / nc > Tp{0.1} && iv.m_abs_error < Tp{2} * h * nc)
	iv.m_abs_error = Tp{2} * h * nc;
      ws.push(iv);

      // Main loop...
      auto igral = iv.m_result;
      auto igral_final = AreaTp{0};
      auto err = iv.m_abs_error;
      auto err_final = Tp{0};
      while (ws.size() > 0 && err > Tp{0} &&
	     !(err <= std::abs(igral) * epsrel || err <= epsabs)
	     && !(err_final > std::abs(igral) * epsrel
		  && err - err_final < std::abs(igral) * epsrel)
	     && !(err_final > epsabs && err - err_final < epsabs))
	{
	  // Put our finger on the interval with the largest error.
	  auto& iv = ws.top();
	  m = (iv.m_lower_lim + iv.m_upper_lim) / Tp{2};
	  h = (iv.m_upper_lim - iv.m_lower_lim) / Tp{2};

	  // Should we try to increase the degree?
	  if (iv.depth < 3)
	    {
	      // Keep tabs on some variables.
	      auto depth = ++iv.depth;

	      // Get the new (missing) function values.
	      for (std::ptrdiff_t i = skip[depth];
			i <= 32; i += 2 * skip[depth])
		iv.fx[i] = func(m + Tp(xi[i]) * h);
	      num_NaNs = 0;
	      for (std::ptrdiff_t i = 0; i <= 32; i += skip[depth])
		if (std::isinf(iv.fx[i]) || std::isnan(iv.fx[i]))
		  {
		    NaN[num_NaNs++] = i;
		    iv.fx[i] = RetTp{0};
		  }

	      // Compute the new coefficients.
	      Vinvfx(iv.fx, &(iv.m_coeff[idx[depth]]), depth);

	      // Downdate any NaNs.
	      if (num_NaNs > 0)
		{
		  downdate(&(iv.m_coeff[idx[depth]]),
			   n[depth], depth, NaN, num_NaNs);
		  for (std::ptrdiff_t i = 0; i < num_NaNs; ++i)
		    iv.fx[NaN[i]] = s_NaN;
		}

	      // Compute the error estimate.
	      nc = AbsRetTp{0};
	      for (std::ptrdiff_t i = n[depth - 1] + 1;
		   i <= n[depth]; ++i)
		{
		  const auto temp = std::abs(iv.m_coeff[idx[depth] + i]);
		  nc += temp * temp;
		}
	      ncdiff = nc;
	      for (std::ptrdiff_t i = 0; i <= n[depth - 1]; ++i)
		{
		  const auto temp
		    = std::abs(iv.m_coeff[idx[depth - 1] + i]
			     - iv.m_coeff[idx[depth] + i]);
		  ncdiff += temp * temp;
		  nc += abs(iv.m_coeff[idx[depth] + i])
			* abs(iv.m_coeff[idx[depth] + i]);
		}
	      ncdiff = std::sqrt(ncdiff);
	      nc = std::sqrt(nc);
	      iv.m_abs_error = ncdiff * Tp{2} * h;

	      // Compute the local integral.
	      iv.m_result = Tp{2} * h * w
			     * iv.m_coeff[idx[depth]];

	      // Split the interval prematurely?
	      split = (nc > Tp{0} && ncdiff / nc > Tp{0.1});
	    }
	  else // Maximum degree reached, just split.
	    split = true;

	  // Should we drop this interval?
	  if ((m + h * Tp(xi[0])) >= (m + h * Tp(xi[1]))
	      || (m + h * Tp(xi[31])) >= (m + h * Tp(xi[32]))
	      || iv.m_abs_error < std::abs(iv.m_result) * s_eps * 10)
	    {
	      // Keep this interval's contribution.
	      err_final += iv.m_abs_error;
	      igral_final += iv.m_result;

	      ws.pop();
	    }
	  else if (split) // Do we need to split this interval?
	    {
	      // Some values we will need often...
	      auto depth = iv.depth;

	      // Generate the interval on the left.
	      cquad_interval<Tp, RetTp> ivl;
	      ivl.m_lower_lim = iv.m_lower_lim;
	      ivl.m_upper_lim = m;
	      ivl.depth = 0;
	      ivl.rdepth = iv.rdepth + 1;
	      ivl.fx[0] = iv.fx[0];
	      ivl.fx[32] = iv.fx[16];
	      for (std::ptrdiff_t i = skip[0]; i < 32; i += skip[0])
		ivl.fx[i] = func((ivl.m_lower_lim + ivl.m_upper_lim)
			      / Tp{2} + Tp(xi[i]) * h / Tp{2});
	      num_NaNs = 0;
	      for (std::ptrdiff_t i = 0; i <= 32; i += skip[0])
		{
		  if (std::isinf(ivl.fx[i]) || std::isnan(ivl.fx[i]))
		    {
		      NaN[num_NaNs++] = i;
		      ivl.fx[i] = RetTp{0};
		    }
		}
	      Vinvfx(ivl.fx, ivl.m_coeff, 0);
	      if (num_NaNs > 0)
		{
		  downdate(ivl.m_coeff, n[0], 0, NaN, num_NaNs);
		  for (std::ptrdiff_t i = 0; i < num_NaNs; ++i)
		    ivl.fx[NaN[i]] = s_NaN;
		}
	      for (std::ptrdiff_t i = 0; i <= n[depth]; ++i)
		{
		  ivl.m_coeff[idx[depth] + i] = Tp{0};
		  for (std::ptrdiff_t j = i; j <= n[depth]; ++j)
		    ivl.m_coeff[idx[depth] + i]
			+= Tp(Tleft[i * 33 + j])
			   * iv.m_coeff[idx[depth] + j];
		}
	      ncdiff = AbsRetTp{0};
	      for (std::ptrdiff_t i = 0; i <= n[0]; ++i)
		{
		  const auto temp
		    = std::abs(ivl.m_coeff[i]
			     - ivl.m_coeff[idx[depth] + i]);
		  ncdiff += temp * temp;
		}
	      for (std::ptrdiff_t i = n[0] + 1; i <= n[depth]; ++i)
		{
		  const auto temp
		    = std::abs(ivl.m_coeff[idx[depth] + i]);
		  ncdiff += temp * temp;
		}
	      ncdiff = std::sqrt(ncdiff);
	      ivl.m_abs_error = ncdiff * h;

	      // Check for divergence.
	      ivl.ndiv = iv.ndiv
			 + (std::abs(iv.m_coeff[0]) > 0
			 && std::abs(ivl.m_coeff[0] / iv.m_coeff[0]) > 2);
	      if (ivl.ndiv > ndiv_max && 2 * ivl.ndiv > ivl.rdepth)
		{
		  //const auto result = std::copysign(s_inf, igral);
		  return {AreaTp(s_inf), s_inf};
		}

	      // Compute the local integral.
	      ivl.m_result = h * w * ivl.m_coeff[0];

	      // Generate the interval on the right.
	      cquad_interval<Tp, RetTp> ivr;
	      ivr.m_lower_lim = m;
	      ivr.m_upper_lim = iv.m_upper_lim;
	      ivr.depth = 0;
	      ivr.rdepth = iv.rdepth + 1;
	      ivr.fx[0] = iv.fx[16];
	      ivr.fx[32] = iv.fx[32];
	      for (std::ptrdiff_t i = skip[0]; i < 32; i += skip[0])
		ivr.fx[i] = func((ivr.m_lower_lim + ivr.m_upper_lim)
			      / Tp{2} + Tp(xi[i]) * h / Tp{2});
	      num_NaNs = 0;
	      for (std::ptrdiff_t i = 0; i <= 32; i += skip[0])
		{
		  if (std::isinf(ivr.fx[i]) || std::isnan(ivr.fx[i]))
		    {
		      NaN[num_NaNs++] = i;
		      ivr.fx[i] = RetTp{0};
		    }
		}
	      Vinvfx (ivr.fx, ivr.m_coeff, 0);
	      if (num_NaNs > 0)
		{
		  downdate(ivr.m_coeff, n[0], 0, NaN, num_NaNs);
		  for (std::ptrdiff_t i = 0; i < num_NaNs; ++i)
		    ivr.fx[NaN[i]] = s_NaN;
		}
	      for (std::ptrdiff_t i = 0; i <= n[depth]; ++i)
		{
		  ivr.m_coeff[idx[depth] + i] = Tp{0};
		  for (std::ptrdiff_t j = i; j <= n[depth]; ++j)
		    ivr.m_coeff[idx[depth] + i]
			+= Tp(Tright[i * 33 + j])
			 * iv.m_coeff[idx[depth] + j];
		}
	      ncdiff = AbsRetTp{0};
	      for (std::ptrdiff_t i = 0; i <= n[0]; ++i)
		{
		  const auto temp
		    = std::abs(ivr.m_coeff[i]
			    - ivr.m_coeff[idx[depth] + i]);
		  ncdiff += temp * temp;
		}
	      for (std::ptrdiff_t i = n[0] + 1; i <= n[depth]; ++i)
		{
		  const auto temp
		    = std::abs(ivr.m_coeff[idx[depth] + i]);
		  ncdiff += temp * temp;
		}
	      ncdiff = std::sqrt(ncdiff);
	      ivr.m_abs_error = ncdiff * h;

	      // Check for divergence.
	      ivr.ndiv = iv.ndiv
			 + (std::abs(iv.m_coeff[0]) > 0
			 && std::abs(ivr.m_coeff[0] / iv.m_coeff[0]) > 2);
	      if (ivr.ndiv > ndiv_max && 2 * ivr.ndiv > ivr.rdepth)
		{
		  //const auto result = std::copysign(s_inf, igral);
		  return {AreaTp(s_inf), s_inf};
		}

	      // Compute the local integral.
	      ivr.m_result = h * w * ivr.m_coeff[0];

	      ws.pop();
	      ws.push(ivl);
	      ws.push(ivr);
	    }
	  else // Otherwise, just fix-up the heap.
	    ws.update();

	  // Collect the value of the integral and error.
	  igral = igral_final + ws.total_integral();
	  err = err_final + ws.total_error();
	}

      return {igral, err};
    }

} // namespace emsr

#endif // CQUAD_INTEGRATE_TCC
