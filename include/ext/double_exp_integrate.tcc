/* quadrature/double_exp_integrate.tcc
 *
 * Copyright (C) 2017-2020 Free Software Foundation, Inc.
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

#ifndef DOUBLE_EXP_INTEGRATE_TCC
#define DOUBLE_EXP_INTEGRATE_TCC 1

#include <type_traits>
#include <cmath>

namespace emsr
{

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
    integrate_tanh_sinh(FuncTp func, Tp lower, Tp upper,
			Tp max_abs_err, Tp max_rel_err,
			int max_iter)
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      const auto s_pi_4 = Tp{3.141592653589793238462643383279502884195L} / 4;

      if (std::isnan(lower) || std::isnan(upper)
          || std::isnan(max_abs_err) || std::isnan(max_rel_err))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}, absarea_t{}};
      else
	{
          int n = 16;
          n /= 2;

          // Find K = ln(ln(max_number))
          const auto k_max
		      = std::log(std::log(std::numeric_limits<Tp>::max()))
		      - Tp{1};
          auto h = k_max / n;

          auto sum = func((lower + upper) / Tp{2}) / Tp{2};
          decltype(sum) sum1{}, sum2{};
          for (int k = -n; k < 0; ++k)
	    {
	      const auto u = h * Tp(k);
	      const auto eu = std::exp(u);
	      const auto cosh = eu + Tp{1} / eu;
	      const auto sinh = eu - Tp{1} / eu;
              const auto esh = std::exp(s_pi_4 * sinh);
	      const auto w = esh + Tp{1} / esh;
	      const auto dxdu = cosh / (w * w);
	      const auto x1 = (upper * esh + lower / esh) / w;
	      if (x1 != lower && x1 != upper) 
	        sum1 += dxdu * func(x1);
	      const auto x2 = (lower * esh + upper / esh) / w;
	      if (x2 != lower && x2 != upper)
	        sum2 += dxdu * func(x2);
	    }

          // Interlace values; don't go past the rightmost point.
          auto prev_sum = sum + sum1 + sum2;
          for (int iter = 0; iter < max_iter; ++iter)
	    {
	      for (int k  = -n; k < 0; ++k)
	        {
	          const auto u = h * Tp(k + 0.5);
	          const auto eu = std::exp(u);
	          // A standard sinhcosh would be a nice idea along with sincos.
	          const auto cosh = eu + Tp{1} / eu;
	          const auto sinh = eu - Tp{1} / eu;
                  const auto esh = std::exp(s_pi_4 * sinh);
	          const auto w = esh + Tp{1} / esh;
	          const auto dxdu = cosh / (w * w);
	          // natural: x1 = (s - 1/s) / (s + 1/s)
	          const auto x1 = (upper * esh + lower / esh) / w;
	          if (x1 != lower && x1 != upper) 
		    sum1 += dxdu * func(x1);
	          // natural: x2 = (-s + 1/s) / (s + 1/s)
	          const auto x2 = (lower * esh + upper / esh) / w;
	          if (x2 != lower && x2 != upper)
		    sum2 += dxdu * func(x2);
	        }

	      n *= 2;
	      h /= Tp{2};

	      const auto curr_sum = sum + sum1 + sum2;
	      if (auto abs_del = std::abs(curr_sum - prev_sum);
                  abs_del < max_abs_err
		  || abs_del < std::abs(max_rel_err * curr_sum)
		  || iter + 1 == max_iter) // Keep prev_sum even if at max_iters.
	        break;

	      prev_sum = curr_sum;
	    }

          const auto fact = Tp{2} * (upper - lower) * s_pi_4 * h;
          const auto tot_sum = sum + sum1 + sum2;
          return {fact * tot_sum,
		  fact * std::abs(tot_sum - Tp{2} * prev_sum)};
	}
    }

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
			int max_iter)
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      const auto s_pi_4 = Tp{3.141592653589793238462643383279502884195L} / 4;

      if (std::isnan(max_abs_err) || std::isnan(max_rel_err))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else
	{
          int n = 16;
          n /= 2;

          // Find K = ln(ln(max_number))
          const auto k_max
		      = std::log(std::log(std::numeric_limits<Tp>::max()))
		      - Tp{1};
          auto h = k_max / n;

          auto sum = func(Tp{0});
          decltype(sum) sum1{}, sum2{};
          for (int k = -n; k < 0; ++k)
	    {
	      const auto u = h * Tp(k);
	      const auto eu = std::exp(u);
	      const auto cosh = eu + Tp{1} / eu;
	      const auto sinh = eu - Tp{1} / eu;
              const auto esh = std::exp(s_pi_4 * sinh);
	      const auto x = (esh - Tp{1} / esh) / Tp{2};
	      const auto w = esh + Tp{1} / esh;
	      const auto dxdu = cosh * w / Tp{4};
	      sum1 += dxdu * func(+x);
	      sum2 += dxdu * func(-x);
	    }

          auto prev_sum = sum + sum1 + sum2;
          for (int iter = 0; iter < max_iter; ++iter)
	    {
	      for (int k  = -n; k < 0; ++k)
	        {
	          const auto u = h * Tp(k + 0.5);
	          const auto eu = std::exp(u);
	          const auto cosh = eu + Tp{1} / eu;
	          const auto sinh = eu - Tp{1} / eu;
                  const auto esh = std::exp(s_pi_4 * sinh);
	          const auto x = (esh - Tp{1} / esh) / Tp{2};
	          const auto w = esh + Tp{1} / esh;
	          const auto dxdu = cosh * w / Tp{4};
	          sum1 += dxdu * func(+x);
	          sum2 += dxdu * func(-x);
	        }

	      n *= 2;
	      h /= Tp{2};

	      const auto curr_sum = sum + sum1 + sum2;
	      if (auto abs_del = std::abs(curr_sum - prev_sum);
                  abs_del < max_abs_err
		  || abs_del < std::abs(max_rel_err * curr_sum)
		  || iter + 1 == max_iter) // Keep prev_sum even if at max_iters.
	        break;

	      prev_sum = curr_sum;
	    }

          const auto fact = Tp{2} * s_pi_4 * h;
          const auto tot_sum = sum + sum1 + sum2;
          return {fact * tot_sum,
		  fact * std::abs(tot_sum - Tp{2} * prev_sum)};
	}
    }

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
    integrate_exp_sinh(FuncTp func, Tp lower,
			Tp max_abs_err, Tp max_rel_err,
			int max_iter)
    {
      using integ_t = adaptive_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;
      using absarea_t = typename integ_t::AbsAreaTp;

      using RetTp = std::invoke_result_t<FuncTp, Tp>;
      using AreaTp = decltype(RetTp{} * Tp{});

      const auto s_pi_4 = Tp{3.141592653589793238462643383279502884195L} / 4;

      if (std::isnan(lower)
          || std::isnan(max_abs_err) || std::isnan(max_rel_err))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN, absarea_t{} * s_NaN};
	}
      else
	{
          int n = 16;

          // Find K = ln(ln(max_number))
          const auto k_max
		      = std::log(std::log(std::numeric_limits<Tp>::max()))
		      - Tp{1};
          auto h = k_max / n;

          auto sum = AreaTp{0};
          for (int k = -n; k <= n; ++k)
	    {
	      const auto u = h * Tp(k);
	      const auto eu = std::exp(u);
	      const auto cosh = eu + Tp{1} / eu;
	      const auto sinh = eu - Tp{1} / eu;
              const auto esh = std::exp(s_pi_4 * sinh);
	      const auto dxdu = cosh * esh;
	      sum += dxdu * func(lower + esh);
	    }

          // Interlace values (don't go past the rightmost point).
          auto prev_sum = sum;
          for (int iter = 0; iter < max_iter; ++iter)
	    {
	      for (int k  = -n; k < n; ++k)
	        {
	          const auto u = h * Tp(k + 0.5);
	          const auto eu = std::exp(u);
	          const auto cosh = eu + Tp{1} / eu;
	          const auto sinh = eu - Tp{1} / eu;
                  const auto esh = std::exp(s_pi_4 * sinh);
	          const auto dxdu = cosh * esh;
	          sum += dxdu * func(lower + esh);
	        }

	      n *= 2;
	      h /= Tp{2};

	      if (auto abs_del = std::abs(sum - prev_sum);
                  abs_del < max_abs_err
		  || abs_del < std::abs(max_rel_err * sum)
		  || iter + 1 == max_iter) // Keep prev_sum even if at max_iters.
	        break;

	      prev_sum = sum;
	    }

          const auto fact = s_pi_4 * h;
          return {fact * sum, fact * std::abs(sum - Tp{2} * prev_sum)};
	}
    }

} // namespace emsr

#endif // DOUBLE_EXP_INTEGRATE_TCC
