/* quadrature/double_exp_integrate.tcc
 *
 * Copyright (C) 2017-2019 Free Software Foundation, Inc.
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

#include <cmath>

namespace __gnu_cxx
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
  template<typename _Tp, typename _FuncTp>
    adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_tanh_sinh(_FuncTp __func, _Tp __a, _Tp __b,
			_Tp __max_abs_err, _Tp __max_rel_err,
			int __max_iter)
    {
      const auto _S_pi_4 = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L} / 4;

      int __n = 16;
      __n /= 2;

      // Find K = ln(ln(max_number))
      const auto __k_max = std::log(std::log(std::numeric_limits<_Tp>::max()))
		- _Tp{1};
      auto __h = __k_max / __n;

      auto __sum = __func((__a + __b) / _Tp{2}) / _Tp{2};
      decltype(__sum) __sum1{}, __sum2{};
      for (int __k = -__n; __k < 0; ++__k)
	{
	  const auto __u = __h * _Tp(__k);
	  const auto __eu = std::exp(__u);
	  const auto __cosh = __eu + _Tp{1} / __eu;
	  const auto __sinh = __eu - _Tp{1} / __eu;
          const auto __esh = std::exp(_S_pi_4 * __sinh);
	  const auto __w = __esh + _Tp{1} / __esh;
	  const auto __dxdu = __cosh / (__w * __w);
	  const auto __x1 = (__b * __esh + __a / __esh) / __w;
	  if (__x1 != __a && __x1 != __b) 
	    __sum1 += __dxdu * __func(__x1);
	  const auto __x2 = (__a * __esh + __b / __esh) / __w;
	  if (__x2 != __a && __x2 != __b)
	    __sum2 += __dxdu * __func(__x2);
	}

      // Interlace values; don't go past the rightmost point.
      decltype(__sum) __prev_sum{};
      for (int __iter = 0; __iter < __max_iter; ++__iter)
	{
	  __prev_sum = __sum + __sum1 + __sum2;
	  for (int __k  = -__n; __k < 0; ++__k)
	    {
	      const auto __u = __h * _Tp(__k + 0.5);
	      const auto __eu = std::exp(__u);
	      // A standard sinhcosh would be a nice idea along with sincos.
	      const auto __cosh = __eu + _Tp{1} / __eu;
	      const auto __sinh = __eu - _Tp{1} / __eu;
              const auto __esh = std::exp(_S_pi_4 * __sinh);
	      const auto __w = __esh + _Tp{1} / __esh;
	      const auto __dxdu = __cosh / (__w * __w);
	      // natural: x1 = (s - 1/s) / (s + 1/s)
	      const auto __x1 = (__b * __esh + __a / __esh) / __w;
	      if (__x1 != __a && __x1 != __b) 
		__sum1 += __dxdu * __func(__x1);
	      // natural: x2 = (-s + 1/s) / (s + 1/s)
	      const auto __x2 = (__a * __esh + __b / __esh) / __w;
	      if (__x2 != __a && __x2 != __b)
		__sum2 += __dxdu * __func(__x2);
	    }

	  const auto __curr_sum = __sum + __sum1 + __sum2;
	  if (auto __abs_del = std::abs(__curr_sum - __prev_sum);
              __abs_del < __max_abs_err
	      || __abs_del < std::abs(__max_rel_err * __curr_sum))
	    break;

	  __n *= 2;
	  __h /= _Tp{2};
	}

      const auto __fact = _Tp{2} * (__b - __a) * _S_pi_4 * __h;
      const auto __tot_sum = __sum + __sum1 + __sum2;
      return {__fact * __tot_sum, __fact * std::abs(__tot_sum - __prev_sum)};
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
  template<typename _Tp, typename _FuncTp>
    adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_sinh_sinh(_FuncTp __func,
			_Tp __max_abs_err, _Tp __max_rel_err,
			int __max_iter)
    {
      const auto _S_pi_4 = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L} / 4;

      int __n = 16;
      __n /= 2;

      // Find K = ln(ln(max_number))
      const auto __k_max = std::log(std::log(std::numeric_limits<_Tp>::max()))
		- _Tp{1};
      auto __h = __k_max / __n;

      auto __sum = __func(_Tp{0});
      decltype(__sum) __sum1{}, __sum2{};
      for (int __k = -__n; __k < 0; ++__k)
	{
	  const auto __u = __h * _Tp(__k);
	  const auto __eu = std::exp(__u);
	  const auto __cosh = __eu + _Tp{1} / __eu;
	  const auto __sinh = __eu - _Tp{1} / __eu;
          const auto __esh = std::exp(_S_pi_4 * __sinh);
	  const auto __x = (__esh - _Tp{1} / __esh) / _Tp{2};
	  const auto __w = __esh + _Tp{1} / __esh;
	  const auto __dxdu = __cosh * __w / _Tp{4};
	  __sum1 += __dxdu * __func(+__x);
	  __sum2 += __dxdu * __func(-__x);
	}

      decltype(__sum) __prev_sum{};
      for (int __iter = 0; __iter < __max_iter; ++__iter)
	{
	  __prev_sum = __sum + __sum1 + __sum2;
	  for (int __k  = -__n; __k < 0; ++__k)
	    {
	      const auto __u = __h * _Tp(__k + 0.5);
	      const auto __eu = std::exp(__u);
	      const auto __cosh = __eu + _Tp{1} / __eu;
	      const auto __sinh = __eu - _Tp{1} / __eu;
              const auto __esh = std::exp(_S_pi_4 * __sinh);
	      const auto __x = (__esh - _Tp{1} / __esh) / _Tp{2};
	      const auto __w = __esh + _Tp{1} / __esh;
	      const auto __dxdu = __cosh * __w / _Tp{4};
	      __sum1 += __dxdu * __func(+__x);
	      __sum2 += __dxdu * __func(-__x);
	    }

	  const auto __curr_sum = __sum + __sum1 + __sum2;
	  if (auto __abs_del = std::abs(__curr_sum - __prev_sum);
              __abs_del < __max_abs_err
	      || __abs_del < std::abs(__max_rel_err * __curr_sum))
	    break;

	  __n *= 2;
	  __h /= _Tp{2};
	}

      const auto __fact = _Tp{2} * _S_pi_4 * __h;
      const auto __tot_sum = __sum + __sum1 + __sum2;
      return {__fact * __tot_sum, __fact * std::abs(__tot_sum - __prev_sum)};
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
  template<typename _Tp, typename _FuncTp>
    adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_exp_sinh(_FuncTp __func, _Tp __a,
			_Tp __max_abs_err, _Tp __max_rel_err,
			int __max_iter)
    {
      using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
      using _AreaTp = decltype(_RetTp{} * _Tp{});

      int __n = 16;

      const auto _S_pi_4 = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L} / 4;

      // Find K = ln(ln(max_number))
      const auto __k_max = std::log(std::log(std::numeric_limits<_Tp>::max()))
		- _Tp{1};
      auto __h = __k_max / __n;

      auto __sum = _AreaTp{0};
      for (int __k = -__n; __k <= __n; ++__k)
	{
	  const auto __u = __h * _Tp(__k);
	  const auto __eu = std::exp(__u);
	  const auto __cosh = __eu + _Tp{1} / __eu;
	  const auto __sinh = __eu - _Tp{1} / __eu;
          const auto __esh = std::exp(_S_pi_4 * __sinh);
	  const auto __dxdu = __cosh * __esh;
	  __sum += __dxdu * __func(__a + __esh);
	}

      // Interlace values (don't go past the rightmost point).
      decltype(__sum) __prev_sum{};
      for (int __iter = 0; __iter < __max_iter; ++__iter)
	{
	  __prev_sum = __sum;
	  for (int __k  = -__n; __k < __n; ++__k)
	    {
	      const auto __u = __h * _Tp(__k + 0.5);
	      const auto __eu = std::exp(__u);
	      const auto __cosh = __eu + _Tp{1} / __eu;
	      const auto __sinh = __eu - _Tp{1} / __eu;
              const auto __esh = std::exp(_S_pi_4 * __sinh);
	      const auto __dxdu = __cosh * __esh;
	      __sum += __dxdu * __func(__a + __esh);
	    }

	  const auto __curr_sum = __sum;
	  if (auto __abs_del = std::abs(__curr_sum - __prev_sum);
              __abs_del < __max_abs_err
	      || __abs_del < std::abs(__max_rel_err * __curr_sum))
	    break;

	  __n *= 2;
	  __h /= _Tp{2};
	}

      const auto __fact = _S_pi_4 * __h;
      const auto __tot_sum = __sum;
      return {__fact * __tot_sum, __fact * std::abs(__tot_sum - __prev_sum)};
    }

} // namespace __gnu_cxx

#endif // DOUBLE_EXP_INTEGRATE_TCC
