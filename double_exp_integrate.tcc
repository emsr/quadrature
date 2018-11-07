/* quadrature/double_exp_integrate.tcc
 *
 * Copyright (C) 2017-2018 Free Software Foundation, Inc.
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

#include <ext/cmath>

namespace __gnu_cxx
{

  /**
   * @f[
   *    \int_{-1}^{+1}f(x)dx
   * @f]
   * Making the change of variables:
   * @f[
   *    x = tanh\left[\frac{\pi}{2}sinh(u)\right],
   *   dx = \frac{\pi}{2}\frac{cosh(u)}{cosh^2\left[\frac{\pi}{2}sinh(u)\right]}du
   * @f]
   * gives the following integral:
   * @f[
   *    \int_{-\infty}^{+\infty}f(tanh\left[\frac{\pi}{2}sinh(u)\right])
   *     
   * @f]
   * 
   * @f[
   *    
   * @f]
   */
  template<typename _Tp, typename _FuncTp>
    _Tp
    tanh_sinh_integrate1(_FuncTp __func, int __n, _Tp __a, _Tp __b, _Tp /*__max_rel_tol*/)
    {
      __n /= 2;
      const auto __h = _Tp{5} / __n; // 5.0 is rough limit of K in exp(exp(K)).

      const auto _S_pi_4 = __gnu_cxx::__const_pi_quarter<_Tp>();
      auto __sum = _Tp{0};
      for (int __k = -__n; __k <= __n; ++__k)
	{
	  const auto __u = __h * _Tp(__k);
	  const auto __eu = std::exp(__u);
	  const auto __cosh = __eu + _Tp{1} / __eu;
	  const auto __sinh = __eu - _Tp{1} / __eu;
          const auto __s = std::exp(_S_pi_4 * __sinh);
	  const auto __w = __s + _Tp{1} / __s;
	  const auto __x = (__b * __s + __a / __s) / __w;
	  if (__x != __a && __x != __b)
	    {
	      const auto __dxdu = _Tp{2} * (__b - __a) * _S_pi_4 * __cosh
	    	       / (__w * __w);
	      __sum += __h * __func(__x) * __dxdu;
	    }
	}
      return __sum;
    }

  /**
   * This uses the symmetry of cosh.
   */
  template<typename _Tp, typename _FuncTp>
    _Tp
    tanh_sinh_integrate2(_FuncTp __func, int __n, _Tp __a, _Tp __b, _Tp /*__max_rel_tol*/)
    {
      __n /= 2;
      const auto __h = _Tp{5} / __n; // 5.0 is rough limit of K in exp(exp(K)).

      const auto _S_pi_4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;
      auto __sum = __func((__a + __b) / _Tp{2}) / _Tp{2};
      for (int __k = -__n; __k < 0; ++__k)
	{
	  const auto __u = __h * _Tp(__k);
	  const auto __eu = std::exp(__u);
	  const auto __cosh = __eu + _Tp{1} / __eu;
	  const auto __sinh = __eu - _Tp{1} / __eu;
          const auto __s = std::exp(_S_pi_4 * __sinh);
	  const auto __w = __s + _Tp{1} / __s;
	  const auto __dxdu = __cosh / (__w * __w);
	  const auto __x1 = (__b * __s + __a / __s) / __w;
	  const auto __x2 = (__a * __s + __b / __s) / __w;
	  if (__x1 != __a && __x1 != __b) 
	    __sum += __dxdu * __func(__x1);
	  if (__x2 != __a && __x2 != __b)
	    __sum += __dxdu * __func(__x2);
	}
      return _Tp{2} * (__b - __a) * _S_pi_4 * __h * __sum;
    }

  /**
   * This is like the first except it interlaces another set of points
   * for a higher-order integral..
   */
  template<typename _Tp, typename _FuncTp>
    adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    tanh_sinh_integrate3(_FuncTp __func, int __n, _Tp __a, _Tp __b, _Tp /*__max_rel_tol*/)
    {
      __n /= 2;
      const auto __h = _Tp{5} / __n; // 5.0 is rough limit of K in exp(exp(K)).
      // @todo Find K = ln(ln(max()))

      const auto _S_pi_4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;
      auto __sum = _Tp{0};
      for (int __k = -__n; __k <= __n; ++__k)
	{
	  const auto __u = __h * _Tp(__k);
	  const auto __eu = std::exp(__u);
	  const auto __cosh = __eu + _Tp{1} / __eu;
	  const auto __sinh = __eu - _Tp{1} / __eu;
          const auto __s = std::exp(_S_pi_4 * __sinh);
	  const auto __w = __s + _Tp{1} / __s;
	  const auto __x = (__b * __s + __a / __s) / __w;
	  if (__x != __a && __x != __b)
	    {
	      const auto __dxdu = __cosh / (__w * __w);
	      __sum += __func(__x) * __dxdu; 
	    }
	}

      // Interlace values (don't go past the rightmost point).
      const auto __sum_prev = __sum;
      for (int __k = -__n; __k < __n; ++__k)
	{
	  const auto __u = __h * (_Tp(__k) + 0.5);
	  const auto __eu = std::exp(__u);
	  const auto __cosh = __eu + _Tp{1} / __eu;
	  const auto __sinh = __eu - _Tp{1} / __eu;
          const auto __s = std::exp(_S_pi_4 * __sinh);
	  const auto __w = __s + _Tp{1} / __s;
	  const auto __x = (__b * __s + __a / __s) / __w;
	  if (__x != __a && __x != __b)
	    {
	      const auto __dxdu = __cosh / (__w * __w);
	      __sum += __func(__x) * __dxdu; 
	    }
	}

      const auto __fact = __h * (__b - __a) * _S_pi_4;
      return {__fact * __sum, __fact * std::abs(__sum - __sum_prev)};
  }

  /**
   * @f[
   *    \int_{-1}^{+1}f(x)dx
   * @f]
   * Making the change of variables:
   * @f[
   *    x = tanh\left[\frac{\pi}{2}sinh(u)\right],
   *   dx = \frac{\pi}{2}\frac{cosh(u)}{cosh^2\left[\frac{\pi}{2}sinh(u)\right]}du
   * @f]
   * gives the following integral:
   * @f[
   *    \int_{-\infty}^{+\infty}f(tanh\left[\frac{\pi}{2}sinh(u)\right])
   *     = \sum_{k=-n}^{+n} 
   * @f]
   */
  template<typename _Tp, typename _FuncTp>
    adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    tanh_sinh_integrate(_FuncTp __func, int __n, _Tp __a, _Tp __b,
			_Tp __max_rel_tol, int __max_iter = 3)
    {
      const auto _S_pi_4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;
      __n /= 2;
      // 5.0 is rough limit of K in exp(exp(K)).
      // Find K = ln(ln(max()))
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
          const auto __s = std::exp(_S_pi_4 * __sinh);
	  const auto __w = __s + _Tp{1} / __s;
	  const auto __dxdu = __cosh / (__w * __w);
	  const auto __x1 = (__b * __s + __a / __s) / __w;
	  if (__x1 != __a && __x1 != __b) 
	    __sum1 += __dxdu * __func(__x1);
	  const auto __x2 = (__a * __s + __b / __s) / __w;
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
              const auto __s = std::exp(_S_pi_4 * __sinh);
	      const auto __w = __s + _Tp{1} / __s;
	      const auto __dxdu = __cosh / (__w * __w);
	      const auto __x1 = (__b * __s + __a / __s) / __w;
	      if (__x1 != __a && __x1 != __b) 
		__sum1 += __dxdu * __func(__x1);
	      const auto __x2 = (__a * __s + __b / __s) / __w;
	      if (__x2 != __a && __x2 != __b)
		__sum2 += __dxdu * __func(__x2);
	    }

	  const auto __curr_sum = __sum + __sum1 + __sum2;
	  if (std::abs(__curr_sum - __prev_sum)
	    < std::abs(__max_rel_tol * __curr_sum))
	    break;

	  __n *= 2;
	  __h /= _Tp{2};
	}

      const auto __fact = _Tp{2} * (__b - __a) * _S_pi_4 * __h;
      const auto __tot_sum = __sum + __sum1 + __sum2;
      return {__fact * __tot_sum, __fact * std::abs(__tot_sum - __prev_sum)};
    }

} // namespace __gnu_cxx

#endif // DOUBLE_EXP_INTEGRATE_TCC
