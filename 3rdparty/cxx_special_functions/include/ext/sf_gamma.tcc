#ifndef SF_GAMMA_TCC
#define SF_GAMMA_TCC 1

#pragma GCC system_header

#include <limits>
#include <cmath>
#include <complex>

#include <ext/math_util.h>
#include <ext/numeric_limits.h>

namespace std
{

// Implementation-space details.
namespace __detail
{

  /**
   * @brief Return the incomplete gamma function by series summation.
   * @f[
   *   \gamma(a,x) = x^a e^{-z}\sum_{k=1}^{\infty} \frac{x^k}{(a)_k}
   * @f]
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    __gamma_series(_Tp __a, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = __gnu_cxx::__num_traits_t<_Val>;
      const auto _S_eps = _Real{3} * __gnu_cxx::__epsilon(__a);
      const auto _S_itmax = 10 * int(10 + std::sqrt(std::abs(__a)));

      auto __lngam = __log_gamma(__a);
      auto __sign = __log_gamma_sign(__a);
      auto __ia = __gnu_cxx::__fp_is_integer(__a);

      if (__ia && __ia() <= 0)
	std::__throw_domain_error(__N("__gamma_series: "
				      "non-positive integer argument a"));
      else if (__x == _Real{0})
	return std::make_pair(_Val{0}, __lngam);
      else if (std::real(__x) < _Real{0})
	std::__throw_domain_error(__N("__gamma_series: negative argument x"));
      else
	{
	  auto __aa = __a;
	  _Val __term, __sum;
	  __term = __sum = _Tp{1} / __a;
	  for (unsigned int __n = 1; __n <= _S_itmax; ++__n)
	    {
	      __aa += _Real{1};
	      __term *= __x / __aa;
	      __sum += __term;
	      if (std::abs(__term) < _S_eps * std::abs(__sum))
		{
		  auto __gamser = std::exp(-__x + __a * std::log(__x) - __lngam)
				* __sum * __sign;
		  return std::make_pair(__gamser, __lngam);
		}
	    }
	  std::__throw_logic_error(__N("__gamma_series: "
	  			"a too large, itmax too small in routine."));
	}
    }

  /**
   * @brief Return the incomplete gamma function by continued fraction.
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    __gamma_cont_frac(_Tp __a, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = __gnu_cxx::__num_traits_t<_Val>;
      const auto _S_fpmin = _Real{3} * __gnu_cxx::__lim_min(__a);
      const auto _S_eps = _Real{3} * __gnu_cxx::__epsilon(__a);
      const auto _S_itmax = 10 * int(10 + std::sqrt(std::abs(__a)));

      auto __lngam = __log_gamma(__a);
      auto __sign = __log_gamma_sign(__a);

      auto __b = __x + _Real{1} - __a;
      auto __c = _Real{1} / _S_fpmin;
      auto __d = _Real{1} / __b;
      auto __h = __d;
      for (unsigned int __n = 1; __n <= _S_itmax; ++__n)
	{
	  auto __an = -_Real{__n} * (_Real{__n} - __a);
	  __b += _Real{2};
	  __d = __an * __d + __b;
	  if (std::abs(__d) < _S_fpmin)
	    __d = _S_fpmin;
	  __c = __b + __an / __c;
	  if (std::abs(__c) < _S_fpmin)
	    __c = _S_fpmin;
	  __d = _Real{1} / __d;
	  auto __del = __d * __c;
	  __h *= __del;
	  if (std::abs(__del - _Real{1}) < _S_eps)
	    {
	      auto __gamcf = std::exp(-__x + __a * std::log(__x) - __lngam)
			  * __h * __sign;
	      return std::make_pair(__gamcf, __lngam);
	    }
	}
      std::__throw_logic_error(__N("__gamma_cont_fraction: "
				   "a too large, itmax too small in routine."));
    }

  /**
   * @brief  Return the lower incomplete gamma function.
   * The lower incomplete gamma function is defined by
   * @f[
   *   \gamma(a,x) = \int_0^x e^{-t}t^{a-1}dt  (a > 0)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __tgamma_lower(_Tp __a, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = __gnu_cxx::__num_traits_t<_Val>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__a);

      if (std::isnan(__a) || std::isnan(__x))
	return _S_NaN;

      auto __ia = __gnu_cxx::__fp_is_integer(__a);
      if (__ia && __ia() <= 0)
	std::__throw_domain_error(__N("__tgamma_lower: "
				      "non-positive integer argument a"));
      else if (std::real(__x) < std::real(__a + _Real{1}))
	{
	  std::pair<_Tp, _Tp> __gp = __gamma_series(__a, __x);
	  return std::exp(__gp.second) * __gp.first;
	}
      else
	{
	  std::pair<_Tp, _Tp> __gp = __gamma_cont_frac(__a, __x);
	  return std::exp(__gp.second) * (_Tp{1} - __gp.first);
	}
    }

  /**
   * @brief  Return the upper incomplete gamma function.
   * The lower incomplete gamma function is defined by
   * @f[
   *   \Gamma(a,x) = \int_x^\infty e^{-t}t^{a-1}dt  (a > 0)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __tgamma(_Tp __a, _Tp __x)
    {
      using _Val = _Tp;
      using _Real = __gnu_cxx::__num_traits_t<_Val>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__a);

      if (std::isnan(__a) || std::isnan(__x))
	return _S_NaN;

      auto __ia = __gnu_cxx::__fp_is_integer(__a);
      if (__ia && __ia() <= 0)
	std::__throw_domain_error(__N("__tgamma: "
				      "non-positive integer argument a"));
      else if (std::real(__x) < std::real(__a + _Real{1}))
	{
	  auto __gp = __gamma_series(__a, __x);
	  return std::exp(__gp.second) * (_Tp{1} - __gp.first);
	}
      else
	{
	  auto __gp = __gamma_cont_frac(__a, __x);
	  return std::exp(__gp.second) * __gp.first;
	}
    }

} // namespace __detail

} // namespace std

#endif // SF_GAMMA_TCC
