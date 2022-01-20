#ifndef SF_GAMMA_TCC
#define SF_GAMMA_TCC 1

#include <limits>
#include <cmath>
#include <complex>

#include <emsr/math_util.h>
#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>

namespace emsr
{

// Implementation-space details.
namespace detail
{

  /**
   * @brief Return @f$ log(|\Gamma(a)|) @f$.
   * 	    This will return values even for @f$ a < 0 @f$.
   * 	    To recover the sign of @f$ \Gamma(a) @f$ for
   * 	    any argument use @a log_gamma_sign.
   *
   * @param a The argument of the log of the gamma function.
   * @return  The logarithm of the gamma function.
   */
  template<typename Tp>
    Tp
    log_gamma(Tp a)
    {
      return std::lgamma(std::abs(a));
    }

  /**
   * @brief Return the sign of @f$ \Gamma(x) @f$.
   * At nonpositive integers zero is returned indicating @f$ \Gamma(x) @f$
   * is undefined.
   *
   * @param __a The argument of the gamma function.
   * @return  The sign of the gamma function.
   */
  template<typename Tp>
    Tp
    log_gamma_sign(Tp a)
    {
      if (a >= Tp{0})
	return Tp{1};
      else if (a == std::nearbyint(a))
	return Tp{0};
      else
	return (int(-a) % 2 == 0) ? -Tp{1} : Tp{1};
    }

  /**
   * @brief Return the incomplete gamma function by series summation.
   * @f[
   *   \gamma(a,x) = x^a e^{-z}\sum_{k=1}^{\infty} \frac{x^k}{(a)_k}
   * @f]
   */
  template<typename Tp>
    std::pair<Tp, Tp>
    gamma_series(Tp a, Tp x)
    {
      using _Val = Tp;
      using _Real = emsr::num_traits_t<_Val>;
      const auto _S_eps = _Real{3} * emsr::epsilon(a);
      const auto _S_itmax = 10 * int(10 + std::sqrt(std::abs(a)));

      auto lngam = log_gamma(a);
      auto sign = log_gamma_sign(a);
      auto ia = emsr::fp_is_integer(a);

      if (ia && ia() <= 0)
	throw std::domain_error("gamma_series: non-positive integer argument a");
      else if (x == _Real{0})
	return std::make_pair(_Val{0}, lngam);
      else if (std::real(x) < _Real{0})
	throw std::domain_error("gamma_series: negative argument x");
      else
	{
	  auto aa = a;
	  _Val term, sum;
	  term = sum = Tp{1} / a;
	  for (unsigned int n = 1; n <= _S_itmax; ++n)
	    {
	      aa += _Real{1};
	      term *= x / aa;
	      sum += term;
	      if (std::abs(term) < _S_eps * std::abs(sum))
		{
		  auto gamser = std::exp(-x + a * std::log(x) - lngam)
				* sum * sign;
		  return std::make_pair(gamser, lngam);
		}
	    }
	  throw std::logic_error("gamma_series: a too large, itmax too small in routine.");
	}
    }

  /**
   * @brief Return the incomplete gamma function by continued fraction.
   */
  template<typename Tp>
    std::pair<Tp, Tp>
    gamma_cont_frac(Tp a, Tp x)
    {
      using _Val = Tp;
      using _Real = emsr::num_traits_t<_Val>;
      const auto _S_fpmin = _Real{3} * emsr::lim_min(a);
      const auto _S_eps = _Real{3} * emsr::epsilon(a);
      const auto _S_itmax = 10 * int(10 + std::sqrt(std::abs(a)));

      auto lngam = log_gamma(a);
      auto sign = log_gamma_sign(a);

      auto b = x + _Real{1} - a;
      auto c = _Real{1} / _S_fpmin;
      auto d = _Real{1} / b;
      auto h = d;
      for (unsigned int n = 1; n <= _S_itmax; ++n)
	{
	  auto an = -_Real{n} * (_Real{n} - a);
	  b += _Real{2};
	  d = an * d + b;
	  if (std::abs(d) < _S_fpmin)
	    d = _S_fpmin;
	  c = b + an / c;
	  if (std::abs(c) < _S_fpmin)
	    c = _S_fpmin;
	  d = _Real{1} / d;
	  auto del = d * c;
	  h *= del;
	  if (std::abs(del - _Real{1}) < _S_eps)
	    {
	      auto gamcf = std::exp(-x + a * std::log(x) - lngam)
			  * h * sign;
	      return std::make_pair(gamcf, lngam);
	    }
	}
      throw std::logic_error("gamma_cont_fraction: a too large, itmax too small in routine.");
    }

  /**
   * @brief  Return the lower incomplete gamma function.
   * The lower incomplete gamma function is defined by
   * @f[
   *   \gamma(a,x) = \int_0^x e^{-t}t^{a-1}dt  (a > 0)
   * @f]
   */
  template<typename Tp>
    Tp
    tgamma_lower(Tp a, Tp x)
    {
      using _Val = Tp;
      using _Real = emsr::num_traits_t<_Val>;
      const auto _S_NaN = emsr::quiet_NaN(a);

      if (std::isnan(a) || std::isnan(x))
	return _S_NaN;

      auto ia = emsr::fp_is_integer(a);
      if (ia && ia() <= 0)
	throw std::domain_error("tgamma_lower: non-positive integer argument a");
      else if (std::real(x) < std::real(a + _Real{1}))
	{
	  std::pair<Tp, Tp> gp = gamma_series(a, x);
	  return std::exp(gp.second) * gp.first;
	}
      else
	{
	  std::pair<Tp, Tp> gp = gamma_cont_frac(a, x);
	  return std::exp(gp.second) * (Tp{1} - gp.first);
	}
    }

  /**
   * @brief  Return the upper incomplete gamma function.
   * The lower incomplete gamma function is defined by
   * @f[
   *   \Gamma(a,x) = \int_x^\infty e^{-t}t^{a-1}dt  (a > 0)
   * @f]
   */
  template<typename Tp>
    Tp
    tgamma(Tp a, Tp x)
    {
      using _Val = Tp;
      using _Real = emsr::num_traits_t<_Val>;
      const auto _S_NaN = emsr::quiet_NaN(a);

      if (std::isnan(a) || std::isnan(x))
	return _S_NaN;

      auto ia = emsr::fp_is_integer(a);
      if (ia && ia() <= 0)
	throw std::domain_error("tgamma: non-positive integer argument a");
      else if (std::real(x) < std::real(a + _Real{1}))
	{
	  auto gp = gamma_series(a, x);
	  return std::exp(gp.second) * (Tp{1} - gp.first);
	}
      else
	{
	  auto gp = gamma_cont_frac(a, x);
	  return std::exp(gp.second) * gp.first;
	}
    }

} // namespace detail

} // namespace emsr

#endif // SF_GAMMA_TCC
