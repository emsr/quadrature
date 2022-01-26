#ifndef SF_GAMMA_TCC
#define SF_GAMMA_TCC 1

#include <limits>
#include <cmath>
#include <complex>

#include <emsr/math_util.h>
#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>
#include <emsr/sf_trig.h>
#include <emsr/sf_bernoulli.h>

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
      const unsigned int _S_itmax = 10 * int(10 + std::sqrt(std::abs(a)));

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
      const unsigned int _S_itmax = 10 * int(10 + std::sqrt(std::abs(a)));

      auto lngam = log_gamma(a);
      auto sign = log_gamma_sign(a);

      auto b = x + _Real{1} - a;
      auto c = _Real{1} / _S_fpmin;
      auto d = _Real{1} / b;
      auto h = d;
      for (unsigned int n = 1; n <= _S_itmax; ++n)
	{
	  auto an = -_Real(n) * (_Real(n) - a);
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

  /**
   * @brief  Return the digamma function of integral argument.
   * The digamma or @f$ \psi(x) @f$ function is defined as the logarithmic
   * derivative of the gamma function:
   * @f[
   *   \psi(x) = \frac{\Gamma'(x)}{\Gamma(x)}
   * @f]
   * The digamma series for integral argument is given by:
   * @f[
   *   \psi(n) = -\gamma_E + \sum_{k=1}^{n-1} \frac{1}{k}
   * @f]
   * The latter sum is called the harmonic number, @f$ H_n @f$.
   */
  template<typename Tp>
    Tp
    digamma(unsigned int n)
    {
      using Val = emsr::num_traits_t<Tp>;
      constexpr auto s_gamma_E = emsr::egamma_v<Val>;
      if (n > 1)
	return -s_gamma_E + harmonic_number<Val>(n - 1);
      else
	return -s_gamma_E;
    }

  /**
   * @brief  Return the digamma function by series expansion.
   * The digamma or @f$ \psi(x) @f$ function is defined by
   * @f[
   *   \psi(x) = \frac{\Gamma'(x)}{\Gamma(x)}
   * @f]
   *
   * The series is given by:
   * @f[
   *   \psi(x) = -\gamma_E - \frac{1}{x}
   *		\sum_{k=1}^{\infty} \frac{x - 1}{(k + 1)(x + k)}
   * @f]
   */
  template<typename Tp>
    Tp
    digamma_series(Tp x)
    {
      using Val = emsr::num_traits_t<Tp>;
      Tp sum = -emsr::egamma_v<Val>;
      const unsigned int s_max_iter = 100000;
      for (unsigned int k = 0; k < s_max_iter; ++k)
	{
	  const auto term = (x - Tp{1})
			    / (Tp(k + 1) * (Tp(k) + x));
	  sum += term;
	  if (std::abs(term) < emsr::epsilon<Val>())
	    break;
	}
      return sum;
    }


  /**
   * @brief  Return the digamma function for large argument.
   * The digamma or @f$ \psi(x) @f$ function is defined by
   * @f[
   *   \psi(x) = \frac{\Gamma'(x)}{\Gamma(x)}
   * @f]
   *
   * The asymptotic series is given by:
   * @f[
   *   \psi(x) = \ln(x) - \frac{1}{2x}
   *	       - \sum_{n=1}^{\infty} \frac{B_{2n}}{2 n x^{2n}}
   * @f]
   */
  template<typename Tp>
    Tp
    digamma_asymp(Tp x)
    {
      using Val = emsr::num_traits_t<Tp>;
      auto sum = std::log(x) - Val{0.5L} / x;
      const auto xx = x * x;
      auto xp = xx;
      const unsigned int max_iter = 100;
      for (unsigned int k = 1; k < max_iter; ++k)
	{
	  const Tp term = bernoulli<Val>(2 * k)
			   / (Val(2 * k) * xp);
	  sum -= term;
	  if (std::abs(term / sum) < emsr::epsilon<Val>())
	    break;
	  xp *= xx;
	}
      return sum;
    }


  /**
   * @brief  Return the digamma function.
   * The digamma or @f$ \psi(x) @f$ function is defined by
   * @f[
   *   \psi(x) = \frac{\Gamma'(x)}{\Gamma(x)}
   * @f]
   * For negative argument the reflection formula is used:
   * @f[
   *   \psi(x) = \psi(1-x) - \pi \cot(\pi x)
   * @f]
   */
  template<typename Tp>
    Tp
    digamma(Tp x)
    {
      using Val = emsr::num_traits_t<Tp>;
      const auto s_x_asymp = Val{20};
      const auto s_gamma_E = emsr::egamma_v<Val>;
      const auto s_2_ln_2 = Tp{2} * emsr::ln2_v<Val>;
      const auto s_pi = emsr::pi_v<Val>;

      const auto n = emsr::fp_is_integer(x);
      const auto m = emsr::fp_is_half_odd_integer(x);
      if (std::real(x) <= Val{0})
	{
	  if (n)
	    return emsr::quiet_NaN(x);
	  else
	    return digamma(Val{1} - x) - s_pi / emsr::tan_pi(x);
	}
      else if (n)
	return digamma<Tp>(n());
      else if (m)
	{
	  Tp sum = -s_gamma_E - s_2_ln_2;
	  for (int k = 1; k <= m(); ++k)
	    sum += Tp{2} / Tp(2 * k - 1);
	  return sum;
	}
      else if (std::real(x) > s_x_asymp)
	return digamma_asymp(x);
      else
	{
	  // The series does not converge quickly enough.
	  // Reflect to larger argument and use asymptotic expansion.
	  auto w = Tp{0};
	  auto y = x;
	  while (std::real(y) <= s_x_asymp)
	    {
	      w += Val{1} / y;
	      y += Val{1};
	    }
	  return digamma_asymp(y) - w;
	}
    }

  constexpr unsigned long long
  s_num_harmonic_numer = 29;
  constexpr unsigned long long
  s_harmonic_numer[s_num_harmonic_numer]
  {
    1ULL,
    3ULL,
    11ULL,
    25ULL,
    137ULL,
    49ULL,
    363ULL,
    761ULL,
    7129ULL,
    7381ULL,
    83711ULL,
    86021ULL,
    1145993ULL,
    1171733ULL,
    1195757ULL,
    2436559ULL,
    42142223ULL,
    14274301ULL,
    275295799ULL,
    55835135ULL,
    18858053ULL,
    19093197ULL,
    444316699ULL,
    1347822955ULL,
    34052522467ULL,
    34395742267ULL,
    312536252003ULL,
    315404588903ULL,
    9227046511387ULL
  };
  constexpr unsigned long long
  s_harmonic_denom[s_num_harmonic_numer]
  {
    1ULL,
    2ULL,
    6ULL,
    12ULL,
    60ULL,
    20ULL,
    140ULL,
    280ULL,
    2520ULL,
    2520ULL,
    27720ULL,
    27720ULL,
    360360ULL,
    360360ULL,
    360360ULL,
    720720ULL,
    12252240ULL,
    4084080ULL,
    77597520ULL,
    15519504ULL,
    5173168ULL,
    5173168ULL,
    118982864ULL,
    356948592ULL,
    8923714800ULL,
    8923714800ULL,
    80313433200ULL,
    80313433200ULL,
    2329089562800ULL
  };

  template<typename Tp>
    Tp
    harmonic_number(unsigned int n)
    {
      if (n <= s_num_harmonic_numer)
	return Tp(s_harmonic_numer[n - 1])
	     / Tp(s_harmonic_denom[n - 1]);
      else
        {
	  unsigned int k = s_num_harmonic_numer - 1;
	  auto _H_k = Tp(s_harmonic_numer[k]) / Tp(s_harmonic_denom[k]);
	  for (k = s_num_harmonic_numer + 1; k <= n; ++k)
	    _H_k += Tp{1} / Tp(k);
	  return _H_k;
	}
    }

} // namespace detail

} // namespace emsr

#endif // SF_GAMMA_TCC
