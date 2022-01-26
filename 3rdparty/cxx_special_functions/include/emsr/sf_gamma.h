#ifndef SF_GAMMA_H
#define SF_GAMMA_H 1

#include <emsr/fp_type_util.h>
#include <emsr/sf_gamma.tcc>

namespace emsr
{

  /**
   * Return the upper incomplete gamma function @f$ \Gamma(a,x) @f$.
   * The (upper) incomplete gamma function is defined by
   * @f[
   *   \Gamma(a,x) = \int_x^\infty t^{a-1}e^{-t}dt
   * @f]
   */
  template<typename _Ta, typename _Tp>
    inline emsr::fp_promote_t<_Ta, _Tp>
    tgamma(_Ta a, _Tp x)
    {
      using type = emsr::fp_promote_t<_Ta, _Tp>;
      return emsr::detail::tgamma<type>(a, x);
    }

  /**
   * Return the lower incomplete gamma function @f$ \gamma(a,x) @f$.
   * The lower incomplete gamma function is defined by
   * @f[
   *   \gamma(a,x) = \int_0^x t^{a-1}e^{-t}dt
   * @f]
   */
  template<typename _Ta, typename _Tp>
    inline emsr::fp_promote_t<_Ta, _Tp>
    tgamma_lower(_Ta a, _Tp x)
    {
      using type = emsr::fp_promote_t<_Ta, _Tp>;
      return emsr::detail::tgamma_lower<type>(a, x);
    }

  /**
   * Return the digamma or psi function of argument @c x.
   *
   * The the digamma or psi function is defined by
   * @f[
   *    \psi(x) = \frac{d}{dx}log\left(\Gamma(x)\right)
   *            = \frac{\Gamma'(x)}{\Gamma(x)},
   * @f]
   * the logarithmic derivative of the gamma function.
   *
   * @param x The parameter
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    digamma(_Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::digamma<type>(x);
    }

  /**
   * Return the harmonic number @f$ H_n @f$.
   *
   * The the harmonic number is defined by
   * @f[
   *    H_n = \sum_{k=1}^{n}\frac{1}{k}
   * @f]
   *
   * @param n The parameter
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    harmonic(unsigned int n)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::harmonic_number<type>(n);
    }

} // namespace emsr

#endif // SF_GAMMA_H
