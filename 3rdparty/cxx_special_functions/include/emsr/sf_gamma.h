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

} // namespace emsr

#endif // SF_GAMMA_H
