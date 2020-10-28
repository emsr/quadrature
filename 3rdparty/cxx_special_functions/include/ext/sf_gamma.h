#ifndef SF_GAMMA_H
#define SF_GAMMA_H 1

#pragma GCC system_header

#include <ext/fp_type_util.h>
#include <ext/sf_gamma.tcc>

namespace __gnu_cxx
{

  /**
   * Return the upper incomplete gamma function @f$ \Gamma(a,x) @f$.
   * The (upper) incomplete gamma function is defined by
   * @f[
   *   \Gamma(a,x) = \int_x^\infty t^{a-1}e^{-t}dt
   * @f]
   */
  template<typename _Ta, typename _Tp>
    inline __gnu_cxx::fp_promote_t<_Ta, _Tp>
    tgamma(_Ta __a, _Tp __x)
    {
      using __type = __gnu_cxx::fp_promote_t<_Ta, _Tp>;
      return std::__detail::__tgamma<__type>(__a, __x);
    }

  /**
   * Return the lower incomplete gamma function @f$ \gamma(a,x) @f$.
   * The lower incomplete gamma function is defined by
   * @f[
   *   \gamma(a,x) = \int_0^x t^{a-1}e^{-t}dt
   * @f]
   */
  template<typename _Ta, typename _Tp>
    inline __gnu_cxx::fp_promote_t<_Ta, _Tp>
    tgamma_lower(_Ta __a, _Tp __x)
    {
      using __type = __gnu_cxx::fp_promote_t<_Ta, _Tp>;
      return std::__detail::__tgamma_lower<__type>(__a, __x);
    }

}

#endif // SF_GAMMA_H
