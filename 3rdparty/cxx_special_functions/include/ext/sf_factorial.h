#ifndef SF_FACTORIAL_H
#define SF_FACTORIAL_H 1

#pragma GCC system_header

#include <ext/fp_type_util.h>
#include <ext/sf_factorial.tcc>

namespace __gnu_cxx
{

  /**
   * @brief Return the factorial @f$ n! @f$ of the argument as a real number.
   * @f[
   *   n! = 1 \times 2 \times ... \times n, 0! = 1
   * @f]
   */
  template<typename _Tp>
    inline __gnu_cxx::fp_promote_t<_Tp>
    factorial(unsigned int __n)
    {
      using __type = __gnu_cxx::fp_promote_t<_Tp>;
      return std::__detail::__factorial<__type>(__n);
    }

}

#endif // SF_FACTORIAL_H
