#ifndef SF_FACTORIAL_H
#define SF_FACTORIAL_H 1

#include <emsr/fp_type_util.h>
#include <emsr/sf_factorial.tcc>

namespace emsr
{

  /**
   * @brief Return the factorial @f$ n! @f$ of the argument as a real number.
   * @f[
   *   n! = 1 \times 2 \times ... \times n, 0! = 1
   * @f]
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    factorial(unsigned int n)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::factorial<type>(n);
    }

} // namespace emsr

#endif // SF_FACTORIAL_H
