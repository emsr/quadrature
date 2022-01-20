
#ifndef ULP_H
#define ULP_H 1

#include <limits>
#include <cmath>
#include <algorithm>

namespace emsr
{

  /**
   * Returns the "Unit in the Last Place" or ulp.
   * This is radix independent.
   *
   * For radix @f$ \beta @f$, exponent e, and precision p value
   * @f[
   *    x \in \left[ \beta^{e}, \beta^{e+1} \right)
   * @f]
   * @f[
   *    ulp(x) = \beta^{max(e,e_{min}) - p + 1}
   * @f[
   *
   * @see Handbook of Floating-Point Arithmetic, Muller, et. al.
   *      Birkhauser, 2010, Chapter 2.6, pp. 32-39.
   */
  template<typename _Tp>
    constexpr _Tp
    ulp(_Tp __x)
    {
      if (std::isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (std::isinf(__x))
	return std::numeric_limits<_Tp>::infinity();
      else
	{
	  const auto __exp = std::ilogb(std::abs(__x));
	  const auto __exp_min = std::numeric_limits<_Tp>::min_exponent;
	  const auto _prec = std::numeric_limits<_Tp>::digits;

	  return std::scalbn(_Tp{1}, std::max(__exp, __exp_min) - _prec + 1);
	}
    }

} // namespace emsr

#endif // ULP_H
