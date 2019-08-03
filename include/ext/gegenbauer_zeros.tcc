#ifndef GEGENBAUER_ZEROS_TCC
#define GEGENBAUER_ZEROS_TCC 1

namespace __gnu_cxx
{

  /**
   * Return a vector containing the zeros of the Gegenbauer or ultraspherical
   * polynomial @f$ C_n^{(\lambda)}@f$.
   * This works for @f$ \lambda > -1/2 @f$
   *
   * @tparam  _Tp  The real type of the order
   * @param[in]  __n  The degree of the Gegenbauer polynomial
   * @param[in]  __lambda  The order of the Gegenbauer polynomial
   */
  template<typename _Tp>
    std::vector<__quadrature_point_t<_Tp>>
    __gegenbauer_zeros(unsigned int __n, _Tp __lambda)
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const unsigned int _S_maxit = 1000u;
      std::vector<__quadrature_point_t<_Tp>> __pt(__n);

      _Tp __z;
      _Tp __w;
      for (auto __i = 1u; __i <= __n; ++__i)
	{
	  if (__i == 1)
	    {
	      auto __an = __lambda / __n;
	      auto __an2 = __an * __an;
	      auto __r1 = (1.0 + __lambda) * (2.78 / (4.0 + __n * __n)
			+ 0.768 * __an / __n);
	      auto __r2 = 1.0 + 1.48 * __an + 0.96 * __an + 1.282 * __an2;
	      __z = 1.0 - __r1 / __r2;
	    }
	  else if (__i == 2)
	    {
	      auto __r1 = (4.1 + __lambda)
			/ ((1.0 + __lambda) * (1.0 + 0.156 * __lambda));
	      auto __r2 = 1.0
			+ 0.06 * (__n - 8.0) * (1.0 + 0.12 * __lambda) / __n;
	      auto __r3 = 1.0
		   + 0.012 * __lambda * (1.0 + 0.25 * std::abs(__lambda)) / __n;
	      __z -= (1.0 - __z) * __r1 * __r2 * __r3;
	    }
	  else if (__i == 3)
	    {
	      auto __r1 = (1.67 + 0.28 * __lambda) / (1.0 + 0.37 * __lambda);
	      auto __r2 = 1.0 + 0.22 * (__n - 8.0) / __n;
	      auto __r3 = 1.0 + 8.0 *__lambda / ((6.28 + __lambda) * __n * __n);
	      __z -= (__pt[0].__point - __z) * __r1 * __r2 * __r3;
	    }
	  else if (__i == __n - 1)
	    {
	      auto __r1 = (1.0 + 0.235 * __lambda) / (0.766 + 0.119 * __lambda);
	      auto __r2 = 1.0 / (1.0 + 0.639 * (__n - 4.0)
						/ (1.0 + 0.71 * (__n - 4.0)));
	      auto __r3 = 1.0 / (1.0 + 20.0 * __lambda
				/ ((7.5 + __lambda) * __n * __n));
	      __z += (__z - __pt[__n - 4].__point) * __r1 * __r2 * __r3;
	    }
	  else if (__i == __n)
	    {
	      auto __r1 = (1.0 + 0.37 * __lambda) / (1.67 + 0.28 * __lambda);
	      auto __r2 = 1.0 / (1.0 + 0.22 * (__n - 8.0) / __n);
	      auto __r3 = 1.0 / (1.0 + 8.0 * __lambda
				 / ((6.28 + __lambda) * __n * __n));
	      __z += (__z - __pt[__n - 3].__point) * __r1 * __r2 * __r3;
	    }
	  else
	    __z = 3.0 * __pt[__i - 2].__point
		- 3.0 * __pt[__i - 3].__point + __pt[__i - 4].__point;

	  auto __2lambda = _Tp{2} * __lambda;
	  for (auto __its = 1u; __its <= _S_maxit; ++__its)
	    {
	      auto __temp = _Tp{2} + __2lambda;
	      auto __C1 = (__temp * __z) / _Tp{2};
	      auto __C2 = _Tp{1};
	      for (auto __j = 2u; __j <= __n; ++__j)
		{
		  auto __C3 = __C2;
		  __C2 = __C1;
		  __temp = _Tp(2 * __j) + __2lambda;
		  auto __a = _Tp(2 * __j) * (__j + __2lambda)
			   * (__temp - _Tp{2});
		  auto __b = (__temp - _Tp{1})
			   * __temp * (__temp - _Tp{2}) * __z;
		  auto __c = _Tp{2} * (__j - 1 + __lambda)
			   * (__j - 1 + __lambda) * __temp;
		  __C1 = (__b * __C2 - __c * __C3) / __a;
		}
	      auto __Cp = (__n * (-__temp * __z) * __C1
			+ _Tp{2} * (__n + __lambda) * (__n + __lambda) * __C2)
			/ (__temp * (_Tp{1} - __z * __z));
	      auto __z1 = __z;
	      __z = __z1 - __C1 / __Cp;
	      if (std::abs(__z - __z1) <= _S_eps)
		{
		  __w = std::exp(std::lgamma(__lambda + _Tp(__n))
			       + std::lgamma(__lambda + _Tp(__n))
			       - std::lgamma(_Tp(__n + 1))
			       - std::lgamma(_Tp(__n + 1) + __2lambda))
		      * __temp * std::pow(_Tp{2}, __2lambda) / (__Cp * __C2);
		  break;
		}
	      if (__its > _S_maxit)
		std::__throw_logic_error("__gegenbauer_zeros: Too many iterations");
	    }
	  __pt[__i - 1].__point = __z;
	  __pt[__i - 1].__weight = __w;
	}

      return __pt;
    }

} // namespace __gnu_cxx

#endif // GEGENBAUER_ZEROS_TCC
