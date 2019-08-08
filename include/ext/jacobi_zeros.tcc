#ifndef JACOBI_ZEROS_TCC
#define JACOBI_ZEROS_TCC 1

namespace __gnu_cxx
{

  /**
   * Return a vector containing the zeros of the Jacobi polynomial
   * @f$ P_n^{(\alpha,\beta)}(x) @f$.
   * Thias works for @f$ \alpha, \beta > -1 @f$.
   *
   * @tparam  _Tp  The real type of the parameters
   * @param[in]  __n  The degree of the Jacobi polynomial
   * @param[in]  __alpha1  The first order parameter of the Jacobi polynomial
   * @param[in]  __beta1  The second order parameter of the Jacobi polynomial
   */
  template<typename _Tp>
    std::vector<__quadrature_point_t<_Tp>>
    __jacobi_zeros(unsigned int __n, _Tp __alpha1, _Tp __beta1)
    {
      const auto _S_eps = __epsilon(__alpha1);
      const unsigned int _S_maxit = 1000u;

      std::vector<__quadrature_point_t<_Tp>> __pt(__n);

      _Tp __z;
      _Tp __w;
      for (auto __i = 1u; __i <= __n; ++__i)
	{
	  if (__i == 1)
	    {
	      const auto __an = __alpha1 / __n;
	      const auto __bn = __beta1 / __n;
	      const auto __r1 = (1.0 + __alpha1) * (2.78 / (4.0 + __n * __n)
			+ 0.768 * __an / __n);
	      const auto __r2 = 1.0 + 1.48 * __an + 0.96 * __bn
			      + 0.452 * __an * __an + 0.83 * __an * __bn;
	      __z = 1.0 - __r1 / __r2;
	    }
	  else if (__i == 2)
	    {
	      const auto __r1 = (4.1 + __alpha1)
			/ ((1.0 + __alpha1) * (1.0 + 0.156 * __alpha1));
	      const auto __r2 = 1.0
			+ 0.06 * (__n - 8.0) * (1.0 + 0.12 * __alpha1) / __n;
	      const auto __r3 = 1.0
			+ 0.012 * __beta1 * (1.0 + 0.25 * std::abs(__alpha1))
			/ __n;
	      __z -= (1.0 - __z) * __r1 * __r2 * __r3;
	    }
	  else if (__i == 3)
	    {
	      const auto __r1 = (1.67 + 0.28 * __alpha1)
			      / (1.0 + 0.37 * __alpha1);
	      const auto __r2 = 1.0 + 0.22 * (__n - 8.0) / __n;
	      const auto __r3 = 1.0 + 8.0 * __beta1
				/ ((6.28 + __beta1) * __n * __n);
	      __z -= (__pt[0].__point - __z) * __r1 * __r2 * __r3;
	    }
	  else if (__i == __n - 1)
	    {
	      const auto __r1 = (1.0 + 0.235 * __beta1)
			      / (0.766 + 0.119 * __beta1);
	      const auto __r2 = 1.0 / (1.0 + 0.639 * (__n - 4.0)
						/ (1.0 + 0.71 * (__n - 4.0)));
	      const auto __r3 = 1.0 / (1.0 + 20.0 * __alpha1
				/ ((7.5 + __alpha1) * __n * __n));
	      __z += (__z - __pt[__n - 4].__point) * __r1 * __r2 * __r3;
	    }
	  else if (__i == __n)
	    {
	      const auto __r1 = (1.0 + 0.37 * __beta1)
			      / (1.67 + 0.28 * __beta1);
	      const auto __r2 = 1.0 / (1.0 + 0.22 * (__n - 8.0) / __n);
	      const auto __r3 = 1.0 / (1.0 + 8.0 * __alpha1
				/ ((6.28 + __alpha1) * __n * __n));
	      __z += (__z - __pt[__n - 3].__point) * __r1 * __r2 * __r3;
	    }
	  else
	    {
	      __z = 3.0 * __pt[__i - 2].__point
		  - 3.0 * __pt[__i - 3].__point + __pt[__i - 4].__point;
	    }

	  auto __alphabeta = __alpha1 + __beta1;
	  for (auto __its = 1u; __its <= _S_maxit; ++__its)
	    {
	      auto __temp = _Tp{2} + __alphabeta;
	      auto __P1 = (__alpha1 - __beta1 + __temp * __z) / _Tp{2};
	      auto __P2 = _Tp{1};
	      for (auto __j = 2u; __j <= __n; ++__j)
		{
		  auto __P3 = __P2;
		  __P2 = __P1;
		  __temp = _Tp{2} * __j + __alphabeta;
		  auto __a = _Tp{2} * __j * (__j + __alphabeta)
			   * (__temp - _Tp{2});
		  auto __b = (__temp - _Tp{1})
			   * (__alpha1 * __alpha1 - __beta1 * __beta1
				+ __temp * (__temp - _Tp{2}) * __z);
		  auto __c = _Tp{2} * (__j - 1 + __alpha1)
			   * (__j - 1 + __beta1) * __temp;
		  __P1 = (__b * __P2 - __c * __P3) / __a;
		}
	      auto __Pp = (__n * (__alpha1 - __beta1 - __temp * __z) * __P1
			   + _Tp{2} * (__n + __alpha1) * (__n + __beta1) * __P2)
			/ (__temp * (_Tp{1} - __z * __z));
	      auto __z1 = __z;
	      __z = __z1 - __P1 / __Pp;
	      if (std::abs(__z - __z1) <= _S_eps)
		{
		  __w = std::exp(std::lgamma(__alpha1 + _Tp(__n))
			       + std::lgamma(__beta1 + _Tp(__n))
			       - std::lgamma(_Tp(__n + 1))
			       - std::lgamma(_Tp(__n + 1) + __alphabeta))
		      * __temp * std::pow(_Tp{2}, __alphabeta) / (__Pp * __P2);
		  break;
		}
	      if (__its > _S_maxit)
		std::__throw_logic_error("__jacobi_zeros: Too many iterations");
	    }
	  __pt[__i - 1].__point = __z;
	  __pt[__i - 1].__weight = __w;
	}

      return __pt;
    }

} // namespace __gnu_cxx

#endif // JACOBI_ZEROS_TCC
