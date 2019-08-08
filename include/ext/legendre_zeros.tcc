#ifndef LEGENDRE_ZEROS_TCC
#define LEGENDRE_ZEROS_TCC 1

namespace __gnu_cxx
{

  /**
   * Build a list of zeros and weights for the Gauss-Legendre integration rule
   * for the Legendre polynomial of degree @c l.
   */
  template<typename _Tp>
    std::vector<__quadrature_point_t<_Tp>>
    __legendre_zeros(unsigned int __l)
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto _S_pi = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
      const unsigned int _S_maxit = 1000u;

      std::vector<__quadrature_point_t<_Tp>> __pt(__l);

      auto __m = __l / 2;

      // Treat the central zero for odd degree specially.
      // Be careful to avoid overflow of the factorials.
      // An alternative would be to proceed with the recursion
      // for large degree.
      if (__l & 1)
	{
	  const auto __lm = __l - 1;
	  const auto __mm = __lm / 2;
	  auto _Am = _Tp{1};
	  for (auto __m = 1u; __m <= __mm; ++__m)
	    _Am *= -_Tp(2 * __m - 1) / _Tp(2 * __m);
	  auto __Plm1 = _Am;
	  auto __Ppl = __l * __Plm1;
	  __pt[__m].__point = _Tp{0};
	  __pt[__m].__weight = _Tp{2} / __Ppl / __Ppl;
	}

      for (auto __i = 1u; __i <= __m; ++__i)
	{
	  // Clever approximation of root.
	  auto __z = std::cos(_S_pi * (__i - _Tp{1} / _Tp{4})
				    / (__l + _Tp{1} / _Tp{2}));
	  auto __z1 = __z;
	  auto __w = _Tp{0};
	  for (auto __its = 0u; __its < _S_maxit; ++__its)
	    {
	      // Compute __P, __P1, and __P2 the Legendre polynomials of degree
	      // l, l-1, l-2 respectively by iterating through the recursion
	      // relation for the Legendre polynomials.
	      // Compute __Pp the derivative of the Legendre polynomial
	      // of degree l.
	      auto __P1 = _Tp{0};
	      auto __P = _Tp{1};
	      for  (auto __k = 1u; __k <= __l; ++__k)
		{
		  auto __P2 = __P1;
		  __P1 = __P;
		  // Recursion for Legendre polynomials.
		  __P = ((_Tp{2} * __k - _Tp{1}) * __z * __P1
		      - (__k - _Tp{1}) * __P2) / __k;
		}
	      // Recursion for the derivative of The Legendre polynomial.
	      auto __Pp = __l * (__z * __P - __P1) / (__z * __z - _Tp{1});
	      __z1 = __z;
	      // Converge on root by Newton's method.
	      __z = __z1 - __P / __Pp;
	      if (std::abs(__z - __z1) < _S_eps)
		{
		  __w = _Tp{2} / ((_Tp{1} - __z * __z) * __Pp * __Pp);
		  break;
		}
	      if (__its > _S_maxit)
		std::__throw_logic_error("__legendre_zeros: "
					 "Too many iterations");
	    }

	  __pt[__i - 1].__point = -__z;
	  __pt[__l - __i].__point = __z;
	  __pt[__i - 1].__weight = __w;
	  __pt[__l - __i].__weight = __w;
	}

      return __pt;
    }

} // namespace __gnu_cxx

#endif // LEGENDRE_ZEROS_TCC
