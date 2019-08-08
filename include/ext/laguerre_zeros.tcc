#ifndef LAGUERRE_ZEROS_TCC
#define LAGUERRE_ZEROS_TCC 1

namespace __gnu_cxx
{

  /**
   * Return an array of abscissae and weights for the Gauss-Laguerre rule.
   */
  template<typename _Tp>
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
    __laguerre_zeros(unsigned int __n, _Tp __alpha1)
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const unsigned int _S_maxit = 1000;

      std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __pt(__n);

      for (auto __i = 1u; __i <= __n; ++__i)
	{
	  auto __z = _Tp{0};
	  auto __w = _Tp{0};
	  // Clever approximations for roots.
	  if (__i == 1)
	    __z += (1.0 + __alpha1)
		 * (3.0 + 0.92 * __alpha1) / (1.0 + 2.4 * __n + 1.8 * __alpha1);
	  else if (__i == 2)
	    __z += (15.0 + 6.25 * __alpha1) / (1.0 + 2.5 * __n + 0.9 * __alpha1);
	  else
	    {
	      auto __ai = __i - 2;
	      __z += ((1.0 + 2.55 * __ai) / (1.9 * __ai)
		     + 1.26 * __ai * __alpha1 / (1.0 + 3.5 * __ai))
		   * (__z - __pt[__i - 3].__point) / (1.0 + 0.3 * __alpha1);
	    }
	  // Iterate TTRR for polynomial values
	  for (auto __its = 1u; __its <= _S_maxit; ++__its)
	    {
	      auto __L2 = _Tp{0};
	      auto __L1 = _Tp{1};
	      for (auto __j = 1u; __j <= __n; ++__j)
		{
		  auto __L3 = __L2;
		  __L2 = __L1;
		  __L1 = ((_Tp(2 * __j - 1 + __alpha1) - __z) * __L2
			- (_Tp(__j - 1 + __alpha1)) * __L3) / _Tp(__j);
		}
	      // Derivative.
	      auto __Lp = (_Tp(__n) * __L1 - _Tp(__n + __alpha1) * __L2) / __z;
	      // Newton's rule for root.
	      auto __z1 = __z;
	      __z = __z1 - __L1 / __Lp;
	      if (std::abs(__z - __z1) <= _S_eps)
		{
		  auto __exparg = std::lgamma(_Tp(__alpha1 + __n))
				- std::lgamma(_Tp(__n));
		  __w = -std::exp(__exparg) / (__Lp * __n * __L2);
		  break;
		}
	      if (__its > _S_maxit)
		std::__throw_logic_error("__laguerre_zeros: "
					 "Too many iterations");
	   }
	  __pt[__i - 1].__point = __z;
	  __pt[__i - 1].__weight = __w;
	}
    return __pt;
  }

} // namespace __gnu_cxx

#endif // LAGUERRE_ZEROS_TCC
