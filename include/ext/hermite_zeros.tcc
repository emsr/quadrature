#ifndef HERMITE_ZEROS_TCC
#define HERMITE_ZEROS_TCC 1

#include <cmath>

namespace __gnu_cxx
{

  /**
   * Build a vector of the Gauss-Hermite integration rule abscissae and weights.
   */
  template<typename _Tp>
    std::vector<__quadrature_point_t<_Tp>>
    __hermite_zeros(unsigned int __n)
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const unsigned int _S_maxit = 1000u;
      const auto _S_pim4 = _Tp{0.7511255444649424828587030047762276930510L};
      const auto _S_sqrt_pi = _Tp{1} / _Tp{5.6418'95835'47756'28694'80794'51560'77258'58438e-1L};

      std::vector<__quadrature_point_t<_Tp>> __pt(__n);

      const auto __m = __n / 2;

      // Treat the central zero for odd order specially.
      // Be careful to avoid overflow of the factorials.
      // An alternative would be to proceed with the recursion
      // for large order.
      if (__n & 1)
	{
	  auto __nm = __n - 1;
	  auto __nmfact = std::lgamma(_Tp(__nm - 1));
	  auto __mm = __nm / 2;
	  auto __mmfact = std::lgamma(_Tp(__mm - 1));
	  __pt[__m].__point = _Tp{0};
	  __pt[__m].__weight = _S_sqrt_pi * std::pow(_Tp{2}, _Tp(__n - 1))
			     *std::exp(-(__nmfact - 2 * __mmfact)) / __n;
	}

      for (auto __i = 1u; __i <= __m; ++__i)
	{
	  _Tp __z;
	  _Tp __w = _Tp{0};
	  if (__i == 1)
	    __z = std::sqrt(_Tp(2 * __n + 1))
		- 1.85575 * std::pow(_Tp(2 * __n + 1), -0.166667);
	  else if (__i == 2)
	    __z -= 1.14 * std::pow(_Tp(__n), 0.426) / __z;
	  else if (__i == 3)
	    __z = 1.86 * __z - 0.86 * __pt[0].__point;
	  else if (__i == 4)
	    __z = 1.91 * __z - 0.91 * __pt[1].__point;
	  else
	    __z = 2.0 * __z - __pt[__i - 3].__point;
	  for (auto __its = 1u; __its <= _S_maxit; ++__its)
	    {
	      auto __H = _S_pim4;
	      auto __H1 = _Tp{0};
	      for (auto __k = 1u; __k <= __n; ++__k)
		{
		  auto __H2 = __H1;
		  __H1 = __H;
		  __H = __z * std::sqrt(_Tp{2} / __k) * __H1
		       - std::sqrt(_Tp(__k - 1) / _Tp(__k)) * __H2;
		}
	      auto __Hp = std::sqrt(_Tp(2 * __n)) * __H1;
	      auto __z1 = __z;
	      __z = __z1 - __H / __Hp;
	      if (std::abs(__z - __z1) <= _S_eps)
		{
		  __w = _Tp{2} / (__Hp * __Hp);
		  break;
		}
	      if (__its > _S_maxit)
		std::__throw_logic_error("__hermite_zeros: "
					 "Too many iterations");
	    }
	  __pt[__n - __i].__point = -__z;
	  __pt[__n - __i].__weight = __w;
	  __pt[__i - 1].__point = __z;
	  __pt[__i - 1].__weight = __w;
	}

      return __pt;
    }

} // namespace __gnu_cxx

#endif // HERMITE_ZEROS_TCC
