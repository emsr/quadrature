
#include <limits>
#include <cassert>
#include <ext/cmath>
#include "integration.h"

// Try to manage the gamma ratio.
template<typename _Tp>
  _Tp
  gamma_ratio(int n, int m)
  {
    auto gaman1 = std::tgamma(_Tp(1) + m);
    auto fact = gaman1;
    for (int k = 1; k <= n; ++k)
      fact *= (_Tp(k) + m) / _Tp(k);
    return fact;
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_assoc_laguerre(int n1, int n2, int m, _Tp x)
  {
    auto norm = gamma_ratio(n1, m);
    return std::pow(x, m) * std::exp(-x)
	 * std::assoc_laguerre(n1, m, x)
	 * std::assoc_laguerre(n2, m, x) / norm;
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_assoc_laguerre(_Tp alpha)
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto integ_prec = _Tp{1000} * eps;
    const auto cmp_prec = _Tp{10} * integ_prec;

    for (int n1 : {0, 2, 5})
      {
	for (int n2 : {0, 2, 5})
	  {
	    auto func = [n1, n2, alpha](_Tp x)
			-> _Tp
			{ return norm_assoc_laguerre<_Tp>(n1, n2, alpha, x); };

	    auto [result, error]
		= __gnu_cxx::integrate_to_infinity(func, _Tp{0},
						   integ_prec, _Tp{0});

	    assert(std::abs(delta<_Tp>(n1, n2) - result) < cmp_prec);
	  }
      }
  }

int
main()
{
  test_assoc_laguerre<float>(0);

  test_assoc_laguerre<double>(0);

  test_assoc_laguerre<long double>(0);
}
