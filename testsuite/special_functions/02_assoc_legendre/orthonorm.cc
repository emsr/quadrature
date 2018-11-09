
#include <limits>
#include <cassert>
#include <ext/cmath>
#include "integration.h"

// Function which should integrate to 1 for l1 == l2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_assoc_legendre(int l1, int m1, int l2, int m2, _Tp x)
  {
    return (_Tp(l1 + l2 + 1) / _Tp{2})
	 * std::assoc_legendre(l1, m1, x)
	 * std::assoc_legendre(l2, m2, x);
  }

template<typename _Tp>
  _Tp
  delta(int l1, int l2)
  { return l1 == l2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_assoc_legendre()
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto integ_prec = _Tp{1000} * eps;
    const auto cmp_prec = _Tp{10} * integ_prec;

    for (int n1 : {0, 2, 5})
      {
	for (int m1 = 0; m1 <= n1; ++m1)
	  {
	    for (int n2 : {0, 2, 5})
	      {
		auto func = [l1, m1, l2, m2](_Tp x)
			    -> _Tp
			    { return norm_assoc_legendre(l1, m1, l2, m2, x); };

		auto [result, error]
		    = __gnu_cxx::integrate(func, _Tp{-1}, _Tp{1},
					   integ_prec, _Tp{0});

		assert(std::abs(delta<_Tp>(n1, n2) - result) < cmp_prec);
	      }
	  }
      }
  }

int
main()
{
  test_assoc_legendre<float>();

  test_assoc_legendre<double>();

  test_assoc_legendre<long double>();
}
