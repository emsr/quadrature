
#include <cmath>
#include "integration.h"

// Function which should integrate to 1 for l1 == l2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_legendre(int l1, int l2, _Tp x)
  {
    return (_Tp(l1 + l2 + 1) / _Tp{2})
	 * std::legendre(l1,x)
	 * std::legendre(l2,x);
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_legendre()
  {
    const _Tp eps = std::numeric_limits<_Tp>::epsilon();
    const _Tp integ_precision = _Tp{100} * eps;
    const _Tp comp_precision = _Tp{10} * integ_precision;

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [l1 = itop, l2](_Tp x)
			-> _Tp
			{ return normalized_legendre(l1, l2, x); };

	    auto [result, error]
		= integrate(func, _Tp{-1}, _Tp{1}, integ_precision, _Tp{0});

	    assert(std::abs(delta<_Tp>(n1, n2) - result) < abs_precision);
	  }
      }
  }

int
main()
{
  test_legendre<float>();

  test_legendre<double>();

  test_legendre<long double>();
}
