
#include <cmath>
#include "integration.h"

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_laguerre(int n1, int n2, _Tp x)
  {
    return std::exp(-x)
	 * std::laguerre(n1, x)
	 * std::laguerre(n2, x);
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_assoc_laguerre()
  {
    const _Tp eps = std::numeric_limits<_Tp>::epsilon();
    const _Tp integ_precision = _Tp{100} * eps;
    const _Tp comp_precision = _Tp{10} * integ_precision;

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [n1 = itop, n2](_Tp x)
			-> _Tp
			{ return normalized_laguerre<_Tp>(n1, n2, x); };

	    auto [result, error]
		= integrate_to_infinity(func, _Tp{0}, integ_precision, _Tp{0});

	    assert(std::abs(delta<_Tp>(n1, n2) - result) < abs_precision);
	  }
      }
  }

int
main()
{
  test_assoc_laguerre<float>();

  test_assoc_laguerre<double>();

  test_assoc_laguerre<long double>();
}
