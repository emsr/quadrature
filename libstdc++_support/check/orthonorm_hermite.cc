
#include <limits>
#include <cassert>
#include <cmath>
#include "integration.h"

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_hermite(int n1, int n2, _Tp x)
  {
    const auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
    _Tp lnorm = _Tp{0.5} * (log(_S_pi) + (n1 + n2) * log(_Tp{2})
			  + __gnu_cxx::lfactorial<_Tp>(n1)
			  + __gnu_cxx::lfactorial<_Tp>(n2));
    return std::hermite(n2, x)
	 * std::exp(-x * x - lnorm)
	 * std::hermite(n1, x);
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_hermite()
  {
    const _Tp eps = std::numeric_limits<_Tp>::epsilon();
    const _Tp infty = std::numeric_limits<_Tp>::infinity();
    const _Tp rel_precision = _Tp{100} * eps;
    const _Tp abs_precision = _Tp{10} * rel_precision;

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [n1, n2](_Tp x)
			-> _Tp
			{ return norm_hermite(n1, n2, x); };
	    using func_t = decltype(func);

	    auto [result, error]
		= __gnu_cxx::integrate_clenshaw_curtis(
				map_minf_pinf<func_t, _Tp>(func),
				_Tp{0}, _Tp{1}, rel_precision, _Tp{0});

	    assert(std::abs(delta<_Tp>(n1, n2) - result) < cmp_prec);
	  }
      }
  }

int
main()
{
  test_hermite<float>();

  test_hermite<double>();

  test_hermite<long double>();
}
