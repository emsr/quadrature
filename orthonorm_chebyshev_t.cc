
#include <limits>
#include <cassert>
#include <ext/cmath>
#include "integration.h"

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_chebyshev_t(int n1, int n2, _Tp x)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_inf = std::numeric_limits<_Tp>::infinity();
    if (std::abs(x - _Tp{1}) < _S_eps)
      return _S_inf;
    else if (std::abs(x + _Tp{1}) < _S_eps)
      return (n1 + n2) & 1 ? -_S_inf : _S_inf;
    else
      return __gnu_cxx::chebyshev_t(n2, x)
	   * __gnu_cxx::chebyshev_t(n1, x)
	   / std::sqrt(_Tp{1} - x * x);
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_chebyshev_t()
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto integ_prec = _Tp{100} * eps;
    const auto cmp_prec = _Tp{10} * integ_prec;

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [n1, n2](_Tp x)->_Tp{return norm_chebyshev_t(n1, n2, x);};

	    auto [result, error]
		= __gnu_cxx::integrate_singular_endpoints(func,
				 _Tp{-1}, _Tp{1},
				 _Tp{-0.5}, _Tp{-0.5}, 0, 0,
				 integ_prec, _Tp{0});

	    assert(std::abs(delta<_Tp>(n1, n2) - result) < cmp_prec);
	  }
      }
  }

int
main()
{
  test_chebyshev_t<float>();

  test_chebyshev_t<double>();

  test_chebyshev_t<long double>();
}