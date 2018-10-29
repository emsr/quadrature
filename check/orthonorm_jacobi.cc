
#include <limits>
#include <cassert>
#include <ext/cmath>
#include "integration.h"

// Try to manage the four-gamma ratio.
// alpha > -1, beta > -1.
template<typename _Tp>
  _Tp
  gamma_ratio(int n, _Tp alpha, _Tp beta)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    if (std::abs(_Tp(1) + alpha + beta) < _S_eps)
      return _Tp(0);
    else
      {
	auto gaman1 = std::tgamma(_Tp(1) + alpha);
	auto gambn1 = std::tgamma(_Tp(1) + beta);
	auto gamabn1 = std::tgamma(_Tp(1) + alpha + beta);
	auto fact = gaman1 * gambn1 / gamabn1;
	for (int k = 1; k <= n; ++k)
	  fact *= (_Tp(k) + alpha) * (_Tp(k) + beta)
		/ (_Tp(k) + alpha + beta) / _Tp(k);
	return fact;
      }
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_jacobi(int n1, int n2, _Tp alpha, _Tp beta, _Tp x)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    if (std::abs(x - _Tp{1}) < _S_eps)
      return _Tp{0};
    else if (std::abs(x + _Tp{1}) < _S_eps)
      return _Tp{0};
    else
      {
	auto gam = gamma_ratio(n1, alpha, beta);
	auto norm = std::pow(_Tp{2}, _Tp{1} + alpha + beta)
		  * gam / (_Tp(2 * n1 + 1) + alpha + beta);
	return std::pow(_Tp{1} - x, alpha) * std::pow(_Tp{1} + x, beta)
	     * __gnu_cxx::jacobi(n1, alpha, beta, x)
	     * __gnu_cxx::jacobi(n2, alpha, beta, x) / norm;
      }
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_jacobi(_Tp alpha, _Tp beta)
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto integ_prec = _Tp{1000} * eps;
    const auto cmp_prec = _Tp{10} * integ_prec;

    const bool singular = (alpha < _Tp{0} || beta < _Tp{0});

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [n1, n2, alpha, beta](_Tp x)
			-> _Tp
		     { return normalized_jacobi<_Tp>(n1, n2, alpha, beta, x); };

	    auto [result, error]
		= singular
		? integrate_singular_endpoints(func,
					       _Tp{-1}, _Tp{1},
					       alpha, beta, 0, 0,
					       integ_prec, _Tp{0})
		: integrate(func, _Tp{-1}, _Tp{1}, integ_prec, _Tp{0});

	    assert(std::abs(delta<_Tp>(n1, n2) - result) < cmp_prec);
	  }
      }
  }

int
main()
{
  test_jacobi<float>();

  test_jacobi<double>();

  test_jacobi<long double>();
}
