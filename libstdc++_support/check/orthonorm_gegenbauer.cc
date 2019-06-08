
#include <limits>
#include <cassert>
#include <cmath>
#include "integration.h"

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_gegenbauer(int n1, int n2, _Tp alpha, _Tp x)
  {
    const auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
    auto gama = std::tgamma(alpha);
    auto gamn2a = std::tgamma(n1 + _Tp{2} * alpha);
    auto norm = _S_pi * std::pow(_Tp{2}, _Tp{1} - _Tp{2} * alpha) * gamn2a
	      / __gnu_cxx::factorial<_Tp>(n1) / (_Tp(n1) + alpha) / gama / gama;
    return std::pow(_Tp{1} - x * x, alpha - _Tp{0.5})
	 * __gnu_cxx::gegenbauer(n1, alpha, x)
	 * __gnu_cxx::gegenbauer(n2, alpha, x) / norm;
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_gegenbauer(_Tp alpha)
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto integ_prec = _Tp{1000} * eps;
    const auto cmp_prec = _Tp{10} * integ_prec;

    const bool singular = (alpha < _Tp{0.5});

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [n1, n2, alpha](_Tp x)
			-> _Tp
			{ return norm_gegenbauer<_Tp>(n1, n2, alpha, x); };

	    // Using integrate_singular works pretty well.
	    auto [result, error]
		= singular
		? __gnu_cxx::integrate_singular_endpoints(func,
					       _Tp{-1}, _Tp{1},
					       alpha - _Tp{0.5},
					       alpha - _Tp{0.5}, 0, 0,
					       integ_prec, _Tp{0})
		: integrate(func, _Tp{-1}, _Tp{1}, integ_prec, _Tp{0});

	    assert(std::abs(delta<_Tp>(n1, n2) - result) < cmp_prec);
	  }
      }
  }

int
main()
{
  test_gegenbauer<float>();

  test_gegenbauer<double>();

  test_gegenbauer<long double>();
}
