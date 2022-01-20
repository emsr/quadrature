
#include <limits>
#include <cassert>
#include <cmath>
#include <numbers>

#include <emsr/integration.h>
#include <emsr/math_constants.h>

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  norm_gegenbauer(int n1, int n2, Tp alpha, Tp x)
  {
    constexpr auto pi = emsr::pi_v<Tp>;
    auto gama = std::tgamma(alpha);
    auto gamn2a = std::tgamma(n1 + Tp{2} * alpha);
    auto norm = pi * std::pow(Tp{2}, Tp{1} - Tp{2} * alpha) * gamn2a
	      / __gnu_cxx::factorial<Tp>(n1) / (Tp(n1) + alpha) / gama / gama;
    return std::pow(Tp{1} - x * x, alpha - Tp{0.5})
	 * __gnu_cxx::gegenbauer(n1, alpha, x)
	 * __gnu_cxx::gegenbauer(n2, alpha, x) / norm;
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  void
  test_gegenbauer(Tp alpha)
  {
    constexpr auto eps = std::numeric_limits<Tp>::epsilon();
    constexpr auto integ_prec = Tp{1000} * eps;
    constexpr auto cmp_prec = Tp{10} * integ_prec;

    const bool singular = (alpha < Tp{0.5});

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [n1, n2, alpha](Tp x)
			-> Tp
			{ return norm_gegenbauer<Tp>(n1, n2, alpha, x); };

	    // Using integrate_singular works pretty well.
	    auto [result, error]
		= singular
		? emsr::integrate_singular_endpoints(func,
					       Tp{-1}, Tp{1},
					       alpha - Tp{0.5},
					       alpha - Tp{0.5}, 0, 0,
					       integ_prec, Tp{0})
		: emsr::integrate(func, Tp{-1}, Tp{1}, integ_prec, Tp{0});

	    assert(std::abs(delta<Tp>(n1, n2) - result) < cmp_prec);
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
