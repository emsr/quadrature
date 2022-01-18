
#include <limits>
#include <cassert>
#include <cmath>

#include <ext/integration.h>

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  norm_hermite(int n1, int n2, Tp x)
  {
    constexpr auto pi = std::numbers::pi_v<Tp>;
    Tp lnorm = Tp{0.5} * (std::log(pi) + (n1 + n2) * std::log(Tp{2})
			  + std::lgamma(Tp(1 + n1))
			  + std::lgamma(Tp(1 + n2)));
    return std::hermite(n2, x)
	 * std::exp(-x * x - lnorm)
	 * std::hermite(n1, x);
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  void
  test_hermite()
  {
    constexpr auto eps = std::numeric_limits<Tp>::epsilon();
    constexpr auto rel_precision = Tp{100} * eps;
    constexpr auto abs_precision = Tp{10} * rel_precision;

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [n1, n2](Tp x)
			-> Tp
			{ return norm_hermite(n1, n2, x); };
	    using func_t = decltype(func);

	    auto [result, error]
		= emsr::integrate_clenshaw_curtis(
				emsr::map_minf_pinf<Tp, func_t>(func),
				Tp{0}, Tp{1}, rel_precision, Tp{0});

	    assert(std::abs(delta<Tp>(n1, n2) - result) < abs_precision);
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
