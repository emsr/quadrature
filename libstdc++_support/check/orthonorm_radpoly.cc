
#include <limits>
#include <cassert>
#include <ext/cmath>
#include "integration.h"

/**
 * Neumann's number
 */
template<typename _Tp>
  _Tp
  epsilon(int m)
  { return m == 0 ? _Tp{2} : _Tp{1}; }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_radpoly(int n1, int m1, int n2, int m2, _Tp rho)
  {
    auto norm = _Tp{1} / std::sqrt(_Tp(2 * n1 + 2) * _Tp(2 * n2 + 2));
    return rho
	 * __gnu_cxx::radpoly(n1, m1, rho)
	 * __gnu_cxx::radpoly(n2, m2, rho) / norm;
  }

template<typename _Tp>
  _Tp
  delta(int n1, int m1, int n2, int m2)
  { return (n1 == n2 && m1 == m2) ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_radpoly()
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto integ_precision = _Tp{1000} * eps;
    const auto comp_precision = _Tp{10} * integ_precision;

    for (int n1 : {0, 2, 5})
      {
	for (int m1 = 0; m1 <= n1; ++m1)
	  {
	    if ((n1 - m1) & 1)
	      continue;
	    for (int n2 : {0, 2, 5})
	      {
		// The orthonormality only works for m2 == m1.
		int m2 = m1;
		if (m2 > n2)
		  continue;
		//for (int m2 = 0; m2 <= n2; ++m2)
		  {
		    if ((n2 - m2) & 1)
		      continue;

		    auto func = [n1, m1, n2, m2](_Tp x)
				-> _Tp
			      { return normalized_radpoly(n1, m1, n2, m2, x); };

		    auto [result, error]
			= integrate_singular(func, _Tp{0}, _Tp{1},
					     integ_precision, _Tp{0});

		    assert(std::abs(delta<_Tp>(n1, n2) - result) < cmp_prec);
		  }
	      }
	  }
      }
  }

int
main()
{
  test_radpoly<float>();

  test_radpoly<double>();

  test_radpoly<long double>();
}
