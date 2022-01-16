#ifndef LEGENDRE_ZEROS_TCC
#define LEGENDRE_ZEROS_TCC 1

#include <cmath>

namespace emsr
{

  /**
   * Build a list of zeros and weights for the Gauss-Legendre integration rule
   * for the Legendre polynomial of degree @c l.
   */
  template<typename Tp>
    std::vector<quadrature_point_t<Tp>>
    legendre_zeros(unsigned int l)
    {
      const auto s_eps = std::numeric_limits<Tp>::epsilon();
      const auto s_pi = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
      const unsigned int s_maxit = 1000u;

      std::vector<quadrature_point_t<Tp>> pt(l);

      auto m = l / 2;

      // Treat the central zero for odd degree specially.
      // Be careful to avoid overflow of the factorials.
      // An alternative would be to proceed with the recursion
      // for large degree.
      if (l & 1)
	{
	  const auto lm = l - 1;
	  const auto mm = lm / 2;
	  auto Am = Tp{1};
	  for (auto m = 1u; m <= mm; ++m)
	    Am *= -Tp(2 * m - 1) / Tp(2 * m);
	  auto Plm1 = Am;
	  auto Ppl = l * Plm1;
	  pt[m].point = Tp{0};
	  pt[m].weight = Tp{2} / Ppl / Ppl;
	}

      for (auto i = 1u; i <= m; ++i)
	{
	  // Clever approximation of root.
	  auto z = std::cos(s_pi * (i - Tp{1} / Tp{4})
				    / (l + Tp{1} / Tp{2}));
	  auto z1 = z;
	  auto w = Tp{0};
	  for (auto its = 0u; its < s_maxit; ++its)
	    {
	      // Compute P, P1, and P2 the Legendre polynomials of degree
	      // l, l-1, l-2 respectively by iterating through the recursion
	      // relation for the Legendre polynomials.
	      // Compute Pp the derivative of the Legendre polynomial
	      // of degree l.
	      auto P1 = Tp{0};
	      auto P = Tp{1};
	      for  (auto k = 1u; k <= l; ++k)
		{
		  auto P2 = P1;
		  P1 = P;
		  // Recursion for Legendre polynomials.
		  P = ((Tp{2} * k - Tp{1}) * z * P1
		      - (k - Tp{1}) * P2) / k;
		}
	      // Recursion for the derivative of The Legendre polynomial.
	      auto Pp = l * (z * P - P1) / (z * z - Tp{1});
	      z1 = z;
	      // Converge on root by Newton's method.
	      z = z1 - P / Pp;
	      if (std::abs(z - z1) < s_eps)
		{
		  w = Tp{2} / ((Tp{1} - z * z) * Pp * Pp);
		  break;
		}
	      if (its > s_maxit)
		throw std::logic_error("legendre_zeros: Too many iterations");
	    }

	  pt[i - 1].point = -z;
	  pt[l - i].point = z;
	  pt[i - 1].weight = w;
	  pt[l - i].weight = w;
	}

      return pt;
    }

} // namespace emsr

#endif // LEGENDRE_ZEROS_TCC
