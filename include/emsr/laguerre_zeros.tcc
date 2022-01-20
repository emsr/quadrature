#ifndef LAGUERRE_ZEROS_TCC
#define LAGUERRE_ZEROS_TCC 1

namespace emsr
{

  /**
   * Return an array of abscissae and weights for the Gauss-Laguerre rule.
   */
  template<typename Tp>
    std::vector<QuadraturePoint<Tp>>
    laguerre_zeros(unsigned int n, Tp alpha1)
    {
      const auto s_eps = std::numeric_limits<Tp>::epsilon();
      const unsigned int s_maxit = 1000;

      std::vector<QuadraturePoint<Tp>> pt(n);

      for (auto i = 1u; i <= n; ++i)
	{
	  auto z = Tp{0};
	  auto w = Tp{0};
	  // Clever approximations for roots.
	  if (i == 1)
	    z += (1.0 + alpha1)
		 * (3.0 + 0.92 * alpha1) / (1.0 + 2.4 * n + 1.8 * alpha1);
	  else if (i == 2)
	    z += (15.0 + 6.25 * alpha1) / (1.0 + 2.5 * n + 0.9 * alpha1);
	  else
	    {
	      auto ai = i - 2;
	      z += ((1.0 + 2.55 * ai) / (1.9 * ai)
		     + 1.26 * ai * alpha1 / (1.0 + 3.5 * ai))
		   * (z - pt[i - 3].point) / (1.0 + 0.3 * alpha1);
	    }
	  // Iterate TTRR for polynomial values
	  for (auto its = 1u; its <= s_maxit; ++its)
	    {
	      auto L2 = Tp{0};
	      auto L1 = Tp{1};
	      for (auto j = 1u; j <= n; ++j)
		{
		  auto L3 = L2;
		  L2 = L1;
		  L1 = ((Tp(2 * j - 1 + alpha1) - z) * L2
		     - (Tp(j - 1 + alpha1)) * L3) / Tp(j);
		}
	      // Derivative.
	      auto Lp = (Tp(n) * L1 - Tp(n + alpha1) * L2) / z;
	      // Newton's rule for root.
	      auto z1 = z;
	      z = z1 - L1 / Lp;
	      if (std::abs(z - z1) <= s_eps)
		{
		  auto exparg = std::lgamma(Tp(alpha1 + n))
				- std::lgamma(Tp(n));
		  w = -std::exp(exparg) / (Lp * n * L2);
		  break;
		}
	      if (its > s_maxit)
		throw std::logic_error("laguerre_zeros: Too many iterations");
	   }
	  pt[i - 1].point = z;
	  pt[i - 1].weight = w;
	}
    return pt;
  }

} // namespace emsr

#endif // LAGUERRE_ZEROS_TCC
