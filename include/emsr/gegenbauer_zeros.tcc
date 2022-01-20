#ifndef GEGENBAUER_ZEROS_TCC
#define GEGENBAUER_ZEROS_TCC 1

namespace emsr
{

  /**
   * Return a vector containing the zeros of the Gegenbauer or ultraspherical
   * polynomial @f$ C_n^{(\lambda)}@f$.
   * This works for @f$ \lambda > -1/2 @f$
   *
   * @tparam  Tp  The real type of the order
   * @param[in]  n  The degree of the Gegenbauer polynomial
   * @param[in]  lambda  The order of the Gegenbauer polynomial
   */
  template<typename Tp>
    std::vector<QuadraturePoint<Tp>>
    gegenbauer_zeros(unsigned int n, Tp lambda)
    {
      const auto s_eps = std::numeric_limits<Tp>::epsilon();
      const unsigned int s_maxit = 1000u;
      std::vector<QuadraturePoint<Tp>> pt(n);

      Tp z;
      Tp w;
      for (auto i = 1u; i <= n; ++i)
	{
	  if (i == 1)
	    {
	      auto an = lambda / n;
	      auto an2 = an * an;
	      auto r1 = (1.0 + lambda) * (2.78 / (4.0 + n * n)
			+ 0.768 * an / n);
	      auto r2 = 1.0 + 1.48 * an + 0.96 * an + 1.282 * an2;
	      z = 1.0 - r1 / r2;
	    }
	  else if (i == 2)
	    {
	      auto r1 = (4.1 + lambda)
			/ ((1.0 + lambda) * (1.0 + 0.156 * lambda));
	      auto r2 = 1.0
			+ 0.06 * (n - 8.0) * (1.0 + 0.12 * lambda) / n;
	      auto r3 = 1.0
		   + 0.012 * lambda * (1.0 + 0.25 * std::abs(lambda)) / n;
	      z -= (1.0 - z) * r1 * r2 * r3;
	    }
	  else if (i == 3)
	    {
	      auto r1 = (1.67 + 0.28 * lambda) / (1.0 + 0.37 * lambda);
	      auto r2 = 1.0 + 0.22 * (n - 8.0) / n;
	      auto r3 = 1.0 + 8.0 *lambda / ((6.28 + lambda) * n * n);
	      z -= (pt[0].point - z) * r1 * r2 * r3;
	    }
	  else if (i == n - 1)
	    {
	      auto r1 = (1.0 + 0.235 * lambda) / (0.766 + 0.119 * lambda);
	      auto r2 = 1.0 / (1.0 + 0.639 * (n - 4.0)
						/ (1.0 + 0.71 * (n - 4.0)));
	      auto r3 = 1.0 / (1.0 + 20.0 * lambda
				/ ((7.5 + lambda) * n * n));
	      z += (z - pt[n - 4].point) * r1 * r2 * r3;
	    }
	  else if (i == n)
	    {
	      auto r1 = (1.0 + 0.37 * lambda) / (1.67 + 0.28 * lambda);
	      auto r2 = 1.0 / (1.0 + 0.22 * (n - 8.0) / n);
	      auto r3 = 1.0 / (1.0 + 8.0 * lambda
				 / ((6.28 + lambda) * n * n));
	      z += (z - pt[n - 3].point) * r1 * r2 * r3;
	    }
	  else
	    z = 3.0 * pt[i - 2].point
		- 3.0 * pt[i - 3].point + pt[i - 4].point;

	  auto __2lambda = Tp{2} * lambda;
	  for (auto its = 1u; its <= s_maxit; ++its)
	    {
	      auto temp = Tp{2} + __2lambda;
	      auto C1 = (temp * z) / Tp{2};
	      auto C2 = Tp{1};
	      for (auto j = 2u; j <= n; ++j)
		{
		  auto C3 = C2;
		  C2 = C1;
		  temp = Tp(2 * j) + __2lambda;
		  auto a = Tp(2 * j) * (j + __2lambda)
			   * (temp - Tp{2});
		  auto b = (temp - Tp{1})
			   * temp * (temp - Tp{2}) * z;
		  auto c = Tp{2} * (j - 1 + lambda)
			   * (j - 1 + lambda) * temp;
		  C1 = (b * C2 - c * C3) / a;
		}
	      auto Cp = (n * (-temp * z) * C1
			+ Tp{2} * (n + lambda) * (n + lambda) * C2)
			/ (temp * (Tp{1} - z * z));
	      auto z1 = z;
	      z = z1 - C1 / Cp;
	      if (std::abs(z - z1) <= s_eps)
		{
		  w = std::exp(std::lgamma(lambda + Tp(n))
			       + std::lgamma(lambda + Tp(n))
			       - std::lgamma(Tp(n + 1))
			       - std::lgamma(Tp(n + 1) + __2lambda))
		      * temp * std::pow(Tp{2}, __2lambda) / (Cp * C2);
		  break;
		}
	      if (its > s_maxit)
		throw std::logic_error("gegenbauer_zeros: Too many iterations");
	    }
	  pt[i - 1].point = z;
	  pt[i - 1].weight = w;
	}

      return pt;
    }

} // namespace emsr

#endif // GEGENBAUER_ZEROS_TCC
