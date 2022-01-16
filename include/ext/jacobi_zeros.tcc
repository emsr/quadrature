#ifndef JACOBI_ZEROS_TCC
#define JACOBI_ZEROS_TCC 1

namespace emsr
{

  /**
   * Return a vector containing the zeros of the Jacobi polynomial
   * @f$ P_n^{(\alpha,\beta)}(x) @f$.
   * Thias works for @f$ \alpha, \beta > -1 @f$.
   *
   * @tparam  Tp  The real type of the parameters
   * @param[in]  n  The degree of the Jacobi polynomial
   * @param[in]  alpha1  The first order parameter of the Jacobi polynomial
   * @param[in]  beta1  The second order parameter of the Jacobi polynomial
   */
  template<typename Tp>
    std::vector<quadrature_point_t<Tp>>
    jacobi_zeros(unsigned int n, Tp alpha1, Tp beta1)
    {
      const auto s_eps = epsilon(alpha1);
      const unsigned int s_maxit = 1000u;

      std::vector<quadrature_point_t<Tp>> pt(n);

      Tp z;
      Tp w;
      for (auto i = 1u; i <= n; ++i)
	{
	  if (i == 1)
	    {
	      const auto an = alpha1 / n;
	      const auto bn = beta1 / n;
	      const auto r1 = (1.0 + alpha1) * (2.78 / (4.0 + n * n)
			+ 0.768 * an / n);
	      const auto r2 = 1.0 + 1.48 * an + 0.96 * bn
			      + 0.452 * an * an + 0.83 * an * bn;
	      z = 1.0 - r1 / r2;
	    }
	  else if (i == 2)
	    {
	      const auto r1 = (4.1 + alpha1)
			/ ((1.0 + alpha1) * (1.0 + 0.156 * alpha1));
	      const auto r2 = 1.0
			+ 0.06 * (n - 8.0) * (1.0 + 0.12 * alpha1) / n;
	      const auto r3 = 1.0
			+ 0.012 * beta1 * (1.0 + 0.25 * std::abs(alpha1))
			/ n;
	      z -= (1.0 - z) * r1 * r2 * r3;
	    }
	  else if (i == 3)
	    {
	      const auto r1 = (1.67 + 0.28 * alpha1)
			      / (1.0 + 0.37 * alpha1);
	      const auto r2 = 1.0 + 0.22 * (n - 8.0) / n;
	      const auto r3 = 1.0 + 8.0 * beta1
				/ ((6.28 + beta1) * n * n);
	      z -= (pt[0].point - z) * r1 * r2 * r3;
	    }
	  else if (i == n - 1)
	    {
	      const auto r1 = (1.0 + 0.235 * beta1)
			      / (0.766 + 0.119 * beta1);
	      const auto r2 = 1.0 / (1.0 + 0.639 * (n - 4.0)
						/ (1.0 + 0.71 * (n - 4.0)));
	      const auto r3 = 1.0 / (1.0 + 20.0 * alpha1
				/ ((7.5 + alpha1) * n * n));
	      z += (z - pt[n - 4].point) * r1 * r2 * r3;
	    }
	  else if (i == n)
	    {
	      const auto r1 = (1.0 + 0.37 * beta1)
			      / (1.67 + 0.28 * beta1);
	      const auto r2 = 1.0 / (1.0 + 0.22 * (n - 8.0) / n);
	      const auto r3 = 1.0 / (1.0 + 8.0 * alpha1
				/ ((6.28 + alpha1) * n * n));
	      z += (z - pt[n - 3].point) * r1 * r2 * r3;
	    }
	  else
	    {
	      z = 3.0 * pt[i - 2].point
		  - 3.0 * pt[i - 3].point + pt[i - 4].point;
	    }

	  auto alphabeta = alpha1 + beta1;
	  for (auto its = 1u; its <= s_maxit; ++its)
	    {
	      auto temp = Tp{2} + alphabeta;
	      auto P1 = (alpha1 - beta1 + temp * z) / Tp{2};
	      auto P2 = Tp{1};
	      for (auto j = 2u; j <= n; ++j)
		{
		  auto P3 = P2;
		  P2 = P1;
		  temp = Tp{2} * j + alphabeta;
		  auto a = Tp{2} * j * (j + alphabeta)
			   * (temp - Tp{2});
		  auto b = (temp - Tp{1})
			   * (alpha1 * alpha1 - beta1 * beta1
				+ temp * (temp - Tp{2}) * z);
		  auto c = Tp{2} * (j - 1 + alpha1)
			   * (j - 1 + beta1) * temp;
		  P1 = (b * P2 - c * P3) / a;
		}
	      auto Pp = (n * (alpha1 - beta1 - temp * z) * P1
			   + Tp{2} * (n + alpha1) * (n + beta1) * P2)
			/ (temp * (Tp{1} - z * z));
	      auto z1 = z;
	      z = z1 - P1 / Pp;
	      if (std::abs(z - z1) <= s_eps)
		{
		  w = std::exp(std::lgamma(alpha1 + Tp(n))
			       + std::lgamma(beta1 + Tp(n))
			       - std::lgamma(Tp(n + 1))
			       - std::lgamma(Tp(n + 1) + alphabeta))
		      * temp * std::pow(Tp{2}, alphabeta) / (Pp * P2);
		  break;
		}
	      if (its > s_maxit)
		throw std::logic_error("jacobi_zeros: Too many iterations");
	    }
	  pt[i - 1].point = z;
	  pt[i - 1].weight = w;
	}

      return pt;
    }

} // namespace emsr

#endif // JACOBI_ZEROS_TCC
