#ifndef HERMITE_ZEROS_TCC
#define HERMITE_ZEROS_TCC 1

#include <cmath>

namespace __gnu_cxx
{

  /**
   * Build a vector of the Gauss-Hermite integration rule abscissae and weights.
   */
  template<typename Tp>
    std::vector<quadrature_point_t<Tp>>
    hermite_zeros(unsigned int n)
    {
      const auto s_eps = std::numeric_limits<Tp>::epsilon();
      const unsigned int s_maxit = 1000u;
      const auto s_pim4 = Tp{0.7511255444649424828587030047762276930510L};
      const auto s_sqrt_pi = Tp{1} / Tp{5.6418'95835'47756'28694'80794'51560'77258'58438e-1L};

      std::vector<quadrature_point_t<Tp>> pt(n);

      const auto m = n / 2;

      // Treat the central zero for odd order specially.
      // Be careful to avoid overflow of the factorials.
      // An alternative would be to proceed with the recursion
      // for large order.
      if (n & 1)
	{
	  auto nm = n - 1;
	  auto nmfact = std::lgamma(Tp(nm - 1));
	  auto mm = nm / 2;
	  auto mmfact = std::lgamma(Tp(mm - 1));
	  pt[m].point = Tp{0};
	  pt[m].weight = s_sqrt_pi * std::pow(Tp{2}, Tp(n - 1))
			     *std::exp(-(nmfact - 2 * mmfact)) / n;
	}

      for (auto i = 1u; i <= m; ++i)
	{
	  Tp z;
	  Tp w = Tp{0};
	  if (i == 1)
	    z = std::sqrt(Tp(2 * n + 1))
		- 1.85575 * std::pow(Tp(2 * n + 1), -0.166667);
	  else if (i == 2)
	    z -= 1.14 * std::pow(Tp(n), 0.426) / z;
	  else if (i == 3)
	    z = 1.86 * z - 0.86 * pt[0].point;
	  else if (i == 4)
	    z = 1.91 * z - 0.91 * pt[1].point;
	  else
	    z = 2.0 * z - pt[i - 3].point;
	  for (auto its = 1u; its <= s_maxit; ++its)
	    {
	      auto H = s_pim4;
	      auto H1 = Tp{0};
	      for (auto k = 1u; k <= n; ++k)
		{
		  auto H2 = H1;
		  H1 = H;
		  H = z * std::sqrt(Tp{2} / k) * H1
		    - std::sqrt(Tp(k - 1) / Tp(k)) * H2;
		}
	      auto Hp = std::sqrt(Tp(2 * n)) * H1;
	      auto z1 = z;
	      z = z1 - H / Hp;
	      if (std::abs(z - z1) <= s_eps)
		{
		  w = Tp{2} / (Hp * Hp);
		  break;
		}
	      if (its > s_maxit)
		throw std::logic_error("hermite_zeros: Too many iterations");
	    }
	  pt[n - i].point = -z;
	  pt[n - i].weight = w;
	  pt[i - 1].point = z;
	  pt[i - 1].weight = w;
	}

      return pt;
    }

} // namespace __gnu_cxx

#endif // HERMITE_ZEROS_TCC
