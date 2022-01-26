//
// Copyright (C) 2018-2020 Free Software Foundation, Inc.
// Copyright (C) 2021-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.

#ifndef GAUSS_KRONROD_RULE_TCC
#define GAUSS_KRONROD_RULE_TCC 1

#include <stdexcept>

namespace std
{
namespace detail
{

  /**
   * @brief Calculates a Gauss-Kronrod abscissa and weight.
   *
   * @see Robert Piessens, Maria Branders,
   * 	A Note on the Optimal Addition of Abscissas to Quadrature Formulas
   * 	of Gauss and Lobatto, thematics of Computation,
   * 	Volume 28, Number 125, January 1974, pages 135-139.
   *
   * @param[in] n       The order of the Gauss rule.
   * @param[in] eps     The requested absolute accuracy of the abscissas.
   * @param[in] coef2   A value needed to compute weights.
   * @param[in] b[m+1]  The Chebyshev coefficients.
   *
   * @param[in,out] x   On input, an estimate for the abscissa,
   * 		     and on output, the computed abscissa.
   * @param[out]   wk  The Kronrod weight.
   */
  template<typename Tp>
    void
    get_kronrod(int n, Tp eps, Tp coef2,
		  const std::vector<Tp>& b,
		  Tp& x, Tp& wk)
    {
      const int max_iter = 100;

      const int m = (n + 1) / 2;
      const bool even = (2 * m == n);

      int ka = x == Tp{0} ? 1 : 0;

      // Iterative process for the computation of a Kronrod abscissa.
      auto delta = Tp{0};
      auto fd = Tp{0};
      for (int iter = 1; iter <= max_iter; ++iter)
      {

	Tp ai, d2, dif;
	if (even)
	  {
	    ai = m + m + 1;
	    d2 = ai * b[m];
	    dif = Tp{2};
	  }
	else
	  {
	    ai = m + 1;
	    d2 = Tp{0};
	    dif = Tp{1};
	  }

	auto d1 = Tp{0};
	auto b0 = Tp{0};
	auto b1 = Tp{0};
	auto b2 = b[m];
	const auto yy = Tp{4} * x * x - Tp{2};
	for (int k = 1; k <= m; ++k)
	  {
	    ai -= dif;
	    int i = m - k + 1;
	    b0 = b1;
	    b1 = b2;
	    auto d0 = d1;
	    d1 = d2;
	    b2 = yy * b1 - b0 + b[i-1];
	    if (!even)
	      ++i;
	    d2 = yy * d1 - d0 + ai * b[i - 1];
	  }

	Tp f;
	if (even)
	  {
	    f = x * (b2 - b1);
	    fd = d2 + d1;
	  }
	else
	  {
	    f = Tp{0.5L} * (b2 - b0);
	    fd = Tp{4} * x * d2;
	  }

	// Newton correction.
	delta = f / fd;
	x -= delta;

	if (ka == 1)
	  break;

	if (std::fabs(delta) <= eps)
	  ka = 1;
      }

      // Catch non-convergence.
      if (ka != 1)
	{
	  std::ostringstream msg;
	  msg << "\nget_kronrod: Iteration limit reached. eps = " << eps
		<< "; Last delta was " << delta << ".\n";
	  throw std::runtime_error(msg.str().c_str());
	}

      // Computation of the weight.
      auto d0 = Tp{1};
      auto d1 = x;
      auto d2 = Tp{0};
      auto ai = Tp{0};
      for (int k = 2; k <= n; ++k)
	{
	  ai += Tp{1};
	  d2 = ((ai + ai + Tp{1}) * x * d1 - ai * d0)
	       / (ai + Tp{1});
	  d0 = d1;
	  d1 = d2;
	}

      wk = coef2 / (fd * d2);

      return;
    }

  /**
   *  @brief Calculates a Gaussian abscissa and two weights.
   *
   *  @see Robert Piessens, Maria Branders,
   *       A Note on the Optimal Addition of Abscissas to Quadrature Formulas
   *       of Gauss and Lobatto, thematics of Computation,
   *       Volume 28, Number 125, January 1974, pages 135-139.
   *
   *  @param[in] n      The order of the Gauss rule.
   *  @param[in] eps    The requested absolute accuracy of the abscissas.
   *  @param[in] coef2  A value needed to compute weights.
   *  @param[in] b[m+1] The Chebyshev coefficients.
   *
   *  @param[in,out] x   On input, an estimate for the abscissa,
   *                    and on output, the computed abscissa.
   *  @param[out]   wk  The Gauss-Kronrod weight.
   *  @param[out]   wg  The Gauss weight.
   */
  template<typename Tp>
    void
    get_gauss(int n, Tp eps, Tp coef2,
		const std::vector<Tp>& b,
		Tp& x, Tp& wk, Tp& wg)
    {
      const int max_iter = 100;

      const int m = (n + 1) / 2;
      const bool even = (2 * m == n);

      int ka = x == Tp{0} ? 1 : 0;

      // Iterative process for the computation of a Gaussian abscissa.
      auto delta = Tp{0};
      Tp p0, p1, p2, pd2;
      for (int iter = 1; iter <= max_iter; ++iter)
	{
	  p0 = Tp{1};
	  p1 = x;
	  auto pd0 = Tp{0};
	  auto pd1 = Tp{1};

	  // When n is 1, we need to initialize p2 and pd2
	  // to avoid problems with delta.
	  if (n <= 1)
	    {
	      if (std::numeric_limits<Tp>::epsilon() < std::fabs(x))
		{
		  p2 = Tp{0.5L} * (Tp{3} * x * x - Tp{1});
		  pd2 = Tp{3} * x;
		}
	      else
		{
		  p2 = Tp{3} * x;
		  pd2 = Tp{3};
		}
	    }

	  auto ai = Tp{0};
	  for (int k = 2; k <= n; ++k)
	    {
	      ai += Tp{1};
	      p2 = ((ai + ai + Tp{1}) * x * p1 - ai * p0)
		    / (ai + Tp{1});
	      pd2 = ((ai + ai + Tp{1})
			     * (p1 + x * pd1) - ai * pd0)
		    / (ai + Tp{1});
	      p0 = p1;
	      p1 = p2;
	      pd0 = pd1;
	      pd1 = pd2;
	    }

	  // Newton correction.
	  delta = p2 / pd2;
	  x -= delta;

	  if (ka == 1)
	    break;

	  if (std::fabs(delta) <= eps)
	    ka = 1;
	}

      // Catch non-convergence.
      if (ka != 1)
	{
	  std::ostringstream msg;
	  msg << "get_gauss: Iteration limit reached."<< " eps = " << eps
		<< "; Last delta was " << delta << ".\n";
	  throw std::runtime_error(msg.str().c_str());
	}

      // Computation of the weights.
      wg = Tp{2} / (Tp(n) * pd2 * p0);

      p1 = Tp{0};
      p2 = b[m];
      const auto yy = Tp{4} * x * x - Tp{2};
      for (int k = 1; k <= m; ++k)
	{
	  const auto i = m - k + 1;
	  p0 = p1;
	  p1 = p2;
	  p2 = yy * p1 - p0 + b[i - 1];
	}

      if (even)
	wk = wg + coef2 / (pd2 * x * (p2 - p1));
      else
	wk = wg + Tp{2} * coef2 / (pd2 * (p2 - p0));

      return;
    }

} // namespace detail
} // namespace std

namespace emsr
{

  /**
   * @brief Add n+1 points to an n-point Gaussian rule.
   *
   * Discussion:
   *
   *   This subroutine calculates the abscissas and weights of the 2n+1
   *   point Gauss Kronrod quadrature formula which is obtained from the
   *   n-point Gauss quadrature formula by the optimal addition of n+1 points.
   *
   *   The optimally added points are called Kronrod abscissas.  The
   *   abscissas and weights for both the Gauss and Gauss-Kronrod rules
   *   are calculated for integration over the interval [-1,+1].
   *
   *   Since the quadrature formula is symmetric with respect to the origin,
   *   only the nonnegative abscissas are calculated.
   *
   * Storage:
   *
   *   Given n, let m = (n + 1) / 2.
   *
   *   The Gauss-Kronrod rule will include 2n+1 points.  However, by symmetry,
   *   only n+1 Kronrod points and only m Gauss weights need to be listed.
   *
   *   The arrays x, wk and wg contain the nonnegative abscissas in decreasing
   *   order, and the weights of each abscissa in the Gauss-Kronrod and
   *   Gauss rules respectively.
   *
   *   For instance, if n = 3, the output is:
   *
   *   i      x 	      wk	      wg
   *   -    --------        --------        --------
   *   0    0.960491	    0.104656	    0.555556
   *   1    0.774597	    0.268488	    0.888889
   *   2    0.434244	    0.401397
   *   3    0.000000	    0.450917
   *
   *   and if n = 4, (notice that 0 is now a Kronrod abscissa) the output is:
   *
   *   i      x 	      wk	      wg
   *   -    --------        --------        --------
   *   0    0.976560	    0.062977	    0.347855
   *   1    0.861136	    0.170054	    0.652145
   *   2    0.640286	    0.266798
   *   3    0.339981	    0.326949
   *   4    0.000000	    0.346443
   *
   * @see Robert Piessens, Maria Branders,
   * 	  A Note on the Optimal Addition of Abscissas to Quadrature Formulas
   * 	  of Gauss and Lobatto, thematics of Computation,
   * 	  Volume 28, Number 125, January 1974, pages 135-139.
   *
   * @param[in] n  The order of the Gauss rule.
   * @param[in] eps  The requested absolute accuracy of the abscissas.
   *
   * @param[out] x[n+1]   the abscissas.
   * @param[out] wk[n+1]  the weights for the Gauss-Kronrod rule.
   * @param[out] wg[n+1]  the weights for the Gauss rule.
   */
  template<typename Tp>
    void
    build_gauss_kronrod(int n, Tp eps,
			  std::vector<Tp>& x,
			  std::vector<Tp>& wg,
			  std::vector<Tp>& wk)
    {
      const int m = (n + 1) / 2;
      const bool even = (2 * m == n);

      x.resize(n + 1);
      wk.resize(n + 1);
      wg.resize(m);

      std::vector<Tp> b(m + 1);
      std::vector<Tp> tau(m);

      auto d = Tp{2};
      auto an = Tp{0};
      for (int k = 1; k <= n; ++k)
	{
	  an += Tp{1};
	  d *= an / (an + Tp{0.5L});
	}

      // Calculation of the Chebyshev coefficients of the orthogonal polynomial.
      tau[0] = (an + Tp{2}) / (an + an + Tp{3});
      b[m - 1] = tau[0] - Tp{1};
      auto ak = an;
      for (int l = 1; l < m; ++l)
	{
	  ak += Tp{2};
	  tau[l] = ((ak - Tp{1}) * ak
		    - an * (an + Tp{1})) * (ak + Tp{2}) * tau[l - 1]
		    / (ak * ((ak + Tp{3}) * (ak + Tp{2})
		    - an * (an + Tp{1})));
	  b[m - l - 1] = tau[l];
	  for (int ll = 1; ll <= l; ++ll)
	    b[m - l - 1] += tau[ll - 1] * b[m - l + ll - 1];
	}
      b[m] = Tp{1};

      // Calculation of approximate values for the abscissas.
      auto bb = std::sin((Tp{0.5L * M_PI}) / (an + an + Tp{1}));
      auto x1 = std::sqrt(Tp{1} - bb * bb);
      auto s = Tp{2} * bb * x1;
      auto c = std::sqrt(Tp{1} - s * s);
      auto coef = Tp{1} - (Tp{1} - Tp{1} / an) / (Tp{8} * an * an);
      auto xx = coef * x1;

      // Coefficient needed for weights.
      // coef2 = 2^{2n+1} n! n! / (2n+1)!
      //       = 2 4^n n! / product((n+1)...(2*n+1))
      auto coef2 = Tp{2} / Tp(2 * n + 1);
      for (int i = 1; i <= n; ++i)
	coef2 *= Tp{4} * Tp(i) / Tp(n + i);

      // Calculation of the k-th abscissa (a Kronrod abscissa) and the
      // corresponding weight.
      for (int k = 1; k <= n; k += 2)
	{
	  std::detail::get_kronrod(n, eps, coef2, b,
					xx, wk[k - 1]);
	  wg[(k - 1) / 2] = Tp{0};
	  x[k - 1] = xx;
	  auto aa = x1;
	  x1 = aa * c - bb * s;
	  bb = aa * s + bb * c;

	  xx = k == n ? Tp{0} : coef * x1;

	  // Calculation of the k+1 abscissa (a Gaussian abscissa) and the
	  // corresponding weights.
	  std::detail::get_gauss(n, eps, coef2, b,
				     xx, wk[k], wg[k / 2]);

	  x[k] = xx;
	  aa = x1;
	  x1 = aa * c - bb * s;
	  bb = aa * s + bb * c;
	  xx = coef * x1;
	}

      // If n is even we must compute the Kronrod abscissa for the origin.
      if (even)
	{
	  xx = Tp{0};
	  std::detail::get_kronrod(n, eps, coef2, b,
					xx, wk[n]);
	  x[n] = xx;
	}

      return;
    }

} // namespace emsr

#endif // GAUSS_KRONROD_RULE_TCC
