/*
$HOME/bin/bin/g++ -std=c++14 -Wall -Wextra -o build_gauss_kronrod build_gauss_kronrod.cpp
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <limits>

/**
 *  @brief Calculates a Gauss-Kronrod abscissa and weight.
 *
 *  @see Robert Piessens, Maria Branders,
 *       A Note on the Optimal Addition of Abscissas to Quadrature Formulas
 *       of Gauss and Lobatto, thematics of Computation,
 *       Volume 28, Number 125, January 1974, pages 135-139.
 *
 *  @param[in] n       The order of the Gauss rule.
 *  @param[in] m       The value of (n + 1) / 2.
 *  @param[in] eps     The requested absolute accuracy of the abscissas.
 *  @param[in] coef2   A value needed to compute weights.
 *  @param[in] even    true if N is even.
 *  @param[in] b[m+1]  The Chebyshev coefficients.
 *
 *  @param[inout] x   On input, an estimate for the abscissa,
 *                    and on output, the computed abscissa.
 *  @param[out]   w   The weight.
 */
template<typename _Tp>
  void
  get_kronrod(int n, int m, _Tp eps, _Tp coef2, bool even,
	      const std::vector<_Tp>& b,
              _Tp& x, _Tp& w)
  {
    const int max_iter = 100;

    int ka;
    if (x == _Tp{0})
      ka = 1;
    else
      ka = 0;

    // Iterative process for the computation of a Kronrod abscissa.
    auto delta = _Tp{0};
    auto fd = _Tp{0};
    for (int iter = 1; iter <= max_iter; ++iter)
    {

      _Tp ai;
      _Tp d2;
      _Tp dif;
      if (even)
	{
	  ai = m + m + 1;
	  d2 = ai * b[m];
	  dif = _Tp{2};
	}
      else
	{
	  ai = m + 1;
	  d2 = _Tp{0};
	  dif = _Tp{1};
	}

      auto d1 = _Tp{0};
      auto b0 = _Tp{0};
      auto b1 = _Tp{0};
      auto b2 = b[m];
      auto yy = _Tp{4} * x * x - _Tp{2};
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

      _Tp f;
      if (even)
	{
	  f = x * (b2 - b1);
	  fd = d2 + d1;
	}
      else
	{
	  f = 0.5 * (b2 - b0);
	  fd = _Tp{4} * x * d2;
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
	std::cout << "\nget_kronrod:";
	std::cout << " Iteration limit reached.";
	std::cout << " eps = " << eps;
	std::cout << "; Last delta was " << delta << ".\n";
	exit(1);
      }

    // Computation of the weight.
    auto d0 = _Tp{1};
    auto d1 = x;
    auto d2 = _Tp{0};
    auto ai = _Tp{0};
    for (int k = 2; k <= n; ++k)
      {
	ai += _Tp{1};
	d2 = ((ai + ai + _Tp{1}) * x * d1 - ai * d0) / (ai + _Tp{1});
	d0 = d1;
	d1 = d2;
      }

    w = coef2 / (fd * d2);

    return;
  }


/**
 *  @brief get_gauss calculates a Gaussian abscissa and two weights.
 *
 *  @see Robert Piessens, Maria Branders,
 *       A Note on the Optimal Addition of Abscissas to Quadrature Formulas
 *       of Gauss and Lobatto, thematics of Computation,
 *       Volume 28, Number 125, January 1974, pages 135-139.
 *
 *  @param[in] n      The order of the Gauss rule.
 *  @param[in] m      The value of (N + 1) / 2.
 *  @param[in] eps    The requested absolute accuracy of the abscissas.
 *  @param[in] coef2  A value needed to compute weights.
 *  @param[in] even   Is TRUE if N is even.
 *  @param[in] b[m+1] The Chebyshev coefficients.
 *
 *  @param[inout] x   On input, an estimate for the abscissa,
 *                    and on output, the computed abscissa.
 *  @param[out]   w1  The Gauss-Kronrod weight.
 *  @param[out]   w2  The Gauss weight.
 */
template<typename _Tp>
  void
  get_gauss(int n, int m, _Tp eps, _Tp coef2, bool even,
	    const std::vector<_Tp>& b,
            _Tp& x, _Tp& w1, _Tp& w2)
  {
    const int max_iter = 100;

    int ka;
    if (x == _Tp{0})
      ka = 1;
    else
      ka = 0;

    // Iterative process for the computation of a Gaussian abscissa.
    auto delta = _Tp{0};
    _Tp p0, p1, p2, pd2;
    for (int iter = 1; iter <= max_iter; ++iter)
    {
      p0 = _Tp{1};
      p1 = x;
      auto pd0 = _Tp{0};
      auto pd1 = _Tp{1};

      // When n is 1, we need to initialize p2 and pd2 to avoid problems with delta.
      if (n <= 1)
	{
	  if (std::numeric_limits<_Tp>::epsilon() < std::fabs(x))
	    {
              p2 = (_Tp{3} * x * x - _Tp{1}) / _Tp{2};
              pd2 = _Tp{3} * x;
	    }
	  else
	    {
              p2 = _Tp{3} * x;
              pd2 = _Tp{3};
	    }
	}

      auto ai = _Tp{0};
      for (int k = 2; k <= n; ++k)
	{
	  ai += _Tp{1};
	  p2 = ((ai + ai + _Tp{1}) * x * p1 - ai * p0) / (ai + _Tp{1});
	  pd2 = ((ai + ai + _Tp{1}) * (p1 + x * pd1) - ai * pd0) 
            / (ai + _Tp{1});
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
	std::cout << "\n ** get_gauss:";
	std::cout << " Iteration limit reached.";
	std::cout << " eps = " << eps;
	std::cout << "; Last delta was " << delta << ".\n";
	exit(1);
      }

    // Computation of the weights.
    w2 = _Tp{2} / (_Tp(n) * pd2 * p0);

    p1 = _Tp{0};
    p2 = b[m];
    auto yy = _Tp{4} * x * x - _Tp{2};
    for (int k = 1; k <= m; ++k)
      {
	auto i = m - k + 1;
	p0 = p1;
	p1 = p2;
	p2 = yy * p1 - p0 + b[i-1];
      }

    if (even)
      w1 = w2 + coef2 / (pd2 * x * (p2 - p1));
    else
      w1 = w2 + _Tp{2} * coef2 / (pd2 * (p2 - p0));

    return;
  }


/**
 *
 *  @brief  adds n+1 points to an n-point Gaussian rule.
 *
 *  Discussion:
 *
 *    This subroutine calculates the abscissas and weights of the 2n+1
 *    point Gauss Kronrod quadrature formula which is obtained from the 
 *    n-point Gauss quadrature formula by the optimal addition of n+1 points.
 *
 *    The optimally added points are called Kronrod abscissas.  The 
 *    abscissas and weights for both the Gauss and Gauss-Kronrod rules
 *    are calculated for integration over the interval [-1,+1].
 *
 *    Since the quadrature formula is symmetric with respect to the origin,
 *    only the nonnegative abscissas are calculated.
 *
 *    Note that the code published in Mathematics of Computation 
 *    omitted the definition of the variable which is here called coef2.
 *
 *  Storage:
 *
 *    Given n, let m = (n + 1) / 2.
 *
 *    The Gauss-Kronrod rule will include 2*n+1 points.  However, by symmetry,
 *    only n + 1 of them need to be listed.
 *
 *    The arrays x, w1 and w2 contain the nonnegative abscissas in decreasing
 *    order, and the weights of each abscissa in the Gauss-Kronrod and
 *    Gauss rules respectively.  This means that about half the entries
 *    in w2 are zero.
 *
 *    For instance, if n = 3, the output is:
 *
 *    I      X               W1              W2
 *
 *    1    0.960491        0.104656         0.000000   
 *    2    0.774597        0.268488         0.555556    
 *    3    0.434244        0.401397         0.000000
 *    4    0.000000        0.450917         0.888889
 *
 *    and if n = 4, (notice that 0 is now a Kronrod abscissa) the output is:
 *
 *    I      X               W1              W2
 *
 *    1    0.976560        0.062977        0.000000   
 *    2    0.861136        0.170054        0.347855    
 *    3    0.640286        0.266798        0.000000   
 *    4    0.339981        0.326949        0.652145    
 *    5    0.000000        0.346443        0.000000
 *
 *  @see Robert Piessens, Maria Branders,
 *       A Note on the Optimal Addition of Abscissas to Quadrature Formulas
 *       of Gauss and Lobatto, thematics of Computation,
 *       Volume 28, Number 125, January 1974, pages 135-139.
 *
 *  @param[in] n, the order of the Gauss rule.
 *  @param[in] eps, the requested absolute accuracy of the abscissas.
 *
 *  @param[out] x[n+1]   the abscissas.
 *  @param[out] w1[n+1]  the weights for the Gauss-Kronrod rule.
 *  @param[out] w2[n+1]  the weights for the Gauss rule.
 */
template<typename _Tp>
  void
  build_gauss_kronrod(int n, _Tp eps,
		      std::vector<_Tp>& x,
		      std::vector<_Tp>& w1, std::vector<_Tp>& w2)
  {
    x.resize(n + 1);
    w1.resize(n + 1);
    w2.resize(n + 1);

    std::vector<_Tp> b((n + 1) / 2 + 1);
    std::vector<_Tp> tau((n + 1) / 2);

    int m = (n + 1) / 2;
    bool even = (2 * m == n);

    auto d = _Tp{2};
    auto an = _Tp{0};
    for (int k = 1; k <= n; ++k)
      {
	an += _Tp{1};
	d *= an / (an + 0.5);
      }

    // Calculation of the Chebyshev coefficients of the orthogonal polynomial.
    tau[0] = (an + _Tp{2}) / (an + an + _Tp{3});
    b[m-1] = tau[0] - _Tp{1};
    auto ak = an;
    for (int l = 1; l < m; ++l)
      {
	ak += _Tp{2};
	tau[l] = ((ak - _Tp{1}) * ak 
	  - an * (an + _Tp{1})) * (ak + _Tp{2}) * tau[l-1] 
	  / (ak * ((ak + _Tp{3}) * (ak + _Tp{2}) 
	  - an * (an + _Tp{1})));
	b[m-l-1] = tau[l];

	for (int ll = 1; ll <= l; ++ll)
	  b[m-l-1] = b[m-l-1] + tau[ll-1] * b[m-l+ll-1];
      }
    b[m] = _Tp{1};

    // Calculation of approximate values for the abscissas.
    auto bb = std::sin(1.570796 / (an + an + _Tp{1}));
    auto x1 = std::sqrt(_Tp{1} - bb * bb);
    auto s = _Tp{2} * bb * x1;
    auto c = std::sqrt(_Tp{1} - s * s);
    auto coef = _Tp{1} - (_Tp{1} - _Tp{1} / an) / (_Tp{8} * an * an);
    auto xx = coef * x1;

    // Coefficient needed for weights.
    // coef2 = 2^{2*n+1} * n! * n! / (2n+1)! 
    //       = 2 * 4^n * n! / product((n+1)*...*(2*n+1))
    auto coef2 = _Tp{2} / _Tp(2 * n + 1);
    for (int i = 1; i <= n; ++i)
      coef2 *= _Tp{4} * _Tp(i) / _Tp(n + i);

    // Calculation of the k-th abscissa (a Kronrod abscissa) and the
    // corresponding weight.
    for (int k = 1; k <= n; k += 2)
      {
	get_kronrod(n, m, eps, coef2, even, b, xx, w1[k-1]);
	w2[k-1] = _Tp{0};
	x[k-1] = xx;
	auto y = x1;
	x1 = y * c - bb * s;
	bb = y * s + bb * c;

	if (k == n)
	  xx = _Tp{0};
	else
	  xx = coef * x1;

	// Calculation of the k+1 abscissa (a Gaussian abscissa) and the
	// corresponding weights.
	get_gauss(n, m, eps, coef2, even, b, xx, w1[k], w2[k]);

	x[k] = xx;
	y = x1;
	x1 = y * c - bb * s;
	bb = y * s + bb * c;
	xx = coef * x1;
      }

    // If N is even, we have one more Kronrod abscissa to compute,
    // namely the origin.
    if (even)
      {
	xx = _Tp{0};
	get_kronrod(n, m, eps, coef2, even, b, xx, w1[n]);
	w2[n] = _Tp{0};
	x[n] = xx;
      }

    return;
  }

int
main()
{
  std::cout.precision(std::numeric_limits<double>::max_digits10);
  auto w = 6 + std::cout.precision();

  const auto eps = 4 * std::numeric_limits<double>::epsilon();

  for (int n : {15, 21, 31, 41, 51, 61})
    {
      std::cout << '\n' << n << "-point Gauss-Kronrod rule\n";
      std::vector<double> x, w1, w2;
      build_gauss_kronrod(n, eps, x, w1, w2);
      for (int i = 0; i < n + 1; ++i)
	{
	  std::cout << ' ' << std::setw(w) << x[i]
		    << ' ' << std::setw(w) << w1[i]
		    << ' ' << std::setw(w) << w2[i] << '\n';
	}
    }
}

