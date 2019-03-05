// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2018-2019 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
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

namespace std
{
namespace __detail
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
  template<typename _Tp>
    void
    __get_kronrod(int __n, _Tp __eps, _Tp __coef2,
		  const std::vector<_Tp>& __b,
		  _Tp& __x, _Tp& __wk)
    {
      const int __max_iter = 100;

      const int __m = (__n + 1) / 2;
      const bool __even = (2 * __m == __n);

      int __ka = __x == _Tp{0} ? 1 : 0;

      // Iterative process for the computation of a Kronrod abscissa.
      auto __delta = _Tp{0};
      auto __fd = _Tp{0};
      for (int __iter = 1; __iter <= __max_iter; ++__iter)
      {

	_Tp __ai, __d2, __dif;
	if (__even)
	  {
	    __ai = __m + __m + 1;
	    __d2 = __ai * __b[__m];
	    __dif = _Tp{2};
	  }
	else
	  {
	    __ai = __m + 1;
	    __d2 = _Tp{0};
	    __dif = _Tp{1};
	  }

	auto __d1 = _Tp{0};
	auto __b0 = _Tp{0};
	auto __b1 = _Tp{0};
	auto __b2 = __b[__m];
	const auto __yy = _Tp{4} * __x * __x - _Tp{2};
	for (int __k = 1; __k <= __m; ++__k)
	  {
	    __ai -= __dif;
	    int __i = __m - __k + 1;
	    __b0 = __b1;
	    __b1 = __b2;
	    auto __d0 = __d1;
	    __d1 = __d2;
	    __b2 = __yy * __b1 - __b0 + __b[__i-1];
	    if (!__even)
	      ++__i;
	    __d2 = __yy * __d1 - __d0 + __ai * __b[__i - 1];
	  }

	_Tp __f;
	if (__even)
	  {
	    __f = __x * (__b2 - __b1);
	    __fd = __d2 + __d1;
	  }
	else
	  {
	    __f = _Tp{0.5L} * (__b2 - __b0);
	    __fd = _Tp{4} * __x * __d2;
	  }

	// Newton correction.
	__delta = __f / __fd;
	__x -= __delta;

	if (__ka == 1)
	  break;

	if (std::fabs(__delta) <= __eps)
	  __ka = 1;
      }

      // Catch non-convergence.
      if (__ka != 1)
	{
	  std::ostringstream __msg;
	  __msg << "\nget_kronrod: Iteration limit reached. eps = " << __eps
		<< "; Last delta was " << __delta << ".\n";
	  std::__throw_runtime_error(__msg.str().c_str());
	}

      // Computation of the weight.
      auto __d0 = _Tp{1};
      auto __d1 = __x;
      auto __d2 = _Tp{0};
      auto __ai = _Tp{0};
      for (int __k = 2; __k <= __n; ++__k)
	{
	  __ai += _Tp{1};
	  __d2 = ((__ai + __ai + _Tp{1}) * __x * __d1 - __ai * __d0)
	       / (__ai + _Tp{1});
	  __d0 = __d1;
	  __d1 = __d2;
	}

      __wk = __coef2 / (__fd * __d2);

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
  template<typename _Tp>
    void
    __get_gauss(int __n, _Tp __eps, _Tp __coef2,
		const std::vector<_Tp>& __b,
		_Tp& __x, _Tp& __wk, _Tp& __wg)
    {
      const int __max_iter = 100;

      const int __m = (__n + 1) / 2;
      const bool __even = (2 * __m == __n);

      int __ka = __x == _Tp{0} ? 1 : 0;

      // Iterative process for the computation of a Gaussian abscissa.
      auto __delta = _Tp{0};
      _Tp __p0, __p1, __p2, __pd2;
      for (int __iter = 1; __iter <= __max_iter; ++__iter)
	{
	  __p0 = _Tp{1};
	  __p1 = __x;
	  auto __pd0 = _Tp{0};
	  auto __pd1 = _Tp{1};

	  // When n is 1, we need to initialize p2 and pd2
	  // to avoid problems with delta.
	  if (__n <= 1)
	    {
	      if (std::numeric_limits<_Tp>::epsilon() < std::fabs(__x))
		{
		  __p2 = _Tp{0.5L} * (_Tp{3} * __x * __x - _Tp{1});
		  __pd2 = _Tp{3} * __x;
		}
	      else
		{
		  __p2 = _Tp{3} * __x;
		  __pd2 = _Tp{3};
		}
	    }

	  auto __ai = _Tp{0};
	  for (int __k = 2; __k <= __n; ++__k)
	    {
	      __ai += _Tp{1};
	      __p2 = ((__ai + __ai + _Tp{1}) * __x * __p1 - __ai * __p0)
		    / (__ai + _Tp{1});
	      __pd2 = ((__ai + __ai + _Tp{1})
			     * (__p1 + __x * __pd1) - __ai * __pd0)
		    / (__ai + _Tp{1});
	      __p0 = __p1;
	      __p1 = __p2;
	      __pd0 = __pd1;
	      __pd1 = __pd2;
	    }

	  // Newton correction.
	  __delta = __p2 / __pd2;
	  __x -= __delta;

	  if (__ka == 1)
	    break;

	  if (std::fabs(__delta) <= __eps)
	    __ka = 1;
	}

      // Catch non-convergence.
      if (__ka != 1)
	{
	  std::ostringstream __msg;
	  __msg << "__get_gauss: Iteration limit reached."<< " eps = " << __eps
		<< "; Last delta was " << __delta << ".\n";
	  std::__throw_runtime_error(__msg.str().c_str());
	}

      // Computation of the weights.
      __wg = _Tp{2} / (_Tp(__n) * __pd2 * __p0);

      __p1 = _Tp{0};
      __p2 = __b[__m];
      const auto __yy = _Tp{4} * __x * __x - _Tp{2};
      for (int __k = 1; __k <= __m; ++__k)
	{
	  const auto __i = __m - __k + 1;
	  __p0 = __p1;
	  __p1 = __p2;
	  __p2 = __yy * __p1 - __p0 + __b[__i - 1];
	}

      if (__even)
	__wk = __wg + __coef2 / (__pd2 * __x * (__p2 - __p1));
      else
	__wk = __wg + _Tp{2} * __coef2 / (__pd2 * (__p2 - __p0));

      return;
    }

} // namespace __detail
} // namespace std

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

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
  template<typename _Tp>
    void
    __build_gauss_kronrod(int __n, _Tp __eps,
			  std::vector<_Tp>& __x,
			  std::vector<_Tp>& __wg,
			  std::vector<_Tp>& __wk)
    {
      const int __m = (__n + 1) / 2;
      const bool __even = (2 * __m == __n);

      __x.resize(__n + 1);
      __wk.resize(__n + 1);
      __wg.resize(__m);

      std::vector<_Tp> __b(__m + 1);
      std::vector<_Tp> __tau(__m);

      auto __d = _Tp{2};
      auto __an = _Tp{0};
      for (int __k = 1; __k <= __n; ++__k)
	{
	  __an += _Tp{1};
	  __d *= __an / (__an + _Tp{0.5L});
	}

      // Calculation of the Chebyshev coefficients of the orthogonal polynomial.
      __tau[0] = (__an + _Tp{2}) / (__an + __an + _Tp{3});
      __b[__m - 1] = __tau[0] - _Tp{1};
      auto __ak = __an;
      for (int __l = 1; __l < __m; ++__l)
	{
	  __ak += _Tp{2};
	  __tau[__l] = ((__ak - _Tp{1}) * __ak
		    - __an * (__an + _Tp{1})) * (__ak + _Tp{2}) * __tau[__l - 1]
		    / (__ak * ((__ak + _Tp{3}) * (__ak + _Tp{2})
		    - __an * (__an + _Tp{1})));
	  __b[__m - __l - 1] = __tau[__l];
	  for (int __ll = 1; __ll <= __l; ++__ll)
	    __b[__m - __l - 1] += __tau[__ll - 1] * __b[__m - __l + __ll - 1];
	}
      __b[__m] = _Tp{1};

      // Calculation of approximate values for the abscissas.
      auto __bb = std::sin((_Tp{0.5L * M_PI}) / (__an + __an + _Tp{1}));
      auto __x1 = std::sqrt(_Tp{1} - __bb * __bb);
      auto __s = _Tp{2} * __bb * __x1;
      auto __c = std::sqrt(_Tp{1} - __s * __s);
      auto __coef = _Tp{1} - (_Tp{1} - _Tp{1} / __an) / (_Tp{8} * __an * __an);
      auto __xx = __coef * __x1;

      // Coefficient needed for weights.
      // coef2 = 2^{2n+1} n! n! / (2n+1)!
      //       = 2 4^n n! / product((n+1)...(2*n+1))
      auto __coef2 = _Tp{2} / _Tp(2 * __n + 1);
      for (int __i = 1; __i <= __n; ++__i)
	__coef2 *= _Tp{4} * _Tp(__i) / _Tp(__n + __i);

      // Calculation of the k-th abscissa (a Kronrod abscissa) and the
      // corresponding weight.
      for (int __k = 1; __k <= __n; __k += 2)
	{
	  std::__detail::__get_kronrod(__n, __eps, __coef2, __b,
					__xx, __wk[__k - 1]);
	  __wg[(__k - 1) / 2] = _Tp{0};
	  __x[__k - 1] = __xx;
	  auto __aa = __x1;
	  __x1 = __aa * __c - __bb * __s;
	  __bb = __aa * __s + __bb * __c;

	  __xx = __k == __n ? _Tp{0} : __coef * __x1;

	  // Calculation of the k+1 abscissa (a Gaussian abscissa) and the
	  // corresponding weights.
	  std::__detail::__get_gauss(__n, __eps, __coef2, __b,
				     __xx, __wk[__k], __wg[__k / 2]);

	  __x[__k] = __xx;
	  __aa = __x1;
	  __x1 = __aa * __c - __bb * __s;
	  __bb = __aa * __s + __bb * __c;
	  __xx = __coef * __x1;
	}

      // If n is even we must compute the Kronrod abscissa for the origin.
      if (__even)
	{
	  __xx = _Tp{0};
	  std::__detail::__get_kronrod(__n, __eps, __coef2, __b,
					__xx, __wk[__n]);
	  __x[__n] = __xx;
	}

      return;
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // GAUSS_KRONROD_RULE_TCC
