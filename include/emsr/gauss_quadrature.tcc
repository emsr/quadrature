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

#ifndef GAUSS_QUADRATURE_TCC
#define GAUSS_QUADRATURE_TCC 1

#include <stdexcept>
#include <type_traits>
#include <cmath>

namespace emsr
{

namespace detail
{

  /**
   * Build the absiscae and weights of a Gauss quadrature rule from
   * a symmetric tridiagonal Jacobi matrix using the Golub-Welsch technique.
   *
   * Sylvan Elhay, Jaroslav Kautsky,
   * Algorithm 655: IQPACK, Fortran Subroutines for the Weights of 
   * Interpolatory Quadrature,
   * ACM Transactions on Mathematical Software,
   * Volume 13, Number 4, December 1987, pages 399-415.
   *
   * @param[in]  n        The number of knots.
   * @param[in]  diag  The diagonal of the Jacobi matrix in [0..n-1].
   * @param[in,out]  subd  The subdiagonal of the Jacobi matrix in [1 ... n-1].
   *                       On output subd is overwritten.
   * @param[in]  moment0  The zero-th moment of the weight function.
   *
   * @param[out]  pt[n]  The points of the integration rule.
   * @param[out]  wt[n]  The weights of the integration rule.
   */
  template<typename Tp, typename InIter, typename OutIter>
    void
    golub_welsch(Tp moment0, int n, InIter& diag, InIter& subd,
		 OutIter& pt, OutIter& wt)
    {
      // Bail if the zero-th moment is not positive.
      if (moment0 <= Tp{0})
	throw std::domain_error("golub_welsch: moment0 <= 0");

      // Set up vectors for matrix diagonalization.
      for (int i = 0; i < n; ++i)
	pt[i] = diag[i];

      wt[0] = std::sqrt(moment0);
      for (int i = 1; i < n; ++i)
	wt[i] = Tp{0};

      // Diagonalize the Jacobi matrix.
      s_tridiag_symm(std::size_t(n), pt, subd, wt);

      for (int i = 0; i < n; ++i)
	wt[i] *= wt[i];

      return;
    }

} // namespace detail


  /**
   * Construct a Gauss-Legendre rule of order @c n.
   *
   * Weight function: @f$ 1 @f$
   * Constraints: @f$ b > a @f$
   */
  template<typename Tp>
    fixed_gauss_legendre_integral<Tp>::
    fixed_gauss_legendre_integral(int n)
    : order(n),
      point(n),
      weight(n)
    {
      const auto mu_0 = Tp{2};

      std::vector<Tp> diag(this->order, Tp{0});
      std::vector<Tp> subd(this->order);

      for (int i = 1; i <= this->order; ++i)
	subd[i - 1] = Tp(i) / std::sqrt(Tp(4 * i * i - 1));

      detail::golub_welsch(mu_0, this->order, diag, subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Legendre rule.
   *
   * Interval: @f$ (a, b) @f$
   */
  template<typename Tp>
    template<typename FuncTp>
      decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
      fixed_gauss_legendre_integral<Tp>::
      operator()(FuncTp func, Tp a, Tp b) const
      {
	using RetTp = std::invoke_result_t<FuncTp, Tp>;
	using AreaTp = decltype(RetTp{} * Tp{});

	auto sign = Tp{1};
	if (b < a)
	  {
	    sign = Tp{-1};
	    std::swap(a, b);
	  }

	const auto shift = (b + a) / Tp{2};
	const auto slope = (b - a) / Tp{2};
	const auto fact = slope;

	auto sum = AreaTp{0};
	for (int i = 0; i < this->order; ++i)
	  sum += this->weight[i]
		 * func(shift + slope * this->point[i]);

	return sign * fact * sum;
      }

  /**
   * Construct a Gauss-Chebyshev-T rule of order @c n.
   *
   * Interval: @f$ (-1, +1) @f$
   * Weight function: @f$ [(1-x)(1+x)]^{-1/2} @f$
   * Jacobi parameter: @f$ \alpha = 1/2 @f$
   * Jacobi parameter: @f$ \beta = 1/2 @f$
   */
  template<typename Tp>
    fixed_gauss_chebyshev_t_integral<Tp>::
    fixed_gauss_chebyshev_t_integral(int n)
    : order(n),
      point(n),
      weight(n)
    {
      const auto mu_0 = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};

      std::vector<Tp> diag(this->order, Tp{0});
      std::vector<Tp> subd(this->order, Tp{0.5L});

      subd[0] = std::sqrt(Tp{0.5L});

      detail::golub_welsch(mu_0, this->order, diag, subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Chebyshev-T rule.
   *
   * Weight function: @f$ ((b-x)(x-a))^{-1/2} @f$
   * Constraints: @f$ b > a @f$
   */
  template<typename Tp>
    template<typename FuncTp>
      decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
      fixed_gauss_chebyshev_t_integral<Tp>::
      operator()(FuncTp func, Tp a, Tp b) const
      {
	using RetTp = std::invoke_result_t<FuncTp, Tp>;
	using AreaTp = decltype(RetTp{} * Tp{});

	auto sign = Tp{1};
	if (b < a)
	  {
	    sign = Tp{-1};
	    std::swap(a, b);
	  }

	const auto shift = (b + a) / Tp{2};
	const auto slope = (b - a) / Tp{2};
	const auto fact = Tp{1};

	auto sum = AreaTp{0};
	for (int i = 0; i < this->order; ++i)
	  sum += this->weight[i]
		 * func(shift + slope * this->point[i]);

	return sign * fact * sum;
      }

  /**
   * Construct a Gauss-Chebyshev-U rule of order @c n.
   *
   * Interval: @f$ (a, b) @f$
   * Weight function: @f$ [(b-x)(x-a)]^{-1/2} @f$
   * Jacobi parameter: @f$ \alpha = -1/2 @f$
   * Jacobi parameter: @f$ \beta = -1/2 @f$
   */
  template<typename Tp>
    fixed_gauss_chebyshev_u_integral<Tp>::
    fixed_gauss_chebyshev_u_integral(int n)
    : order(n),
      point(n),
      weight(n)
    {
      const auto mu_0 = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L} / Tp{2};

      std::vector<Tp> diag(this->order, Tp{0});
      std::vector<Tp> subd(this->order, Tp{0.5L});

      detail::golub_welsch(mu_0, this->order, diag, subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Chebyshev-U rule.
   *
   * Interval: @f$ (a, b), b > a @f$
   * Weight function: @f$ [(b-x)(x-a)]^{1/2} @f$
   */
  template<typename Tp>
    template<typename FuncTp>
      decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
      fixed_gauss_chebyshev_u_integral<Tp>::
      operator()(FuncTp func, Tp a, Tp b) const
      {
	using RetTp = std::invoke_result_t<FuncTp, Tp>;
	using AreaTp = decltype(RetTp{} * Tp{});

	auto sign = Tp{1};
	if (b < a)
	  {
	    sign = Tp{-1};
	    std::swap(a, b);
	  }

	const auto shift = (b + a) / Tp{2};
	const auto slope = (b - a) / Tp{2};
	const auto fact = slope * slope;

	auto sum = AreaTp{0};
	for (int i = 0; i < this->order; ++i)
	  sum += this->weight[i]
		 * func(shift + slope * this->point[i]);

	return sign * fact * sum;
      }

  /**
   * Construct a Gauss-Chebyshev-V rule of order @c n.
   *
   * Interval: @f$ (-1, +1) @f$
   * Weight function: @f$ [(x-a)/(b-x)]^{1/2} @f$
   * Jacobi parameter: @f$ \alpha = +1/2 @f$
   * Jacobi parameter: @f$ \beta = -1/2 @f$
   */
  template<typename Tp>
    fixed_gauss_chebyshev_v_integral<Tp>::
    fixed_gauss_chebyshev_v_integral(int n)
    : order(n),
      point(n),
      weight(n)
    {
      const auto mu_0 = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};

      std::vector<Tp> diag(this->order, Tp{0});
      std::vector<Tp> subd(this->order, Tp{0.5L});

      diag[0] = Tp{-0.5L};
      subd[0] = std::sqrt(Tp{2});

      detail::golub_welsch(mu_0, this->order, diag, subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Chebyshev-V rule.
   *
   * Interval: @f$ (a, b), b > a @f$
   * Weight function: @f$ [(x-a)/(b-x)]^{1/2} @f$
   * 
   */
  template<typename Tp>
    template<typename FuncTp>
      decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
      fixed_gauss_chebyshev_v_integral<Tp>::
      operator()(FuncTp func, Tp a, Tp b) const
      {
	using RetTp = std::invoke_result_t<FuncTp, Tp>;
	using AreaTp = decltype(RetTp{} * Tp{});

	auto sign = Tp{1};
	if (b < a)
	  {
	    sign = Tp{-1};
	    std::swap(a, b);
	  }

	const auto shift = (b + a) / Tp{2};
	const auto slope = (b - a) / Tp{2};
	const auto fact = slope;

	auto sum = AreaTp{0};
	for (int i = 0; i < this->order; ++i)
	  sum += this->weight[i]
		 * func(shift + slope * this->point[i]);

	return sign * fact * sum;
      }

  /**
   * Construct a Gauss-Chebyshev-W rule of order @c n.
   *
   * Interval: @f$ (-1, +1) @f$
   * Weight function: @f$ [(1-x)/(1+x)]^{1/2} @f$
   * Jacobi parameter: @f$ \alpha = -1/2 @f$
   * Jacobi parameter: @f$ \beta = +1/2 @f$
   */
  template<typename Tp>
    fixed_gauss_chebyshev_w_integral<Tp>::
    fixed_gauss_chebyshev_w_integral(int n)
    : order(n),
      point(n),
      weight(n)
    {
      const auto mu_0 = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};

      std::vector<Tp> diag(this->order, Tp{0});
      std::vector<Tp> subd(this->order, Tp{0.5L});

      diag[0] = Tp{0.5L};
      subd[0] = std::sqrt(Tp{2});

      detail::golub_welsch(mu_0, this->order, diag, subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Chebyshev-W rule.
   *
   * Interval: @f$ (a, b), b > a @f$
   * Weight function: @f$ [(b-x)/(x-a)]^{1/2} @f$
   */
  template<typename Tp>
    template<typename FuncTp>
      decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
      fixed_gauss_chebyshev_w_integral<Tp>::
      operator()(FuncTp func, Tp a, Tp b) const
      {
	using RetTp = std::invoke_result_t<FuncTp, Tp>;
	using AreaTp = decltype(RetTp{} * Tp{});

	auto sign = Tp{1};
	if (b < a)
	  {
	    sign = Tp{-1};
	    std::swap(a, b);
	  }

	const auto shift = (b + a) / Tp{2};
	const auto slope = (b - a) / Tp{2};
	const auto fact = slope;

	auto sum = AreaTp{0};
	for (int i = 0; i < this->order; ++i)
	  sum += this->weight[i]
		 * func(shift + slope * this->point[i]);

	return sign * fact * sum;
      }

  /**
   * Construct a Gauss-Gegenbauer rule of order @c n.
   *
   * Weight function: @f$ [(b-x)(x-a)]^\lambda @f$
   * Constraints: @f$ \lambda > -1 @f$
   */
  template<typename Tp>
    fixed_gauss_gegenbauer_integral<Tp>::
    fixed_gauss_gegenbauer_integral(int n, Tp lam)
    : order(n),
      lambda(lam),
      point(n),
      weight(n)
    {
      const auto ab = Tp{2} * this->lambda;
      const auto gam = std::tgamma(this->lambda + Tp{1});
      const auto mu_0 = std::pow(Tp{2}, ab + Tp{1})
			* gam * gam / std::tgamma(ab + Tp{2});

      std::vector<Tp> diag(this->order, Tp{0});
      std::vector<Tp> subd(this->order);

      subd[0] = std::sqrt(Tp{1} / (Tp{2} * this->lambda + Tp{3}));

      for (int i = 2; i <= this->order; ++i)
	subd[i - 1] = std::sqrt(Tp(i) * (ab + Tp(i))
		        / (Tp{4}
			  * std::pow(this->lambda + Tp(i), Tp{2}) - Tp{1}));

      detail::golub_welsch(mu_0, this->order, diag, subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Gegenbauer rule.
   *
   * Interval: @f$ (a, b) @f$
   * Constraints: @f$ b > a @f$
   */
  template<typename Tp>
    template<typename FuncTp>
      decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
      fixed_gauss_gegenbauer_integral<Tp>::
      operator()(FuncTp func, Tp a, Tp b) const
      {
	using RetTp = std::invoke_result_t<FuncTp, Tp>;
	using AreaTp = decltype(RetTp{} * Tp{});

	auto sign = Tp{1};
	if (b < a)
	  {
	    sign = Tp{-1};
	    std::swap(a, b);
	  }

	const auto shift = (b + a) / Tp{2};
	const auto slope = (b - a) / Tp{2};
	const auto fact = std::pow(slope, Tp{2} * this->lambda + Tp{1});

	auto sum = AreaTp{0};
	for (int i = 0; i < this->order; ++i)
	  sum += this->weight[i]
		 * func(shift + slope * this->point[i]);

	return sign * fact * sum;
      }

  /**
   * Construct a Gauss-Jacobi rule of order @c n.
   *
   * Interval: @f$ (-1, +1) @f$
   * Weight function: @f$ (1-x)^\alpha (1+x)^\beta @f$
   * Constraints: @f$ \alpha, \beta > -1 @f$
   */
  template<typename Tp>
    fixed_gauss_jacobi_integral<Tp>::
    fixed_gauss_jacobi_integral(int n, Tp alf, Tp bet)
    : order(n),
      alpha(alf),
      beta(bet),
      point(n),
      weight(n)
    {
      const auto ab = this->alpha + this->beta;
      auto abp2i = ab + Tp{2};
      const auto mu_0 = std::pow(Tp{2}, ab + Tp{1})
			* std::tgamma(this->alpha + Tp{1})
			* std::tgamma(this->beta + Tp{1})
			/ std::tgamma(abp2i);

      std::vector<Tp> diag(this->order);
      std::vector<Tp> subd(this->order);

      diag[0] = (this->beta - this->alpha) / abp2i;
      subd[0] = Tp{2} * std::sqrt((this->alpha + Tp{1})
				   * (this->beta + Tp{1})
				   / (abp2i + Tp{1})) / abp2i;
      const auto a2mb2 = (this->beta - this->alpha)
			* (this->beta + this->alpha);
      for (int i = 1; i < this->order; ++i)
	{
	  const auto abp2ip2 = abp2i + Tp{2};
	  diag[i] = a2mb2 / abp2i / abp2ip2;
	  const auto ip1 = Tp(i + 1);
	  subd[i] = std::sqrt(Tp(4 * ip1) * (this->alpha + ip1)
				  * (this->beta + ip1) * (ab + ip1)
				  / (abp2ip2 * abp2ip2 - Tp{1})) / abp2ip2;
	  abp2i += Tp{2};
	}

      detail::golub_welsch(mu_0, this->order, diag, subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Jacobi rule.
   *
   * Interval: @f$ (a, b) @f$
   * Weight function: @f$ (b-x)^\alpha (x-a)^\beta @f$
   * Constraints: @f$ b > a @f$
   */
  template<typename Tp>
    template<typename FuncTp>
      decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
      fixed_gauss_jacobi_integral<Tp>::
      operator()(FuncTp func, Tp a, Tp b) const
      {
	using RetTp = std::invoke_result_t<FuncTp, Tp>;
	using AreaTp = decltype(RetTp{} * Tp{});

	auto sign = Tp{1};
	if (b < a)
	  {
	    sign = Tp{-1};
	    std::swap(a, b);
	  }

	const auto shift = (b + a) / Tp{2};
	const auto slope = (b - a) / Tp{2};
	const auto fact = std::pow(slope, this->alpha + this->beta + Tp{1});

	auto sum = AreaTp{0};
	for (int i = 0; i < this->order; ++i)
	  sum += this->weight[i]
		 * func(shift + slope * this->point[i]);

	return sign * fact * sum;
      }

  /**
   * Construct a Gauss-Laguerre rule of order @c n.
   *
   * Interval: @f$ (0, \infty) @f$
   * Weight function: @f$ x^\alpha \exp(-b(x-a)) @f$
   * Constraints: @f$ \alpha > -1 @f$
   */
  template<typename Tp>
    fixed_gauss_laguerre_integral<Tp>::
    fixed_gauss_laguerre_integral(int n, Tp alf)
    : order(n),
      alpha(alf),
      point(n),
      weight(n)
    {
      const auto mu_0 = std::tgamma(this->alpha + Tp{1});

      std::vector<Tp> diag(this->order);
      std::vector<Tp> subd(this->order);

      for (int i = 0; i < this->order; ++i)
	{
	  diag[i] = Tp(2 * i + 1) + this->alpha;
	  subd[i] = std::sqrt(Tp(i + 1) * (this->alpha + Tp(i + 1)));
	}

      detail::golub_welsch(mu_0, this->order, diag, subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Laguerre rule.
   *
   * Interval: @f$ (a, \infty) @f$
   * Weight function: @f$ (x-a)^\alpha \exp(-b(x-a)) @f$
   * Constraints: @f$ b > 0 @f$
   */
  template<typename Tp>
    template<typename FuncTp>
      decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
      fixed_gauss_laguerre_integral<Tp>::
      operator()(FuncTp func, Tp a, Tp b) const
      {
	using RetTp = std::invoke_result_t<FuncTp, Tp>;
	using AreaTp = decltype(RetTp{} * Tp{});

	const auto shift = a;
	const auto slope = Tp{1} / b;
	const auto fact = std::pow(slope, this->alpha + Tp{1});

	auto sum = AreaTp{0};
	for (int i = 0; i < this->order; ++i)
	  sum += this->weight[i]
		 * func(shift + slope * this->point[i]);

	return fact * sum;
      }

  /**
   * Construct a Generalized Gauss-Hermite rule of order @c n.
   *
   * Interval: @f$ (-\infty, \infty) @f$
   * Weight function: @f$ |x-a|^\alpha \exp(-b(x-a)^2) @f$
   * Constraints: @f$ \alpha > -1 @f$
   */
  template<typename Tp>
    fixed_gauss_hermite_integral<Tp>::
    fixed_gauss_hermite_integral(int n, Tp alf)
    : order(n),
      alpha(alf),
      point(n),
      weight(n)
    {
      const auto mu_0 = std::tgamma((this->alpha + Tp{1}) / Tp{2});

      std::vector<Tp> diag(this->order, Tp{0});
      std::vector<Tp> subd(this->order);

      for (int i = 1; i <= this->order; ++i)
	subd[i - 1] = std::sqrt((Tp(i) + Tp(i % 2) * this->alpha)
				   / Tp{2});

      detail::golub_welsch(mu_0, this->order, diag, subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Generalized Gauss-Hermite rule.
   *
   * Interval: @f$ (-\infty, \infty) @f$
   * Weight function: @f$ |x-a|^\alpha \exp(-b(x-a)^2) @f$
   * Constraints: @f$ b > 0 @f$
   */
  template<typename Tp>
    template<typename FuncTp>
      decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
      fixed_gauss_hermite_integral<Tp>::
      operator()(FuncTp func, Tp a, Tp b) const
      {
	using RetTp = std::invoke_result_t<FuncTp, Tp>;
	using AreaTp = decltype(RetTp{} * Tp{});

	const auto shift = a;
	const auto slope = Tp{1} / std::sqrt(b);
	const auto fact = std::pow(slope, this->alpha + Tp{1});

	auto sum = AreaTp{0};
	for (int i = 0; i < this->order; ++i)
	  sum += this->weight[i]
		 * func(shift + slope * this->point[i]);

	return fact * sum;
      }

  /**
   * Construct a Generalized Gauss-Exponential rule of order @c n.
   *
   * Interval: @f$ (a, b) @f$
   * Weight function: @f$ |x - (a + b) / 2|^\alpha @f$
   * Constraints: @f$ \alpha > -1 @f$
   */
  template<typename Tp>
    fixed_gauss_exponential_integral<Tp>::
    fixed_gauss_exponential_integral(int n, Tp alf)
    : order(n),
      alpha(alf),
      point(n),
      weight(n)
    {
      const auto mu_0 = Tp{2} / (this->alpha + Tp{1});

      std::vector<Tp> diag(this->order, Tp{0});
      std::vector<Tp> subd(this->order);

      auto ap2i = this->alpha;
      for (int i = 1; i <= this->order; ++i)
	{
	  ap2i += Tp{2};
	  subd[i - 1] = (Tp(i) + Tp(i % 2) * this->alpha)
			  / std::sqrt((ap2i * ap2i - Tp{1}));
	}

      detail::golub_welsch(mu_0, this->order, diag, subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Generalized Gauss-Exponential rule.
   *
   * Constraints: @f$ b > a @f$
   */
  template<typename Tp>
    template<typename FuncTp>
      decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
      fixed_gauss_exponential_integral<Tp>::
      operator()(FuncTp func, Tp a, Tp b) const
      {
	using RetTp = std::invoke_result_t<FuncTp, Tp>;
	using AreaTp = decltype(RetTp{} * Tp{});

	auto sign = Tp{1};
	if (b < a)
	  {
	    sign = Tp{-1};
	    std::swap(a, b);
	  }

	const auto shift = (b + a) / Tp{2};
	const auto slope = (b - a) / Tp{2};
	const auto fact = std::pow(slope, this->alpha + Tp{1});

	auto sum = AreaTp{0};
	for (int i = 0; i < this->order; ++i)
	  sum += this->weight[i]
		 * func(shift + slope * this->point[i]);

	return sign * fact * sum;
      }

  /**
   * Construct a Gauss-Rational rule of order @c n.
   *
   * Weight function: @f$ (x - a)^\alpha (b + x)^\beta @f$
   * Constraints: @f$ \alpha > -1, \alpha + \beta + 2n < 0
   */
  template<typename Tp>
    fixed_gauss_rational_integral<Tp>::
    fixed_gauss_rational_integral(int n, Tp alf, Tp bet)
    : order(n),
      alpha(alf),
      beta(bet),
      point(n),
      weight(n)
    {
      const auto ab = this->alpha + this->beta;
      const auto mu_0 = std::tgamma(this->alpha + Tp{1})
			* std::tgamma(-(ab + Tp{1}))
			/ std::tgamma(-this->beta);
      const auto ap1 = this->alpha + Tp{1};
      const auto aba = ab * ap1;

      std::vector<Tp> diag(this->order);
      std::vector<Tp> subd(this->order);

      diag[0] = -ap1 / (ab + Tp{2});
      subd[0] = -diag[0] * (this->beta + Tp{1})
		/ (ab + Tp{2}) / (ab + Tp{3});
      for (int i = 2; i <= this->order; ++i)
	{
	  const auto abp2i = ab + Tp(2 * i);
	  diag[i - 1] = aba + Tp{2} * (ab + Tp(i)) * Tp(i - 1);
	  diag[i - 1] = -diag[i-1] / abp2i / (abp2i - Tp{2});
	}

      for (int i = 2; i <= this->order - 1; ++i)
	{
	  const auto abp2i = ab + Tp(2 * i);
	  subd[i - 1] = Tp(i) * (this->alpha + Tp(i))
			  / (abp2i - Tp{1}) * (this->beta + Tp(i))
			  / (abp2i * abp2i) * (ab + Tp(i))
			  / (abp2i + Tp{1});
	}
      subd[this->order - 1] = Tp{0};
      for (int i = 0; i < this->order; ++i)
	subd[i] = std::sqrt(subd[i]);

      detail::golub_welsch(mu_0, this->order, diag, subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Rational rule.
   *
   * Interval: @f$ (a, \infty) @f$
   * Constraints: @f$ a + b > 0 @f$
   */
  template<typename Tp>
    template<typename FuncTp>
      decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
      fixed_gauss_rational_integral<Tp>::
      operator()(FuncTp func, Tp a, Tp b) const
      {
	using RetTp = std::invoke_result_t<FuncTp, Tp>;
	using AreaTp = decltype(RetTp{} * Tp{});

	auto sign = Tp{1};
	if (b < a)
	  {
	    sign = Tp{-1};
	    std::swap(a, b);
	  }

	const auto shift = a;
	const auto slope = b + a;
	const auto fact = std::pow(slope, this->alpha + this->beta + Tp{1});

	auto sum = AreaTp{0};
	for (int i = 0; i < this->order; ++i)
	  sum += this->weight[i]
		 * func(shift + slope * this->point[i]);

	return sign * fact * sum;
      }

} // namespace emsr

#endif // GAUSS_QUADRATURE_TCC
