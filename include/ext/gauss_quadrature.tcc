// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2018-2020 Free Software Foundation, Inc.
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

#ifndef GAUSS_QUADRATURE_TCC
#define GAUSS_QUADRATURE_TCC 1

#include <type_traits>
#include <cmath>

namespace __gnu_cxx
{

namespace __detail
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
  template<typename _Tp, typename _InIter, typename _OutIter>
    void
    __golub_welsch(_Tp __moment0, int __n, _InIter& __diag, _InIter& __subd,
		   _OutIter& __pt, _OutIter& __wt)
    {
      // Bail if the zero-th moment is not positive.
      if (__moment0 <= _Tp{0})
	std::__throw_domain_error("__golub_welsch: moment0 <= 0");

      // Set up vectors for matrix diagonalization.
      for (int i = 0; i < __n; ++i)
	__pt[i] = __diag[i];

      __wt[0] = std::sqrt(__moment0);
      for (int i = 1; i < __n; ++i)
	__wt[i] = _Tp{0};

      // Diagonalize the Jacobi matrix.
      _S_tridiag_symm(std::size_t(__n), __pt, __subd, __wt);

      for (int i = 0; i < __n; ++i)
	__wt[i] *= __wt[i];

      return;
    }

} // namespace __detail


  /**
   * Construct a Gauss-Legendre rule of order @c n.
   *
   * Weight function: @f$ 1 @f$
   * Constraints: @f$ b > a @f$
   */
  template<typename _Tp>
    fixed_gauss_legendre_integral<_Tp>::
    fixed_gauss_legendre_integral(int __n)
    : order(__n),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = _Tp{2};

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order);

      for (int __i = 1; __i <= this->order; ++__i)
	__subd[__i - 1] = _Tp(__i) / std::sqrt(_Tp(4 * __i * __i - 1));

      __detail::__golub_welsch(__mu_0, this->order, __diag, __subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Legendre rule.
   *
   * Interval: @f$ (a, b) @f$
   */
  template<typename _Tp>
    template<typename _FuncTp>
      decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
      fixed_gauss_legendre_integral<_Tp>::
      operator()(_FuncTp __func, _Tp __a, _Tp __b) const
      {
	using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
	using _AreaTp = decltype(_RetTp{} * _Tp{});

	auto __sign = _Tp{1};
	if (__b < __a)
	  {
	    __sign = _Tp{-1};
	    std::swap(__a, __b);
	  }

	const auto __shift = (__b + __a) / _Tp{2};
	const auto __slope = (__b - __a) / _Tp{2};
	const auto __fact = __slope;

	auto __sum = _AreaTp{0};
	for (int __i = 0; __i < this->order; ++__i)
	  __sum += this->weight[__i]
		 * __func(__shift + __slope * this->point[__i]);

	return __sign * __fact * __sum;
      }

  /**
   * Construct a Gauss-Chebyshev-T rule of order @c n.
   *
   * Interval: @f$ (-1, +1) @f$
   * Weight function: @f$ [(1-x)(1+x)]^{-1/2} @f$
   * Jacobi parameter: @f$ \alpha = 1/2 @f$
   * Jacobi parameter: @f$ \beta = 1/2 @f$
   */
  template<typename _Tp>
    fixed_gauss_chebyshev_t_integral<_Tp>::
    fixed_gauss_chebyshev_t_integral(int __n)
    : order(__n),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order, _Tp{0.5L});

      __subd[0] = std::sqrt(_Tp{0.5L});

      __detail::__golub_welsch(__mu_0, this->order, __diag, __subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Chebyshev-T rule.
   *
   * Weight function: @f$ ((b-x)(x-a))^{-1/2} @f$
   * Constraints: @f$ b > a @f$
   */
  template<typename _Tp>
    template<typename _FuncTp>
      decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
      fixed_gauss_chebyshev_t_integral<_Tp>::
      operator()(_FuncTp __func, _Tp __a, _Tp __b) const
      {
	using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
	using _AreaTp = decltype(_RetTp{} * _Tp{});

	auto __sign = _Tp{1};
	if (__b < __a)
	  {
	    __sign = _Tp{-1};
	    std::swap(__a, __b);
	  }

	const auto __shift = (__b + __a) / _Tp{2};
	const auto __slope = (__b - __a) / _Tp{2};
	const auto __fact = _Tp{1};

	auto __sum = _AreaTp{0};
	for (int __i = 0; __i < this->order; ++__i)
	  __sum += this->weight[__i]
		 * __func(__shift + __slope * this->point[__i]);

	return __sign * __fact * __sum;
      }

  /**
   * Construct a Gauss-Chebyshev-U rule of order @c n.
   *
   * Interval: @f$ (a, b) @f$
   * Weight function: @f$ [(b-x)(x-a)]^{-1/2} @f$
   * Jacobi parameter: @f$ \alpha = -1/2 @f$
   * Jacobi parameter: @f$ \beta = -1/2 @f$
   */
  template<typename _Tp>
    fixed_gauss_chebyshev_u_integral<_Tp>::
    fixed_gauss_chebyshev_u_integral(int __n)
    : order(__n),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L} / _Tp{2};

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order, _Tp{0.5L});

      __detail::__golub_welsch(__mu_0, this->order, __diag, __subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Chebyshev-U rule.
   *
   * Interval: @f$ (a, b), b > a @f$
   * Weight function: @f$ [(b-x)(x-a)]^{1/2} @f$
   */
  template<typename _Tp>
    template<typename _FuncTp>
      decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
      fixed_gauss_chebyshev_u_integral<_Tp>::
      operator()(_FuncTp __func, _Tp __a, _Tp __b) const
      {
	using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
	using _AreaTp = decltype(_RetTp{} * _Tp{});

	auto __sign = _Tp{1};
	if (__b < __a)
	  {
	    __sign = _Tp{-1};
	    std::swap(__a, __b);
	  }

	const auto __shift = (__b + __a) / _Tp{2};
	const auto __slope = (__b - __a) / _Tp{2};
	const auto __fact = __slope * __slope;

	auto __sum = _AreaTp{0};
	for (int __i = 0; __i < this->order; ++__i)
	  __sum += this->weight[__i]
		 * __func(__shift + __slope * this->point[__i]);

	return __sign * __fact * __sum;
      }

  /**
   * Construct a Gauss-Chebyshev-V rule of order @c n.
   *
   * Interval: @f$ (-1, +1) @f$
   * Weight function: @f$ [(x-a)/(b-x)]^{1/2} @f$
   * Jacobi parameter: @f$ \alpha = +1/2 @f$
   * Jacobi parameter: @f$ \beta = -1/2 @f$
   */
  template<typename _Tp>
    fixed_gauss_chebyshev_v_integral<_Tp>::
    fixed_gauss_chebyshev_v_integral(int __n)
    : order(__n),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order, _Tp{0.5L});

      __diag[0] = _Tp{-0.5L};
      __subd[0] = std::sqrt(_Tp{2});

      __detail::__golub_welsch(__mu_0, this->order, __diag, __subd,
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
  template<typename _Tp>
    template<typename _FuncTp>
      decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
      fixed_gauss_chebyshev_v_integral<_Tp>::
      operator()(_FuncTp __func, _Tp __a, _Tp __b) const
      {
	using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
	using _AreaTp = decltype(_RetTp{} * _Tp{});

	auto __sign = _Tp{1};
	if (__b < __a)
	  {
	    __sign = _Tp{-1};
	    std::swap(__a, __b);
	  }

	const auto __shift = (__b + __a) / _Tp{2};
	const auto __slope = (__b - __a) / _Tp{2};
	const auto __fact = __slope;

	auto __sum = _AreaTp{0};
	for (int __i = 0; __i < this->order; ++__i)
	  __sum += this->weight[__i]
		 * __func(__shift + __slope * this->point[__i]);

	return __sign * __fact * __sum;
      }

  /**
   * Construct a Gauss-Chebyshev-W rule of order @c n.
   *
   * Interval: @f$ (-1, +1) @f$
   * Weight function: @f$ [(1-x)/(1+x)]^{1/2} @f$
   * Jacobi parameter: @f$ \alpha = -1/2 @f$
   * Jacobi parameter: @f$ \beta = +1/2 @f$
   */
  template<typename _Tp>
    fixed_gauss_chebyshev_w_integral<_Tp>::
    fixed_gauss_chebyshev_w_integral(int __n)
    : order(__n),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order, _Tp{0.5L});

      __diag[0] = _Tp{0.5L};
      __subd[0] = std::sqrt(_Tp{2});

      __detail::__golub_welsch(__mu_0, this->order, __diag, __subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Chebyshev-W rule.
   *
   * Interval: @f$ (a, b), b > a @f$
   * Weight function: @f$ [(b-x)/(x-a)]^{1/2} @f$
   */
  template<typename _Tp>
    template<typename _FuncTp>
      decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
      fixed_gauss_chebyshev_w_integral<_Tp>::
      operator()(_FuncTp __func, _Tp __a, _Tp __b) const
      {
	using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
	using _AreaTp = decltype(_RetTp{} * _Tp{});

	auto __sign = _Tp{1};
	if (__b < __a)
	  {
	    __sign = _Tp{-1};
	    std::swap(__a, __b);
	  }

	const auto __shift = (__b + __a) / _Tp{2};
	const auto __slope = (__b - __a) / _Tp{2};
	const auto __fact = __slope;

	auto __sum = _AreaTp{0};
	for (int __i = 0; __i < this->order; ++__i)
	  __sum += this->weight[__i]
		 * __func(__shift + __slope * this->point[__i]);

	return __sign * __fact * __sum;
      }

  /**
   * Construct a Gauss-Gegenbauer rule of order @c n.
   *
   * Weight function: @f$ [(b-x)(x-a)]^\lambda @f$
   * Constraints: @f$ \lambda > -1 @f$
   */
  template<typename _Tp>
    fixed_gauss_gegenbauer_integral<_Tp>::
    fixed_gauss_gegenbauer_integral(int __n, _Tp __lam)
    : order(__n),
      lambda(__lam),
      point(__n),
      weight(__n)
    {
      const auto __ab = _Tp{2} * this->lambda;
      const auto __gam = std::tgamma(this->lambda + _Tp{1});
      const auto __mu_0 = std::pow(_Tp{2}, __ab + _Tp{1})
			* __gam * __gam / std::tgamma(__ab + _Tp{2});

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order);

      __subd[0] = std::sqrt(_Tp{1} / (_Tp{2} * this->lambda + _Tp{3}));

      for (int __i = 2; __i <= this->order; ++__i)
	__subd[__i - 1] = std::sqrt(_Tp(__i) * (__ab + _Tp(__i))
		        / (_Tp{4}
			  * std::pow(this->lambda + _Tp(__i), _Tp{2}) - _Tp{1}));

      __detail::__golub_welsch(__mu_0, this->order, __diag, __subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Gegenbauer rule.
   *
   * Interval: @f$ (a, b) @f$
   * Constraints: @f$ b > a @f$
   */
  template<typename _Tp>
    template<typename _FuncTp>
      decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
      fixed_gauss_gegenbauer_integral<_Tp>::
      operator()(_FuncTp __func, _Tp __a, _Tp __b) const
      {
	using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
	using _AreaTp = decltype(_RetTp{} * _Tp{});

	auto __sign = _Tp{1};
	if (__b < __a)
	  {
	    __sign = _Tp{-1};
	    std::swap(__a, __b);
	  }

	const auto __shift = (__b + __a) / _Tp{2};
	const auto __slope = (__b - __a) / _Tp{2};
	const auto __fact = std::pow(__slope, _Tp{2} * this->lambda + _Tp{1});

	auto __sum = _AreaTp{0};
	for (int __i = 0; __i < this->order; ++__i)
	  __sum += this->weight[__i]
		 * __func(__shift + __slope * this->point[__i]);

	return __sign * __fact * __sum;
      }

  /**
   * Construct a Gauss-Jacobi rule of order @c n.
   *
   * Interval: @f$ (-1, +1) @f$
   * Weight function: @f$ (1-x)^\alpha (1+x)^\beta @f$
   * Constraints: @f$ \alpha, \beta > -1 @f$
   */
  template<typename _Tp>
    fixed_gauss_jacobi_integral<_Tp>::
    fixed_gauss_jacobi_integral(int __n, _Tp __alf, _Tp __bet)
    : order(__n),
      alpha(__alf),
      beta(__bet),
      point(__n),
      weight(__n)
    {
      const auto __ab = this->alpha + this->beta;
      auto __abp2i = __ab + _Tp{2};
      const auto __mu_0 = std::pow(_Tp{2}, __ab + _Tp{1})
			* std::tgamma(this->alpha + _Tp{1})
			* std::tgamma(this->beta + _Tp{1})
			/ std::tgamma(__abp2i);

      std::vector<_Tp> __diag(this->order);
      std::vector<_Tp> __subd(this->order);

      __diag[0] = (this->beta - this->alpha) / __abp2i;
      __subd[0] = _Tp{2} * std::sqrt((this->alpha + _Tp{1})
				   * (this->beta + _Tp{1})
				   / (__abp2i + _Tp{1})) / __abp2i;
      const auto a2mb2 = (this->beta - this->alpha)
			* (this->beta + this->alpha);
      for (int __i = 1; __i < this->order; ++__i)
	{
	  const auto abp2ip2 = __abp2i + _Tp{2};
	  __diag[__i] = a2mb2 / __abp2i / abp2ip2;
	  const auto __ip1 = _Tp(__i + 1);
	  __subd[__i] = std::sqrt(_Tp(4 * __ip1) * (this->alpha + __ip1)
				  * (this->beta + __ip1) * (__ab + __ip1)
				  / (abp2ip2 * abp2ip2 - _Tp{1})) / abp2ip2;
	  __abp2i += _Tp{2};
	}

      __detail::__golub_welsch(__mu_0, this->order, __diag, __subd,
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
  template<typename _Tp>
    template<typename _FuncTp>
      decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
      fixed_gauss_jacobi_integral<_Tp>::
      operator()(_FuncTp __func, _Tp __a, _Tp __b) const
      {
	using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
	using _AreaTp = decltype(_RetTp{} * _Tp{});

	auto __sign = _Tp{1};
	if (__b < __a)
	  {
	    __sign = _Tp{-1};
	    std::swap(__a, __b);
	  }

	const auto __shift = (__b + __a) / _Tp{2};
	const auto __slope = (__b - __a) / _Tp{2};
	const auto __fact = std::pow(__slope, this->alpha + this->beta + _Tp{1});

	auto __sum = _AreaTp{0};
	for (int __i = 0; __i < this->order; ++__i)
	  __sum += this->weight[__i]
		 * __func(__shift + __slope * this->point[__i]);

	return __sign * __fact * __sum;
      }

  /**
   * Construct a Gauss-Laguerre rule of order @c n.
   *
   * Interval: @f$ (0, \infty) @f$
   * Weight function: @f$ x^\alpha \exp(-b(x-a)) @f$
   * Constraints: @f$ \alpha > -1 @f$
   */
  template<typename _Tp>
    fixed_gauss_laguerre_integral<_Tp>::
    fixed_gauss_laguerre_integral(int __n, _Tp __alf)
    : order(__n),
      alpha(__alf),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = std::tgamma(this->alpha + _Tp{1});

      std::vector<_Tp> __diag(this->order);
      std::vector<_Tp> __subd(this->order);

      for (int __i = 0; __i < this->order; ++__i)
	{
	  __diag[__i] = _Tp(2 * __i + 1) + this->alpha;
	  __subd[__i] = std::sqrt(_Tp(__i + 1) * (this->alpha + _Tp(__i + 1)));
	}

      __detail::__golub_welsch(__mu_0, this->order, __diag, __subd,
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
  template<typename _Tp>
    template<typename _FuncTp>
      decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
      fixed_gauss_laguerre_integral<_Tp>::
      operator()(_FuncTp __func, _Tp __a, _Tp __b) const
      {
	using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
	using _AreaTp = decltype(_RetTp{} * _Tp{});

	const auto __shift = __a;
	const auto __slope = _Tp{1} / __b;
	const auto __fact = std::pow(__slope, this->alpha + _Tp{1});

	auto __sum = _AreaTp{0};
	for (int __i = 0; __i < this->order; ++__i)
	  __sum += this->weight[__i]
		 * __func(__shift + __slope * this->point[__i]);

	return __fact * __sum;
      }

  /**
   * Construct a Generalized Gauss-Hermite rule of order @c n.
   *
   * Interval: @f$ (-\infty, \infty) @f$
   * Weight function: @f$ |x-a|^\alpha \exp(-b(x-a)^2) @f$
   * Constraints: @f$ \alpha > -1 @f$
   */
  template<typename _Tp>
    fixed_gauss_hermite_integral<_Tp>::
    fixed_gauss_hermite_integral(int __n, _Tp __alf)
    : order(__n),
      alpha(__alf),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = std::tgamma((this->alpha + _Tp{1}) / _Tp{2});

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order);

      for (int __i = 1; __i <= this->order; ++__i)
	__subd[__i - 1] = std::sqrt((_Tp(__i) + _Tp(__i % 2) * this->alpha)
				   / _Tp{2});

      __detail::__golub_welsch(__mu_0, this->order, __diag, __subd,
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
  template<typename _Tp>
    template<typename _FuncTp>
      decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
      fixed_gauss_hermite_integral<_Tp>::
      operator()(_FuncTp __func, _Tp __a, _Tp __b) const
      {
	using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
	using _AreaTp = decltype(_RetTp{} * _Tp{});

	const auto __shift = __a;
	const auto __slope = _Tp{1} / std::sqrt(__b);
	const auto __fact = std::pow(__slope, this->alpha + _Tp{1});

	auto __sum = _AreaTp{0};
	for (int __i = 0; __i < this->order; ++__i)
	  __sum += this->weight[__i]
		 * __func(__shift + __slope * this->point[__i]);

	return __fact * __sum;
      }

  /**
   * Construct a Generalized Gauss-Exponential rule of order @c n.
   *
   * Interval: @f$ (a, b) @f$
   * Weight function: @f$ |x - (a + b) / 2|^\alpha @f$
   * Constraints: @f$ \alpha > -1 @f$
   */
  template<typename _Tp>
    fixed_gauss_exponential_integral<_Tp>::
    fixed_gauss_exponential_integral(int __n, _Tp __alf)
    : order(__n),
      alpha(__alf),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = _Tp{2} / (this->alpha + _Tp{1});

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order);

      auto __ap2i = this->alpha;
      for (int __i = 1; __i <= this->order; ++__i)
	{
	  __ap2i += _Tp{2};
	  __subd[__i - 1] = (_Tp(__i) + _Tp(__i % 2) * this->alpha)
			  / std::sqrt((__ap2i * __ap2i - _Tp{1}));
	}

      __detail::__golub_welsch(__mu_0, this->order, __diag, __subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Generalized Gauss-Exponential rule.
   *
   * Constraints: @f$ b > a @f$
   */
  template<typename _Tp>
    template<typename _FuncTp>
      decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
      fixed_gauss_exponential_integral<_Tp>::
      operator()(_FuncTp __func, _Tp __a, _Tp __b) const
      {
	using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
	using _AreaTp = decltype(_RetTp{} * _Tp{});

	auto __sign = _Tp{1};
	if (__b < __a)
	  {
	    __sign = _Tp{-1};
	    std::swap(__a, __b);
	  }

	const auto __shift = (__b + __a) / _Tp{2};
	const auto __slope = (__b - __a) / _Tp{2};
	const auto __fact = std::pow(__slope, this->alpha + _Tp{1});

	auto __sum = _AreaTp{0};
	for (int __i = 0; __i < this->order; ++__i)
	  __sum += this->weight[__i]
		 * __func(__shift + __slope * this->point[__i]);

	return __sign * __fact * __sum;
      }

  /**
   * Construct a Gauss-Rational rule of order @c n.
   *
   * Weight function: @f$ (x - a)^\alpha (b + x)^\beta @f$
   * Constraints: @f$ \alpha > -1, \alpha + \beta + 2n < 0
   */
  template<typename _Tp>
    fixed_gauss_rational_integral<_Tp>::
    fixed_gauss_rational_integral(int __n, _Tp __alf, _Tp __bet)
    : order(__n),
      alpha(__alf),
      beta(__bet),
      point(__n),
      weight(__n)
    {
      const auto __ab = this->alpha + this->beta;
      const auto __mu_0 = std::tgamma(this->alpha + _Tp{1})
			* std::tgamma(-(__ab + _Tp{1}))
			/ std::tgamma(-this->beta);
      const auto __ap1 = this->alpha + _Tp{1};
      const auto __aba = __ab * __ap1;

      std::vector<_Tp> __diag(this->order);
      std::vector<_Tp> __subd(this->order);

      __diag[0] = -__ap1 / (__ab + _Tp{2});
      __subd[0] = -__diag[0] * (this->beta + _Tp{1})
		/ (__ab + _Tp{2}) / (__ab + _Tp{3});
      for (int __i = 2; __i <= this->order; ++__i)
	{
	  const auto __abp2i = __ab + _Tp(2 * __i);
	  __diag[__i - 1] = __aba + _Tp{2} * (__ab + _Tp(__i)) * _Tp(__i - 1);
	  __diag[__i - 1] = -__diag[__i-1] / __abp2i / (__abp2i - _Tp{2});
	}

      for (int __i = 2; __i <= this->order - 1; ++__i)
	{
	  const auto __abp2i = __ab + _Tp(2 * __i);
	  __subd[__i - 1] = _Tp(__i) * (this->alpha + _Tp(__i))
			  / (__abp2i - _Tp{1}) * (this->beta + _Tp(__i))
			  / (__abp2i * __abp2i) * (__ab + _Tp(__i))
			  / (__abp2i + _Tp{1});
	}
      __subd[this->order - 1] = _Tp{0};
      for (int __i = 0; __i < this->order; ++__i)
	__subd[__i] = std::sqrt(__subd[__i]);

      __detail::__golub_welsch(__mu_0, this->order, __diag, __subd,
			       this->point, this->weight);
    }

  /**
   * Evaluate the integral of a function @c func from limits @c a to @c b
   * using the Gauss-Rational rule.
   *
   * Interval: @f$ (a, \infty) @f$
   * Constraints: @f$ a + b > 0 @f$
   */
  template<typename _Tp>
    template<typename _FuncTp>
      decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
      fixed_gauss_rational_integral<_Tp>::
      operator()(_FuncTp __func, _Tp __a, _Tp __b) const
      {
	using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
	using _AreaTp = decltype(_RetTp{} * _Tp{});

	auto __sign = _Tp{1};
	if (__b < __a)
	  {
	    __sign = _Tp{-1};
	    std::swap(__a, __b);
	  }

	const auto __shift = __a;
	const auto __slope = __b + __a;
	const auto __fact = std::pow(__slope, this->alpha + this->beta + _Tp{1});

	auto __sum = _AreaTp{0};
	for (int __i = 0; __i < this->order; ++__i)
	  __sum += this->weight[__i]
		 * __func(__shift + __slope * this->point[__i]);

	return __sign * __fact * __sum;
      }

} // namespace __gnu_cxx

#endif // GAUSS_QUADRATURE_TCC
