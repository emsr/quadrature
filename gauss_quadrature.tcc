// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2018 Free Software Foundation, Inc.
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

#include <ext/cmath>


namespace __gnu_cxx
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
    golub_welsch(_Tp __moment0, int __n, _InIter __diag, _InIter __subd,
		_OutIter __pt, _OutIter __wt)
    {
      // Bail if the zero-th moment is not positive.
      if (__moment0 <= _Tp{0})
	std::__throw_domain_error("golub_welsch: moment0 <= 0");

      // Set up vectors for matrix diagonalization.
      for (int i = 0; i < __n; ++i)
	__pt[i] = __diag[i];

      __wt[0] = std::sqrt(__moment0);
      for (int i = 1; i < __n; ++i)
	__wt[i] = _Tp{0};

      // Diagonalize the Jacobi matrix.
      _S_tridiag_symm(std::size_t(__n), __diag, __subd, __wt);

      for (int i = 0; i < __n; ++i)
	__wt[i] *= __wt[i];

      return;
    }

  template<typename _Tp, typename _FuncTp, typename _OutIter>
    _Tp
    gauss_eval(_FuncTp __func, _Tp __a, _Tp __b,
	       int __n, const _OutIter __pt, const _OutIter __wt)
    {
      auto __sign = _Tp{1};
      if (__b < __a)
	{
	  __sign = _Tp{-1};
	  std::swap(__a, __b);
	}

      auto __sum = _Tp{0};
      for (int __i = 0; __i < __n; ++__i)
	__sum += __wt[__i] * __func(__pt[__i]); // FIXME: Map to a: b!!!

      return __sign * __sum;
    }


  /**
   * Interval: @f$ (a, b) @f$
   * Weight function: @f$ 1 @f$
   * Constraints: @f$ b > a @f$
   */
  template<typename _Tp>
    gauss_legendre_rule<_Tp>::gauss_legendre_rule(int __n)
    : order(__n),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = _Tp{2};

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order);

      for (int __i = 1; __i <= this->order; ++__i)
	__subd[__i - 1] = _Tp(__i) / std::sqrt(_Tp(4 * __i * __i - 1));

      golub_welsch(__mu_0, this->order, __diag, __subd,
		   this->point, this->weight);
for (int __i = 0; __i < this->order; ++__i) std::cout << ' ' << this->point[__i] << ' ' << this->point[__i] << '\n';
    }

  template<typename _Tp>
    template<typename _FuncTp>
      _Tp
      gauss_legendre_rule<_Tp>::operator()(_FuncTp __func,
					   _Tp __a, _Tp __b) const
      {
	return gauss_eval(__func, __a, __b,
			  this->order, this->point, this->weight);
      }

  /**
   * Interval: @f$ (a, b) @f$
   * Weight function: @f$ ((b-x)(x-a))^{-1/2} @f$
   * Constraints: @f$ b > a @f$
   * Jacobi parameter: @f$ \alpha = 1/2 @f$
   * Jacobi parameter: @f$ \beta = 1/2 @f$
   */
  template<typename _Tp>
    gauss_chebyshev_t_rule<_Tp>::gauss_chebyshev_t_rule(int __n)
    : order(__n),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = __gnu_cxx::__const_pi<_Tp>();

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order, _Tp{0.5L});

      __subd[0] = std::sqrt(_Tp{0.5L});

      golub_welsch(__mu_0, this->order, __diag, __subd,
		   this->point, this->weight);
    }

  template<typename _Tp>
    template<typename _FuncTp>
      _Tp
      gauss_chebyshev_t_rule<_Tp>::operator()(_FuncTp __func,
					      _Tp __a, _Tp __b) const
      {
	return gauss_eval(__func, __a, __b,
			  this->order, this->point, this->weight);
      }

  /**
   * Interval: @f$ (a, b), b > a @f$
   * Weight function: @f$ ((b-x)(x-a))^{1/2} @f$
   * Constraints: @f$ \alpha = -1/2 @f$
   * Constraints: @f$ \beta = -1/2 @f$
   */
  template<typename _Tp>
    gauss_chebyshev_u_rule<_Tp>::gauss_chebyshev_u_rule(int __n)
    : order(__n),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = __gnu_cxx::__const_pi<_Tp>() / _Tp{2};

      this->point.resize(this->order);
      this->weight.resize(this->order);

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order, _Tp{0.5L});

      golub_welsch(__mu_0, this->order, __diag, __subd,
		   this->point, this->weight);
    }

  template<typename _Tp>
    template<typename _FuncTp>
      _Tp
      gauss_chebyshev_u_rule<_Tp>::operator()(_FuncTp __func,
					      _Tp __a, _Tp __b) const
      {
	return gauss_eval(__func, __a, __b,
			  this->order, this->point, this->weight);
      }

  /**
   * Interval: @f$ (a, b), b > a @f$
   * Weight function: @f$ ((x-a)/(b-x))^{1/2} @f$
   * Constraints: @f$ \alpha = +1/2 @f$
   * Constraints: @f$ \beta = -1/2 @f$
   */
  template<typename _Tp>
    gauss_chebyshev_v_rule<_Tp>::gauss_chebyshev_v_rule(int __n)
    : order(__n),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = __gnu_cxx::__const_pi<_Tp>();

      this->point.resize(this->order);
      this->weight.resize(this->order);

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order, _Tp{0.5L});

      __diag[0] = _Tp{-0.5L};
      __subd[0] = std::sqrt(_Tp{2});

      golub_welsch(__mu_0, this->order, __diag, __subd,
		   this->point, this->weight);
    }

  template<typename _Tp>
    template<typename _FuncTp>
      _Tp
      gauss_chebyshev_v_rule<_Tp>::operator()(_FuncTp __func,
					      _Tp __a, _Tp __b) const
      {
	return gauss_eval(__func, __a, __b,
			  this->order, this->point, this->weight);
      }

  /**
   * Interval: @f$ (a, b), b > a @f$
   * Weight function: @f$ ((b-x)/(x-a))^{1/2} @f$
   * Constraints: @f$ \alpha = -1/2 @f$
   * Constraints: @f$ \beta = +1/2 @f$
   */
  template<typename _Tp>
    gauss_chebyshev_w_rule<_Tp>::gauss_chebyshev_w_rule(int __n)
    : order(__n),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = __gnu_cxx::__const_pi<_Tp>();

      this->point.resize(this->order);
      this->weight.resize(this->order);

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order, _Tp{0.5L});

      __diag[0] = _Tp{0.5L};
      __subd[0] = std::sqrt(_Tp{2});

      golub_welsch(__mu_0, this->order, __diag, __subd,
		   this->point, this->weight);
    }

  template<typename _Tp>
    template<typename _FuncTp>
      _Tp
      gauss_chebyshev_w_rule<_Tp>::operator()(_FuncTp __func,
					      _Tp __a, _Tp __b) const
      {
	return gauss_eval(__func, __a, __b,
			  this->order, this->point, this->weight);
      }

  /**
   * Interval: @f$ (a, b) @f$
   * Weight function: @f$ ((b-x)(x-a))^\alpha @f$
   * Constraints: @f$ \alpha > -1 @f$
   * Constraints: @f$ b > a @f$
   */
  template<typename _Tp>
    gauss_gegenbauer_rule<_Tp>::gauss_gegenbauer_rule(int __n, _Tp __alf)
    : order(__n),
      alpha(__alf),
      point(__n),
      weight(__n)
    {
      const auto __ab = _Tp{2} * this->alpha;
      const auto __gam = std::tgamma(this->alpha + _Tp{1});
      const auto __mu_0 = std::pow(_Tp{2}, __ab + _Tp{1})
			* __gam * __gam / std::tgamma(__ab + _Tp{2});

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order);

      __subd[0] = std::sqrt(_Tp{1} / (_Tp{2} * this->alpha + _Tp{3}));

      for (int __i = 2; __i <= this->order; ++__i)
	__subd[__i-1] = std::sqrt(_Tp(__i) * (__ab + _Tp(__i))
		      / (_Tp{4}
			* std::pow(this->alpha + _Tp(__i), _Tp{2}) - _Tp{1}));

      golub_welsch(__mu_0, this->order, __diag, __subd,
		   this->point, this->weight);
    }

  template<typename _Tp>
    template<typename _FuncTp>
      _Tp
      gauss_gegenbauer_rule<_Tp>::operator()(_FuncTp __func,
					     _Tp __a, _Tp __b) const
      {
	return gauss_eval(__func, __a, __b,
			  this->order, this->point, this->weight);
      }

  /**
   * Interval: @f$ (a, b) @f$
   * Weight function: @f$ (b-x)^\alpha (x-a)^\beta @f$
   * Constraints: @f$ \alpha, \beta > -1 @f$
   * Constraints: @f$ b > a @f$
   */
  template<typename _Tp>
    gauss_jacobi_rule<_Tp>::gauss_jacobi_rule(int __n, _Tp __alf, _Tp __bet)
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
	     * std::tgamma(this->beta + _Tp{1}) / std::tgamma(__abp2i);

      std::vector<_Tp> __diag(this->order);
      std::vector<_Tp> __subd(this->order);

      __diag[0] = (this->beta - this->alpha) / __abp2i;
      __subd[0] = _Tp{2} * std::sqrt((this->alpha + _Tp{1})
				   * (this->beta + _Tp{1})
				   / (__abp2i + _Tp{1}) / __abp2i);
      const auto a2mb2 = (this->beta - this->alpha)
			* (this->beta + this->alpha);
      for (int __i = 1; __i <= this->order; ++__i)
	{
	  __abp2i += _Tp{2};
	  const auto abp2ip2 = __abp2i + _Tp{2};
	  __diag[__i] = a2mb2 / __abp2i / abp2ip2;
	  const auto ip1 = _Tp(__i + 1);
	  __subd[__i] = std::sqrt(_Tp(4 * __i) * (this->alpha + ip1)
				  * (this->beta + ip1) * (__ab + ip1)
				  / (abp2ip2 * abp2ip2 - _Tp{1})) / abp2ip2;
	}

      golub_welsch(__mu_0, this->order, __diag, __subd,
		   this->point, this->weight);
    }

  template<typename _Tp>
    template<typename _FuncTp>
      _Tp
      gauss_jacobi_rule<_Tp>::operator()(_FuncTp __func,
					 _Tp __a, _Tp __b) const
      {
	return gauss_eval(__func, __a, __b,
			  this->order, this->point, this->weight);
      }

  /**
   * Interval: @f$ (a, \infty) @f$
   * Weight function: @f$ (x-a)^\alpha \exp(-b(x-a)) @f$
   * Constraints: @f$ \alpha > -1 @f$
   * Constraints: @f$ b > 0 @f$
   */
  template<typename _Tp>
    gauss_laguerre_rule<_Tp>::gauss_laguerre_rule(int __n, _Tp __alf)
    : order(__n),
      alpha(__alf),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = std::tgamma(this->alpha + _Tp{1});

      std::vector<_Tp> __diag(this->order);
      std::vector<_Tp> __subd(this->order);

      for (int __i = 1; __i <= this->order; ++__i)
	{
	  __diag[__i - 1] = _Tp(2 * __i - 1) + this->alpha;
	  __subd[__i - 1] = std::sqrt(_Tp(__i) * (_Tp(__i) + this->alpha));
	}

      golub_welsch(__mu_0, this->order, __diag, __subd,
		   this->point, this->weight);
    }

  template<typename _Tp>
    template<typename _FuncTp>
      _Tp
      gauss_laguerre_rule<_Tp>::operator()(_FuncTp __func,
					   _Tp __a, _Tp __b) const
      {
	return gauss_eval(__func, __a, __b,
			  this->order, this->point, this->weight);
      }

  /**
   * Interval: @f$ (-\infty, \infty) @f$
   * Weight function: @f$ |x-a|^\alpha \exp(-b(x-a)^2) @f$
   * Constraints: @f$ \alpha > -1 @f$
   * Constraints: @f$ b > 0 @f$
   */
  template<typename _Tp>
    gauss_hermite_rule<_Tp>::gauss_hermite_rule(int __n, _Tp __alf)
    : order(__n),
      alpha(__alf),
      point(__n),
      weight(__n)
    {
      const auto __mu_0 = std::tgamma((this->alpha + 1) / _Tp{2});

      std::vector<_Tp> __diag(this->order, _Tp{0});
      std::vector<_Tp> __subd(this->order);

      for (int __i = 1; __i <= this->order; ++__i)
	__subd[__i - 1] = std::sqrt((_Tp(__i) + this->alpha * _Tp(__i % 2))
				   / _Tp{2});

      golub_welsch(__mu_0, this->order, __diag, __subd,
		   this->point, this->weight);
    }

  template<typename _Tp>
    template<typename _FuncTp>
      _Tp
      gauss_hermite_rule<_Tp>::operator()(_FuncTp __func,
					  _Tp __a, _Tp __b) const
      {
	return gauss_eval(__func, __a, __b,
			  this->order, this->point, this->weight);
      }

  /**
   * Interval: @f$ (a, b) @f$
   * Weight function: @f$ |x-(a+b)/2|^\alpha @f$
   * Constraints: @f$ \alpha > -1 @f$
   * Constraints: @f$ b > a @f$
   */
  template<typename _Tp>
    gauss_exponential_rule<_Tp>::gauss_exponential_rule(int __n, _Tp __alf)
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

      golub_welsch(__mu_0, this->order, __diag, __subd,
		   this->point, this->weight);
    }

  template<typename _Tp>
    template<typename _FuncTp>
      _Tp
      gauss_exponential_rule<_Tp>::operator()(_FuncTp __func,
					      _Tp __a, _Tp __b) const
      {
	return gauss_eval(__func, __a, __b,
			  this->order, this->point, this->weight);
      }

  /**
   * Interval: @f$ (a, \infty) @f$
   * Weight function: @f$ (x-a)^\alpha*(b+x)^\beta @f$
   * Constraints: @f$ \alpha > -1, \alpha + \beta + 2n < 0, a + b > 0 @f$
   */
  template<typename _Tp>
    gauss_rational_rule<_Tp>::gauss_rational_rule(int __n, _Tp __alf, _Tp __bet)
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
	  __subd[__i-1] = _Tp(__i) * (this->alpha + _Tp(__i))
			/ (__abp2i - _Tp{1}) * (this->beta + _Tp(__i))
			/ (__abp2i * __abp2i) * (__ab + _Tp(__i))
			/ (__abp2i + _Tp{1});
	}
      __subd[this->order - 1] = _Tp{0};
      for (int __i = 0; __i < this->order; ++__i)
	__subd[__i] = std::sqrt(__subd[__i]);

      golub_welsch(__mu_0, this->order, __diag, __subd,
		   this->point, this->weight);
    }

  template<typename _Tp>
    template<typename _FuncTp>
      _Tp
      gauss_rational_rule<_Tp>::operator()(_FuncTp __func,
					   _Tp __a, _Tp __b) const
      {
	return gauss_eval(__func, __a, __b,
			  this->order, this->point, this->weight);
      }

} // namespace __gnu_cxx

#endif // GAUSS_QUADRATURE_TCC
