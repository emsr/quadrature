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

#ifndef GAUSS_QUADRATURE_H
#define GAUSS_QUADRATURE_H 1

#include <vector>

namespace __gnu_cxx
{

  /**
   * Build the absiscae and weights of a Gauss quadrature rule from
   * a symmetric tridiagonal Jacobi matrix.
   */
  template<typename _Tp>
    struct fixed_gauss_legendre_integral
    {
      int order;

      explicit fixed_gauss_legendre_integral(int __n);

      template<typename _FuncTp>
	decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
	operator()(_FuncTp __func, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct fixed_gauss_chebyshev_t_integral
    {
      int order;

      explicit fixed_gauss_chebyshev_t_integral(int __n);

      template<typename _FuncTp>
	decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
	operator()(_FuncTp __func, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct fixed_gauss_chebyshev_u_integral
    {
      int order;

      explicit fixed_gauss_chebyshev_u_integral(int __n);

      template<typename _FuncTp>
	decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
	operator()(_FuncTp __func, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * Majid Tavassoli Kajani, Adem KJlJçman, and Mohammad Maleki
   * The Rational Third-Kind Chebyshev Pseudospectral Method
   * for the Solution of the Thomas-Fermi Equation over Infinite Interval
   *
   */
  template<typename _Tp>
    struct fixed_gauss_chebyshev_v_integral
    {
      int order;

      explicit fixed_gauss_chebyshev_v_integral(int __n);

      template<typename _FuncTp>
	decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
	operator()(_FuncTp __func, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct fixed_gauss_chebyshev_w_integral
    {
      int order;

      explicit fixed_gauss_chebyshev_w_integral(int __n);

      template<typename _FuncTp>
	decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
	operator()(_FuncTp __func, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct fixed_gauss_gegenbauer_integral
    {
      int order;
      _Tp lambda;

      explicit fixed_gauss_gegenbauer_integral(int __n, _Tp __lam);

      template<typename _FuncTp>
	decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
	operator()(_FuncTp __func, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct fixed_gauss_jacobi_integral
    {
      int order;
      _Tp alpha;
      _Tp beta;

      explicit fixed_gauss_jacobi_integral(int __n, _Tp __alf, _Tp __bet);

      template<typename _FuncTp>
	decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
	operator()(_FuncTp __func, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct fixed_gauss_laguerre_integral
    {
      int order;
      _Tp alpha;

      explicit fixed_gauss_laguerre_integral(int n, _Tp __alf);

      template<typename _FuncTp>
	decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
	operator()(_FuncTp __func, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct fixed_gauss_hermite_integral
    {
      int order;
      _Tp alpha;

      explicit fixed_gauss_hermite_integral(int __n, _Tp __alf);

      template<typename _FuncTp>
	decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
	operator()(_FuncTp __func, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct fixed_gauss_exponential_integral
    {
      int order;
      _Tp alpha;

      explicit fixed_gauss_exponential_integral(int __n, _Tp __alf);

      template<typename _FuncTp>
	decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
	operator()(_FuncTp __func, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct fixed_gauss_rational_integral
    {
      int order;
      _Tp alpha;
      _Tp beta;

      explicit fixed_gauss_rational_integral(int __n, _Tp __alf, _Tp __bet);

      template<typename _FuncTp>
	decltype(std::invoke_result_t<_FuncTp, _Tp>{} * _Tp{})
	operator()(_FuncTp __func, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_legendre(int __n,
				   _FuncTp __func, _Tp __a, _Tp __b)
    {
      fixed_gauss_legendre_integral<_Tp> __integ(__n);
      return { __integ(__func, __a, __b) };
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_chebyshev_t(int __n,
				      _FuncTp __func, _Tp __a, _Tp __b)
    {
      fixed_gauss_chebyshev_t_integral<_Tp> __integ(__n);
      return { __integ(__func, __a, __b) };
    }


  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_chebyshev_u(int __n,
				      _FuncTp __func, _Tp __a, _Tp __b)
    {
      fixed_gauss_chebyshev_u_integral<_Tp> __integ(__n);
      return { __integ(__func, __a, __b) };
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_chebyshev_v(int __n,
				      _FuncTp __func, _Tp __a, _Tp __b)
    {
      fixed_gauss_chebyshev_v_integral<_Tp> __integ(__n);
      return { __integ(__func, __a, __b) };
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_chebyshev_w(int __n,
				      _FuncTp __func, _Tp __a, _Tp __b)
    {
      fixed_gauss_chebyshev_w_integral<_Tp> __integ(__n);
      return { __integ(__func, __a, __b) };
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_gegenbauer(int __n, _Tp __lambda,
				     _FuncTp __func, _Tp __a, _Tp __b)
    {
      fixed_gauss_gegenbauer_integral<_Tp> __integ(__n, __lambda);
      return { __integ(__func, __a, __b) };
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_jacobi(int __n, _Tp __alf, _Tp __bet,
				 _FuncTp __func, _Tp __a, _Tp __b)
    {
      fixed_gauss_jacobi_integral<_Tp> __integ(__n, __alf, __bet);
      return { __integ(__func, __a, __b) };
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_laguerre(int __n, _Tp __alf,
				   _FuncTp __func, _Tp __a, _Tp __b)
    {
      fixed_gauss_laguerre_integral<_Tp> __integ(__n, __alf);
      return { __integ(__func, __a, __b) };
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_hermite(int __n, _Tp __alf,
				  _FuncTp __func, _Tp __a, _Tp __b)
    {
      fixed_gauss_hermite_integral<_Tp> __integ(__n, __alf);
      return { __integ(__func, __a, __b) };
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_exponential(int __n, _Tp __alf,
				      _FuncTp __func, _Tp __a, _Tp __b)
    {
      fixed_gauss_exponential_integral<_Tp> __integ(__n, __alf);
      return { __integ(__func, __a, __b) };
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_rational(int __n, _Tp __alf, _Tp __bet,
				   _FuncTp __func, _Tp __a, _Tp __b)
    {
      fixed_gauss_rational_integral<_Tp> __integ(__n, __alf, __bet);
      return { __integ(__func, __a, __b) };
    }

} // namespace __gnu_cxx

#endif // GAUSS_QUADRATURE_H
