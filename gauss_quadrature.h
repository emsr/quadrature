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
    struct gauss_legendre_rule
    {
      int order;

      explicit gauss_legendre_rule(int __n);

      template<typename _FuncTp>
	_Tp
	operator()(_FuncTp, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct gauss_chebyshev_t_rule
    {
      int order;

      explicit gauss_chebyshev_t_rule(int __n);

      template<typename _FuncTp>
	_Tp
	operator()(_FuncTp, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct gauss_chebyshev_u_rule
    {
      int order;

      explicit gauss_chebyshev_u_rule(int __n);

      template<typename _FuncTp>
	_Tp
	operator()(_FuncTp, _Tp __a, _Tp __b) const;

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
    struct gauss_chebyshev_v_rule
    {
      int order;

      explicit gauss_chebyshev_v_rule(int __n);

      template<typename _FuncTp>
	_Tp
	operator()(_FuncTp, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct gauss_chebyshev_w_rule
    {
      int order;

      explicit gauss_chebyshev_w_rule(int __n);

      template<typename _FuncTp>
	_Tp
	operator()(_FuncTp, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct gauss_gegenbauer_rule
    {
      int order;
      _Tp lambda;

      explicit gauss_gegenbauer_rule(int __n, _Tp __lam);

      template<typename _FuncTp>
	_Tp
	operator()(_FuncTp, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct gauss_jacobi_rule
    {
      int order;
      _Tp alpha;
      _Tp beta;

      explicit gauss_jacobi_rule(int __n, _Tp __alf, _Tp __bet);

      template<typename _FuncTp>
	_Tp
	operator()(_FuncTp, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct gauss_laguerre_rule
    {
      int order;
      _Tp alpha;

      explicit gauss_laguerre_rule(int n, _Tp __alf);

      template<typename _FuncTp>
	_Tp
	operator()(_FuncTp, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct gauss_hermite_rule
    {
      int order;
      _Tp alpha;

      explicit gauss_hermite_rule(int __n, _Tp __alf);

      template<typename _FuncTp>
	_Tp
	operator()(_FuncTp, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct gauss_exponential_rule
    {
      int order;
      _Tp alpha;

      explicit gauss_exponential_rule(int __n, _Tp __alf);

      template<typename _FuncTp>
	_Tp
	operator()(_FuncTp, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  /**
   * 
   */
  template<typename _Tp>
    struct gauss_rational_rule
    {
      int order;
      _Tp alpha;
      _Tp beta;

      explicit gauss_rational_rule(int __n, _Tp __alf, _Tp __bet);

      template<typename _FuncTp>
	_Tp
	operator()(_FuncTp, _Tp __a, _Tp __b) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

} // namespace __gnu_cxx

#include "gauss_quadrature.tcc"

#endif // GAUSS_QUADRATURE_H
