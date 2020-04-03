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
	operator()(_FuncTp __func, _Tp __lower, _Tp __upper) const;

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
	operator()(_FuncTp __func, _Tp __lower, _Tp __upper) const;

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
	operator()(_FuncTp __func, _Tp __lower, _Tp __upper) const;

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
	operator()(_FuncTp __func, _Tp __lower, _Tp __upper) const;

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
	operator()(_FuncTp __func, _Tp __lower, _Tp __upper) const;

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
	operator()(_FuncTp __func, _Tp __lower, _Tp __upper) const;

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
	operator()(_FuncTp __func, _Tp __lower, _Tp __upper) const;

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
	operator()(_FuncTp __func, _Tp __lower, _Tp __upper) const;

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
	operator()(_FuncTp __func, _Tp __lower, _Tp __upper) const;

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
	operator()(_FuncTp __func, _Tp __lower, _Tp __upper) const;

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
	operator()(_FuncTp __func, _Tp __lower, _Tp __upper) const;

    private:
      std::vector<_Tp> point;
      std::vector<_Tp> weight;
    };

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_legendre(int __n,
				   _FuncTp __func, _Tp __lower, _Tp __upper)
    {
      using __integ_t = fixed_integral_t<_Tp,
				     std::invoke_result_t<_FuncTp, _Tp>>;
      using __area_t = typename __integ_t::_AreaTp;

      if (std::isnan(__lower) || std::isnan(__upper))
	{
	  const auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
	  return {__area_t{} * _S_NaN};
	}
      else if (__lower == __upper)
	return {__area_t{}};
      else
	{
          fixed_gauss_legendre_integral<_Tp> __integ(__n);
          return { __integ(__func, __lower, __upper) };
	}
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_chebyshev_t(int __n,
				      _FuncTp __func, _Tp __lower, _Tp __upper)
    {
      using __integ_t = fixed_integral_t<_Tp,
				     std::invoke_result_t<_FuncTp, _Tp>>;
      using __area_t = typename __integ_t::_AreaTp;

      if (std::isnan(__lower) || std::isnan(__upper))
	{
	  const auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
	  return {__area_t{} * _S_NaN};
	}
      else if (__lower == __upper)
	return {__area_t{}};
      else
	{
          fixed_gauss_chebyshev_t_integral<_Tp> __integ(__n);
          return { __integ(__func, __lower, __upper) };
	}
    }


  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_chebyshev_u(int __n,
				      _FuncTp __func, _Tp __lower, _Tp __upper)
    {
      using __integ_t = fixed_integral_t<_Tp,
				     std::invoke_result_t<_FuncTp, _Tp>>;
      using __area_t = typename __integ_t::_AreaTp;

      if (std::isnan(__lower) || std::isnan(__upper))
	{
	  const auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
	  return {__area_t{} * _S_NaN};
	}
      else if (__lower == __upper)
	return {__area_t{}};
      else
	{
          fixed_gauss_chebyshev_u_integral<_Tp> __integ(__n);
          return { __integ(__func, __lower, __upper) };
	}
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_chebyshev_v(int __n,
				      _FuncTp __func, _Tp __lower, _Tp __upper)
    {
      using __integ_t = fixed_integral_t<_Tp,
				     std::invoke_result_t<_FuncTp, _Tp>>;
      using __area_t = typename __integ_t::_AreaTp;

      if (std::isnan(__lower) || std::isnan(__upper))
	{
	  const auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
	  return {__area_t{} * _S_NaN};
	}
      else if (__lower == __upper)
	return {__area_t{}};
      else
	{
          fixed_gauss_chebyshev_v_integral<_Tp> __integ(__n);
          return { __integ(__func, __lower, __upper) };
	}
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_chebyshev_w(int __n,
				      _FuncTp __func, _Tp __lower, _Tp __upper)
    {
      using __integ_t = fixed_integral_t<_Tp,
				     std::invoke_result_t<_FuncTp, _Tp>>;
      using __area_t = typename __integ_t::_AreaTp;

      if (std::isnan(__lower) || std::isnan(__upper))
	{
	  const auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
	  return {__area_t{} * _S_NaN};
	}
      else if (__lower == __upper)
	return {__area_t{}};
      else
	{
          fixed_gauss_chebyshev_w_integral<_Tp> __integ(__n);
          return { __integ(__func, __lower, __upper) };
	}
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_gegenbauer(int __n, _Tp __lambda,
				     _FuncTp __func, _Tp __lower, _Tp __upper)
    {
      using __integ_t = fixed_integral_t<_Tp,
				     std::invoke_result_t<_FuncTp, _Tp>>;
      using __area_t = typename __integ_t::_AreaTp;

      if (std::isnan(__lower) || std::isnan(__upper) || std::isnan(__lambda))
	{
	  const auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
	  return {__area_t{} * _S_NaN};
	}
      else if (__lower == __upper)
	return {__area_t{}};
      else
	{
          fixed_gauss_gegenbauer_integral<_Tp> __integ(__n, __lambda);
          return { __integ(__func, __lower, __upper) };
	}
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_jacobi(int __n, _Tp __alf, _Tp __bet,
				 _FuncTp __func, _Tp __lower, _Tp __upper)
    {
      using __integ_t = fixed_integral_t<_Tp,
				     std::invoke_result_t<_FuncTp, _Tp>>;
      using __area_t = typename __integ_t::_AreaTp;

      if (std::isnan(__lower) || std::isnan(__upper)
          || std::isnan(__alf) || std::isnan(__bet))
	{
	  const auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
	  return {__area_t{} * _S_NaN};
	}
      else if (__lower == __upper)
	return {__area_t{}};
      else
	{
          fixed_gauss_jacobi_integral<_Tp> __integ(__n, __alf, __bet);
          return { __integ(__func, __lower, __upper) };
	}
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_laguerre(int __n, _Tp __alf,
				   _FuncTp __func, _Tp __lower, _Tp __upper)
    {
      using __integ_t = fixed_integral_t<_Tp,
				     std::invoke_result_t<_FuncTp, _Tp>>;
      using __area_t = typename __integ_t::_AreaTp;

      if (std::isnan(__lower) || std::isnan(__upper) || std::isnan(__alf))
	{
	  const auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
	  return {__area_t{} * _S_NaN};
	}
      else if (__lower == __upper)
	return {__area_t{}};
      else
	{
          fixed_gauss_laguerre_integral<_Tp> __integ(__n, __alf);
          return { __integ(__func, __lower, __upper) };
	}
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_hermite(int __n, _Tp __alf,
				  _FuncTp __func, _Tp __lower, _Tp __upper)
    {
      using __integ_t = fixed_integral_t<_Tp,
				     std::invoke_result_t<_FuncTp, _Tp>>;
      using __area_t = typename __integ_t::_AreaTp;

      if (std::isnan(__lower) || std::isnan(__upper) || std::isnan(__alf))
	{
	  const auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
	  return {__area_t{} * _S_NaN};
	}
      else if (__lower == __upper)
	return {__area_t{}};
      else
	{
          fixed_gauss_hermite_integral<_Tp> __integ(__n, __alf);
          return { __integ(__func, __lower, __upper) };
	}
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_exponential(int __n, _Tp __alf,
				      _FuncTp __func, _Tp __lower, _Tp __upper)
    {
      using __integ_t = fixed_integral_t<_Tp,
				     std::invoke_result_t<_FuncTp, _Tp>>;
      using __area_t = typename __integ_t::_AreaTp;

      if (std::isnan(__lower) || std::isnan(__upper) || std::isnan(__alf))
	{
	  const auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
	  return {__area_t{} * _S_NaN};
	}
      else if (__lower == __upper)
	return {__area_t{}};
      else
	{
          fixed_gauss_exponential_integral<_Tp> __integ(__n, __alf);
          return { __integ(__func, __lower, __upper) };
	}
    }

  template<typename _Tp, typename _FuncTp>
    fixed_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    integrate_fixed_gauss_rational(int __n, _Tp __alf, _Tp __bet,
				   _FuncTp __func, _Tp __lower, _Tp __upper)
    {
      using __integ_t = fixed_integral_t<_Tp,
				     std::invoke_result_t<_FuncTp, _Tp>>;
      using __area_t = typename __integ_t::_AreaTp;

      if (std::isnan(__lower) || std::isnan(__upper)
          || std::isnan(__alf) || std::isnan(__bet))
	{
	  const auto _S_NaN = std::numeric_limits<_Tp>::quiet_NaN();
	  return {__area_t{} * _S_NaN};
	}
      else if (__lower == __upper)
	return {__area_t{}};
      else
	{
          fixed_gauss_rational_integral<_Tp> __integ(__n, __alf, __bet);
          return { __integ(__func, __lower, __upper) };
	}
    }

} // namespace __gnu_cxx

#endif // GAUSS_QUADRATURE_H
