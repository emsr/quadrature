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

#ifndef GAUSS_QUADRATURE_H
#define GAUSS_QUADRATURE_H 1

#include <type_traits>
#include <vector>

namespace emsr
{

  /**
   * Build the absiscae and weights of a Gauss quadrature rule from
   * a symmetric tridiagonal Jacobi matrix.
   */
  template<typename Tp>
    struct fixed_gauss_legendre_integral
    {
      int order;

      explicit fixed_gauss_legendre_integral(int n);

      template<typename FuncTp>
	decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
	operator()(FuncTp func, Tp lower, Tp upper) const;

    private:
      std::vector<Tp> point;
      std::vector<Tp> weight;
    };

  /**
   * 
   */
  template<typename Tp>
    struct fixed_gauss_chebyshev_t_integral
    {
      int order;

      explicit fixed_gauss_chebyshev_t_integral(int n);

      template<typename FuncTp>
	decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
	operator()(FuncTp func, Tp lower, Tp upper) const;

    private:
      std::vector<Tp> point;
      std::vector<Tp> weight;
    };

  /**
   * 
   */
  template<typename Tp>
    struct fixed_gauss_chebyshev_u_integral
    {
      int order;

      explicit fixed_gauss_chebyshev_u_integral(int n);

      template<typename FuncTp>
	decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
	operator()(FuncTp func, Tp lower, Tp upper) const;

    private:
      std::vector<Tp> point;
      std::vector<Tp> weight;
    };

  /**
   * Majid Tavassoli Kajani, Adem KJlJçman, and Mohammad Maleki
   * The Rational Third-Kind Chebyshev Pseudospectral Method
   * for the Solution of the Thomas-Fermi Equation over Infinite Interval
   *
   */
  template<typename Tp>
    struct fixed_gauss_chebyshev_v_integral
    {
      int order;

      explicit fixed_gauss_chebyshev_v_integral(int n);

      template<typename FuncTp>
	decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
	operator()(FuncTp func, Tp lower, Tp upper) const;

    private:
      std::vector<Tp> point;
      std::vector<Tp> weight;
    };

  /**
   * 
   */
  template<typename Tp>
    struct fixed_gauss_chebyshev_w_integral
    {
      int order;

      explicit fixed_gauss_chebyshev_w_integral(int n);

      template<typename FuncTp>
	decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
	operator()(FuncTp func, Tp lower, Tp upper) const;

    private:
      std::vector<Tp> point;
      std::vector<Tp> weight;
    };

  /**
   * 
   */
  template<typename Tp>
    struct fixed_gauss_gegenbauer_integral
    {
      int order;
      Tp lambda;

      explicit fixed_gauss_gegenbauer_integral(int n, Tp lam);

      template<typename FuncTp>
	decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
	operator()(FuncTp func, Tp lower, Tp upper) const;

    private:
      std::vector<Tp> point;
      std::vector<Tp> weight;
    };

  /**
   * 
   */
  template<typename Tp>
    struct fixed_gauss_jacobi_integral
    {
      int order;
      Tp alpha;
      Tp beta;

      explicit fixed_gauss_jacobi_integral(int n, Tp alf, Tp bet);

      template<typename FuncTp>
	decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
	operator()(FuncTp func, Tp lower, Tp upper) const;

    private:
      std::vector<Tp> point;
      std::vector<Tp> weight;
    };

  /**
   * 
   */
  template<typename Tp>
    struct fixed_gauss_laguerre_integral
    {
      int order;
      Tp alpha;

      explicit fixed_gauss_laguerre_integral(int n, Tp alf);

      template<typename FuncTp>
	decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
	operator()(FuncTp func, Tp lower, Tp upper) const;

    private:
      std::vector<Tp> point;
      std::vector<Tp> weight;
    };

  /**
   * 
   */
  template<typename Tp>
    struct fixed_gauss_hermite_integral
    {
      int order;
      Tp alpha;

      explicit fixed_gauss_hermite_integral(int n, Tp alf);

      template<typename FuncTp>
	decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
	operator()(FuncTp func, Tp lower, Tp upper) const;

    private:
      std::vector<Tp> point;
      std::vector<Tp> weight;
    };

  /**
   * 
   */
  template<typename Tp>
    struct fixed_gauss_exponential_integral
    {
      int order;
      Tp alpha;

      explicit fixed_gauss_exponential_integral(int n, Tp alf);

      template<typename FuncTp>
	decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
	operator()(FuncTp func, Tp lower, Tp upper) const;

    private:
      std::vector<Tp> point;
      std::vector<Tp> weight;
    };

  /**
   * 
   */
  template<typename Tp>
    struct fixed_gauss_rational_integral
    {
      int order;
      Tp alpha;
      Tp beta;

      explicit fixed_gauss_rational_integral(int n, Tp alf, Tp bet);

      template<typename FuncTp>
	decltype(std::invoke_result_t<FuncTp, Tp>{} * Tp{})
	operator()(FuncTp func, Tp lower, Tp upper) const;

    private:
      std::vector<Tp> point;
      std::vector<Tp> weight;
    };

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_legendre(int n,
				   FuncTp func, Tp lower, Tp upper)
    {
      using integ_t = fixed_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;

      if (std::isnan(lower) || std::isnan(upper))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}};
      else
	{
          fixed_gauss_legendre_integral<Tp> integ(n);
          return { integ(func, lower, upper) };
	}
    }

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_chebyshev_t(int n,
				      FuncTp func, Tp lower, Tp upper)
    {
      using integ_t = fixed_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;

      if (std::isnan(lower) || std::isnan(upper))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}};
      else
	{
          fixed_gauss_chebyshev_t_integral<Tp> integ(n);
          return { integ(func, lower, upper) };
	}
    }


  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_chebyshev_u(int n,
				      FuncTp func, Tp lower, Tp upper)
    {
      using integ_t = fixed_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;

      if (std::isnan(lower) || std::isnan(upper))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}};
      else
	{
          fixed_gauss_chebyshev_u_integral<Tp> integ(n);
          return { integ(func, lower, upper) };
	}
    }

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_chebyshev_v(int n,
				      FuncTp func, Tp lower, Tp upper)
    {
      using integ_t = fixed_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;

      if (std::isnan(lower) || std::isnan(upper))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}};
      else
	{
          fixed_gauss_chebyshev_v_integral<Tp> integ(n);
          return { integ(func, lower, upper) };
	}
    }

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_chebyshev_w(int n,
				      FuncTp func, Tp lower, Tp upper)
    {
      using integ_t = fixed_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;

      if (std::isnan(lower) || std::isnan(upper))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}};
      else
	{
          fixed_gauss_chebyshev_w_integral<Tp> integ(n);
          return { integ(func, lower, upper) };
	}
    }

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_gegenbauer(int n, Tp lambda,
				     FuncTp func, Tp lower, Tp upper)
    {
      using integ_t = fixed_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;

      if (std::isnan(lower) || std::isnan(upper) || std::isnan(lambda))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}};
      else
	{
          fixed_gauss_gegenbauer_integral<Tp> integ(n, lambda);
          return { integ(func, lower, upper) };
	}
    }

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_jacobi(int n, Tp alf, Tp bet,
				 FuncTp func, Tp lower, Tp upper)
    {
      using integ_t = fixed_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;

      if (std::isnan(lower) || std::isnan(upper)
          || std::isnan(alf) || std::isnan(bet))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}};
      else
	{
          fixed_gauss_jacobi_integral<Tp> integ(n, alf, bet);
          return { integ(func, lower, upper) };
	}
    }

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_laguerre(int n, Tp alf,
				   FuncTp func, Tp lower, Tp upper)
    {
      using integ_t = fixed_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;

      if (std::isnan(lower) || std::isnan(upper) || std::isnan(alf))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}};
      else
	{
          fixed_gauss_laguerre_integral<Tp> integ(n, alf);
          return { integ(func, lower, upper) };
	}
    }

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_hermite(int n, Tp alf,
				  FuncTp func, Tp lower, Tp upper)
    {
      using integ_t = fixed_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;

      if (std::isnan(lower) || std::isnan(upper) || std::isnan(alf))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}};
      else
	{
          fixed_gauss_hermite_integral<Tp> integ(n, alf);
          return { integ(func, lower, upper) };
	}
    }

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_exponential(int n, Tp alf,
				      FuncTp func, Tp lower, Tp upper)
    {
      using integ_t = fixed_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;

      if (std::isnan(lower) || std::isnan(upper) || std::isnan(alf))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}};
      else
	{
          fixed_gauss_exponential_integral<Tp> integ(n, alf);
          return { integ(func, lower, upper) };
	}
    }

  template<typename Tp, typename FuncTp>
    fixed_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    integrate_fixed_gauss_rational(int n, Tp alf, Tp bet,
				   FuncTp func, Tp lower, Tp upper)
    {
      using integ_t = fixed_integral_t<Tp,
				     std::invoke_result_t<FuncTp, Tp>>;
      using area_t = typename integ_t::AreaTp;

      if (std::isnan(lower) || std::isnan(upper)
          || std::isnan(alf) || std::isnan(bet))
	{
	  const auto s_NaN = std::numeric_limits<Tp>::quiet_NaN();
	  return {area_t{} * s_NaN};
	}
      else if (lower == upper)
	return {area_t{}};
      else
	{
          fixed_gauss_rational_integral<Tp> integ(n, alf, bet);
          return { integ(func, lower, upper) };
	}
    }

} // namespace emsr

#endif // GAUSS_QUADRATURE_H
