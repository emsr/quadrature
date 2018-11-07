// qaws_integrate.tcc
//
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
// Copyright (C) 2016-2018 Free Software Foundation, Inc.
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
//
// Ported from GSL by Ed Smith-Rowland
// Originally written by Brian Gaugh
//
// Implements QAWS integration algorithm
// Based on gsl/integration/qaws.c

#ifndef QAWS_INTEGRATE_TCC
#define QAWS_INTEGRATE_TCC 1

#include <array>

#include "integration_workspace.h"
#include "qaws_integration_table.h"

namespace __gnu_cxx
{

  template<typename _Tp, typename _FuncTp>
    struct fn_qaws;

  template<typename _AreaTp>
    struct compute_result_t
    {
      _AreaTp __res12;
      _AreaTp __res24;
    };

  template<typename _Tp, typename _RetTp>
    auto
    compute_result(const std::array<_Tp, 25>& __r,
		   const std::array<_RetTp, 13>& __cheb12,
		   const std::array<_RetTp, 25>& __cheb24)
    -> compute_result_t<decltype(_Tp{} * _RetTp{})>;

 template<typename _Tp, typename _FuncTp>
    std::tuple<_Tp, _Tp, bool>
    qc25s(qaws_integration_table<_Tp, _FuncTp>& __t,
	  _FuncTp __func, _Tp __lower, _Tp __upper, _Tp __a1, _Tp __b1);

  /**
   * The singular weight function is defined by:
   * @f[
   *    W(x) = (x-a)^\alpha (b-x)^\beta log^\mu (x-a) log^\nu (b-x)
   * @f]
   * where @f$ \alpha > -1 @f$, @f$ \beta > -1 @f$,
   * and @f$ \mu = 0 @f$, 1, @f$ \nu = 0, 1 @f$.
   *
   * The weight function can take four different forms depending
   * on the values of \mu and \nu,
   * @f[
   *    W(x) = (x-a)^\alpha (b-x)^\beta                   (\mu = 0, \nu = 0)
   *    W(x) = (x-a)^\alpha (b-x)^\beta log(x-a)          (\mu = 1, \nu = 0)
   *    W(x) = (x-a)^\alpha (b-x)^\beta log(b-x)          (\mu = 0, \nu = 1)
   *    W(x) = (x-a)^\alpha (b-x)^\beta log(x-a) log(b-x) (\mu = 1, \nu = 1)
   * @f]
   *
   * The QAWS algorithm is designed for integrands with algebraic-logarithmic
   * singularities at the end-points of an integration region.
   *
   * In order to work efficiently the algorithm requires a precomputed table
   * of Chebyshev moments.
   */
  template<typename _Tp, typename _FuncTp>
    auto
    qaws_integrate(integration_workspace<_Tp,
			std::invoke_result_t<_FuncTp, _Tp>>& __workspace,
		   qaws_integration_table<_Tp, _FuncTp>& __table,
		   _FuncTp __func,
		   _Tp __lower, _Tp __upper,
		   _Tp __max_abs_err, _Tp __max_rel_err)
    -> adaptive_integral_t<_Tp, std::invoke_result_t<_FuncTp, _Tp>>
    {
      // Try to adjust tests for varing precision.
      const auto _M_rel_err = std::pow(_Tp{10},
				 -std::numeric_limits<_Tp>::digits / _Tp{10});

      if (__upper <= __lower)
	std::__throw_runtime_error("qaws_integrate: "
				   "Limits must form an ascending sequence");
      if (!valid_tolerances(__max_abs_err, __max_rel_err))
	{
	  std::ostringstream __msg;
	  __msg << "qaws_integrate: Tolerance cannot be achieved with given "
		   "absolute (" << __max_abs_err << ") and relative ("
		<< __max_rel_err << ") error limits.";
	  std::__throw_runtime_error(__msg.str().c_str());
	}

      const auto __limit = __workspace.capacity();
      __workspace.clear();

      // Perform the first integration.
      _Tp __result0, __abserr0;
      {
	const auto __a1 = __lower;
	const auto __mid = (__lower + __upper) / _Tp{2};
	const auto __a2 = __mid;
	const auto __b2 = __upper;

	_Tp __area1, __error1;
	bool __err_reliable1;
	std::tie(__area1, __error1, __err_reliable1)
	  = qc25s(__table, __func, __lower, __upper, __a1, __mid);
	__workspace.append(__a1, __mid, __area1, __error1);

	_Tp __area2, __error2;
	bool __err_reliable2;
	std::tie(__area2, __error2, __err_reliable2)
	  = qc25s(__table, __func, __lower, __upper, __a2, __b2);
	__workspace.append(__a2, __b2, __area2, __error2);

	__result0 = __area1 + __area2;
	__abserr0 = __error1 + __error2;
      }

      // Test on accuracy; Use 0.01 relative error as an extra safety
      // margin on the first iteration (ignored for subsequent iterations).
      auto __tolerance = std::max(__max_abs_err, __max_rel_err * std::abs(__result0));
      if (__abserr0 < __tolerance && __abserr0 < 0.01 * std::abs(__result0))
	return {__result0, __abserr0};
      else if (__limit == 1)
	__throw_integration_error("qaws_integrate: "
				  "a maximum of one iteration was insufficient",
				  MAX_ITER_ERROR, __result0, __abserr0);

      auto __area = __result0;
      auto __errsum = __abserr0;
      auto __iteration = 2u;
      int __error_type = NO_ERROR;
      int __roundoff_type1 = 0, __roundoff_type2 = 0;
      do
	{
	  // Bisect the subinterval with the largest error estimate.
	  const auto& __curr = __workspace.retrieve();

	  const auto __a1 = __curr.__lower_lim;
	  const auto __mid = (__curr.__lower_lim + __curr.__upper_lim) / _Tp{2};
	  const auto __a2 = __mid;
	  const auto __b2 = __curr.__upper_lim;

	  _Tp __area1, __error1;
	  bool __err_reliable1;
	  std::tie(__area1, __error1, __err_reliable1)
	    = qc25s(__table, __func, __lower, __upper, __a1, __mid);

	  _Tp __area2, __error2;
	  bool __err_reliable2;
	  std::tie(__area2, __error2, __err_reliable2)
	    = qc25s(__table, __func, __lower, __upper, __a2, __b2);

	  const auto __area12 = __area1 + __area2;
	  const auto __error12 = __error1 + __error2;

	  __errsum += __error12 - __curr.__abs_error;
	  __area += __area12 - __curr.__result;

	  if (__err_reliable1 && __err_reliable2)
	    {
	      const auto __delta = __curr.__result - __area12;

	      if (std::abs (__delta) <= _M_rel_err * std::abs(__area12)
		 && __error12 >= 0.99 * __curr.__abs_error)
		++__roundoff_type1;
	      if (__iteration >= 10 && __error12 > __curr.__abs_error)
		++__roundoff_type2;
	    }

	  __tolerance = std::max(__max_abs_err, __max_rel_err * std::abs(__area));
	  if (__errsum > __tolerance)
	    {
	      if (__roundoff_type1 >= 6 || __roundoff_type2 >= 20)
		__error_type = ROUNDOFF_ERROR;

	      // set error flag in the case of bad integrand behaviour at
	      // a point of the integration range
	      if (__workspace.subinterval_too_small(__a1, __a2, __b2))
		__error_type = SINGULAR_ERROR;
	    }

	  __workspace.split(__mid, __area1, __error1, __area2, __error2);

	  ++__iteration;
	}
      while (__iteration < __limit && !__error_type && __errsum > __tolerance);

      const auto __result = __workspace.total_integral();
      const auto __abserr = __errsum;

      if (__iteration == __limit)
	__error_type = MAX_SUBDIV_ERROR;

      if (__errsum <= __tolerance)
	return {__result, __abserr};

      if (__error_type == NO_ERROR)
	return {__result, __abserr};

      __check_error<_Tp>(__func__, __error_type, __result, __abserr);
      __throw_integration_error("qaws_integrate: Unknown error.",
				UNKNOWN_ERROR, __result, __abserr);
    }

  /**
   *
   */
  template<typename _Tp, typename _FuncTp>
    std::tuple<_Tp, _Tp, bool>
    qc25s(qaws_integration_table<_Tp, _FuncTp>& __t,
	  _FuncTp __func, _Tp __lower, _Tp __upper, _Tp __a1, _Tp __b1)
    {
      fn_qaws<_Tp, _FuncTp> __fqaws(&__t, __func, __lower, __upper);

      if (__a1 == __lower && (__t.alpha != _Tp{0} || __t.mu != 0))
	{
	  const auto __factor
	    = std::pow(0.5 * (__b1 - __a1), __t.alpha + _Tp{1});

	  auto __f = [__fqaws](_Tp __x)
		     -> _Tp { return __fqaws.eval_right(__x); };
	  auto __chout = qcheb_integrate(__f, __a1, __b1);
	  const auto& __cheb12 = __chout.__cheb12;
	  const auto& __cheb24 = __chout.__cheb24;

	  if (__t.mu == 0)
	    {
	      const auto __u = __factor;

	      auto [__res12, __res24]
		= compute_result(__t.ri, __cheb12, __cheb24);

	      const auto __result = __u * __res24;
	      const auto __abserr = std::abs(__u * (__res24 - __res12));
	      return std::make_tuple(__result, __abserr, false);
	    }
	  else
	    {
	      const auto __u = __factor * std::log(__b1 - __a1);
	      const auto __v = __factor;

	      auto [__res12a, __res24a]
		= compute_result(__t.ri, __cheb12, __cheb24);

	      auto [__res12b, __res24b]
		= compute_result(__t.rg, __cheb12, __cheb24);

	      const auto __result = __u * __res24a + __v * __res24b;
	      const auto __abserr = std::abs(__u * (__res24a - __res12a))
				  + std::abs(__v * (__res24b - __res12b));
	      return std::make_tuple(__result, __abserr, false);
	    }
	}
      else if (__b1 == __upper && (__t.beta != _Tp{0} || __t.nu != 0))
	{
	  auto __factor = std::pow(0.5 * (__b1 - __a1), __t.beta + _Tp{1});

	  auto __f = [__fqaws](_Tp __x)
		     -> _Tp { return __fqaws.eval_left(__x); };
	  auto __chout = qcheb_integrate(__f, __a1, __b1);
	  const auto& __cheb12 = __chout.__cheb12;
	  const auto& __cheb24 = __chout.__cheb24;

	  if (__t.nu == 0)
	    {
	      const auto __u = __factor;

	      auto [__res12, __res24]
		= compute_result(__t.rj, __cheb12, __cheb24);

	      const auto __result = __u * __res24;
	      const auto __abserr = std::abs(__u * (__res24 - __res12));
	      return std::make_tuple(__result, __abserr, false);
	    }
	  else
	    {
	      const auto __u = __factor * std::log(__b1 - __a1);
	      const auto __v = __factor;

	      auto [__res12a, __res24a]
		= compute_result(__t.rj, __cheb12, __cheb24);

	      auto [__res12b, __res24b]
		= compute_result(__t.rh, __cheb12, __cheb24);

	      const auto __result = __u * __res24a + __v * __res24b;
	      const auto __abserr = std::abs(__u * (__res24a - __res12a))
				  + std::abs(__v * (__res24b - __res12b));
	      return std::make_tuple(__result, __abserr, false);
	    }
	}
      else
	{
	  auto __f = [__fqaws](_Tp __x)
		     ->_Tp
		     { return __fqaws.eval_middle(__x); };

	  auto [__result, __abserr, __resabs, __resasc]
	    = qk_integrate(__f, __a1, __b1, Kronrod_15);

	  bool __err_reliable;
	  if (__abserr == __resasc)
	    __err_reliable = false;
	  else
	    __err_reliable = true;

	  return std::make_tuple(__result, __abserr, __err_reliable);
	}
    }

  /*
   *
   */
  template<typename _Tp, typename _FuncTp>
    struct fn_qaws
    {
      using _RetTp = std::invoke_result_t<_FuncTp, _Tp>;
      using _AreaTp = decltype(_RetTp{} * _Tp{});

      const qaws_integration_table<_Tp, _FuncTp>* table;

      _FuncTp func;
      _Tp a;
      _Tp b;

      fn_qaws(const qaws_integration_table<_Tp, _FuncTp>* __tab,
	      _FuncTp __func, _Tp __a_in, _Tp __b_in)
      : table(__tab),
	func(__func), a(__a_in), b(__b_in)
      { }

      _RetTp eval_middle(_Tp) const;
      _RetTp eval_left(_Tp) const;
      _RetTp eval_right(_Tp) const;
    };

  template<typename _Tp, typename _FuncTp>
    std::invoke_result_t<_FuncTp, _Tp>
    fn_qaws<_Tp, _FuncTp>::eval_middle(_Tp __x) const
    {
      auto __factor = _Tp{1};

      if (this->table->alpha != _Tp{0})
	__factor *= std::pow(__x - this->a, this->table->alpha);

      if (table->mu == 1)
	__factor *= std::log(__x - this->a);

      if (this->table->beta != _Tp{0})
	__factor *= std::pow(this->b - __x, this->table->beta);

      if (table->nu == 1)
	__factor *= std::log(this->b - __x);

      return __factor * this->func(__x);
    }

  template<typename _Tp, typename _FuncTp>
    std::invoke_result_t<_FuncTp, _Tp>
    fn_qaws<_Tp, _FuncTp>::eval_left(_Tp __x) const
    {
      auto __factor = _Tp{1};

      if (this->table->alpha != _Tp{0})
	__factor *= std::pow(__x - this->a, this->table->alpha);

      if (this->table->mu == 1)
	__factor *= std::log(__x - this->a);

      return __factor * this->func(__x);
    }

  template<typename _Tp, typename _FuncTp>
    std::invoke_result_t<_FuncTp, _Tp>
    fn_qaws<_Tp, _FuncTp>::eval_right(_Tp __x) const
    {
      auto __factor = _Tp{1};

      if (this->table->beta != _Tp{0})
	__factor *= std::pow(this->b - __x, this->table->beta);

      if (this->table->nu == 1)
	__factor *= std::log(this->b - __x);

      return __factor * this->func(__x);
    }

  template<typename _Tp, typename _RetTp>
    auto
    compute_result(const std::array<_Tp, 25>& __r,
		   const std::array<_RetTp, 13>& __cheb12,
		   const std::array<_RetTp, 25>& __cheb24)
    -> compute_result_t<decltype(_Tp{} * _RetTp{})>
    {
      using _AreaTp = decltype(_RetTp{} * _Tp{});

      auto __res12 = _AreaTp{};
      for (size_t __i = 0; __i < __cheb12.size(); ++__i)
	__res12 += __r[__i] * __cheb12[__i];

      auto __res24 = _AreaTp{};
      for (size_t __i = 0; __i < __cheb24.size(); ++__i)
	__res24 += __r[__i] * __cheb24[__i];

      return compute_result_t<_AreaTp>{__res12, __res24};
    }

} // namespace __gnu_cxx

#endif // QAWS_INTEGRATE_TCC
