// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2011-2020 Free Software Foundation, Inc.
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
// This file implements an extrapolation table for use in integration schemes
// Based on gsl/integration/qelg.c

#ifndef EXTRAPOLATION_TABLE_TCC
#define EXTRAPOLATION_TABLE_TCC 1

namespace __gnu_cxx
{

  template<typename AreaTp, typename AbsAreaTp>
    std::tuple<AreaTp, AbsAreaTp>
    extrapolation_table<AreaTp, AbsAreaTp>::qelg()
    {
      const auto cur_n = this->m_nn - 1;
      const auto current = this->m_rlist2[cur_n];

      const auto s_eps = std::numeric_limits<AbsAreaTp>::epsilon();
      // Less than max to prevent overflow.
      const auto s_max = std::numeric_limits<AbsAreaTp>::max() / AbsAreaTp{100};

      auto absolute = s_max;
      auto relative = AbsAreaTp{5} * s_eps * std::abs(current);

      const auto newelm = cur_n / 2;
      const auto n_orig = cur_n;
      auto n_final = cur_n;

      const auto nres_orig = this->m_nres;

      auto result = current;
      auto abserr = s_max;

      if (cur_n < 2)
	{
	  result = current;
	  abserr = std::max(absolute, relative);
	  return std::make_tuple(result, abserr);
	}

      this->m_rlist2[cur_n + 2] = this->m_rlist2[cur_n];
      this->m_rlist2[cur_n] = s_max;

      for (size_t ii = 0; ii < newelm; ++ii)
	{
	  auto res = this->m_rlist2[cur_n - 2 * ii + 2];
	  const auto e0 = this->m_rlist2[cur_n - 2 * ii - 2];
	  const auto e1 = this->m_rlist2[cur_n - 2 * ii - 1];
	  const auto e2 = res;

	  const auto e1abs = std::abs(e1);
	  const auto delta2 = e2 - e1;
	  const auto err2 = std::abs(delta2);
	  const auto tol2 = std::max(std::abs(e2), e1abs) * s_eps;
	  const auto delta3 = e1 - e0;
	  const auto err3 = std::abs(delta3);
	  const auto tol3 = std::max(e1abs, std::abs(e0)) * s_eps;

	  if (err2 <= tol2 && err3 <= tol3)
	    {
	      // If e0, e1 and e2 are equal to within machine accuracy,
	      // convergence is assumed.
	      result = res;
	      absolute = err2 + err3;
	      relative = 5 * s_eps * std::abs(res);
	      abserr = std::max(absolute, relative);
	      return std::make_tuple(result, abserr);
	    }

	  const auto e3 = this->m_rlist2[cur_n - 2 * ii];
	  this->m_rlist2[cur_n - 2 * ii] = e1;
	  const auto delta1 = e1 - e3;
	  const auto err1 = std::abs(delta1);
	  const auto tol1 = std::max(e1abs, std::abs(e3)) * s_eps;

	  // If two elements are very close to each other, omit a part of
	  // the table by adjusting the value of n.
	  if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3)
	    {
	      n_final = 2 * ii;
	      break;
	    }

	  const auto ss = AbsAreaTp{1} / delta1 + AbsAreaTp{1} / delta2
			  - AbsAreaTp{1} / delta3;

	  // Test to detect irregular behaviour in the table,
	  // and eventually omit a part of the table by adjusting
	  // the value of n.
	  if (std::abs(ss * e1) <= this->m_irreg_test)
	    {
	      n_final = 2 * ii;
	      break;
	    }

	  // Compute a new element and eventually adjust the value of result.
	  res = e1 + AbsAreaTp{1} / ss;
	  this->m_rlist2[cur_n - 2 * ii] = res;
	  {
	    const auto error = err2 + std::abs(res - e2) + err3;

	    if (error <= abserr)
	      {
		abserr = error;
		result = res;
	      }
	  }
	}

      // Shift the table.

      {
	const size_t limexp = 50 - 1;

	if (n_final == limexp)
	  n_final = 2 * (limexp / 2);
      }

      if (n_orig % 2 == 1)
	{
	  for (size_t ii = 0; ii <= newelm; ++ii)
	    this->m_rlist2[ii * 2 + 1] = this->m_rlist2[ii * 2 + 3];
	}
      else
	{
	  for (size_t ii = 0; ii <= newelm; ++ii)
	    this->m_rlist2[ii * 2] = this->m_rlist2[ii * 2 + 2];
	}

      if (n_orig != n_final)
	{
	  for (size_t ii = 0; ii <= n_final; ++ii)
	    this->m_rlist2[ii] = this->m_rlist2[ii + n_orig - n_final];
	}

      this->m_nn = n_final + 1;

      if (nres_orig < 3)
	{
	  this->m_res3la[nres_orig] = result;
	  abserr = s_max;
	}
      else
	{ // Compute error estimate.
	  abserr = (std::abs(result - this->m_res3la[2])
		    + std::abs(result - this->m_res3la[1])
		    + std::abs(result - this->m_res3la[0]));

	  this->m_res3la[0] = this->m_res3la[1];
	  this->m_res3la[1] = this->m_res3la[2];
	  this->m_res3la[2] = result;
	}

      /* In QUADPACK the variable table->nres is incremented at the top of
	qelg, so it increases on every call. This leads to the array
	res3la being accessed when its elements are still undefined, so I
	have moved the update to this point so that its value more
	useful. */

      this->m_nres = nres_orig + 1;

      abserr = std::max(abserr, 5 * s_eps * std::abs(result));

      return std::make_tuple(result, abserr);
    }

} // namespace __gnu_cxx

#endif // EXTRAPOLATION_TABLE_TCC
