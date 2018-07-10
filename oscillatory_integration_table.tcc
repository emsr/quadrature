// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2011-2018 Free Software Foundation, Inc.
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
// This file implements an oscillatory integrand table for use
// in integration schemes.
// Based on gsl/integration/qmomof.c

#ifndef OSCILLATORY_INTEGRATION_TABLE_TCC
#define OSCILLATORY_INTEGRATION_TABLE_TCC 1

#include <type_traits> // For decay.

namespace __gnu_cxx
{
  template<typename _RandAccIter, typename _RandAccIterRHS>
    int
    _S_tridiag(size_t __n,
	       _RandAccIter __dsub, _RandAccIter __diag, _RandAccIter __dsup,
	       _RandAccIterRHS __b);

  /**
   * Compute Chebyshev moments at level @c level.
   */
  template<typename _Tp>
    void
    oscillatory_integration_table<_Tp>::
    compute_moments(_Tp __par, std::size_t __level)
    {
      _Tp __v[28];
      std::array<_Tp, 25> __diag, __dsub, __dsup;

      const size_t __noeq = 25;

      const auto __par2 = __par * __par;
      const auto __par4 = __par2 * __par2;
      const auto __par22 = __par2 + _Tp{2};

      const auto __sinpar = std::sin(__par);
      const auto __cospar = std::cos(__par);

      //
      // Compute the Chebyschev moments with respect to cosine.
      //

      auto __ac = _Tp{8} * __cospar;
      auto __as = _Tp{24} * __par * __sinpar;

      __v[0] = _Tp{2} * __sinpar / __par;
      __v[1] = (_Tp{8} * __cospar + (_Tp{2} * __par2 - _Tp{8})
			 * __sinpar / __par) / __par2;
      __v[2] = (_Tp{32} * (__par2 - _Tp{12}) * __cospar
	   + (_Tp{2} * ((__par2 - _Tp{80}) * __par2 + _Tp{192}) * __sinpar)
		 / __par) / __par4;

      if (std::abs(__par) <= _Tp{24})
	{
	  // Compute the moments as the solution of a boundary value
	  // problem using the asyptotic expansion as an endpoint.
	  auto __an = _Tp{6};
	  for (auto __k = 0u; __k < __noeq - 1; ++__k)
	    {
	      auto __an2 = __an * __an;
	      __diag[__k] = _Tp{-2} * (__an2 - _Tp{4})
			 * (__par22 - _Tp{2} * __an2);
	      __dsup[__k] = (__an - 1) * (__an - _Tp{2}) * __par2;
	      __dsub[__k + 1] = (__an + _Tp{3}) * (__an + _Tp{4}) * __par2;
	      __v[__k + 3] = __as - (__an2 - _Tp{4}) * __ac;
	      __an += _Tp{2};
	    }

	  const auto __an2 = __an * __an;

	  __diag[__noeq - 1] = _Tp{-2} * (__an2 - _Tp{4})
			  * (__par22 - _Tp{2} * __an2);
	  __v[__noeq + 2] = __as - (__an2 - _Tp{4}) * __ac;
	  __v[3] = __v[3] - _Tp{56} * __par2 * __v[2];

	  const auto __ass = __par * __sinpar;
	  const auto __asap = (((((_Tp{210} * __par2 - 1) * __cospar
			    - (_Tp{105} * __par2 - _Tp{63}) * __ass) / __an2
			   - (_Tp{1} - _Tp{15} * __par2) * __cospar
				 + _Tp{15} * __ass) / __an2
			  - __cospar + _Tp{3} * __ass) / __an2
			 - __cospar) / __an2;
	  __v[__noeq + 2] -= _Tp{2} * __asap * __par2
			   * (__an - _Tp{1}) * (__an - _Tp{2});

	  _S_tridiag(__noeq, __dsub, __diag, __dsup, __v + 3);
	}
      else
	{
	  // Compute the moments by forward recursion.
	  auto __an = _Tp{4};
	  for (auto __k = 3u; __k < 13u; ++__k)
	    {
	      auto __an2 = __an * __an;
	      __v[__k] = ((__an2 - _Tp{4})
		       * (_Tp{2} * (__par22 - _Tp{2} * __an2)
				 * __v[__k - 1] - __ac)
		  + __as
		  - __par2 * (__an + _Tp{1}) * (__an + _Tp{2}) * __v[__k - 2])
		  / (__par2 * (__an - _Tp{1}) * (__an - _Tp{2}));
	      __an += _Tp{2};
	    }
	}


      for (auto __i = 0u; __i < 13u; ++__i)
	this->chebmo[25 * __level + 2 * __i] = __v[__i];

      //
      // Compute the Chebyschev moments with respect to sine.
      //

      __v[0] = _Tp{2} * (__sinpar - __par * __cospar) / __par2;
      __v[1] = (_Tp{18} - _Tp{48} / __par2) * __sinpar / __par2
	   + (_Tp{-2} + _Tp{48} / __par2) * __cospar / __par;

      __ac = _Tp{-24} * __par * __cospar;
      __as = _Tp{-8} * __sinpar;

      if (std::abs(__par) <= 24)
	{
	  // Compute the moments as the solution of a boundary value
	  // problem using the asyptotic expansion as an endpoint.
	  auto __an = _Tp{5};
	  for (auto __k = 0u; __k < __noeq - 1; ++__k)
	    {
	      auto __an2 = __an * __an;
	      __diag[__k] = -_Tp{2} * (__an2 - _Tp{4})
			 * (__par22 - _Tp{2} * __an2);
	      __dsup[__k] = (__an - 1) * (__an - _Tp{2}) * __par2;
	      __dsub[__k + 1] = (__an + 3) * (__an + _Tp{4}) * __par2;
	      __v[__k + 2] = __ac + (__an2 - _Tp{4}) * __as;
	      __an += _Tp{2};
	    }
	  const auto __an2 = __an * __an;

	  __diag[__noeq - 1] = -_Tp{2} * (__an2 - _Tp{4})
			     * (__par22 - 2 * __an2);
	  __v[__noeq + 1] = __ac + (__an2 - _Tp{4}) * __as;
	  __v[2] = __v[2] - _Tp{42} * __par2 * __v[1];

	  const auto __ass = __par * __cospar;
	  const auto __asap = (((((_Tp{105} * __par2 - _Tp{63}) * __ass
			 - (_Tp{210} * __par2 - _Tp{1}) * __sinpar) / __an2
		    + (_Tp{15} * __par2 - 1) * __sinpar
		    - _Tp{15} * __ass) / __an2 - __sinpar - _Tp{3} * __ass)
			 / __an2 - __sinpar) / __an2;
	  __v[__noeq + 1] -= _Tp{2} * __asap * __par2
			   * (__an - _Tp{1}) * (__an - _Tp{2});

	  _S_tridiag(__noeq, __dsub, __diag, __dsup, __v + 2);
	}
      else
	{
	  // Compute the moments by forward recursion.
	  auto __an = _Tp{3};
	  for (auto __k = 2u; __k < 12u; ++__k)
	    {
	      const auto __an2 = __an * __an;
	      __v[__k] = ((__an2 - _Tp{4})
		       * (_Tp{2} * (__par22 - _Tp{2} * __an2)
				 * __v[__k - 1] + __as)
		   + __ac
		   - __par2 * (__an + _Tp{1}) * (__an + _Tp{2}) * __v[__k - 2])
		   / (__par2 * (__an - _Tp{1}) * (__an - _Tp{2}));
	      __an += _Tp{2};
	    }
	}

      for (auto __i = 0u; __i < 12u; ++__i)
	this->chebmo[25 * __level + 2 * __i + 1] = __v[__i];
    }

  /**
   * Solve a tridiagonal system A x = b:
   *
   *   dsub[1 ... n - 1]   subdiagonal of the matrix A
   *   diag[0 ... n - 1]   diagonal of the matrix A
   *   dsup[0 ... n - 2]   superdiagonal of the matrix A
   *
   *   b[0 ... n - 1]   right hand side, replaced by the solution vector x
   */
  template<typename _RandAccIter, typename _RandAccIterRHS>
    int
    _S_tridiag(size_t __n,
	       _RandAccIter __dsub, _RandAccIter __diag, _RandAccIter __dsup,
	       _RandAccIterRHS __b)
    {
      using _Tp = std::decay_t<decltype(__dsub[0])>;

      __dsub[0] = __diag[0];

      if (__n == 0)
	return 0;

      if (__n == 1)
	{
	  __b[0] = __b[0] / __diag[0];
	  return 0;
	}

      __diag[0] = __dsup[0];
      __dsup[0] = _Tp{0};
      __dsup[__n - 1] = 0;

      for (auto __k = 0u; __k < __n - 1; ++__k)
	{
	  const auto __k1 = __k + 1;

	  if (std::abs(__dsub[__k1]) >= std::abs(__dsub[__k]))
	    {
	      std::swap(__dsub[__k1], __dsub[__k]);
	      std::swap(__diag[__k1], __diag[__k]);
	      std::swap(__dsup[__k1], __dsup[__k]);
	      std::swap(__b[__k1], __b[__k]);
	    }

	  if (__dsub[__k] == 0)
	    return 1;

	  {
	    const auto __t = -__dsub[__k1] / __dsub[__k];

	    __dsub[__k1] = __diag[__k1] + __t * __diag[__k];
	    __diag[__k1] = __dsup[__k1] + __t * __dsup[__k];
	    __dsup[__k1] = _Tp{0};
	    __b[__k1] = __b[__k1] + __t * __b[__k];
	  }
	}

      if (__dsub[__n - 1] == 0)
	return 1;

      __b[__n - 1] = __b[__n - 1] / __dsub[__n - 1];

      __b[__n - 2] = (__b[__n - 2] - __diag[__n - 2] * __b[__n - 1])
		   / __dsub[__n - 2];

      for (std::ptrdiff_t __k = __n ; __k > 2; --__k)
	{
	  const auto __kb = __k - 3;
	  __b[__kb] = (__b[__kb] - __diag[__kb] * __b[__kb + 1]
		     - __dsup[__kb] * __b[__kb + 2]) / __dsub[__kb];
	}

      return 0;
    }

} // namespace __gnu_cxx

#endif // OSCILLATORY_INTEGRATION_TABLE_TCC
