// -*- C++ -*-
// Integration utilities for C++.
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

#include <ext/matrix.h>

namespace __gnu_cxx
{

  /**
   * Compute Chebyshev moments at level @c level.
   */
  template<typename _Tp>
    void
    oscillatory_integration_table<_Tp>::
    compute_moments(_Tp __par, std::size_t __level)
    {
      std::array<_Tp, 28> __v;
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

	  _S_tridiag(__noeq, __dsup, __diag, __dsub, __v.begin() + 3);
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

	  _S_tridiag(__noeq, __dsup, __diag, __dsub, __v.begin() + 2);
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

} // namespace __gnu_cxx

#endif // OSCILLATORY_INTEGRATION_TABLE_TCC
