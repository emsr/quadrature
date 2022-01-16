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
// This file implements an oscillatory integrand table for use
// in integration schemes.
// Based on gsl/integration/qmomof.c

#ifndef OSCILLATORY_INTEGRATION_TABLE_TCC
#define OSCILLATORY_INTEGRATION_TABLE_TCC 1

#include <ext/matrix.h>

namespace emsr
{

  /**
   * Compute Chebyshev moments at level @c level.
   */
  template<typename Tp>
    void
    oscillatory_integration_table<Tp>::
    compute_moments(Tp par, std::size_t level)
    {
      std::array<Tp, 28> v;
      std::array<Tp, 25> diag, dsub, dsup;

      const size_t noeq = 25;

      const auto par2 = par * par;
      const auto par4 = par2 * par2;
      const auto par22 = par2 + Tp{2};

      const auto sinpar = std::sin(par);
      const auto cospar = std::cos(par);

      //
      // Compute the Chebyschev moments with respect to cosine.
      //

      auto ac = Tp{8} * cospar;
      auto as = Tp{24} * par * sinpar;

      v[0] = Tp{2} * sinpar / par;
      v[1] = (Tp{8} * cospar + (Tp{2} * par2 - Tp{8})
			 * sinpar / par) / par2;
      v[2] = (Tp{32} * (par2 - Tp{12}) * cospar
	   + (Tp{2} * ((par2 - Tp{80}) * par2 + Tp{192}) * sinpar)
		 / par) / par4;

      if (std::abs(par) <= Tp{24})
	{
	  // Compute the moments as the solution of a boundary value
	  // problem using the asyptotic expansion as an endpoint.
	  auto an = Tp{6};
	  for (auto k = 0u; k < noeq - 1; ++k)
	    {
	      auto an2 = an * an;
	      diag[k] = Tp{-2} * (an2 - Tp{4})
			 * (par22 - Tp{2} * an2);
	      dsup[k] = (an - 1) * (an - Tp{2}) * par2;
	      dsub[k + 1] = (an + Tp{3}) * (an + Tp{4}) * par2;
	      v[k + 3] = as - (an2 - Tp{4}) * ac;
	      an += Tp{2};
	    }

	  const auto an2 = an * an;

	  diag[noeq - 1] = Tp{-2} * (an2 - Tp{4})
			  * (par22 - Tp{2} * an2);
	  v[noeq + 2] = as - (an2 - Tp{4}) * ac;
	  v[3] = v[3] - Tp{56} * par2 * v[2];

	  const auto ass = par * sinpar;
	  const auto asap = (((((Tp{210} * par2 - 1) * cospar
			    - (Tp{105} * par2 - Tp{63}) * ass) / an2
			   - (Tp{1} - Tp{15} * par2) * cospar
				 + Tp{15} * ass) / an2
			  - cospar + Tp{3} * ass) / an2
			 - cospar) / an2;
	  v[noeq + 2] -= Tp{2} * asap * par2
			   * (an - Tp{1}) * (an - Tp{2});

	  s_tridiag(noeq, dsup, diag, dsub, v.begin() + 3);
	}
      else
	{
	  // Compute the moments by forward recursion.
	  auto an = Tp{4};
	  for (auto k = 3u; k < 13u; ++k)
	    {
	      auto an2 = an * an;
	      v[k] = ((an2 - Tp{4})
		       * (Tp{2} * (par22 - Tp{2} * an2)
				 * v[k - 1] - ac)
		  + as
		  - par2 * (an + Tp{1}) * (an + Tp{2}) * v[k - 2])
		  / (par2 * (an - Tp{1}) * (an - Tp{2}));
	      an += Tp{2};
	    }
	}


      for (auto i = 0u; i < 13u; ++i)
	this->chebmo[25 * level + 2 * i] = v[i];

      //
      // Compute the Chebyschev moments with respect to sine.
      //

      v[0] = Tp{2} * (sinpar - par * cospar) / par2;
      v[1] = (Tp{18} - Tp{48} / par2) * sinpar / par2
	   + (Tp{-2} + Tp{48} / par2) * cospar / par;

      ac = Tp{-24} * par * cospar;
      as = Tp{-8} * sinpar;

      if (std::abs(par) <= 24)
	{
	  // Compute the moments as the solution of a boundary value
	  // problem using the asyptotic expansion as an endpoint.
	  auto an = Tp{5};
	  for (auto k = 0u; k < noeq - 1; ++k)
	    {
	      auto an2 = an * an;
	      diag[k] = -Tp{2} * (an2 - Tp{4})
			 * (par22 - Tp{2} * an2);
	      dsup[k] = (an - 1) * (an - Tp{2}) * par2;
	      dsub[k + 1] = (an + 3) * (an + Tp{4}) * par2;
	      v[k + 2] = ac + (an2 - Tp{4}) * as;
	      an += Tp{2};
	    }
	  const auto an2 = an * an;

	  diag[noeq - 1] = -Tp{2} * (an2 - Tp{4})
			     * (par22 - 2 * an2);
	  v[noeq + 1] = ac + (an2 - Tp{4}) * as;
	  v[2] = v[2] - Tp{42} * par2 * v[1];

	  const auto ass = par * cospar;
	  const auto asap = (((((Tp{105} * par2 - Tp{63}) * ass
			 - (Tp{210} * par2 - Tp{1}) * sinpar) / an2
		    + (Tp{15} * par2 - 1) * sinpar
		    - Tp{15} * ass) / an2 - sinpar - Tp{3} * ass)
			 / an2 - sinpar) / an2;
	  v[noeq + 1] -= Tp{2} * asap * par2
			   * (an - Tp{1}) * (an - Tp{2});

	  s_tridiag(noeq, dsup, diag, dsub, v.begin() + 2);
	}
      else
	{
	  // Compute the moments by forward recursion.
	  auto an = Tp{3};
	  for (auto k = 2u; k < 12u; ++k)
	    {
	      const auto an2 = an * an;
	      v[k] = ((an2 - Tp{4})
		       * (Tp{2} * (par22 - Tp{2} * an2)
				 * v[k - 1] + as)
		   + ac
		   - par2 * (an + Tp{1}) * (an + Tp{2}) * v[k - 2])
		   / (par2 * (an - Tp{1}) * (an - Tp{2}));
	      an += Tp{2};
	    }
	}

      for (auto i = 0u; i < 12u; ++i)
	this->chebmo[25 * level + 2 * i + 1] = v[i];
    }

} // namespace emsr

#endif // OSCILLATORY_INTEGRATION_TABLE_TCC
