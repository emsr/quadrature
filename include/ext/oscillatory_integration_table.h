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
// in integration schemes
// Based on gsl/integration/qmomof.c

#ifndef OSCILLATORY_INTEGRATION_TABLE_H
#define OSCILLATORY_INTEGRATION_TABLE_H 1

namespace __gnu_cxx
{

  template<typename Tp>
    struct oscillatory_integration_table
    {
      enum circular_function
      {
	INTEG_COSINE,
	INTEG_SINE
      };

      std::size_t n;
      Tp omega;
      Tp length;
      Tp par;
      enum circular_function circfun;
      std::vector<Tp> chebmo;

      oscillatory_integration_table(Tp omega_in, Tp length_in,
				    circular_function circfun_in,
				    std::size_t n_in)
      : n(n_in),
	omega(omega_in),
	length(length_in),
	par(0.5 * omega_in * length_in),
	circfun(circfun_in),
	chebmo(25 * n_in)
      {
	auto scale = Tp{1};
	for (auto i = 0u; i < this->n; ++i)
	  {
	    this->compute_moments(this->par * scale, i);
	    scale *= 0.5;
	    // Prevent divide by zero.
	    if (const auto scale2 = scale * scale;
		scale2 * scale2 == Tp{0})
	      {
		this->n = i;
		break;
	      }
	  }
      }

      Tp
      get_length() const
      { return this->length; }

      oscillatory_integration_table<Tp>&
      set_length(Tp length_in)
      {
	this->length = length_in;
	this->par = 0.5 * this->omega * this->length;
	this->reset(this->omega, this->length, this->circfun);
	return *this;
      }

      void
      reset(Tp omega_in, Tp length_in,
	    circular_function circfun_in)
      {
	this->omega = omega_in;
	this->length = length_in;
	this->par = 0.5 * omega_in * length_in;
	this->circfun = circfun_in;
	auto scale = Tp{1};
	for (auto i = 0u; i < this->n; ++i)
	  {
	    this->compute_moments(this->par * scale, i);
	    scale *= 0.5;
	  }
      }

      inline const Tp*
      get_moments(std::size_t level) const
      { return this->chebmo.data() + 25 * level; }

    private:

      void compute_moments(Tp par, std::size_t level);

    };

} // namespace __gnu_cxx

#include <ext/oscillatory_integration_table.tcc>

#endif // OSCILLATORY_INTEGRATION_TABLE_H
