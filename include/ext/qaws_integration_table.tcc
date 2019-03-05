// quadrature/qaws_integration_table.tcc
//
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Brian Gough
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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

#ifndef QAWS_INTEGRATION_TABLE_TCC
#define QAWS_INTEGRATION_TABLE_TCC 1

namespace __gnu_cxx
{

  template<typename _Tp>
    qaws_integration_table<_Tp>::
    qaws_integration_table(_Tp __alpha, _Tp __beta,
			   int __mu, int __nu)
    : alpha(__alpha),
      beta(__beta),
      mu(__mu),
      nu(__nu)
    {
      if (this->alpha < _Tp{-1})
	std::__throw_domain_error("qaws_integration_table: "
				  "alalpha must be greater than -1.0");
      if (this->beta < _Tp{-1})
	std::__throw_domain_error("qaws_integration_table: "
				  "albeta must be greater than -1.0");
      if (this->mu != 0 && this->mu != 1)
	std::__throw_domain_error("qaws_integration_table: "
				  "almu must be 0 or 1");
      if (this->nu != 0 && this->nu != 1)
	std::__throw_domain_error("qaws_integration_table: "
				  "alnu must be 0 or 1");

      this->initialise();
    }

  template<typename _Tp>
    void
    qaws_integration_table<_Tp>::set(_Tp __alpha, _Tp __beta,
				     int __mu, int __nu)
    {
      if (__alpha < _Tp{-1})
	std::__throw_domain_error("qaws_integration_table: "
				  "alpha must be greater than -1.0");
      if (__beta < _Tp{-1})
	std::__throw_domain_error("qaws_integration_table: "
				  "beta must be greater than -1.0");
      if (__mu != 0 && __mu != 1)
	std::__throw_domain_error("qaws_integration_table: "
				  "mu must be 0 or 1");
      if (__nu != 0 && __nu != 1)
	std::__throw_domain_error("qaws_integration_table: "
				  "nu must be 0 or 1");

      this->alpha = __alpha;
      this->beta = __beta;
      this->mu = __mu;
      this->nu = __nu;

      this->initialise();
    }

  template<typename _Tp>
    void
    qaws_integration_table<_Tp>::initialise()
    {
      const auto __alpha_p1 = this->alpha + _Tp{1};
      const auto __beta_p1 = this->beta + _Tp{1};

      const auto __alpha_p2 = this->alpha + _Tp{2};
      const auto __beta_p2 = this->beta + _Tp{2};

      const auto __r_alpha = std::pow(_Tp{2}, __alpha_p1);
      const auto __r_beta = std::pow(_Tp{2}, __beta_p1);

      this->ri[0] = __r_alpha / __alpha_p1;
      this->ri[1] = this->ri[0] * this->alpha / __alpha_p2;
      auto __an = _Tp{2};
      auto __anm1 = _Tp{1};
      for (size_t __i = 2; __i < this->ri.size(); ++__i)
	{
	  this->ri[__i] = -(__r_alpha
			  + __an * (__an - __alpha_p2) * this->ri[__i - 1])
			/ (__anm1 * (__an + __alpha_p1));
	  __anm1 = __an;
	  __an += _Tp{1};
	}

      this->rj[0] = __r_beta / __beta_p1;
      this->rj[1] = this->rj[0] * this->beta / __beta_p2;
      __an = _Tp{2};
      __anm1 = _Tp{1};
      for (size_t __i = 2; __i < this->rj.size(); ++__i)
	{
	  this->rj[__i] = -(__r_beta
			  + __an * (__an - __beta_p2) * this->rj[__i - 1])
			/ (__anm1 * (__an + __beta_p1));
	  __anm1 = __an;
	  __an += _Tp{1};
	}

      this->rg[0] = -this->ri[0] / __alpha_p1;
      this->rg[1] = -this->rg[0]
		  - _Tp{2} * __r_alpha / (__alpha_p2 * __alpha_p2);
      __an = _Tp{2};
      __anm1 = _Tp{1};
      for (size_t __i = 2; __i < this->rg.size(); ++__i)
	{
	  this->rg[__i] = -(__an * (__an - __alpha_p2) * this->rg[__i - 1]
				- __an * this->ri[__i - 1]
				+ __anm1 * this->ri[__i])
			/ (__anm1 * (__an + __alpha_p1));
	  __anm1 = __an;
	  __an += _Tp{1};
	}

      this->rh[0] = -this->rj[0] / __beta_p1;
      this->rh[1] = -this->rh[0] - _Tp{2} * __r_beta / (__beta_p2 * __beta_p2);
      __an = _Tp{2};
      __anm1 = _Tp{1};
      for (size_t __i = 2; __i < this->rh.size(); ++__i)
	{
	  this->rh[__i] = -(__an * (__an - __beta_p2) * this->rh[__i - 1]
				- __an * this->rj[__i - 1]
				+ __anm1 * this->rj[__i])
			/ (__anm1 * (__an + __beta_p1));
	  __anm1 = __an;
	  __an += _Tp{1};
	}

      for (size_t __i = 1; __i < this->rj.size(); __i += 2)
	this->rj[__i] *= _Tp{-1};
      for (size_t __i = 1; __i < this->rh.size(); __i += 2)
	this->rh[__i] *= _Tp{-1};
    }

} // namespace __gnu_cxx

#endif // QAWS_INTEGRATION_TABLE_TCC
