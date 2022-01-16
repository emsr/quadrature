// quadrature/qaws_integration_table.tcc
//
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Brian Gough
// Copyright (C) 2016-2020 Free Software Foundation, Inc.
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

namespace emsr
{

  template<typename Tp>
    qaws_integration_table<Tp>::
    qaws_integration_table(Tp alpha, Tp beta, int mu, int nu)
    : alpha(alpha), beta(beta),
      mu(mu), nu(nu)
    {
      if (this->alpha < Tp{-1})
	throw std::domain_error("qaws_integration_table: alalpha must be greater than -1.0");
      if (this->beta < Tp{-1})
	throw std::domain_error("qaws_integration_table: albeta must be greater than -1.0");
      if (this->mu != 0 && this->mu != 1)
	throw std::domain_error("qaws_integration_table: almu must be 0 or 1");
      if (this->nu != 0 && this->nu != 1)
	throw std::domain_error("qaws_integration_table: alnu must be 0 or 1");

      this->initialise();
    }

  template<typename Tp>
    void
    qaws_integration_table<Tp>::set(Tp alpha, Tp beta, int mu, int nu)
    {
      if (alpha < Tp{-1})
	throw std::domain_error("qaws_integration_table: alpha must be greater than -1.0");
      if (beta < Tp{-1})
	throw std::domain_error("qaws_integration_table: beta must be greater than -1.0");
      if (mu != 0 && mu != 1)
	throw std::domain_error("qaws_integration_table: mu must be 0 or 1");
      if (nu != 0 && nu != 1)
	throw std::domain_error("qaws_integration_table: nu must be 0 or 1");

      this->alpha = alpha;
      this->beta = beta;
      this->mu = mu;
      this->nu = nu;

      this->initialise();
    }

  template<typename Tp>
    void
    qaws_integration_table<Tp>::initialise()
    {
      const auto alpha_p1 = this->alpha + Tp{1};
      const auto beta_p1 = this->beta + Tp{1};

      const auto alpha_p2 = this->alpha + Tp{2};
      const auto beta_p2 = this->beta + Tp{2};

      const auto r_alpha = std::pow(Tp{2}, alpha_p1);
      const auto r_beta = std::pow(Tp{2}, beta_p1);

      this->ri[0] = r_alpha / alpha_p1;
      this->ri[1] = this->ri[0] * this->alpha / alpha_p2;
      auto an = Tp{2};
      auto anm1 = Tp{1};
      for (size_t i = 2; i < this->ri.size(); ++i)
	{
	  this->ri[i] = -(r_alpha
			  + an * (an - alpha_p2) * this->ri[i - 1])
			/ (anm1 * (an + alpha_p1));
	  anm1 = an;
	  an += Tp{1};
	}

      this->rj[0] = r_beta / beta_p1;
      this->rj[1] = this->rj[0] * this->beta / beta_p2;
      an = Tp{2};
      anm1 = Tp{1};
      for (size_t i = 2; i < this->rj.size(); ++i)
	{
	  this->rj[i] = -(r_beta
			  + an * (an - beta_p2) * this->rj[i - 1])
			/ (anm1 * (an + beta_p1));
	  anm1 = an;
	  an += Tp{1};
	}

      this->rg[0] = -this->ri[0] / alpha_p1;
      this->rg[1] = -this->rg[0]
		  - Tp{2} * r_alpha / (alpha_p2 * alpha_p2);
      an = Tp{2};
      anm1 = Tp{1};
      for (size_t i = 2; i < this->rg.size(); ++i)
	{
	  this->rg[i] = -(an * (an - alpha_p2) * this->rg[i - 1]
				- an * this->ri[i - 1]
				+ anm1 * this->ri[i])
			/ (anm1 * (an + alpha_p1));
	  anm1 = an;
	  an += Tp{1};
	}

      this->rh[0] = -this->rj[0] / beta_p1;
      this->rh[1] = -this->rh[0] - Tp{2} * r_beta / (beta_p2 * beta_p2);
      an = Tp{2};
      anm1 = Tp{1};
      for (size_t i = 2; i < this->rh.size(); ++i)
	{
	  this->rh[i] = -(an * (an - beta_p2) * this->rh[i - 1]
				- an * this->rj[i - 1]
				+ anm1 * this->rj[i])
			/ (anm1 * (an + beta_p1));
	  anm1 = an;
	  an += Tp{1};
	}

      for (size_t i = 1; i < this->rj.size(); i += 2)
	this->rj[i] *= Tp{-1};
      for (size_t i = 1; i < this->rh.size(); i += 2)
	this->rh[i] *= Tp{-1};
    }

} // namespace emsr

#endif // QAWS_INTEGRATION_TABLE_TCC
