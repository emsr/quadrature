/* cxx_integration/testcase.h
 *
 * Copyright (C) 2016-2018 Free Software Foundation, Inc.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef QUADRATURE_TESTCASE_H
#define QUADRATURE_TESTCASE_H 1

template<typename Tp>
  struct monomial
  {
    int degree;
    Tp constant;

    monomial(int deg, Tp con)
    : degree(deg),
      constant(con)
    { }

    Tp
    operator()(Tp x) const
    { return constant * std::pow(x, degree); }

    monomial
    integral() const
    { return monomial(degree + 1, constant / Tp(degree + 1)); }

    monomial
    derivative() const
    {
      auto deg = std::max(0, degree - 1);
      return monomial(deg, Tp(deg) * constant);
    }
  };

template<typename Tp>
  Tp
  integrate(const monomial<Tp>& mon, Tp a, Tp b)
  {
    auto integ = mon.integral();
    return integ(b) - integ(a);
  }

template<typename Tp>
  struct func_test
  {
    typedef Tp (*fptr) (Tp);

    fptr fun;
    Tp a;
    Tp b;
    Tp exact;
  };

#include "testcase.tcc"

#endif // QUADRATURE_TESTCASE_H
