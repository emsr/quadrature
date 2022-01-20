//
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Brian Gough
// Copyright (C) 2016-2020 Free Software Foundation, Inc.
// Copyright (C) 2021-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//

#ifndef QAWS_INTEGRATION_TABLE_H
#define QAWS_INTEGRATION_TABLE_H 1

#include <array>

namespace emsr
{

  /**
   * This structure manages integration of functions
   * with optional singular factors
   * @f[
   *   log(x - \alpha) log(\beta - x)
   * @f]
   */
  template<typename Tp>
    struct qaws_integration_table
    {
      Tp alpha;
      Tp beta;
      int mu;
      int nu;
      std::array<Tp, 25> ri;
      std::array<Tp, 25> rj;
      std::array<Tp, 25> rg;
      std::array<Tp, 25> rh;

      qaws_integration_table(Tp alpha, Tp beta, int mu, int nu);
      void set(Tp alpha, Tp beta, int mu, int nu);
      void initialise();
    };

} // namespace emsr

#include <emsr/qaws_integration_table.tcc>

#endif // QAWS_INTEGRATION_TABLE_H
