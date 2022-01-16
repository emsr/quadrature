// quadrature/qcheb_integrate.tcc
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

#ifndef QCHEB_INTEGRATE_TCC
#define QCHEB_INTEGRATE_TCC 1

#include <type_traits>
#include <array>

namespace __gnu_cxx
{

  template<typename RetTp>
    struct chebyshev_integral_t
    {
      std::array<RetTp, 13> cheb12;
      std::array<RetTp, 25> cheb24;
    };

  template<typename Tp, typename FuncTp>
    auto
    qcheb_integrate(FuncTp func, Tp lower, Tp upper)
    -> chebyshev_integral_t<std::invoke_result_t<FuncTp, Tp>>
    {
      using RetTp = std::invoke_result_t<FuncTp, Tp>;

      chebyshev_integral_t<RetTp> out;
      auto& cheb12 = out.cheb12;
      auto& cheb24 = out.cheb24;
      RetTp fval[25], v[12];

      // These are the values of cos(pi*k/24) for k=1..11 needed for the
      // Chebyshev expansion of f(x).  These are the zeros of the Chebyshev
      // function of the second kind of order 23: U_23(x).

      constexpr Tp
      x[11]
      {
	9.914448613738104111442846968605486e-01L,
	9.659258262890682867486612158530536e-01L,
	9.238795325112867561257834975394469e-01L,
	8.660254037844386467595427060757126e-01L,
	7.933533402912351645734146973742314e-01L,
	7.071067811865475243919762573395221e-01L,
	6.087614290087206394044894932434070e-01L,
	4.999999999999999999855184455596035e-01L,
	3.826834323650897717110798781478690e-01L,
	2.588190451025207623287087436359508e-01L,
	1.305261922200515915256103766723547e-01L,
      };

      const auto center = (upper + lower) / Tp{2};
      const auto half_length = (upper - lower) / Tp{2};

      fval[0] = func(upper) / Tp{2};
      fval[12] = func(center);
      fval[24] = func(lower) / Tp{2};

      for (int i = 1; i < 12; ++i)
	{
	  const std::size_t j = 24 - i;
	  const auto u = half_length * x[i - 1];
	  fval[i] = func(center + u);
	  fval[j] = func(center - u);
	}

      for (int i = 0; i < 12; ++i)
	{
	  const std::size_t j = 24 - i;
	  v[i] = fval[i] - fval[j];
	  fval[i] = fval[i] + fval[j];
	}

      {
	const auto alam1 = v[0] - v[8];
	const auto alam2 = x[5] * (v[2] - v[6] - v[10]);
	cheb12[3] = alam1 + alam2;
	cheb12[9] = alam1 - alam2;
      }

      {
	const auto alam1 = v[1] - v[7] - v[9];
	const auto alam2 = v[3] - v[5] - v[11];
	{
	  const auto alam = x[2] * alam1 + x[8] * alam2;
	  cheb24[3] = cheb12[3] + alam;
	  cheb24[21] = cheb12[3] - alam;
	}

	{
	  const auto alam = x[8] * alam1 - x[2] * alam2;
	  cheb24[9] = cheb12[9] + alam;
	  cheb24[15] = cheb12[9] - alam;
	}
      }

      {
	const auto part1 = x[3] * v[4];
	const auto part2 = x[7] * v[8];
	const auto part3 = x[5] * v[6];

	{
	  const auto alam1 = v[0] + part1 + part2;
	  const auto alam2 = x[1] * v[2] + part3 + x[9] * v[10];

	  cheb12[1] = alam1 + alam2;
	  cheb12[11] = alam1 - alam2;
	}

	{
	  const Tp alam1 = v[0] - part1 + part2;
	  const Tp alam2 = x[9] * v[2] - part3 + x[1] * v[10];
	  cheb12[5] = alam1 + alam2;
	  cheb12[7] = alam1 - alam2;
	}
      }

      {
	const auto alam = (x[0] * v[1] + x[2] * v[3]
			  + x[4] * v[5] + x[6] * v[7]
			  + x[8] * v[9] + x[10] * v[11]);
	cheb24[1] = cheb12[1] + alam;
	cheb24[23] = cheb12[1] - alam;
      }

      {
	const auto alam = (x[10] * v[1] - x[8] * v[3]
			  + x[6] * v[5] - x[4] * v[7]
			  + x[2] * v[9] - x[0] * v[11]);
	cheb24[11] = cheb12[11] + alam;
	cheb24[13] = cheb12[11] - alam;
      }

      {
	const auto alam = (x[4] * v[1] - x[8] * v[3]
			  - x[0] * v[5] - x[10] * v[7]
			  + x[2] * v[9] + x[6] * v[11]);
	cheb24[5] = cheb12[5] + alam;
	cheb24[19] = cheb12[5] - alam;
      }

      {
	const auto alam = (x[6] * v[1] - x[2] * v[3]
			  - x[10] * v[5] + x[0] * v[7]
			  - x[8] * v[9] - x[4] * v[11]);
	cheb24[7] = cheb12[7] + alam;
	cheb24[17] = cheb12[7] - alam;
      }

      for (int i = 0; i < 6; ++i)
	{
	  const std::size_t j = 12 - i;
	  v[i] = fval[i] - fval[j];
	  fval[i] = fval[i] + fval[j];
	}

      {
	const auto alam1 = v[0] + x[7] * v[4];
	const auto alam2 = x[3] * v[2];

	cheb12[2] = alam1 + alam2;
	cheb12[10] = alam1 - alam2;
      }

      cheb12[6] = v[0] - v[4];

      {
	const auto alam = x[1] * v[1] + x[5] * v[3] + x[9] * v[5];
	cheb24[2] = cheb12[2] + alam;
	cheb24[22] = cheb12[2] - alam;
      }

      {
	const auto alam = x[5] * (v[1] - v[3] - v[5]);
	cheb24[6] = cheb12[6] + alam;
	cheb24[18] = cheb12[6] - alam;
      }

      {
	const auto alam = x[9] * v[1] - x[5] * v[3] + x[1] * v[5];
	cheb24[10] = cheb12[10] + alam;
	cheb24[14] = cheb12[10] - alam;
      }

      for (int i = 0; i < 3; ++i)
	{
	  const std::size_t j = 6 - i;
	  v[i] = fval[i] - fval[j];
	  fval[i] = fval[i] + fval[j];
	}

      cheb12[4] = v[0] + x[7] * v[2];
      cheb12[8] = fval[0] - x[7] * fval[2];

      {
	const auto alam = x[3] * v[1];
	cheb24[4] = cheb12[4] + alam;
	cheb24[20] = cheb12[4] - alam;
      }

      {
	const auto alam = x[7] * fval[1] - fval[3];
	cheb24[8] = cheb12[8] + alam;
	cheb24[16] = cheb12[8] - alam;
      }

      cheb12[0] = fval[0] + fval[2];

      {
	const auto alam = fval[1] + fval[3];
	cheb24[0] = cheb12[0] + alam;
	cheb24[24] = cheb12[0] - alam;
      }

      cheb12[12] = v[0] - v[2];
      cheb24[12] = cheb12[12];

      for (int i = 1; i < 12; ++i)
	cheb12[i] *= Tp{1} / Tp{6};

      cheb12[0] *= Tp{1} / Tp{12};
      cheb12[12] *= Tp{1} / Tp{12};

      for (int i = 1; i < 24; ++i)
	cheb24[i] *= Tp{1} / Tp{12};

      cheb24[0] *= Tp{1} / Tp{24};
      cheb24[24] *= Tp{1} / Tp{24};

      return out;
    }

} // namespace __gnu_cxx

#endif // QCHEB_INTEGRATE_TCC
