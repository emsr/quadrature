// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2018-2020 Free Software Foundation, Inc.
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

#ifndef MATRIX_TCC
#define MATRIX_TCC 1

#include <type_traits> // For decay.

namespace emsr
{

  /**
   * Solve a tridiagonal system A x = b:
   *
   * @param[in]  supd  The superdiagonal of the matrix in [0, n - 2]
   * @param[in]  diag  The diagonal of the matrix in [0, n - 1]
   * @param[in]  subd  The subdiagonal of the matrix in [1, n - 1]
   *
   * @param[out]  rhs  The right hand side in [0, n - 1],
   *                   replaced by the solution vector x
   */
  template<typename RandAccIter, typename RandAccIterRHS>
    int
    s_tridiag(size_t n,
	       RandAccIter supd, RandAccIter diag, RandAccIter subd,
	       RandAccIterRHS rhs)
    {
      using Tp = std::decay_t<decltype(subd[0])>;

      subd[0] = diag[0];

      if (n == 0)
	return 0;

      if (n == 1)
	{
	  rhs[0] /= diag[0];
	  return 0;
	}

      diag[0] = supd[0];
      supd[0] = Tp{0};
      supd[n - 1] = Tp{0};

      for (auto k = 0u; k < n - 1; ++k)
	{
	  const auto k1 = k + 1;

	  if (std::abs(subd[k1]) >= std::abs(subd[k]))
	    {
	      std::swap(subd[k1], subd[k]);
	      std::swap(diag[k1], diag[k]);
	      std::swap(supd[k1], supd[k]);
	      std::swap(rhs[k1], rhs[k]);
	    }

	  if (subd[k] == Tp{0})
	    return 1;

	  {
	    const auto t = -subd[k1] / subd[k];

	    subd[k1] = diag[k1] + t * diag[k];
	    diag[k1] = supd[k1] + t * supd[k];
	    supd[k1] = Tp{0};
	    rhs[k1] = rhs[k1] + t * rhs[k];
	  }
	}

      if (subd[n - 1] == Tp{0})
	return 1;

      rhs[n - 1] = rhs[n - 1] / subd[n - 1];

      rhs[n - 2] = (rhs[n - 2] - diag[n - 2] * rhs[n - 1])
		     / subd[n - 2];

      for (std::ptrdiff_t k = n ; k > 2; --k)
	{
	  const auto kb = k - 3;
	  rhs[kb] = (rhs[kb] - diag[kb] * rhs[kb + 1]
		     - supd[kb] * rhs[kb + 2]) / subd[kb];
	}

      return 0;
    }

  /**
   *  @brief Diagonalize a symmetric tridiagonal matrix.
   *
   *  This routine is a slightly modified version of the EISPACK routine to 
   *  perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
   *
   *  It has been modified to produce the product Q' * Z, where Z is an input 
   *  vector and Q is the orthogonal matrix diagonalizing the input matrix.  
   *  The changes consist (essentialy) of applying the orthogonal
   *  transformations directly to Z as they are generated.
   *
   *  Original C++ version by John Burkardt (GPL).
   *
   *  @see Sylvan Elhay, Jaroslav Kautsky,
   *  Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
   *  Interpolatory Quadrature,
   *  ACM Transactions on Mathematical Software,
   *  Volume 13, Number 4, December 1987, pages 399-415.
   *
   *  @see Roger Martin, James Wilkinson,
   *  The Implicit QL Algorithm,
   *  Numerische Mathematik,
   *  Volume 12, Number 5, December 1968, pages 377-383.
   *
   *  @param[in]     n     The order of the matrix.
   *  @param[in,out] diag  The diagonal entries [0, n-1] of the matrix.
   *                       On output, the information in diag is overwritten.
   *  @param[in,out] subd  The subdiagonal entries of the matrix,
   *                       in entries [0, n-2]. On output,
   *                       the information in subd is overwritten.
   *  @param[in,out] rhs   On input, a vector.  On output, the value of Q' * Z,
   *                       where Q is the matrix that diagonalizes the
   *                       input symmetric tridiagonal matrix.
   */
  template<typename RandAccIter, typename RandAccIterRHS>
    int
    s_tridiag_symm(std::size_t n, RandAccIter& diag, RandAccIter& subd,
		   RandAccIterRHS& rhs)
    {
      using Tp = std::decay_t<decltype(subd[0])>;

      const int max_iter = 50;
      int upper = 0;

      const auto prec = std::numeric_limits<Tp>::epsilon();

      if (n == 1)
	return 0;

      subd[n - 1] = Tp{0};

      for (int lower = 1; lower <= int(n); ++lower)
	{
	  int j = 0;
	  while (true)
	    {
	      for (upper = lower; upper <= int(n); ++upper)
		{
        	  if (upper == int(n))
        	    break;
        	  if (std::abs(subd[upper - 1])
			<= prec * (std::abs(diag[upper - 1])
				   + std::abs(diag[upper])))
        	    break;
		}
	      auto p = diag[lower - 1];
	      if (upper == lower)
        	break;
	      if (max_iter <= j)
		throw std::runtime_error("s_tridiag_symm:"
					   " Iteration limit exceeded");
	      ++j;
	      auto g = (diag[lower] - p)
			 / (Tp{2} * subd[lower - 1]);
	      auto r = std::sqrt(g * g + Tp{1});
	      g = diag[upper - 1] - p
		  + subd[lower - 1]
			/ (g + std::copysign(std::abs(r), g));
	      auto s = Tp{1};
	      auto c = Tp{1};
	      p = Tp{0};

	      const auto mml = upper - lower;
	      for (int ii = 1; ii <= mml; ii++)
		{
        	  const auto i = upper - ii;
        	  auto f = s * subd[i - 1];
        	  auto b = c * subd[i - 1];

        	  if (std::abs(g) <= std::abs(f))
        	    {
        	      c = g / f;
        	      r =  std::sqrt(c * c + Tp{1});
        	      subd[i] = f * r;
        	      s = Tp{1} / r;
        	      c *= s;
        	    }
        	  else
        	    {
        	      s = f / g;
        	      r =  std::sqrt(s * s + Tp{1});
        	      subd[i] = g * r;
        	      c = Tp{1} / r;
        	      s *= c;
        	    }

        	  g = diag[i] - p;
        	  r = (diag[i - 1] - g) * s + Tp{2} * c * b;
        	  p = s * r;
        	  diag[i] = g + p;
        	  g = c * r - b;
        	  f = rhs[i];
        	  rhs[i] = s * rhs[i - 1] + c * f;
        	  rhs[i - 1] = c * rhs[i - 1] - s * f;
		}
	      diag[lower - 1] = diag[lower - 1] - p;
	      subd[lower - 1] = g;
	      subd[upper - 1] = Tp{0};
	    }
	}

      // Sorting.
      for (int ii = 2; ii <= upper; ++ii)
	{
	  int i = ii - 1;
	  int k = i;
	  auto p = diag[i - 1];

	  for (int j = ii; j <= int(n); ++j)
	    if (diag[j - 1] < p)
              {
        	k = j;
        	p = diag[j - 1];
              }

	  if (k != i)
	    {
	      diag[k - 1] = diag[i - 1];
	      diag[i - 1] = p;
	      std::swap(rhs[i - 1], rhs[k - 1]);
	    }
	}

      return 0;
    }

} // namespace emsr

#endif // MATRIX_TCC
