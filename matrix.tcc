// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2018 Free Software Foundation, Inc.
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

namespace __gnu_cxx
{

  /**
   * Solve a tridiagonal system A x = b:
   *
   * @param[in]  dsup  The superdiagonal of the matrix in [0, n - 2]
   * @param[in]  diag  The diagonal of the matrix in [0, n - 1]
   * @param[in]  dsub  The subdiagonal of the matrix in [1, n - 1]
   *
   * @param[out]  rhs  The right hand side in [0 ... n - 1],
   *                   replaced by the solution vector x
   */
  template<typename _RandAccIter, typename _RandAccIterRHS>
    int
    _S_tridiag(size_t __n,
	       _RandAccIter __dsup, _RandAccIter __diag, _RandAccIter __dsub,
	       _RandAccIterRHS __rhs)
    {
      using _Tp = std::decay_t<decltype(__dsub[0])>;

      __dsub[0] = __diag[0];

      if (__n == 0)
	return 0;

      if (__n == 1)
	{
	  __rhs[0] = __rhs[0] / __diag[0];
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
	      std::swap(__rhs[__k1], __rhs[__k]);
	    }

	  if (__dsub[__k] == 0)
	    return 1;

	  {
	    const auto __t = -__dsub[__k1] / __dsub[__k];

	    __dsub[__k1] = __diag[__k1] + __t * __diag[__k];
	    __diag[__k1] = __dsup[__k1] + __t * __dsup[__k];
	    __dsup[__k1] = _Tp{0};
	    __rhs[__k1] = __rhs[__k1] + __t * __rhs[__k];
	  }
	}

      if (__dsub[__n - 1] == 0)
	return 1;

      __rhs[__n - 1] = __rhs[__n - 1] / __dsub[__n - 1];

      __rhs[__n - 2] = (__rhs[__n - 2] - __diag[__n - 2] * __rhs[__n - 1])
		   / __dsub[__n - 2];

      for (std::ptrdiff_t __k = __n ; __k > 2; --__k)
	{
	  const auto __kb = __k - 3;
	  __rhs[__kb] = (__rhs[__kb] - __diag[__kb] * __rhs[__kb + 1]
		     - __dsup[__kb] * __rhs[__kb + 2]) / __dsub[__kb];
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
   *  The changes consist (essentialy) of applying the orthogonal transformations
   *  directly to Z as they are generated.
   *
   *  Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
   *  C++ version by John Burkardt.
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
   *  @param[in,out] diag  The diagonal entries of the matrix. On output,
   *                       the information in diag has been overwritten.
   *  @param[in,out] dsub  The subdiagonal entries of the matrix,
   *                       in entries dsub[1] through dsub[n-1]. On output,
   *                       the information in dsub has been overwritten.
   *  @param[in,out] z     On input, a vector.  On output, the value of Q' * Z,
   *                       where Q is the matrix that diagonalizes the
   *                       input symmetric tridiagonal matrix.
   */
  template<typename _RandAccIter, typename _RandAccIterRHS>
    void
    _S_tridiag_symm(std::size_t n, _RandAccIter __diag, _RandAccIter __dsub,
		    _RandAccIterRHS z)
    {
      using _Tp = std::decay_t<decltype(__dsub[0])>;

      const int max_iter = 30;
      _Tp b;
      _Tp c;
      _Tp f;
      _Tp g;
      _Tp p;
      _Tp prec;
      _Tp r;
      _Tp s;
      int __upper = 0;

      prec = std::numeric_limits<_Tp>::epsilon();

      if (n == 1)
	return;

      __dsub[n-1] = _Tp{0};

      for (int __lower = 1; __lower <= n; ++__lower)
	{
	  int __j = 0;
	  while (true)
	    {
	      for (__upper = __lower; __upper <= n; ++__upper)
		{
        	  if (__upper == n)
        	    break;
        	  if (std::abs(__dsub[__upper-1])
			<= prec * (std::abs(__diag[__upper-1]) + std::abs(__diag[__upper])))
        	    break;
		}
	      auto p = __diag[__lower-1];
	      if (__upper == __lower)
        	break;
	      if (max_iter <= __j)
		{
        	  std::cerr << ": Iteration limit exceeded";
        	  exit(1);
		}
	      ++__j;
	      auto g = (__diag[__lower] - p) / (_Tp{2} * __dsub[__lower-1]);
	      auto r = std::sqrt(g * g + _Tp{1});
	      g = __diag[__upper-1] - p
		+ __dsub[__lower-1] / (g + std::copysign(std::abs(r), g));
	      auto s = _Tp{1};
	      auto c = _Tp{1};
	      p = _Tp{0};

	      const auto mml = __upper - __lower;
	      for (int __ii = 1; __ii <= mml; __ii++)
		{
        	  auto __i = __upper - __ii;
        	  auto f = s * __dsub[__i-1];
        	  auto b = c * __dsub[__i-1];

        	  if (std::abs(g) <= std::abs(f))
        	    {
        	      c = g / f;
        	      r =  std::sqrt (c * c + _Tp{1});
        	      __dsub[__i] = f * r;
        	      s = _Tp{1} / r;
        	      c *= s;
        	    }
        	  else
        	    {
        	      s = f / g;
        	      r =  sqrt (s * s + _Tp{1});
        	      __dsub[__i] = g * r;
        	      c = _Tp{1} / r;
        	      s *= c;
        	    }

        	  g = __diag[__i] - p;
        	  r = (__diag[__i-1] - g) * s + _Tp{2} * c * b;
        	  p = s * r;
        	  __diag[__i] = g + p;
        	  g = c * r - b;
        	  f = z[__i];
        	  z[__i] = s * z[__i-1] + c * f;
        	  z[__i-1] = c * z[__i-1] - s * f;
		}
	      __diag[__lower-1] = __diag[__lower-1] - p;
	      __dsub[__lower-1] = g;
	      __dsub[__upper-1] = _Tp{0};
	    }
	}

      // Sorting.
      for (int __ii = 2; __ii <= __upper; ++__ii)
	{
	  int __i = __ii - 1;
	  int __k = __i;
	  auto p = __diag[__i-1];

	  for (int __j = __ii; __j <= n; ++__j)
	    if (__diag[__j-1] < p)
              {
        	__k = __j;
        	p = __diag[__j-1];
              }

	  if (__k != __i)
	    {
	      __diag[__k-1] = __diag[__i-1];
	      __diag[__i-1] = p;
	      std::swap(z[__i-1], z[__k-1]);
	    }
	}

      return;
    }

} // namespace __gnu_cxx

#endif // MATRIX_H
