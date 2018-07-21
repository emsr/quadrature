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
   * @param[in]  subd  The subdiagonal of the matrix in [1, n - 1]
   *
   * @param[out]  rhs  The right hand side in [0, n - 1],
   *                   replaced by the solution vector x
   */
  template<typename _RandAccIter, typename _RandAccIterRHS>
    int
    _S_tridiag(size_t __n,
	       _RandAccIter __dsup, _RandAccIter __diag, _RandAccIter __subd,
	       _RandAccIterRHS __rhs)
    {
      using _Tp = std::decay_t<decltype(__subd[0])>;

      __subd[0] = __diag[0];

      if (__n == 0)
	return 0;

      if (__n == 1)
	{
	  __rhs[0] /= __diag[0];
	  return 0;
	}

      __diag[0] = __dsup[0];
      __dsup[0] = _Tp{0};
      __dsup[__n - 1] = _Tp{0};

      for (auto __k = 0u; __k < __n - 1; ++__k)
	{
	  const auto __k1 = __k + 1;

	  if (std::abs(__subd[__k1]) >= std::abs(__subd[__k]))
	    {
	      std::swap(__subd[__k1], __subd[__k]);
	      std::swap(__diag[__k1], __diag[__k]);
	      std::swap(__dsup[__k1], __dsup[__k]);
	      std::swap(__rhs[__k1], __rhs[__k]);
	    }

	  if (__subd[__k] == _Tp{0})
	    return 1;

	  {
	    const auto __t = -__subd[__k1] / __subd[__k];

	    __subd[__k1] = __diag[__k1] + __t * __diag[__k];
	    __diag[__k1] = __dsup[__k1] + __t * __dsup[__k];
	    __dsup[__k1] = _Tp{0};
	    __rhs[__k1] = __rhs[__k1] + __t * __rhs[__k];
	  }
	}

      if (__subd[__n - 1] == _Tp{0})
	return 1;

      __rhs[__n - 1] = __rhs[__n - 1] / __subd[__n - 1];

      __rhs[__n - 2] = (__rhs[__n - 2] - __diag[__n - 2] * __rhs[__n - 1])
		     / __subd[__n - 2];

      for (std::ptrdiff_t __k = __n ; __k > 2; --__k)
	{
	  const auto __kb = __k - 3;
	  __rhs[__kb] = (__rhs[__kb] - __diag[__kb] * __rhs[__kb + 1]
		     - __dsup[__kb] * __rhs[__kb + 2]) / __subd[__kb];
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
  template<typename _RandAccIter, typename _RandAccIterRHS>
    int
    _S_tridiag_symm(std::size_t __n, _RandAccIter __diag, _RandAccIter __subd,
		    _RandAccIterRHS __rhs)
    {
      using _Tp = std::decay_t<decltype(__subd[0])>;

      const int __max_iter = 50;
      int __upper = 0;

      const auto __prec = std::numeric_limits<_Tp>::epsilon();

      if (__n == 1)
	return 0;

      __subd[__n - 1] = _Tp{0};

      for (int __lower = 1; __lower <= int(__n); ++__lower)
	{
	  int __j = 0;
	  while (true)
	    {
	      for (__upper = __lower; __upper <= int(__n); ++__upper)
		{
        	  if (__upper == int(__n))
        	    break;
        	  if (std::abs(__subd[__upper - 1])
			<= __prec * (std::abs(__diag[__upper - 1])
				   + std::abs(__diag[__upper])))
        	    break;
		}
	      auto __p = __diag[__lower - 1];
	      if (__upper == __lower)
        	break;
	      if (__max_iter <= __j)
		std::__throw_runtime_error("_S_tridiag_symm:"
					   " Iteration limit exceeded");
	      ++__j;
	      auto __g = (__diag[__lower] - __p)
			 / (_Tp{2} * __subd[__lower - 1]);
	      auto __r = std::sqrt(__g * __g + _Tp{1});
	      __g = __diag[__upper - 1] - __p
		  + __subd[__lower - 1]
			/ (__g + std::copysign(std::abs(__r), __g));
	      auto __s = _Tp{1};
	      auto __c = _Tp{1};
	      __p = _Tp{0};

	      const auto __mml = __upper - __lower;
	      for (int __ii = 1; __ii <= __mml; __ii++)
		{
        	  const auto __i = __upper - __ii;
        	  auto __f = __s * __subd[__i - 1];
        	  auto __b = __c * __subd[__i - 1];

        	  if (std::abs(__g) <= std::abs(__f))
        	    {
        	      __c = __g / __f;
        	      __r =  std::sqrt(__c * __c + _Tp{1});
        	      __subd[__i] = __f * __r;
        	      __s = _Tp{1} / __r;
        	      __c *= __s;
        	    }
        	  else
        	    {
        	      __s = __f / __g;
        	      __r =  std::sqrt(__s * __s + _Tp{1});
        	      __subd[__i] = __g * __r;
        	      __c = _Tp{1} / __r;
        	      __s *= __c;
        	    }

        	  __g = __diag[__i] - __p;
        	  __r = (__diag[__i - 1] - __g) * __s + _Tp{2} * __c * __b;
        	  __p = __s * __r;
        	  __diag[__i] = __g + __p;
        	  __g = __c * __r - __b;
        	  __f = __rhs[__i];
        	  __rhs[__i] = __s * __rhs[__i - 1] + __c * __f;
        	  __rhs[__i - 1] = __c * __rhs[__i - 1] - __s * __f;
		}
	      __diag[__lower - 1] = __diag[__lower - 1] - __p;
	      __subd[__lower - 1] = __g;
	      __subd[__upper - 1] = _Tp{0};
	    }
	}

      // Sorting.
      for (int __ii = 2; __ii <= __upper; ++__ii)
	{
	  int __i = __ii - 1;
	  int __k = __i;
	  auto __p = __diag[__i - 1];

	  for (int __j = __ii; __j <= int(__n); ++__j)
	    if (__diag[__j - 1] < __p)
              {
        	__k = __j;
        	__p = __diag[__j - 1];
              }

	  if (__k != __i)
	    {
	      __diag[__k - 1] = __diag[__i - 1];
	      __diag[__i - 1] = __p;
	      std::swap(__rhs[__i - 1], __rhs[__k - 1]);
	    }
	}

      return 0;
    }

} // namespace __gnu_cxx

#endif // MATRIX_TCC
