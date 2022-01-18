//
// Copyright (C) 2018-2020 Free Software Foundation, Inc.
// Copyright (C) 2021-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.

#ifndef MATRIX_H
#define MATRIX_H 1

namespace emsr
{

  template<typename RandAccIter, typename RandAccIterRHS>
    int
    s_tridiag(std::size_t n,
	      RandAccIter supd, RandAccIter diag, RandAccIter subd,
	      RandAccIterRHS rhs);

  template<typename RandAccIter, typename RandAccIterRHS>
    int
    s_tridiag_symm(std::size_t n,
		   RandAccIter& diag, RandAccIter& subd,
		   RandAccIterRHS& rhs);

} // namespace emsr

#include <ext/matrix.tcc>

#endif // MATRIX_H
