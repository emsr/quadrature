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

#ifndef MATRIX_H
#define MATRIX_H 1

namespace __gnu_cxx
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

} // namespace __gnu_cxx

#include <ext/matrix.tcc>

#endif // MATRIX_H
