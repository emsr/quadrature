
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
// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

#ifndef NUMERIC_LIMITS_H
#define NUMERIC_LIMITS_H 1

#define _GLIBCXX_NO_SPECFUN 1
#include <cmath>
#undef _GLIBCXX_NO_SPECFUN

#include <limits>

namespace emsr
{

  /**
   *  @brief Part of std::numeric_limits.
   *  The idea is that types with, say, non-constexpr or even dynamic epsilon()
   *  can participate in this.
   *  I think variable templates could be specialized with non-constexpr types
   *  but I need something to work in C++11 and variable templates won't allow
   *  extraction of variable max from a mp number.
   */

  // Constexpr function template versions of std::numeric_limits.

  template<typename _Tp>
    constexpr bool
    is_specialized(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::is_specialized; }

  template<typename _Tp>
    constexpr _Tp
    lim_min(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::min(); }

  template<typename _Tp>
    constexpr _Tp
    lim_max(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::max(); }

  template<typename _Tp>
    constexpr _Tp
    lowest(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::lowest(); }

  template<typename _Tp>
    constexpr int
    digits(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::digits; }

  template<typename _Tp>
    constexpr int
    digits10(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::digits10; }

  template<typename _Tp>
    constexpr int
    max_digits10(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::max_digits10; }

  template<typename _Tp>
    constexpr bool
    is_signed(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::is_signed; }

  template<typename _Tp>
    constexpr bool
    is_integer(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::is_integer; }

  template<typename _Tp>
    constexpr bool
    is_exact(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::is_exact; }

  template<typename _Tp>
    constexpr int
    radix(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::radix; }

  template<typename _Tp>
    constexpr _Tp
    epsilon(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::epsilon(); }

  template<typename _Tp>
    constexpr _Tp
    round_error(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::round_error(); }

  template<typename _Tp>
    constexpr int
    min_exponent(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::min_exponent; }

  template<typename _Tp>
    constexpr int
    min_exponent10(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::min_exponent10; }

  template<typename _Tp>
    constexpr int
    max_exponent(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::max_exponent; }

  template<typename _Tp>
    constexpr int
    max_exponent10(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::max_exponent10; }

  template<typename _Tp>
    constexpr bool
    has_infinity(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::has_infinity; }

  template<typename _Tp>
    constexpr bool
    has_quiet_NaN(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::has_quiet_NaN; }

  template<typename _Tp>
    constexpr bool
    has_signaling_NaN(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::has_signaling_NaN; }

  template<typename _Tp>
    constexpr std::float_denorm_style
    has_denorm(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::has_denorm; }

  template<typename _Tp>
    constexpr bool
    has_denorm_loss(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::has_denorm_loss; }

  template<typename _Tp>
    constexpr _Tp
    infinity(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::infinity(); }

  template<typename _Tp>
    constexpr _Tp
    quiet_NaN(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::quiet_NaN(); }

  template<typename _Tp>
    constexpr _Tp
    signaling_NaN(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::signaling_NaN(); }

  template<typename _Tp>
    constexpr _Tp
    denorm_min(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::denorm_min(); }

  template<typename _Tp>
    constexpr bool
    is_iec559(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::is_iec559; }

  template<typename _Tp>
    constexpr bool
    is_bounded(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::is_bounded; }

  template<typename _Tp>
    constexpr bool
    is_modulo(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::is_modulo; }

  template<typename _Tp>
    constexpr bool
    traps(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::traps; }

  template<typename _Tp>
    constexpr bool
    tinyness_before(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::tinyness_before; }

  template<typename _Tp>
    constexpr std::float_round_style
    round_style(_Tp = _Tp{}) noexcept
    { return std::numeric_limits<_Tp>::round_style; }

  // Extra bits to help with numerics...
  // These depend on constexpr math functions.

  template<typename _Tp>
    inline _Tp
    max_integer(_Tp = _Tp{}) noexcept
    { return std::ldexp(_Tp{1}, digits(_Tp{})); }

  template<typename _Tp>
    inline _Tp
    sqrt_max(_Tp = _Tp{}) noexcept
    { return std::sqrt(lim_max(_Tp{})); }

  template<typename _Tp>
    inline _Tp
    cbrt_max(_Tp = _Tp{}) noexcept
    { return std::cbrt(lim_max(_Tp{})); }

  template<typename _Tp>
    inline _Tp
    root_max(_Tp root) noexcept
    { return std::pow(lim_max(_Tp{}), 1 / root); }

  template<typename _Tp>
    inline _Tp
    log_max(_Tp = _Tp{}) noexcept
    { return std::log(lim_max(_Tp{})); }

  template<typename _Tp>
    inline _Tp
    log10_max(_Tp = _Tp{}) noexcept
    { return std::log10(lim_max(_Tp{})); }


  template<typename _Tp>
    inline _Tp
    sqrt_min(_Tp = _Tp{}) noexcept
    { return std::sqrt(lim_min(_Tp{})); }

  template<typename _Tp>
    inline _Tp
    cbrt_min(_Tp = _Tp{}) noexcept
    { return std::cbrt(lim_min(_Tp{})); }

  template<typename _Tp>
    inline _Tp
    root_min(_Tp root) noexcept
    { return std::pow(lim_min(_Tp{}), 1 / root); }

  template<typename _Tp>
    inline _Tp
    log_min(_Tp = _Tp{}) noexcept
    { return std::log(lim_min(_Tp{})); }

  template<typename _Tp>
    inline _Tp
    log10_min(_Tp = _Tp{}) noexcept
    { return std::log10(lim_min(_Tp{})); }

  template<typename _Tp>
    inline _Tp
    sqrt_eps(_Tp = _Tp{}) noexcept
    { return std::sqrt(epsilon(_Tp{})); }

  template<typename _Tp>
    inline _Tp
    cbrt_eps(_Tp = _Tp{}) noexcept
    { return std::cbrt(epsilon(_Tp{})); }

  template<typename _Tp>
    inline _Tp
    root_eps(_Tp root) noexcept
    { return std::pow(epsilon(_Tp{}), 1 / root); }

  template<typename _Tp>
    inline _Tp
    log_eps(_Tp = _Tp{}) noexcept
    { return std::log(epsilon(_Tp{})); }

  template<typename _Tp>
    inline _Tp
    log10_eps(_Tp = _Tp{}) noexcept
    { return std::log10(epsilon(_Tp{})); }

} // namespace emsr

#endif // NUMERIC_LIMITS_H
