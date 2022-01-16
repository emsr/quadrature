// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2017-2020 Free Software Foundation, Inc.
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


#ifndef FOURIER_TRANSFORM_H
#define FOURIER_TRANSFORM_H 1

#include <complex>
#include <vector>

namespace emsr
{

  /**
   * 
   */
  template<typename Tp>
    class phase_iterator
    : public std::iterator<std::input_iterator_tag,
                           std::complex<Tp>,
                           std::ptrdiff_t,
                           const std::complex<Tp>*,
                           const std::complex<Tp>&>
    {
    private:

      std::complex<Tp> m_omega_pow_i;
      std::complex<Tp> m_omega_pow_ik;
      std::size_t m_k;

      static Tp
      s_rational_arg(std::size_t i, std::size_t m)
      {
#if REPERIOD
	return Tp(2 * i) / Tp(m);
#else
	const auto s_2pi = Tp{2}
			  * Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
	return s_2pi * Tp(i) / Tp(m);
#endif
      }

    public:

      phase_iterator(Tp sign,
                       std::size_t i,
                       std::size_t len,
                       bool past_end = false)
      : m_omega_pow_i(std::polar(Tp{1},
				  -sign * s_rational_arg(i, len))),
	m_omega_pow_ik(Tp{1}),
	m_k(past_end ? len : 0)
      { }

      phase_iterator(Tp delta)
      : m_omega_pow_i(std::polar(Tp{1}, delta)),
	m_omega_pow_ik(Tp{1}),
	m_k(0)
      { }

      std::complex<Tp>
      operator*() const
      { return m_omega_pow_ik; }

      phase_iterator&
      operator++()
      {
	++this->m_k;
	this->m_omega_pow_ik *= this->m_omega_pow_i;
	return *this;
      }

      phase_iterator
      operator++(int)
      {
	phase_iterator dummy(*this);
	++this->m_k;
	this->m_omega_pow_ik *= this->m_omega_pow_i;
	return dummy;
      }

      Tp
      cos() const
      { return std::real(this->m_omega_pow_ik); }

      Tp
      sin() const
      { return std::imag(this->m_omega_pow_ik); }

      bool
      operator==(const phase_iterator& other) const
      {
	return (this->m_omega_pow_i == other.m_omega_pow_i)
	    && (this->m_k == other.m_k);
      }

      bool
      operator!=(const phase_iterator& other) const
      { return !(*this == other); }

    }; // phase_iterator

/**
 * Fast Fourier Transform
 *
 * Discrete Fourier Transform can be regarded as evaluating a
 * polynomial of degree N-1 on the powers
 * @f$ \omega^0, \omega, \omega^2, ..., \omega^(N-1) @f$
 * where @f$ \omega @f$ is the Nth root of unity.
 *
 * Given a polynomial of even degree
 * @f[
 *    p(t) = a_0 + a_1 t + a_2 t^2 + ...
 * @f]
 * we find:
 * @f[
 *    p(t) = a_0 + a_2 t^2 + a_4 t^4 + ... + t(a_1 + a_3 t^2 + ...)
 *         = q_{even}(t^2) + t q_{odd} (t^2)
 * @f]
 * where q_e and q_o are polynomials formed with the even and odd
 * coefficients of p(t). Thus, we get
 *
 * p(1) = q_{even}(1) + \omega^0 q_{odd}(1)
 * p(\omega) = q_{even}(\omega^2) + \omega^1 q_{odd}(\omega^2)
 * p(\omega^2) = q_{even}(\omega^4) + \omega^2 q_{odd}(\omega^4)
 * ...
 *
 * Note how on the RHS the Fourier transforms of q_e and q_o appear. Thus:
 */

  /**
   * Discrete Fourier transform type.
   */
  template<typename Tp>
    class fourier_transform_t
    {
    private:

      std::vector<std::complex<Tp>> m_xform;

    public:

      fourier_transform_t(std::size_t n)
      : m_xform{}
      { this->m_xform.reserve(n / 2 + 1); }

      fourier_transform_t(std::vector<Tp> data);

      std::size_t
      size() const
      { return 2 * this->m_xform.size() - 2; }

      std::complex<Tp>
      operator[](std::size_t k) const
      {
	if (k < this->m_xform.size())
	  return this->m_xform[k];
	else
	  return std::conj(this->m_xform[this->m_xform.size() - k]);
	// This is real array indexing.
	//if (k < len / 2)
	//	? std::complex(xform[2 * k], xform[2 * k + 1])
	//	: std::complex(xform[2 * len - 2 * i - 2],
	//		      -xform[2 * len - 2 * k - 1]);
      }
    };

  /**
   * Discrete Fourier transform specialization for transform of complex data.
   */
  template<typename Tp>
    class fourier_transform_t<std::complex<Tp>>
    {
    private:

      std::vector<std::complex<Tp>> m_xform;

    public:

      fourier_transform_t(std::size_t n)
      : m_xform{}
      { this->m_xform.reserve(n / 2 + 1); }

      fourier_transform_t(const std::vector<std::complex<Tp>>& data);

      std::size_t
      size() const
      { return 2 * this->m_xform.size() - 2; }

      std::complex<Tp>
      operator[](std::size_t k) const
      { return this->m_xform[k]; }
    };

  /**
   * Discrete Fourier transform on complex data.
   */
  template<typename Tp>
    void
    discrete_fourier_transform(bool do_forward,
				 std::vector<std::complex<Tp>>& z);

  /**
   * Fast Fourier Transform on complex data.
   */
  template<typename Tp>
    void
    fast_fourier_transform(std::vector<std::complex<Tp>>& z);

  /**
   * Inverse Fast Fourier Transform on complex data.
   */
  template<typename Tp>
    void
    inv_fast_fourier_transform(std::vector<std::complex<Tp>>& z);

  /**
   * Fast Sine Transform on real data.
   */
  template <typename Tp>
    void
    fast_sine_transform(std::vector<Tp>& x);

  /**
   * Fast Sine Transform on real data.
   */
  template <typename Tp>
    void
    inv_fast_sine_transform(std::vector<Tp>& x);

  /**
   * Fast Fourier Transform on real data.
   */
  template<typename Tp>
    void
    fast_fourier_transform(std::vector<Tp>& x);

  /**
   * Inverse Fast Fourier Transform on real data.
   */
  template<typename Tp>
    void
    inv_fast_fourier_transform(std::vector<Tp>& x);

  /**
   * Fast Fourier Transform on input range.
   */
  template <typename CmplxIter>
    void
    fast_fourier_transform(const CmplxIter& from, const CmplxIter& to);

  /**
   * Inverse Fast Fourier Transform on input range.
   */
  template <typename CmplxIter>
    void
    inv_fast_fourier_transform(const CmplxIter& from, const CmplxIter& to);

} // namespace emsr

#include <ext/fourier_transform.tcc>

#endif // FOURIER_TRANSFORM_H
