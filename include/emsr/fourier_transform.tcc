//
// Copyright (C) 2017-2020 Free Software Foundation, Inc.
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

#ifndef FOURIER_TRANSFORM_TCC
#define FOURIER_TRANSFORM_TCC 1

#include <stdexcept>
#include <numeric>
#include <vector>
#include <complex>

namespace emsr
{

  /**
   * Discrete Fourier transform on complex data.
   */
  template<typename Tp>
    void
    discrete_fourier_transform(bool do_forward,
				 std::vector<std::complex<Tp>>& z)
    {
      const auto sign(do_forward ? Tp{+1} : Tp{-1});
      const auto len = z.size();

      std::vector<std::complex<Tp>> result;
      result.reserve(len);

      // Do matrix multiplication.
      for (std::size_t i = 0; i < len; ++i)
	{
	  phase_iterator coefficient_iter(sign, i, len);
	  result.push_back(std::inner_product(z.begin(), z.end(),
			     coefficient_iter, std::complex<Tp>{0}));
	}

      // Rescale if forward.
      if (do_forward)
	{
	  const auto norm = Tp{1} / Tp(len);
	  for (std::size_t i = 0; i < len; ++i)
	    result[i] *= norm;
	}

      // Copy data back (swap is fast!).
      z.swap(result);
    }

  /**
   * Fast Fourier Transform on complex data.
   */
  template<typename Tp>
    void
    fast_fourier_transform(std::vector<std::complex<Tp>>& z)
    {
      const auto len = z.size();
      if (len % 2 == 1) // Too bad, we're odd.
	discrete_fourier_transform(true, z);
      else // Good, we're even.
	{
	  // Let's divide and conquer!
	  std::vector<std::complex<Tp>> odd;
	  std::vector<std::complex<Tp>> even;
	  const auto halflen = len / 2;
	  odd.reserve(halflen);
	  even.reserve(halflen);
	  for (auto run = z.cbegin(); run != z.cend(); ++run)
	    {
              even.push_back(*run);
              ++run;
              odd.push_back(*run);
	    }
	  fast_fourier_transform(even);
	  fast_fourier_transform(odd);
	  phase_iterator omega_iter(Tp{1}, 1, len);
	  for (std::size_t i = 0, j = halflen; i < halflen;
		++i, ++j, ++omega_iter)
	    {
              z[i] = (even[i] + *omega_iter * odd[i]) / Tp{2};
              // The next line works because omega^(length/2) = -1.
              z[j] = (even[i] - *omega_iter * odd[i]) / Tp{2};
	    }
	}
    }

  /**
   * Inverse Fast Fourier Transform on complex data.
   */
  template<typename Tp>
    void
    inv_fast_fourier_transform(std::vector<std::complex<Tp>>& z)
    {
      const std::size_t len = z.size();
      if (len % 2 == 1) // Too bad, we're odd.
	discrete_fourier_transform(false, z);
      else // Good, we are even.
	{
	  // Let's divide and conquer!
	  std::vector<std::complex<Tp>> odd;
	  std::vector<std::complex<Tp>> even;
	  const auto halflen = len / 2;
	  odd.reserve(halflen);
	  even.reserve(halflen);
	  for (auto run = z.cbegin(); run != z.cend(); ++run)
	    {
	      even.push_back(*run);
	      ++run;
	      odd.push_back(*run);
	    }
	  inv_fast_fourier_transform(even);
	  inv_fast_fourier_transform(odd);
	  phase_iterator omega_iter(Tp{-1}, 1, len);
	  for (std::size_t i = 0, j = halflen; i < halflen;
		++i, ++j, ++omega_iter)
	    {
	      z[i] = even[i] + *omega_iter * odd[i];
	      // The next line works because omega^(length/2) = -1.
	      z[j] = even[i] - *omega_iter * odd[i];
	    }
	}
    }

  /**
   * Fast Fourier Transform on real data.
   */
  template<typename Tp>
    void
    fast_fourier_transform(std::vector<Tp>& x)
    {
      const auto len = x.size();
      if (len % 2 == 1) // Too bad, we're odd.
	throw std::domain_error("fast_fourier_transform: "
				  "Real data must have even length.");
      else // Good, we're even.
	{
	  // @todo I really want to just make a view to a complex vector.
	  const auto halflen = len / 2;
	  std::vector<std::complex<Tp>> z;
	  z.reserve(halflen + 1);
	  for (std::size_t i = 0; i < halflen; ++i)
	    z.emplace_back(x[2 * i], x[2 * i + 1]);
	  fast_fourier_transform(z);
	  z.emplace_back(z[0]); // Use symmetry.  We need N/2 transform.
	  const auto s_i2 = std::complex<Tp>{0, 2};
	  phase_iterator omega_iter(Tp{+1}, 1, halflen);
	  for (std::size_t i = 0; i < halflen; ++i)
	    {
	      const auto z1 = z[i];
	      const auto z2 = std::conj(z[halflen - i]);
	      const auto f = (z1 + z2) / Tp{2}
			     + *omega_iter * (z1 - z2) / s_i2;
	      x[2 * i] = f.real();
	      x[2 * i + 1] = f.imag();
	    }
	}
    }

  /**
   * Inverse Fast Fourier Transform on real data.
   */
  template<typename Tp>
    void
    inv_fast_fourier_transform(std::vector<Tp>& x)
    {
      const auto len = x.size();
      if (len % 2 == 1) // Too bad, we're odd.
	throw std::domain_error("inv_fast_fourier_transform: "
				  "Real data must have even length.");
      else // Good, we're even.
	{
	  // @todo I really want to just make a view to a complex vector.
	  const auto halflen = len / 2;
	  std::vector<std::complex<Tp>> z;
	  z.reserve(halflen);
	  const auto s_i2 = std::complex<Tp>{0, 2};
	  phase_iterator omega_iter(Tp{-1}, 1, halflen);
	  for (std::size_t i = 0; i < halflen; ++i)
	    {
	      const auto z1 = std::complex<Tp>(x[2 * i], x[2 * i + 1]);
	      const auto z2 = std::complex<Tp>(x[len - 2 * i - 2], -x[len - 2 * i - 1]);
	      const auto ze = (z1 + z2) / Tp{2};
	      const auto zo = -*omega_iter * (z1 - z2) / s_i2;
	      z.emplace_back(ze + zo);
	    }
	  inv_fast_fourier_transform(z);
	  for (std::size_t i = 0; i < halflen; ++i)
	    {
	      x[2 * i] = z[i].real();
	      x[2 * i + 1] = z[i].imag();
	    }
	}
    }

  /**
   * Fast Sine Transform on real data.
   */
  template <typename Tp>
    void
    fast_sine_transform(std::vector<Tp>& x)
    {
      const auto len = x.size();
      const auto halflen = len / 2;
      const auto n2 = len - 1;
      phase_iterator omega_iter(Tp{+1}, 1, len);
      x[0] = Tp{0};
      for (std::size_t k = 1; k <= halflen; ++k, ++omega_iter)
	{
	  const auto y1 = omega_iter.sin()
			  * (x[k] + x[n2 - k]);
	  const auto y2 = (x[k] - x[n2 - k]) / Tp{2};
	  x[k] = y1 + y2;
	  x[n2 - k] = y1 - y2;
	}
      fast_fourier_transform(x);
      x[0] *= Tp{0.5L};
      x[1] = Tp{0};
      auto sum = Tp{0};
      for (std::size_t i = 2; i < halflen; i += 2)
	{
	  sum += std::exchange(x[i], x[i + 1]);
	  x[i + 1] = sum;
	}
    }

  /**
   * Fast Sine Transform on real data.
   */
  template <typename Tp>
    void
    inv_fast_sine_transform(std::vector<Tp>& x)
    {
      fast_sine_transform(x);
      const auto norm = Tp{2} / x.size();
      for (std::size_t i = 0; i < x.size(); ++i)
	x[i] *= norm;
    }

  /**
   * Fast Fourier Transform on input range.
   */
  template <typename CmplxIter>
    void
    fast_fourier_transform(const CmplxIter& from, const CmplxIter& to)
    {
      using Cmplx = typename CmplxIter::value_type;
      using Tp = typename Cmplx::value_type;
      std::vector<std::complex<Tp>> z(from, to);
      fast_fourier_transform(z);
      std::copy(z.begin(), z.end(), from);
    }

  /**
   * Inverse Fast Fourier Transform on input range.
   */
  template <typename CmplxIter>
    void
    inv_fast_fourier_transform(const CmplxIter& from, const CmplxIter& to)
    {
      using Cmplx = typename CmplxIter::value_type;
      using Tp = typename Cmplx::value_type;
      std::vector<std::complex<Tp>> z(from, to);
      inv_fast_fourier_transform(z);
      std::copy(z.begin(), z.end(), from);
    }

} // namespace emsr

#endif // FOURIER_TRANSFORM_TCC
