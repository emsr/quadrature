//
// Copyright 2008-2017
// Alion Science and Technology
// US Govt Retains rights in accordance
// with DoD FAR Supp 252.227 - 7013.
//


// FFT
// ===
/*
 | Discrete Fourier Transform can be regarded as evaluating a
 | polynomial of degree N-1 on the powers
 | omega^0, omega, omega^2, ..., omega^(N-1) where omega is a the
 | Nth root of unity.
 |
 | Given a polynomial of even degree
 |
 | p(t) = a_0 + a_1 t + a_2 t^2 + ...
 |
 | we find:
 |
 | p(t) = a_0 + a_2 t^2 + a_4 t^4 + ... + t(a_1 + a_3 t^2 + ...)
 |      = q_e(t^2) + t q_o (t^2)
 |
 | where q_e and q_o are polynomials formed with the even and odd
 | coefficients of p(t). Thus, we get
 |
 | p(1) = q_e(1) + omega^0 q_o(1)
 | p(omega) = q_e(omega^2) + omega^1 q_o(omega^2)
 | p(omega^2) = q_e(omega^4) + omega^2 q_o(omega^4)
 | ...
 |
 | Note how on the rhs the Fourier transforms of q_e and q_o appear. Thus:
 */

#include <numeric>
#include <vector>
#include <complex>
#include <ext/math_const.h>

namespace __gnu_cxx
{

  template<typename _Tp>
    _Tp
    rational_arg(std::size_t __i, std::size_t __m)
    {
      const auto _S_2pi = __gnu_cxx::__const_2_pi<_Tp>();
      return _S_2pi * _Tp(__i) / _Tp(__m);
    }

  template<typename _Tp>
    class DFT_Pseudoiterator
    : public std::iterator<std::input_iterator_tag,
                           std::complex<_Tp>,
                           std::ptrdiff_t,
                           std::complex<_Tp> const *,
                           std::complex<_Tp> const &>
    {

    private:

      std::complex<_Tp> __omega_pow_i;
      std::complex<_Tp> __omega_pow_ik;
      std::size_t __k;

    public:

      DFT_Pseudoiterator(_Tp __sign,
                	 std::size_t __i,
                	 std::size_t __length,
                	 bool __past_end = false)
      : __omega_pow_i(std::polar(_Tp{1},
		      -__sign * rational_arg<_Tp>(__i, __length))),
	__omega_pow_ik(_Tp{1}),
	__k(__past_end ? __length : 0)
      {}

      std::complex<_Tp>
      operator*() const
      { return __omega_pow_ik; }

      DFT_Pseudoiterator &
      operator++()
      {
	__omega_pow_ik *= __omega_pow_i;
	++__k;
	return *this;
      }

      DFT_Pseudoiterator
      operator++(int)
      {
	DFT_Pseudoiterator __dummy(*this);
	++__k;
	return __dummy;
      }

      bool
      operator==(DFT_Pseudoiterator const & __other) const
      {
	return (this->__omega_pow_i == __other.__omega_pow_i)
	    && (this->__k == __other.__k);
      }

      bool
      operator!=(DFT_Pseudoiterator const & __other) const
      { return ! (*this == __other); }

    }; // DFT_Pseudoiterator


  /**
   * Discreet Fourier transform.
   */
  template<typename _Tp>
    void
    dft(bool __do_forward, std::vector<std::complex<_Tp>> & __z)
    {
      const _Tp __sign(__do_forward ? _Tp{+1} : _Tp{-1});
      const std::size_t __length = __z.size();

      std::vector<std::complex<_Tp>> __result;
      __result.reserve(__length);

      // Do matrix multiplication
      for (std::size_t __i = 0; __i < __length; ++__i)
	{
	  DFT_Pseudoiterator __coefficient_iter(__sign, __i, __length);
	  __result.push_back(std::inner_product(__z.begin(), __z.end(),
			     __coefficient_iter, std::complex<_Tp>(0)));
	}

      // Rescale if forward
      if (__do_forward)
	{
	  const auto __len = _Tp(__length);
	  for (std::size_t __i = 0; __i < __length; ++__i)
	    __result[__i] /= __len;
	}

      // Copy data back (swap is fast!).
      __z.swap(__result);
    }

  /**
   *  Do Fast Fourier Transform.
   */
  template<typename _Tp>
    void
    fft(std::vector<std::complex<_Tp>> & __z)
    {
      std::size_t const __length = __z.size();
      if (__length % 2 == 1) // Too bad, we're odd.
	dft(true, __z);
      else // Good, we're even.
	{
	  // Let's divide and conquer!
	  std::vector<std::complex<_Tp>> __odd;
	  std::vector<std::complex<_Tp>> __even;
	  __odd.reserve(__length / 2);
	  __even.reserve(__length / 2);
	  for (auto __run = __z.cbegin(); __run != __z.cend(); ++__run)
	    {
              __even.push_back(*__run);
              ++__run;
              __odd.push_back(*__run);
	    }
	  fft(__even);
	  fft(__odd);
	  DFT_Pseudoiterator __omega_iter(_Tp{1}, 1, __length);
	  for (std::size_t __i = 0; __i < __length / 2; ++__i, ++__omega_iter)
	    {
              __z[__i] = 0.5 * (__even[__i] + *__omega_iter * __odd[__i]);
              // The next line works because omega^(length/2) = -1 !!!
              __z[__i + __length/2] = 0.5 * (__even[__i] - *__omega_iter * __odd[__i]);
	    }
	}
    }

  /**
   *  Do inverse Fast Fourier Transform.
   */
  template<typename _Tp>
    void
    ifft(std::vector<std::complex<_Tp>> & __z)
    {
      std::size_t const __length = __z.size();
      if (__length % 2 == 1) // Too bad, we're odd.
	dft(false, __z);
      else // Good, we are even.
	{
	  // Let's divide and conquer!
	  std::vector<std::complex<_Tp>> __odd;
	  std::vector<std::complex<_Tp>> __even;
	  __odd.reserve(__length / 2);
	  __even.reserve(__length / 2);
	  for (auto __run = __z.cbegin(); __run != __z.cend(); ++__run)
	    {
	      __even.push_back(*__run);
	      ++__run;
	      __odd.push_back(*__run);
	    }
	  ifft(__even);
	  ifft(__odd);
	  DFT_Pseudoiterator __omega_iter(_Tp{-1}, 1, __length);
	  for (std::size_t __i = 0; __i < __length / 2; ++__i, ++__omega_iter)
	    {
	      __z[__i] = __even[__i] + *__omega_iter * __odd[__i];
	      // The next line works because omega^(length/2) = -1 !!!
	      __z[__i + __length/2] = __even[__i] - *__omega_iter * __odd[__i];
	    }
	}
    }

  /**
   *  Do Fast Fourier Transform on input range.
   */
  template <typename ComplexIterator>
    void
    fft(bool __do_forward,
        ComplexIterator const & __from, ComplexIterator const & __to)
    {
      using _Cmplx = typename ComplexIterator::value_type;
      using _Tp = typename _Cmplx::value_type;
      std::vector<std::complex<_Tp>> __z(__from, __to);
      if (__do_forward)
	fft(__z);
      else
	ifft(__z);
      std::copy(__z.begin(), __z.end(), __from);
    }

} // namespace __gnu_cxx
