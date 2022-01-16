
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <cmath>

#include <fourier_transform.h>

template<typename Tp>
  void
  test_fft_seq()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = 8 + std::cout.precision();
    auto cw = 4 + 2 * w;
    const auto len = 1000u;
    using Cmplx = std::complex<Tp>;

    // Complex doesn't have ++.
    std::vector<Cmplx> vec(len);
    for (auto k = 0u; k < len; ++k)
      vec[k] = k;
    auto xform = vec;
    __gnu_cxx::fast_fourier_transform(xform);

    auto mean_abs_diff = Tp{0};
    std::cout << '\n';
    auto diff = xform[0] - Cmplx(len * (len - 1) / Tp{2});
    auto abs_diff = std::abs(diff);
    mean_abs_diff += abs_diff;
    __gnu_cxx::__phase_iterator omega_k(Tp{-1}, 1, len);
    ++omega_k;
    for (auto k = 1u; k < len; ++k, ++omega_k)
      {
        auto xact = Tp(len) / (Tp{1} - *omega_k);
        diff = xform[k] - xact;
	abs_diff = std::abs(diff);
	mean_abs_diff += abs_diff;
	std::cout << ' ' << std::setw(6) << k
		  << ' ' << std::setw(cw) << vec[k]
		  << ' ' << std::setw(cw) << xact
		  << ' ' << std::setw(cw) << diff
		  << ' ' << std::setw(w) << abs_diff
		  << '\n';
      }
    mean_abs_diff /= Tp(len);
    std::cout << "mean_abs_diff = " << mean_abs_diff << '\n';
  }

template<typename Tp>
  void
  test_fft()
  {
    const auto s_2pi = 2 * Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = 8 + std::cout.precision();
    auto cw = 4 + 2 * w;
    const auto len = 1000u;

    std::default_random_engine re;
    std::uniform_real_distribution<Tp> ud(Tp{0}, s_2pi);
    auto gen = [&ud, &re]()->Tp{ return ud(re); };
    std::vector<std::complex<Tp>> vec;
    vec.reserve(len);
    for (auto i = 0u; i < len; ++i)
      vec.push_back(std::polar(Tp{1}, gen()));

    auto xform = vec;
    __gnu_cxx::fast_fourier_transform(xform);

    auto iform = xform;
    __gnu_cxx::inv_fast_fourier_transform(iform);

    auto mean_abs_diff = Tp{0};
    std::cout << '\n';
    for (auto i = 0u; i < len; ++i)
      {
	auto diff = vec[i] - iform[i];
	auto abs_diff = std::abs(diff);
	mean_abs_diff += abs_diff;
	std::cout << ' ' << std::setw(cw) << vec[i]
		  << ' ' << std::setw(cw) << xform[i]
		  << ' ' << std::setw(cw) << iform[i]
		  << ' ' << std::setw(cw) << diff
		  << ' ' << std::setw(w) << abs_diff
		  << '\n';
      }
    mean_abs_diff /= Tp(len);
    std::cout << "mean_abs_diff = " << mean_abs_diff << '\n';
  }

template<typename Tp>
  void
  test_real_fft()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = 8 + std::cout.precision();
    auto cw = 4 + 2 * w;
    const auto len = 1000u;

    std::default_random_engine re;
    std::uniform_real_distribution<Tp> ud(Tp{-5}, Tp{+5});
    auto gen = [&ud, &re]()->Tp{ return ud(re); };
    std::vector<Tp> vec;
    vec.reserve(len);
    for (auto i = 0u; i < len; ++i)
      vec.push_back(gen());

    auto xform = vec;
    __gnu_cxx::fast_fourier_transform(xform);

    auto iform = xform;
    __gnu_cxx::inv_fast_fourier_transform(iform);

    // This has the correct indexing you would want for fourier_transform_t<real_t>.
    auto mean_abs_diff = Tp{0};
    std::cout << '\n';
    for (auto i = 0u; i < len; ++i)
      {
	auto xf = i < len / 2
		? std::complex(xform[2 * i], xform[2 * i + 1])
		: std::complex(xform[2 * len - 2 * i - 2], -xform[2 * len - 2 * i - 1]);
	auto diff = vec[i] - iform[i];
	auto abs_diff = std::abs(diff);
	mean_abs_diff += abs_diff;
	std::cout << ' ' << std::setw(6) << i
		  << ' ' << std::setw(w) << vec[i]
		  << ' ' << std::setw(cw) << xf
		  << ' ' << std::setw(w) << iform[i]
		  << ' ' << std::setw(w) << diff
		  << ' ' << std::setw(w) << abs_diff
		  << '\n';
      }
    mean_abs_diff /= Tp(len);
    std::cout << "mean_abs_diff = " << mean_abs_diff << '\n';
  }

template<typename Tp>
  void
  test_fst()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = 8 + std::cout.precision();
    const auto len = 1000u;

    std::default_random_engine re;
    std::uniform_real_distribution<Tp> ud(Tp{-5}, Tp{+5});
    auto gen = [&ud, &re]()->Tp{ return ud(re); };
    std::vector<Tp> vec;
    vec.reserve(len);
    for (auto i = 0u; i < len; ++i)
      vec.push_back(gen());

    auto xform = vec;
    __gnu_cxx::fast_sine_transform(xform);

    auto iform = xform;
    __gnu_cxx::inv_fast_sine_transform(iform);

    // This has the correct indexing you would want for fourier_transform_t<real_t>.
    auto mean_abs_diff = Tp{0};
    std::cout << '\n';
    for (auto i = 0u; i < len; ++i)
      {
	auto diff = vec[i] - iform[i];
	auto abs_diff = std::abs(diff);
	mean_abs_diff += abs_diff;
	std::cout << ' ' << std::setw(6) << i
		  << ' ' << std::setw(w) << vec[i]
		  << ' ' << std::setw(w) << xform[i]
		  << ' ' << std::setw(w) << iform[i]
		  << ' ' << std::setw(w) << diff
		  << ' ' << std::setw(w) << abs_diff
		  << '\n';
      }
    mean_abs_diff /= Tp(len);
    std::cout << "mean_abs_diff = " << mean_abs_diff << '\n';
  }

int
main()
{
  test_fft_seq<double>();

  test_fft<float>();

  test_fft<double>();

  test_fft<long double>();

  test_fft<__float128>();

  test_real_fft<double>();

  test_fst<double>();
}
