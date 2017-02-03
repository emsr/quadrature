/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I.. -DREPERIOD=0 -o test_fft test_fft.cpp -lquadmath
./test_fft > test_fft.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I.. -DREPERIOD=1 -o test_fft test_fft.cpp -lquadmath
./test_fft > test_fft_pi.txt
*/

#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <ext/cmath>
#include "fft.h"

template<typename _Tp>
  void
  test_fft()
  {
    const auto _S_2pi = __gnu_cxx::__const_2_pi<_Tp>();
    std::cout.precision(__gnu_cxx::__digits10<_Tp>());
    auto w = 8 + std::cout.precision();
    auto cw = 4 + 2 * w;
    const auto len = 1000u;

    std::default_random_engine re;
    std::uniform_real_distribution<_Tp> ud(_Tp{0}, _S_2pi);
    auto gen = [&ud, &re]()->_Tp{ return ud(re); };
    std::vector<std::complex<_Tp>> vec;
    vec.reserve(len);
    for (auto i = 0u; i < len; ++i)
      vec.push_back(std::polar(_Tp{1}, gen()));

    auto xform = vec;
    __gnu_cxx::fft(xform);

    auto iform = xform;
    __gnu_cxx::ifft(iform);

    auto mean_abs_diff = _Tp{0};
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
    mean_abs_diff /= len;
    std::cout << "mean_abs_diff = " << mean_abs_diff << '\n';
  }

int
main()
{
  test_fft<float>();

  test_fft<double>();

  test_fft<long double>();

  test_fft<__float128>();
}
