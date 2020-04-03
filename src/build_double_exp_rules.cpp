
#include <cmath>
#include <iostream>
#include <iomanip>

  template<typename _Tp>
    void
    build_tanh_sinh_rule(int __n = 128)
    {
      std::cout.precision(std::numeric_limits<_Tp>::digits10);
      auto w = 8 + std::cout.precision();

      const auto _S_pi_4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;

      __n /= 2;

      // Find K = ln(ln(max_number))
      const auto __k_max = std::log(std::log(std::numeric_limits<_Tp>::max()))
		- _Tp{1};
      auto __h = __k_max / __n;

      std::cout << " k_max = " << __k_max << '\n';
      std::cout << " h     = " << __h << '\n';

      for (int __k = -__n; __k <= __n; ++__k)
	{
	  const auto __u = __h * _Tp(__k);
	  const auto __eu = std::exp(__u);
	  const auto __cosh = __eu + _Tp{1} / __eu;
	  const auto __sinh = __eu - _Tp{1} / __eu;
          const auto __esh = std::exp(_S_pi_4 * __sinh);
	  const auto __w = __esh + _Tp{1} / __esh;
	  const auto __dxdu = __cosh / (__w * __w);
	  const auto __x = (__esh - _Tp{1} / __esh) / __w;
	  std::cout << ' ' << std::setw(w) << __x << ' ' << std::setw(w) << __dxdu << '\n';
	}
    }

  template<typename _Tp>
    void
    build_sinh_sinh_rule(int __n = 128)
    {
      std::cout.precision(std::numeric_limits<_Tp>::digits10);
      auto w = 8 + std::cout.precision();

      const auto _S_pi_4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;

      __n /= 2;

      // Find K = ln(ln(max_number))
      const auto __k_max = std::log(std::log(std::numeric_limits<_Tp>::max()))
		- _Tp{1};
      auto __h = __k_max / __n;

      std::cout << " k_max = " << __k_max << '\n';
      std::cout << " h     = " << __h << '\n';

      for (int __k = -__n; __k <= __n; ++__k)
	{
	  const auto __u = __h * _Tp(__k);
	  const auto __eu = std::exp(__u);
	  const auto __cosh = __eu + _Tp{1} / __eu;
	  const auto __sinh = __eu - _Tp{1} / __eu;
          const auto __esh = std::exp(_S_pi_4 * __sinh);
	  const auto __x = (__esh - _Tp{1} / __esh) / _Tp{2};
	  const auto __w = __esh + _Tp{1} / __esh;
	  const auto __dxdu = __cosh * __w / _Tp{4};
	  std::cout << ' ' << std::setw(w) << __x << ' ' << std::setw(w) << __dxdu << '\n';
	}
    }

  template<typename _Tp>
    void
    build_exp_sinh_rule(int __n = 128)
    {
      std::cout.precision(std::numeric_limits<_Tp>::digits10);
      auto w = 8 + std::cout.precision();

      const auto _S_pi_4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;

      __n /= 2;

      // Find K = ln(ln(max_number))
      const auto __k_max = std::log(std::log(std::numeric_limits<_Tp>::max()))
		- _Tp{1};
      auto __h = __k_max / __n;

      std::cout << " k_max = " << __k_max << '\n';
      std::cout << " h     = " << __h << '\n';

      for (int __k = -__n; __k <= __n; ++__k)
	{
	  const auto __u = __h * _Tp(__k);
	  const auto __eu = std::exp(__u);
	  const auto __cosh = __eu + _Tp{1} / __eu;
	  const auto __sinh = __eu - _Tp{1} / __eu;
          const auto __esh = std::exp(_S_pi_4 * __sinh);
	  const auto __dxdu = __cosh * __esh;
	  std::cout << ' ' << std::setw(w) << __esh << ' ' << std::setw(w) << __dxdu << '\n';
	}
    }

int
main()
{
  std::cout << "\n\n tanh_sinh\n===============\n";
  std::cout << " float\n---------------\n";
  build_tanh_sinh_rule<float>();
  std::cout << " double\n---------------\n";
  build_tanh_sinh_rule<double>();
  std::cout << " long double\n---------------\n";
  build_tanh_sinh_rule<long double>();
  std::cout << " __float128\n---------------\n";
  build_tanh_sinh_rule<__float128>();

  std::cout << "\n\n sinh_sinh\n===============\n";
  std::cout << " float\n---------------\n";
  build_sinh_sinh_rule<float>();
  std::cout << " double\n---------------\n";
  build_sinh_sinh_rule<double>();
  std::cout << " long double\n---------------\n";
  build_sinh_sinh_rule<long double>();
  std::cout << " __float128\n---------------\n";
  build_sinh_sinh_rule<__float128>();

  std::cout << "\n\n exp_sinh\n===============\n";
  std::cout << " float\n---------------\n";
  build_exp_sinh_rule<float>();
  std::cout << " double\n---------------\n";
  build_exp_sinh_rule<double>();
  std::cout << " long double\n---------------\n";
  build_exp_sinh_rule<long double>();
  std::cout << " __float128\n---------------\n";
  build_exp_sinh_rule<__float128>();
}
