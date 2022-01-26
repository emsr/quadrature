
#include <cmath>
#include <iostream>
#include <iomanip>

  template<typename Tp>
    void
    build_tanh_sinh_rule(int n = 128)
    {
      std::cout.precision(std::numeric_limits<Tp>::digits10);
      auto w = 8 + std::cout.precision();

      const auto s_pi_4 = emsr::math_constants<Tp>::pi_v / 4;

      n /= 2;

      // Find K = ln(ln(max_number))
      const auto k_max = std::log(std::log(std::numeric_limits<Tp>::max()))
		- Tp{1};
      auto h = k_max / n;

      std::cout << " k_max = " << k_max << '\n';
      std::cout << " h     = " << h << '\n';

      for (int k = -n; k <= n; ++k)
	{
	  const auto u = h * Tp(k);
	  const auto eu = std::exp(u);
	  const auto cosh = eu + Tp{1} / eu;
	  const auto sinh = eu - Tp{1} / eu;
          const auto esh = std::exp(s_pi_4 * sinh);
	  const auto w = esh + Tp{1} / esh;
	  const auto dxdu = cosh / (w * w);
	  const auto x = (esh - Tp{1} / esh) / w;
	  std::cout << ' ' << std::setw(w) << x << ' ' << std::setw(w) << dxdu << '\n';
	}
    }

  template<typename Tp>
    void
    build_sinh_sinh_rule(int n = 128)
    {
      std::cout.precision(std::numeric_limits<Tp>::digits10);
      auto w = 8 + std::cout.precision();

      const auto s_pi_4 = gnu_cxx::math_constants<Tp>::pi_quarter;

      n /= 2;

      // Find K = ln(ln(max_number))
      const auto k_max = std::log(std::log(std::numeric_limits<Tp>::max()))
		- Tp{1};
      auto h = k_max / n;

      std::cout << " k_max = " << k_max << '\n';
      std::cout << " h     = " << h << '\n';

      for (int k = -n; k <= n; ++k)
	{
	  const auto u = h * Tp(k);
	  const auto eu = std::exp(u);
	  const auto cosh = eu + Tp{1} / eu;
	  const auto sinh = eu - Tp{1} / eu;
          const auto esh = std::exp(s_pi_4 * sinh);
	  const auto x = (esh - Tp{1} / esh) / Tp{2};
	  const auto w = esh + Tp{1} / esh;
	  const auto dxdu = cosh * w / Tp{4};
	  std::cout << ' ' << std::setw(w) << x << ' ' << std::setw(w) << dxdu << '\n';
	}
    }

  template<typename Tp>
    void
    build_exp_sinh_rule(int n = 128)
    {
      std::cout.precision(std::numeric_limits<Tp>::digits10);
      auto w = 8 + std::cout.precision();

      const auto s_pi_4 = emsr::math_constants<Tp>::pi / 4;

      n /= 2;

      // Find K = ln(ln(max_number))
      const auto k_max = std::log(std::log(std::numeric_limits<Tp>::max()))
		- Tp{1};
      auto h = k_max / n;

      std::cout << " k_max = " << k_max << '\n';
      std::cout << " h     = " << h << '\n';

      for (int k = -n; k <= n; ++k)
	{
	  const auto u = h * Tp(k);
	  const auto eu = std::exp(u);
	  const auto cosh = eu + Tp{1} / eu;
	  const auto sinh = eu - Tp{1} / eu;
          const auto esh = std::exp(s_pi_4 * sinh);
	  const auto dxdu = cosh * esh;
	  std::cout << ' ' << std::setw(w) << esh << ' ' << std::setw(w) << dxdu << '\n';
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
