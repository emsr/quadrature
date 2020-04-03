
#include <iostream>
#include <iomanip>
#include <ext/integration_transform.h>

template<typename _Tp>
  void
  test_map_minf_pinf()
  {
    auto fun = [](_Tp x)->_Tp{ return x; };
    using func_t = decltype(fun);

    auto x_minf_pinf = [](auto t)->_Tp{ return -_Tp{1} / t + _Tp{1} / (_Tp{1} - t); };

    for (int i = 0; i <= 100; ++i)
      {
	auto x = _Tp{0.01} * i;
        auto y = __gnu_cxx::map_minf_pinf<_Tp, func_t>(fun)(x);
	std::cout << ' ' << std::setw(12) << x
		  << ' ' << std::setw(12) << x_minf_pinf(x)
		  << ' ' << std::setw(12) << y
		  << '\n';
      }
  }

template<typename _Tp>
  void
  test_map_a_pinf()
  {
    auto fun = [](_Tp x)->_Tp{ return x; };
    using func_t = decltype(fun);

    const auto a = _Tp{-3};

    auto x_a_pinf = [a](auto t)->_Tp{ return a + t / (_Tp{1} - t); };

    for (int i = 0; i <= 100; ++i)
      {
	auto x = _Tp{0.01} * i;
        auto y = __gnu_cxx::map_a_pinf<_Tp, func_t>(fun, a)(x);
	std::cout << ' ' << std::setw(12) << x
		  << ' ' << std::setw(12) << x_a_pinf(x)
		  << ' ' << std::setw(12) << y
		  << '\n';
      }
  }

template<typename _Tp>
  void
  test_map_minf_b()
  {
    auto fun = [](_Tp x)->_Tp{ return x; };
    using func_t = decltype(fun);

    const auto b = _Tp{5};

    auto x_inf_b = [b](auto t)->_Tp{ return b - (_Tp{1} - t) / t; };

    for (int i = 0; i <= 100; ++i)
      {
	auto x = _Tp{0.01} * i;
        auto y = __gnu_cxx::map_a_pinf<_Tp, func_t>(fun, b)(x);
	std::cout << ' ' << std::setw(12) << x
		  << ' ' << std::setw(12) << x_inf_b(x)
		  << ' ' << std::setw(12) << y
		  << '\n';
      }
  }

// You have dx/dt in here too.
int
main()
{
  std::cout << "\nminf_pinf\n";
  test_map_minf_pinf<double>();

  std::cout << "\na_pinf\n";
  test_map_a_pinf<double>();

  std::cout << "\nminf_b\n";
  test_map_minf_b<double>();
}
