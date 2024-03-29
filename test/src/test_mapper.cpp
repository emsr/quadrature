
#include <iostream>
#include <iomanip>

#include <emsr/integration_transform.h>

template<typename Tp>
  void
  test_map_minf_pinf()
  {
    auto fun = [](Tp x)->Tp{ return x; };
    using func_t = decltype(fun);

    auto x_minf_pinf = [](auto t)->Tp{ return -Tp{1} / t + Tp{1} / (Tp{1} - t); };

    for (int i = 0; i <= 100; ++i)
      {
	auto x = Tp{0.01} * i;
        auto y = emsr::map_minf_pinf<Tp, func_t>(fun)(x);
	std::cout << ' ' << std::setw(12) << x
		  << ' ' << std::setw(12) << x_minf_pinf(x)
		  << ' ' << std::setw(12) << y
		  << '\n';
      }
  }

template<typename Tp>
  void
  test_map_a_pinf()
  {
    auto fun = [](Tp x)->Tp{ return x; };
    using func_t = decltype(fun);

    const auto a = Tp{-3};

    auto x_a_pinf = [a](auto t)->Tp{ return a + t / (Tp{1} - t); };

    for (int i = 0; i <= 100; ++i)
      {
	auto x = Tp{0.01} * i;
        auto y = emsr::map_a_pinf<Tp, func_t>(fun, a)(x);
	std::cout << ' ' << std::setw(12) << x
		  << ' ' << std::setw(12) << x_a_pinf(x)
		  << ' ' << std::setw(12) << y
		  << '\n';
      }
  }

template<typename Tp>
  void
  test_map_minf_b()
  {
    auto fun = [](Tp x)->Tp{ return x; };
    using func_t = decltype(fun);

    const auto b = Tp{5};

    auto x_inf_b = [b](auto t)->Tp{ return b - (Tp{1} - t) / t; };

    for (int i = 0; i <= 100; ++i)
      {
	auto x = Tp{0.01} * i;
        auto y = emsr::map_a_pinf<Tp, func_t>(fun, b)(x);
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
