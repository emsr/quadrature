
#include <iostream>
#include <iomanip>
#include <ext/cmath>

#include "integration.h"

template<typename Tp>
  void
  test_integral(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout.flags(std::ios::showpoint);
    const auto w = 8 + std::cout.precision();

    const auto PI = __gnu_cxx::__const_pi(proto);

    auto sin2 = [](Tp x) -> Tp { Tp s = std::sin(x); return s * s; };
    auto cos2 = [](Tp x) -> Tp { Tp c = std::cos(x); return c * c; };
    auto j1 = [](Tp x) -> Tp { return std::cyl_bessel_j(Tp{1}, x); };
    auto foo = [](Tp x) -> Tp { return (Tp{1} - x) * std::exp(-x / Tp{2}); };
    auto foonum = [](Tp x) -> Tp { return (Tp{1} - x); };
    auto funk1 = [PI](Tp x) -> Tp { return std::cos(x) / std::sqrt(x * (PI - x)); };
    auto funk1num = [](Tp x) -> Tp { return std::cos(x); };
    auto funk2 = [PI](Tp x) -> Tp { return (Tp{2} + std::sin(x)) / std::sqrt(x * (PI - x)); };
    auto funk2num = [](Tp x) -> Tp { return Tp{2} + std::sin(x); };
    auto one = [](Tp) -> Tp { return Tp{1}; };
    auto ex = [](Tp x) -> Tp { return x; };
    __gnu_cxx::_Polynomial<Tp> poly1({1.0l, -0.5l, -3.5l, 2.0l});

    auto fun = [](Tp x) -> Tp { return std::sin(x); };
    using fun_t = decltype(fun);
    auto mq = __gnu_cxx::double_exp_integrate4<Tp, fun_t>(fun, 20, Tp{0}, PI, Tp{0.0000001});
    std::cout << mq.__result << ' ' << mq.__abserr << '\n';

    Tp a = Tp{0};
    Tp b = PI;
    Tp err = Tp{0.0000000001};

    auto a0 = __gnu_cxx::double_exp_integrate4<Tp, decltype(one)>(one, 20, a, b, err);
    auto e0 = b - a;
    std::cout << "one     : "
	      << ' ' << std::setw(w) << a0.__result
	      << ' ' << std::setw(w) << e0
	      << ' ' << std::setw(w) << a0.__result - e0
	      //<< ' ' << std::setw(w) << t0.abs_error()
	      << '\n';

    auto a1 = __gnu_cxx::double_exp_integrate4<Tp, decltype(ex)>(ex, 20, a, b, err);
    auto e1 = (b * b - a * a) / Tp{2};
    std::cout << "ex      : "
	      << ' ' << std::setw(w) << a1.__result
	      << ' ' << std::setw(w) << e1
	      << ' ' << std::setw(w) << a1.__result - e1
	      //<< ' ' << std::setw(w) << t1.abs_error()
	      << '\n';

    auto a2 = __gnu_cxx::double_exp_integrate4<Tp, decltype(cos2)>(cos2, 20, a, b, err);
    auto e2 = PI / Tp{2};
    std::cout << "cos2    : "
	      << ' ' << std::setw(w) << a2.__result
	      << ' ' << std::setw(w) << e2
	      << ' ' << std::setw(w) << a2.__result - e2
	      //<< ' ' << std::setw(w) << t2.abs_error()
	      << '\n';

    auto a3 = __gnu_cxx::double_exp_integrate4<Tp, decltype(sin2)>(sin2, 20, a, b, err);
    auto e3 = PI / Tp{2};
    std::cout << "sin2    : "
	      << ' ' << std::setw(w) << a3.__result
	      << ' ' << std::setw(w) << e3
	      << ' ' << std::setw(w) << a3.__result - e3
	      //<< ' ' << std::setw(w) << t3.abs_error()
	      << '\n';

    auto a4 = __gnu_cxx::double_exp_integrate4<Tp, decltype(j1)>(j1, 20, a, b, err);
    auto e4 = std::cyl_bessel_j(Tp{0}, Tp{0}) - std::cyl_bessel_j(Tp{0}, PI);
    std::cout << "j1      : "
	      << ' ' << std::setw(w) << a4.__result
	      << ' ' << std::setw(w) << e4
	      << ' ' << std::setw(w) << a4.__result - e4
	      //<< ' ' << std::setw(w) << t4.abs_error()
	      << '\n';

    a = Tp{0};
    b = Tp{10} * PI;
    auto a5 = __gnu_cxx::double_exp_integrate4<Tp, decltype(foo)>(foo, 20, a, b, err);
    auto e5 = Tp{2} * (Tp{1} + b) * std::exp(-b / Tp{2})
	    - Tp{2} * (Tp{1} + a) * std::exp(-a / Tp{2});
    std::cout << "foo     : "
	      << ' ' << std::setw(w) << a5.__result
	      << ' ' << std::setw(w) << e5
	      << ' ' << std::setw(w) << a5.__result - e5
	      //<< ' ' << std::setw(w) << t5.abs_error()
	      << '\n';

    auto a5n = __gnu_cxx::double_exp_integrate4<Tp, decltype(foonum)>(foonum, 20, a, b, err);
    auto e5n = b * (Tp{1} - b / Tp{2})
	     - a * (Tp{1} - a / Tp{2});
    std::cout << "foonum  : "
	      << ' ' << std::setw(w) << a5n.__result
	      << ' ' << std::setw(w) << e5n
	      << ' ' << std::setw(w) << a5n.__result - e5n
	      //<< ' ' << std::setw(w) << t5n.abs_error()
	      << '\n';

    auto a6 = __gnu_cxx::double_exp_integrate4<Tp, __gnu_cxx::_Polynomial<Tp>>(poly1, 20, a, b, err);
    auto e6 = poly1.integral()(b) - poly1.integral()(a);
    std::cout << "poly1   : "
	      << ' ' << std::setw(w) << a6.__result
	      << ' ' << std::setw(w) << e6
	      << ' ' << std::setw(w) << a6.__result - e6
	      //<< ' ' << std::setw(w) << t6.abs_error()
	      << '\n';

    a = Tp{0};
    b = PI;
    auto a7 = __gnu_cxx::double_exp_integrate4<Tp, decltype(funk1)>(funk1, 20, a, b, err);
    auto e7 = Tp{0};
    std::cout << "funk1   : "
	      << ' ' << std::setw(w) << a7.__result
	      << ' ' << std::setw(w) << e7
	      << ' ' << std::setw(w) << a7.__result - e7
	      //<< ' ' << std::setw(w) << t7.abs_error()
	      << '\n';

    auto a7n = __gnu_cxx::double_exp_integrate4<Tp, decltype(funk1num)>(funk1num, 20, a, b, err);
    auto e7n = Tp{0};
    std::cout << "funk1num: "
	      << ' ' << std::setw(w) << a7n.__result
	      << ' ' << std::setw(w) << e7n
	      << ' ' << std::setw(w) << a7n.__result - e7n
	      //<< ' ' << std::setw(w) << t7n.abs_error()
	      << '\n';

    auto a8 = __gnu_cxx::double_exp_integrate4<Tp, decltype(funk2)>(funk2, 20, a, b, err);
    auto e8 = Tp{0};
    std::cout << "funk2   : "
	      << ' ' << std::setw(w) << a8.__result
	      << ' ' << std::setw(w) << e8
	      << ' ' << std::setw(w) << a8.__result - e8
	      //<< ' ' << std::setw(w) << t8.abs_error()
	      << '\n';

    auto a8n = __gnu_cxx::double_exp_integrate4<Tp, decltype(funk2num)>(funk2num, 20, a, b, err);
    auto e8n = Tp{2} * (b - a) - std::cos(b) + std::cos(a);
    std::cout << "funk2num: "
	      << ' ' << std::setw(w) << a8n.__result
	      << ' ' << std::setw(w) << e8n
	      << ' ' << std::setw(w) << a8n.__result - e8n
	      //<< ' ' << std::setw(w) << t8n.abs_error()
	      << '\n';
  }

int
main()
{
  test_integral<double>();
}
