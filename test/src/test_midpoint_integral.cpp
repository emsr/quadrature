
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>

#include <emsr/integration.h>
#include <emsr/polynomial.h>

template<typename Tp>
  std::complex<Tp>
  test_cyl_hankel_2(Tp nu, Tp x)
  { return {std::cyl_bessel_j(nu, x), std::cyl_neumann(nu, x)}; }

template<typename Tp>
  void
  test_integral()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout.flags(std::ios::showpoint);
    const auto w = 8 + std::cout.precision();

    const auto s_pi = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};

    auto sin2 = [](Tp x) ->Tp { Tp s = std::sin(x); return s * s; };
    auto cos2 = [](Tp x) ->Tp { Tp c = std::cos(x); return c * c; };
    auto j1 = [](Tp x) ->Tp { return std::cyl_bessel_j(Tp{1}, x); };
    auto foo = [](Tp x) ->Tp { return (Tp{1} - x) * std::exp(-x / Tp{2}); };
    auto foonum = [](Tp x) ->Tp { return (Tp{1} - x); };
    auto funk1 = [s_pi](Tp x) ->Tp { return std::cos(x) / std::sqrt(x * (s_pi - x)); };
    auto funk1num = [](Tp x) ->Tp { return std::cos(x); };
    auto funk2 = [s_pi](Tp x) ->Tp { return (Tp{2} + std::sin(x)) / std::sqrt(x * (s_pi - x)); };
    auto funk2num = [](Tp x) ->Tp { return Tp{2} + std::sin(x); };
    auto one = [](Tp) ->Tp { return Tp{1}; };
    auto ex = [](Tp x) ->Tp { return x; };
    emsr::Polynomial<Tp> poly1({1.0l, -0.5l, -3.5l, 2.0l});
    auto chank2 = [](Tp x) ->std::complex<Tp> { return test_cyl_hankel_2(Tp{1}, x); };

    auto a = Tp{0};
    auto b = Tp(s_pi);
    auto abs_err = Tp{0};
    auto rel_err = Tp{1.0e-10};

    auto sine = [](Tp x) -> Tp { return std::sin(x); };
    emsr::midpoint_integral<Tp, decltype(sine)> mq(sine, a, b, abs_err, rel_err);
    std::cout << mq() << '\n';

    emsr::midpoint_integral<Tp, decltype(one)> t0(one, a, b, abs_err, rel_err);
    Tp a0 = t0();
    Tp e0 = b - a;
    std::cout << "one     : "
	      << ' ' << std::setw(w) << a0
	      << ' ' << std::setw(w) << e0
	      << ' ' << std::setw(w) << a0 - e0
	      << ' ' << std::setw(w) << t0.abs_error() << '\n';

    emsr::midpoint_integral<Tp, decltype(ex)> t1(ex, a, b, abs_err, rel_err);
    Tp a1 = t1();
    Tp e1 = (b * b - a * a) / Tp{2};
    std::cout << "ex      : "
	      << ' ' << std::setw(w) << a1
	      << ' ' << std::setw(w) << e1
	      << ' ' << std::setw(w) << a1 - e1
	      << ' ' << std::setw(w) << t1.abs_error() << '\n';

    emsr::midpoint_integral<Tp, decltype(cos2)> t2(cos2, a, b, abs_err, rel_err);
    Tp a2 = t2();
    Tp e2 = s_pi / Tp{2};
    std::cout << "cos2    : "
	      << ' ' << std::setw(w) << a2
	      << ' ' << std::setw(w) << e2
	      << ' ' << std::setw(w) << a2 - e2
	      << ' ' << std::setw(w) << t2.abs_error() << '\n';

    emsr::midpoint_integral<Tp, decltype(sin2)> t3(sin2, a, b, abs_err, rel_err);
    Tp a3 = t3();
    Tp e3 = s_pi / Tp{2};
    std::cout << "sin2    : "
	      << ' ' << std::setw(w) << a3
	      << ' ' << std::setw(w) << e3
	      << ' ' << std::setw(w) << a3 - e3
	      << ' ' << std::setw(w) << t3.abs_error() << '\n';

    emsr::midpoint_integral<Tp, decltype(j1)> t4(j1, a, b, abs_err, rel_err);
    Tp a4 = t4();
    Tp e4 = std::cyl_bessel_j(Tp{0}, Tp{0}) - std::cyl_bessel_j(Tp{0}, s_pi);
    std::cout << "j1      : "
	      << ' ' << std::setw(w) << a4
	      << ' ' << std::setw(w) << e4
	      << ' ' << std::setw(w) << a4 - e4
	      << ' ' << std::setw(w) << t4.abs_error() << '\n';

    a = Tp{0};
    b = Tp{10} * s_pi;
    emsr::midpoint_integral<Tp, decltype(foo)> t5(foo, a, b, abs_err, rel_err);
    Tp a5 = t5();
    Tp e5 = Tp{2} * (Tp{1} + b) * std::exp(-b / Tp{2})
	  - Tp{2} * (Tp{1} + a) * std::exp(-a / Tp{2});
    std::cout << "foo     : "
	      << ' ' << std::setw(w) << a5
	      << ' ' << std::setw(w) << e5
	      << ' ' << std::setw(w) << a5 - e5
	      << ' ' << std::setw(w) << t5.abs_error() << '\n';

    emsr::midpoint_integral<Tp, decltype(foonum)> t5n(foonum, a, b, abs_err, rel_err);
    Tp a5n = t5n();
    Tp e5n = b * (Tp{1} - b / Tp{2})
	   - a * (Tp{1} - a / Tp{2});
    std::cout << "foonum  : "
	      << ' ' << std::setw(w) << a5n
	      << ' ' << std::setw(w) << e5n
	      << ' ' << std::setw(w) << a5n - e5n
	      << ' ' << std::setw(w) << t5n.abs_error() << '\n';

    emsr::midpoint_integral<Tp, emsr::Polynomial<Tp>> t6(poly1, a, b, abs_err, rel_err);
    Tp a6 = t6();
    Tp e6 = poly1.integral()(b) - poly1.integral()(a);
    std::cout << "poly1   : "
	      << ' ' << std::setw(w) << a6
	      << ' ' << std::setw(w) << e6
	      << ' ' << std::setw(w) << a6 - e6
	      << ' ' << std::setw(w) << t6.abs_error() << '\n';

    a = Tp{0};
    b = s_pi;
    emsr::midpoint_integral<Tp, decltype(funk1)> t7(funk1, a, b, abs_err, rel_err);
    Tp a7 = t7();
    Tp e7 = Tp{0};
    std::cout << "funk1   : "
	      << ' ' << std::setw(w) << a7
	      << ' ' << std::setw(w) << e7
	      << ' ' << std::setw(w) << a7 - e7
	      << ' ' << std::setw(w) << t7.abs_error() << '\n';

    emsr::midpoint_integral<Tp, decltype(funk1num)> t7n(funk1num, a, b, abs_err, rel_err);
    Tp a7n = t7n();
    Tp e7n = Tp{0};
    std::cout << "funk1num: "
	      << ' ' << std::setw(w) << a7n
	      << ' ' << std::setw(w) << e7n
	      << ' ' << std::setw(w) << a7n - e7n
	      << ' ' << std::setw(w) << t7n.abs_error() << '\n';

    emsr::midpoint_integral<Tp, decltype(funk2)> t8(funk2, a, b, abs_err, rel_err);
    Tp a8 = t8();
    Tp e8 = Tp{0};
    std::cout << "funk2   : "
	      << ' ' << std::setw(w) << a8
	      << ' ' << std::setw(w) << e8
	      << ' ' << std::setw(w) << a8 - e8
	      << ' ' << std::setw(w) << t8.abs_error() << '\n';

    emsr::midpoint_integral<Tp, decltype(funk2num)> t8n(funk2num, a, b, abs_err, rel_err);
    Tp a8n = t8n();
    auto e8n = Tp{2} * (b - a) - std::cos(b) + std::cos(a);
    std::cout << "funk2num: "
	      << ' ' << std::setw(w) << a8n
	      << ' ' << std::setw(w) << e8n
	      << ' ' << std::setw(w) << a8n - e8n
	      << ' ' << std::setw(w) << t8n.abs_error() << '\n';

    emsr::midpoint_integral<Tp, decltype(chank2)> thank2(chank2, b / Tp{2}, b, abs_err, rel_err);
    auto ahank2 = thank2();
    auto reehank2 = std::cyl_bessel_j(Tp{0}, b / Tp{2}) - std::cyl_bessel_j(Tp{0}, b);
    auto imehank2 = std::cyl_neumann(Tp{0}, b / Tp{2}) - std::cyl_neumann(Tp{0}, b);
    auto ehank2 = std::complex(reehank2, -imehank2);
    std::cout << "cyl_hankel_2: "
	      << ' ' << std::setw(w) << ahank2
	      << ' ' << std::setw(w) << ehank2
	      << ' ' << std::setw(w) << ahank2 - ehank2
	      << ' ' << std::setw(w) << thank2.abs_error() << '\n';
  }

int
main()
{
  test_integral<double>();
}
