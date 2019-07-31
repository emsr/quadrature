
#include <limits>
#include <iostream>
#include <iomanip>

#include <cmath>
#include <ext/integration.h>
#include <ext/math_constants.h>

template<typename Tp>
  void
  test_tanh_sinh()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout.flags(std::ios::showpoint);
    const auto w = 8 + std::cout.precision();

    const auto pi = __gnu_cxx::numbers::__pi_v<Tp>;
    const auto sqrt_pi = __gnu_cxx::numbers::__root_pi_v<Tp>;

    auto sin2 = [](Tp x) -> Tp { Tp s = std::sin(x); return s * s; };
    auto cos2 = [](Tp x) -> Tp { Tp c = std::cos(x); return c * c; };
    auto j1 = [](Tp x) -> Tp { return std::cyl_bessel_j(Tp{1}, x); };
    auto foo = [](Tp x) -> Tp { return (Tp{1} - x) * std::exp(-x / Tp{2}); };
    auto foonum = [](Tp x) -> Tp { return (Tp{1} - x); };
    auto funk1 = [pi](Tp x) -> Tp { return std::cos(x) / std::sqrt(x * (pi - x)); };
    auto funk1num = [](Tp x) -> Tp { return std::cos(x); };
    auto funk2 = [pi](Tp x) -> Tp { return (Tp{2} + std::sin(x)) / std::sqrt(x * (pi - x)); };
    auto funk2num = [](Tp x) -> Tp { return Tp{2} + std::sin(x); };
    auto one = [](Tp) -> Tp { return Tp{1}; };
    auto ex = [](Tp x) -> Tp { return x; };
    __gnu_cxx::_Polynomial<Tp> poly1({1.0l, -0.5l, -3.5l, 2.0l});

    auto a = Tp{0};
    auto b = pi;
    auto abs_err = Tp{0};
    auto rel_err = Tp{1.0e-10};

    auto sine = [](Tp x) -> Tp { return std::sin(x); };
    auto mq = __gnu_cxx::integrate_tanh_sinh<Tp, decltype(sine)>(sine, a, b, abs_err, rel_err);
    std::cout << mq.__result << ' ' << mq.__abserr << '\n';

    auto a0 = __gnu_cxx::integrate_tanh_sinh<Tp>(one, a, b, abs_err, rel_err);
    auto e0 = b - a;
    std::cout << "one     : "
	      << ' ' << std::setw(w) << a0.__result
	      << ' ' << std::setw(w) << e0
	      << ' ' << std::setw(w) << a0.__result - e0
	      << ' ' << std::setw(w) << std::abs(a0.__result - e0) / std::abs(e0)
	      << '\n';

    auto a1 = __gnu_cxx::integrate_tanh_sinh<Tp>(ex, a, b, abs_err, rel_err);
    auto e1 = (b * b - a * a) / Tp{2};
    std::cout << "ex      : "
	      << ' ' << std::setw(w) << a1.__result
	      << ' ' << std::setw(w) << e1
	      << ' ' << std::setw(w) << a1.__result - e1
	      << ' ' << std::setw(w) << std::abs(a1.__result - e1) / std::abs(e1)
	      << '\n';

    auto a2 = __gnu_cxx::integrate_tanh_sinh<Tp>(cos2, a, b, abs_err, rel_err);
    auto e2 = pi / Tp{2};
    std::cout << "cos2    : "
	      << ' ' << std::setw(w) << a2.__result
	      << ' ' << std::setw(w) << e2
	      << ' ' << std::setw(w) << a2.__result - e2
	      << ' ' << std::setw(w) << std::abs(a2.__result - e2) / std::abs(e2)
	      << '\n';

    auto a3 = __gnu_cxx::integrate_tanh_sinh<Tp>(sin2, a, b, abs_err, rel_err);
    auto e3 = pi / Tp{2};
    std::cout << "sin2    : "
	      << ' ' << std::setw(w) << a3.__result
	      << ' ' << std::setw(w) << e3
	      << ' ' << std::setw(w) << a3.__result - e3
	      << ' ' << std::setw(w) << std::abs(a3.__result - e3) / std::abs(e3)
	      << '\n';

    auto a4 = __gnu_cxx::integrate_tanh_sinh<Tp>(j1, a, b, abs_err, rel_err);
    auto e4 = std::cyl_bessel_j(Tp{0}, Tp{0}) - std::cyl_bessel_j(Tp{0}, pi);
    std::cout << "j1      : "
	      << ' ' << std::setw(w) << a4.__result
	      << ' ' << std::setw(w) << e4
	      << ' ' << std::setw(w) << a4.__result - e4
	      << ' ' << std::setw(w) << std::abs(a4.__result - e4) / std::abs(e4)
	      << '\n';

    a = Tp{0};
    b = Tp{10} * pi;
    auto a5 = __gnu_cxx::integrate_tanh_sinh<Tp>(foo, a, b, abs_err, rel_err);
    auto e5 = Tp{2} * (Tp{1} + b) * std::exp(-b / Tp{2})
	    - Tp{2} * (Tp{1} + a) * std::exp(-a / Tp{2});
    std::cout << "foo     : "
	      << ' ' << std::setw(w) << a5.__result
	      << ' ' << std::setw(w) << e5
	      << ' ' << std::setw(w) << a5.__result - e5
	      << ' ' << std::setw(w) << std::abs(a5.__result - e5) / std::abs(e5)
	      << '\n';

    auto a5n = __gnu_cxx::integrate_tanh_sinh<Tp>(foonum, a, b, abs_err, rel_err);
    auto e5n = b * (Tp{1} - b / Tp{2})
	     - a * (Tp{1} - a / Tp{2});
    std::cout << "foonum  : "
	      << ' ' << std::setw(w) << a5n.__result
	      << ' ' << std::setw(w) << e5n
	      << ' ' << std::setw(w) << a5n.__result - e5n
	      << ' ' << std::setw(w) << std::abs(a5n.__result - e5n) / std::abs(e5n)
	      << '\n';

    auto a6 = __gnu_cxx::integrate_tanh_sinh<Tp>(poly1, a, b, abs_err, rel_err);
    auto e6 = poly1.integral()(b) - poly1.integral()(a);
    std::cout << "poly1   : "
	      << ' ' << std::setw(w) << a6.__result
	      << ' ' << std::setw(w) << e6
	      << ' ' << std::setw(w) << a6.__result - e6
	      << ' ' << std::setw(w) << std::abs(a6.__result - e6) / std::abs(e6)
	      << '\n';

    a = Tp{0};
    b = pi;
    auto a7 = __gnu_cxx::integrate_tanh_sinh<Tp>(funk1, a, b, abs_err, rel_err);
    auto e7 = Tp{0};
    std::cout << "funk1   : "
	      << ' ' << std::setw(w) << a7.__result
	      << ' ' << std::setw(w) << e7
	      << ' ' << std::setw(w) << a7.__result - e7
	      << ' ' << std::setw(w) << std::abs(a7.__result - e7) / std::abs(e7)
	      << '\n';

    auto a7n = __gnu_cxx::integrate_tanh_sinh<Tp>(funk1num, a, b, abs_err, rel_err);
    auto e7n = Tp{0};
    std::cout << "funk1num: "
	      << ' ' << std::setw(w) << a7n.__result
	      << ' ' << std::setw(w) << e7n
	      << ' ' << std::setw(w) << a7n.__result - e7n
	      << ' ' << std::setw(w) << std::abs(a7n.__result - e7n) / std::abs(e7n)
	      << '\n';

    auto a8 = __gnu_cxx::integrate_tanh_sinh<Tp>(funk2, a, b, abs_err, rel_err);
    auto e8 = Tp{0};
    std::cout << "funk2   : "
	      << ' ' << std::setw(w) << a8.__result
	      << ' ' << std::setw(w) << e8
	      << ' ' << std::setw(w) << a8.__result - e8
	      << ' ' << std::setw(w) << std::abs(a8.__result - e8) / std::abs(e8)
	      << '\n';

    auto a8n = __gnu_cxx::integrate_tanh_sinh<Tp>(funk2num, a, b, abs_err, rel_err);
    auto e8n = Tp{2} * (b - a) - std::cos(b) + std::cos(a);
    std::cout << "funk2num: "
	      << ' ' << std::setw(w) << a8n.__result
	      << ' ' << std::setw(w) << e8n
	      << ' ' << std::setw(w) << a8n.__result - e8n
	      << ' ' << std::setw(w) << std::abs(a8n.__result - e8n) / std::abs(e8n)
	      << '\n';

    // sinh-sinh

    auto mu = Tp{5};
    auto gauss = [sqrt_pi, mu](Tp x)
                 -> Tp
                 { auto t = x - mu; return Tp{2} * std::exp(-t * t) / sqrt_pi; };
    auto ag = __gnu_cxx::integrate_sinh_sinh<Tp>(gauss, abs_err, rel_err);
    auto eg = Tp{2};
    std::cout << "gauss   : "
	      << ' ' << std::setw(w) << ag.__result
	      << ' ' << std::setw(w) << eg
	      << ' ' << std::setw(w) << ag.__result - eg
	      << ' ' << std::setw(w) << std::abs(ag.__result - eg) / std::abs(eg)
	      << '\n';

    // exp-sinh

    auto ag2 = __gnu_cxx::integrate_exp_sinh<Tp>(gauss, mu, abs_err, rel_err);
    auto eg2 = Tp{1};
    std::cout << "gauss   : "
	      << ' ' << std::setw(w) << ag2.__result
	      << ' ' << std::setw(w) << eg2
	      << ' ' << std::setw(w) << ag2.__result - eg2
	      << ' ' << std::setw(w) << std::abs(ag2.__result - eg2) / std::abs(eg2)
	      << '\n';
  }

int
main()
{
  test_tanh_sinh<double>();

  std::cout << "K(float)       = " << std::log(std::log(std::numeric_limits<float>::max())) << '\n';
  std::cout << "K(double)      = " << std::log(std::log(std::numeric_limits<double>::max())) << '\n';
  std::cout << "K(long double) = " << std::log(std::log(std::numeric_limits<long double>::max())) << '\n';
}
