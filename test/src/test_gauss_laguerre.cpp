
#include <cmath>
#include <iostream>
#include <iomanip>

#include <emsr/gauss_laguerre_integrate.h>
#include <emsr/quadrature_point.h>

template<typename Tp>
  struct my_f
  {
    Tp a;
    Tp b;
    Tp c;
    Tp d;
    Tp e;

    Tp
    operator()(Tp x) const
    { return (((a * x + b) * x + c) * x + d) * x + e; }

    Tp
    integ(Tp x) const
    { return x * ((((a / 5 * x + b / 4) * x + c / 3) * x + d / 2) * x + e / 1); }
  };

template<typename Tp>
  struct my_fc
  {
    Tp a;
    Tp b;
    Tp c;

    Tp
    operator()(Tp x) const
    { return a * std::cos(b * x + c); }

    Tp
    integ(Tp x) const
    { return a * std::sin(b * x + c) / b; }
  };

template<typename Tp>
  void
  test_laguerre_quad()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout.flags(std::ios::showpoint);
    auto width = 8 + std::cout.precision();

    std::cout << '\n';

    my_f<Tp> poly{ 0.5, 0.0, 0.0, -1.0, 1.3 };
    my_fc<Tp> cosine{ 1.0, 1.0, 0.0 };

    //pfoo2 = emsr::gauss_laguerre_prob_integrate<Tp>(poly, 3);
    //pfooc = emsr::gauss_laguerre_prob_integrate<Tp>(cosine, 29);

    auto xfoo2 = emsr::gauss_laguerre_integrate<Tp>(poly, 3, 0.0);
    std::cout << "integral = " << std::setw(width) << xfoo2 << '\n';
    std::cout << "delta    = " << std::setw(width) << xfoo2 - 0.0L << '\n';

    auto xfooc = emsr::gauss_laguerre_integrate<Tp>(cosine, 29, 0.0);
    std::cout << "integral = " << std::setw(width) << xfooc << '\n';
    std::cout << "delta    = " << std::setw(width) << xfooc - 0.0L << '\n';
  }

int
main()
{
  std::cout << "\n\ndouble\n";
  test_laguerre_quad<double>();

  std::cout << "\n\nlong double\n";
  test_laguerre_quad<long double>();
}
