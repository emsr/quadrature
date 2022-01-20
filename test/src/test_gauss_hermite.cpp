
#include <cmath>
#include <iostream>
#include <iomanip>

#include <emsr/gauss_hermite_integrate.h>

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
    { return (((a * x + b) * x + c ) * x + d ) * x + e; }

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
  test_hermite_quad()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout.flags(std::ios::showpoint);
    auto width = 8 + std::cout.precision();

    std::cout << '\n';

    my_f<Tp> params{ 0.5, 0.0, 0.0, -1.0, 1.3 };
    my_fc<Tp> cparams{ 1.0, 1.0, 0.0 };

    //pfoo2 = emsr::gauss_hermite_prob_integrate<Tp>(params, 3);
    //std::cout << "integral = "  << std::setw(width)<< pfoo2 << '\n';
    //std::cout << "delta = " << std::setw(width) << pfoo2 - Tp{7.01855916896680140676414279747L} << '\n';
    //pfooc = emsr::gauss_hermite_prob_integrate<Tp>(cparams, 29);
    //std::cout << "integral = " << std::setw(width) << pfooc << '\n';
    //std::cout << "delta = "  << std::setw(width)<< pfooc - Tp{1.52034690106628080561194014675L} << '\n';

    auto xfoo2 = emsr::gauss_hermite_integrate<Tp>(params, 3);
    std::cout << "integral = " << std::setw(width) << xfoo2 << '\n';
    std::cout << "delta    = " << std::setw(width) << xfoo2 - Tp{2.96886020026673934572443053460L} << '\n';
    //std::cout << -params.integ(1) + params.integ(0) << '\n';

    auto xfooc = emsr::gauss_hermite_integrate<Tp>(cparams, 29);
    std::cout << "integral = " << std::setw(width) << xfooc << '\n';
    std::cout << "delta    = " << std::setw(width) << xfooc - Tp{1.38038844704314297477341524673L} << '\n';
    //std::cout << -cparams.integ(1) + cparams.integ(0) << '\n';
  }

int
main()
{
  std::cout << "\n\ndouble\n";
  test_hermite_quad<double>();

  std::cout << "\n\nlong double\n";
  test_hermite_quad<long double>();
}
