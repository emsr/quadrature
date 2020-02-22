
#include <cmath>
#include <iostream>
#include <iomanip>

#include <ext/gauss_hermite_integrate.h>

template<typename _Tp>
  struct my_f
  {
    _Tp a;
    _Tp b;
    _Tp c;
    _Tp d;
    _Tp e;

    _Tp
    operator()(_Tp x) const
    { return (((a * x + b) * x + c ) * x + d ) * x + e; }

    _Tp
    integ(_Tp x) const
    { return x * ((((a / 5 * x + b / 4) * x + c / 3) * x + d / 2) * x + e / 1); }
  };

template<typename _Tp>
  struct my_fc
  {
    _Tp a;
    _Tp b;
    _Tp c;

    _Tp
    operator()(_Tp x) const
    { return a * std::cos(b * x + c); }

    _Tp
    integ(_Tp x) const
    { return a * std::sin(b * x + c) / b; }
  };

template<typename _Tp>
  void
  test_hermite_quad()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout.flags(std::ios::showpoint);
    auto width = 8 + std::cout.precision();

    std::cout << '\n';

    my_f<_Tp> params{ 0.5, 0.0, 0.0, -1.0, 1.3 };
    my_fc<_Tp> cparams{ 1.0, 1.0, 0.0 };

    //pfoo2 = __gnu_cxx::gauss_hermite_prob_integrate<_Tp>(params, 3);
    //std::cout << "integral = "  << std::setw(width)<< pfoo2 << '\n';
    //std::cout << "delta = " << std::setw(width) << pfoo2 - _Tp{7.01855916896680140676414279747L} << '\n';
    //pfooc = __gnu_cxx::gauss_hermite_prob_integrate<_Tp>(cparams, 29);
    //std::cout << "integral = " << std::setw(width) << pfooc << '\n';
    //std::cout << "delta = "  << std::setw(width)<< pfooc - _Tp{1.52034690106628080561194014675L} << '\n';

    auto xfoo2 = __gnu_cxx::gauss_hermite_integrate<_Tp>(params, 3);
    std::cout << "integral = " << std::setw(width) << xfoo2 << '\n';
    std::cout << "delta    = " << std::setw(width) << xfoo2 - _Tp{2.96886020026673934572443053460L} << '\n';
    //std::cout << -params.integ(1) + params.integ(0) << '\n';

    auto xfooc = __gnu_cxx::gauss_hermite_integrate<_Tp>(cparams, 29);
    std::cout << "integral = " << std::setw(width) << xfooc << '\n';
    std::cout << "delta    = " << std::setw(width) << xfooc - _Tp{1.38038844704314297477341524673L} << '\n';
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
