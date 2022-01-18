#include <ext/integration.h>

#include <complex>

template <typename Tp>
std::complex<Tp>
fun(Tp x)
{
  //using namespace std::complex_literals;
  constexpr std::complex<Tp> j{0, 1};
  return 1.0 / (x*x - j * x);
}

// simple test
template <typename Tp>
void test()
{
    auto ans = emsr::integrate_singular(fun<Tp>, Tp{0.001}, Tp{1}, Tp{1e-8}, Tp{1e-8});
    std::cout << "integral = " << ans.result << std::endl;
}

int
main()
{
  test<double>();
}
