
#include <complex>

#include <emsr/integration.h>

std::complex<double>
fun(double x)
{
  using namespace std::complex_literals;
  return 1.0 / (x*x - 1.0i * x);
}

int
main()
{
  const auto I = emsr::integrate_clenshaw_curtis(fun, 1.0, 4.0, 0.000001, 0.0000111, 12);
}
