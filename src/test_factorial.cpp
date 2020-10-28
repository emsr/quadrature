
#include <cmath>
#include <ext/sf_factorial.h>

// Factorial.
template<typename _Tp>
  struct testcase_factorial
  {
    _Tp f0;
    unsigned int n;
    _Tp f;
  };

const testcase_factorial<double>
data001[21] =
{
  { 1.0000000000000000, 0, 0.0 },
  { 1.0000000000000000, 1, 0.0 },
  { 2.0000000000000000, 2, 0.0 },
  { 6.0000000000000000, 3, 0.0 },
  { 24.000000000000000, 4, 0.0 },
  { 120.00000000000000, 5, 0.0 },
  { 720.00000000000000, 6, 0.0 },
  { 5040.0000000000000, 7, 0.0 },
  { 40320.000000000000, 8, 0.0 },
  { 362880.00000000000, 9, 0.0 },
  { 3628800.0000000000, 10, 0.0 },
  { 39916800.000000000, 11, 0.0 },
  { 479001600.00000000, 12, 0.0 },
  { 6227020800.0000000, 13, 0.0 },
  { 87178291200.000000, 14, 0.0 },
  { 1307674368000.0000, 15, 0.0 },
  { 20922789888000.000, 16, 0.0 },
  { 355687428096000.00, 17, 0.0 },
  { 6402373705728000.0, 18, 0.0 },
  { 1.2164510040883200e+17, 19, 0.0 },
  { 2.4329020081766400e+18, 20, 0.0 }
};
const double toler001 = 2.5000000000000020e-13;

template<typename Ret, unsigned int Num>
  int
  test(const testcase_factorial<Ret> (&data)[Num], Ret toler)
  {
    bool test __attribute__((unused)) = true;
    const Ret eps = std::numeric_limits<Ret>::epsilon();
    Ret max_abs_diff = Ret(-1);
    Ret max_abs_frac = Ret(-1);
    bool failure = false;
    unsigned int num_datum = Num;
    for (unsigned int i = 0; i < num_datum; ++i)
      {
	const Ret f = __gnu_cxx::factorial<Ret>(data[i].n);
	const bool failure_f = std::isnan(f);
	if (!failure && failure_f)
	  failure = true;
	if (!failure_f)
	  {
	    const Ret f0 = data[i].f0;
	    const Ret diff = f - f0;
	    if (std::abs(diff) > max_abs_diff)
	      max_abs_diff = std::abs(diff);
	    if (std::abs(f0) > Ret(10) * eps
	     && std::abs(f) > Ret(10) * eps)
	      {
		const Ret frac = diff / f0;
		if (std::abs(frac) > max_abs_frac)
		  max_abs_frac = std::abs(frac);
	      }
	  }
      }
    return failure || max_abs_frac >= toler;
  }

int
main()
{
  test(data001, toler001);
  return 0;
}
