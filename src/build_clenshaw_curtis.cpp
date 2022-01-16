
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <wrap_burkhardt.h>

#include <ext/quadrature_point.h>

  /**
   * Return a vector of zeros of the Chebyshev function of the second kind
   * of order @f$ n @f$, @f$ U_n(x) @f$.
   * The zeros are given by:
   * @f[
   *   x_k = \cos\left(\frac{k\pi}{n+1}\right), k \elem {1, ..., n}
   * @f]
   */
  template<typename Tp>
    std::vector<emsr::quadrature_point_t<Tp>>
    chebyshev_u_zeros(unsigned int n)
    {
      const auto s_pi = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
      std::vector<emsr::quadrature_point_t<Tp>> pt(n);
      for (unsigned int k = 1; k <= n; ++k)
	{
	  auto arg = Tp(k) / Tp(n + 1);
	  auto half = gnu_cxx::fp_is_equal<Tp>(arg, Tp{0.5L});
	  auto z = (half ? Tp{0} : gnu_cxx::cos_pi(arg));
	  auto w = s_pi * (Tp{1} - z * z) / Tp(n + 1);
	  pt[k - 1].point = z;
	  pt[k - 1].weight = w;
	}
      return pt;
    }

/**
 * 
 *
 * @f[
 *    w_k = frac{c_k}{n}\left[1-\sum_{j=1}^{[n/2]}\frac{b_k}{4j^2-1}
 *            \cos\left(\frac{2jk\pi}{n}\right)\right]
 *    \mbox{   } k = 0, 1, ..., n
 * @f]
 * 
 * @see Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules
 */
template<typename Tp>
  std::vector<emsr::quadrature_point_t<Tp>>
  build_clenshaw_curtis_sum(std::size_t __n)
  {
    std::vector<emsr::quadrature_point_t<Tp>> out(n + 1);
    if (n == 1)
      {
	out[0].point = Tp{0};
	out[0].weight = Tp{2};
	return out;
      }
    else
      {
	const auto s_pi = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
	auto uz = chebyshev_u_zeros<Tp>(n - 1);
	out[0].point = Tp{+1};
	out[0].weight = Tp{1} / (n * n - 1 + n % 2);
	for (auto k = 1u; k <= uz.size(); ++k)
	  {
	    out[k].point = uz[k - 1].point;

	    auto sum = Tp{0};
	    for (auto j = 1u; j <= n / 2; ++j)
	      {
		auto b = Tp(j == n / 2 ? 1 : 2);
		sum += b * std::cos(2 * j * k * s_pi / n)
		       / Tp(4 * j * j - 1);
	      }
	    auto w = Tp{2} * (Tp{1} - sum) / Tp(n);
	    out[k].weight = w;
	  }
	out[n].point = Tp{-1};
	out[n].weight = out[0].weight;
	return out;
      }
  }

/**
 * 
 *
 * @f[
 *    w_k = frac{c_k}{n}\left[1-\sum_{j=1}^{[n/2]}\frac{b_k}{4j^2-1}
 *            \cos\left(\frac{2jk\pi}{n}\right)\right]
 *    \mbox{   } k = 0, 1, ..., n
 * @f]
 * 
 * @see Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules
 */
template<typename Tp>
  std::vector<emsr::quadrature_point_t<Tp>>
  build_clenshaw_curtis_fft(std::size_t __n)
  {
    std::vector<emsr::quadrature_point_t<Tp>> __out(__n + 1);
    return __out;
  }

/**
 * 
 *
 * @f[
 *    w_k = frac{2}{n}\left[1-2\sum_{j=1}^{[n/2]}\frac{1}{4j^2-1}
 *            \cos\left(j\frac{(2k+1)\pi}{n}\right)\right]
 *    \mbox{   } k = 0, 1, ..., n-1
 * @f]
 */
template<typename Tp>
  std::vector<emsr::quadrature_point_t<Tp>>
  build_fejer_1_sum(std::size_t n)
  {
    std::vector<emsr::quadrature_point_t<Tp>> out(n + 1);
    if (n == 1)
      {
	out[0].point = Tp{0};
	out[0].weight = Tp{2};
	return out;
      }
    else
      {
	const auto s_pi = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
	auto uz = chebyshev_u_zeros<Tp>(n - 1);
	out[0].point = Tp{+1};
	out[0].weight = Tp{1} / (n * n - 1 + n % 2);
	for (auto k = 1u; k <= uz.size(); ++k)
	  {
	    out[k].point = uz[k - 1].point;

	    auto sum = Tp{0};
	    for (auto j = 1u; j <= n / 2; ++j)
	      sum += Tp{2} * std::cos(j * (2 * k + 1) * s_pi / n)
		     / Tp(4 * j * j - 1);
	    auto w = Tp{2} * (Tp{1} - sum) / Tp(n);
	    out[k].weight = w;
	  }
	out[n].point = Tp{-1};
	out[n].weight = out[0].weight;
	return out;
      }
  }

/**
 * 
 *
 * @f[
 *    w_k = frac{2}{n}\left[1-2\sum_{j=1}^{[n/2]}\frac{1}{4j^2-1}
 *            \cos\left(j\frac{(2k+1)\pi}{n}\right)\right]
 *    \mbox{   } k = 0, 1, ..., n-1
 * @f]
 */
template<typename Tp>
  std::vector<emsr::quadrature_point_t<Tp>>
  build_fejer_1_fft(std::size_t n)
  {
    const auto s_pi = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
    std::vector<emsr::quadrature_point_t<Tp>> out(n);
    std::vector<std::complex<Tp>> v(n + 1);
    for (auto k = 0u; k < n / 2; ++k)
      v[n - k] = v[k] = std::polar(Tp{2}, k * s_pi / Tp{2})
      		     / Tp(1 - 4 * k * k);
    v[n / 2] = Tp(n - 3) / Tp(2 * (n / 2) - 1) - Tp(1);
    emsr::inv_fast_fourier_transform(v);
std::cout << '\n';
for (auto x : v)
std::cout << x << '\n';
    return out;
  }

/**
 * 
 *
 * @f[
 *    w_k = frac{4}{n}\sin\left(\frac{k\pi}{n}\right)
 *        \sum_{j=1}^{[n/2]}\frac{1}{2j-1}\sin left((2j-1)\frac{k\pi}{n}\right)
 *    \mbox{   } k = 0, 1, ..., n
 * @f]
 */
template<typename Tp>
  std::vector<emsr::__quadrature_point_t<Tp>>
  build_fejer_2_sum(std::size_t n)
  {
    std::vector<emsr::quadrature_point_t<Tp>> out(n + 1);
    if (n == 0)
      {
	out[0].point = Tp{0};
	out[0].weight = Tp{2};
	return out;
      }
    else if (n == 1)
      {
	out[0].point = Tp{-0.5L};
	out[0].weight = Tp{0};
	out[1].point = Tp{+0.5L};
	out[1].weight = Tp{0};
	return out;
      }
    else
      {
	const auto s_pi = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
	auto uz = chebyshev_u_zeros<Tp>(n - 1);
	out[0].point = Tp{+1};
	out[0].weight = Tp{0};
	for (auto k = 1u; k <= uz.size(); ++k)
	  {
	    out[k].point = uz[k - 1].point;

	    auto sum = Tp{0};
	    for (auto j = 1u; j <= n / 2; ++j)
	      sum += std::sin((2 * j - 1) * k * s_pi / n)
		     / Tp(2 * j - 1);
	    auto w = Tp{4} * std::sin(k * s_pi / n) * sum
		     / Tp(n);
	    out[k].weight = w;
	  }
	out[n].point = Tp{-1};
	out[n].weight = Tp{0};
	return out;
      }
  }

/**
 * 
 *
 * @f[
 *    w_k = frac{4}{n}\sin\left(j\frac{k\pi}{n}\right)
 *        \sum_{j=1}^{[n/2]}\frac{1}{2j-1}\sin left((2j-1)\frac{k\pi}{n}\right)
 *    \mbox{   } k = 0, 1, ..., n
 * @f]
 */
template<typename Tp>
  std::vector<emsr::__quadrature_point_t<Tp>>
  build_fejer_2_fft(std::size_t n)
  {
    std::vector<emsr::quadrature_point_t<Tp>> out(n + 1);
    std::vector<std::complex<Tp>> v(n + 1);
    for (auto k = 0u; k < n / 2; ++k)
      v[n - k] = v[k] = Tp{2} / Tp(1 - 4 * k * k);
    v[n / 2] = Tp(n - 3) / Tp(2 * (n / 2) - 1) - Tp(1);
    emsr::inv_fast_fourier_transform(v);
std::cout << '\n';
for (auto x : v)
std::cout << x << '\n';
    return out;
  }

int
main()
{
  std::cout.precision(std::numeric_limits<long double>::digits10);
  auto w = 8 + std::cout.precision();

  auto cc24b = burkhardt::clenshaw_curtis_rule(25);

  auto cc24s = build_clenshaw_curtis_sum<long double>(24);
  auto cc24f = build_clenshaw_curtis_fft<long double>(24);
  std::cout << "\nClenshaw-Curtis 24\n";
  int i = 0;
  for (const auto& cc : cc24s)
    {
      std::cout << std::setw(w) << cc.point << ' '
		<< std::setw(w) << cc.weight << ' '
		<< std::setw(w) << cc.point - cc24b[i].point << ' '
		<< std::setw(w) << cc.weight - cc24b[i].weight
		<< '\n';
      ++i;
    }

  auto cc48b = burkhardt::clenshaw_curtis_rule(49);

  std::cout << "\nClenshaw-Curtis 48\n";
  auto cc48 = build_clenshaw_curtis_sum<long double>(48);
  i = 0;
  for (const auto& cc : cc48)
    {
      std::cout << std::setw(w) << cc.point << ' '
		<< std::setw(w) << cc.weight << ' '
		<< std::setw(w) << cc.point - cc48b[i].point << ' '
		<< std::setw(w) << cc.weight - cc48b[i].weight
		<< '\n';
      ++i;
    }

  auto f1_24b = burkhardt::fejer_1_rule(25);

  auto f1_24s = build_fejer_1_sum<long double>(24);
  auto f1_24f = build_fejer_1_fft<long double>(24);

  std::cout << "\nFejer1 24\n";
  i = 0;
  for (const auto& f1 : f1_24s)
    {
      std::cout << std::setw(w) << f1.point << ' '
		<< std::setw(w) << f1.weight << ' '
		<< std::setw(w) << f1.point - f1_24b[i].point << ' '
		<< std::setw(w) << f1.weight - f1_24b[i].weight
		<< '\n';
      ++i;
    }

  auto f2_24b = burkhardt::fejer_2_rule(25);

  auto f2_24s = build_fejer_2_sum<long double>(24);
  auto f2_24f = build_fejer_2_fft<long double>(24);
  std::cout << "\nFejer2 24\n";
  i = 0;
  for (const auto& f2 : f2_24s)
    {
      std::cout << std::setw(w) << f2.point << ' '
		<< std::setw(w) << f2.weight << ' '
		<< std::setw(w) << f2.point - f2_24b[i].point << ' '
		<< std::setw(w) << f2.weight - f2_24b[i].weight
		<< '\n';
      ++i;
    }

  std::cout << "\n\nCQUAD Rules\n";
  for (const auto& n : {4, 8, 16, 32, 64})
    {
      std::cout << "\nClenshaw-Curtis " << n << "\n";
      for (const auto& cc : build_clenshaw_curtis_sum<long double>(n))
	{
	  std::cout << std::setw(w) << cc.point << ' '
		    << std::setw(w) << cc.weight << ' '
		    << '\n';
	}
    }
}
