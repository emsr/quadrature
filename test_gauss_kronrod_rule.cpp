/*
g++ -std=c++14 -Wall -Wextra -Iinclude -o test_gauss_kronrod_rule test_gauss_kronrod_rule.cpp
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "ext/gauss_kronrod_rule.tcc"

template<typename Tp>
  void
  test_gauss_kronrod_rule()
  {
    std::cout.precision(std::numeric_limits<Tp>::max_digits10);
    const auto w = 6 + std::cout.precision();

    const auto eps = 4 * std::numeric_limits<Tp>::epsilon();

    for (int n : {7, 10, 15, 20, 25, 30})
      {
	std::cout << '\n' << (2 * n + 1) << "-point Gauss-Kronrod rule\n";
	std::vector<Tp> x, wk, wg;
	__gnu_cxx::__build_gauss_kronrod(n, eps, x, wg, wk);
	std::cout << "x\n";
	for (int i = 0; i < n + 1; ++i)
	  std::cout << ' ' << std::setw(w) << x[i] << '\n';
	std::cout << "wg\n";
	for (int i = 0; i < (n + 1) / 2; ++i)
	  std::cout << ' ' << std::setw(w) << wg[i] << '\n';
	std::cout << "wk\n";
	for (int i = 0; i < n + 1; ++i)
	  std::cout << ' ' << std::setw(w) << wk[i] << '\n';
      }
  }

int
main()
{
  std::cout << "\n\nTesting float Gauss-Fronrod rule ...\n\n";
  test_gauss_kronrod_rule<float>();

  std::cout << "\n\nTesting double Gauss-Fronrod rule ...\n\n";
  test_gauss_kronrod_rule<double>();

  std::cout << "\n\nTesting long double Gauss-Fronrod rule ...\n\n";
  test_gauss_kronrod_rule<long double>();
}

