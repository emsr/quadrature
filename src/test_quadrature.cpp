/* quadrature/test_quadrature.cpp
 *
 * Copyright (C) 1996-2000, 2007 Brian Gough
 * Copyright (C) 2016-2020 Free Software Foundation, Inc.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <sstream>
#include <memory>
#include <cfenv>

#include <ext/integration.h>
#include "testcase.h"
#include "func_utils.h"
#include <ext/special_functions.h>

template<typename _Tp>
  inline constexpr _Tp
//  eps_ratio = _Tp{1};
  eps_ratio = std::numeric_limits<_Tp>::epsilon()
	    / std::numeric_limits<double>::epsilon();

template<typename _Tp>
  inline constexpr _Tp
  prec_fixed = _Tp{1.0e-12}; // 10^{digits/10}*epsilon()?

template<>
  inline constexpr float
  prec_fixed<float> = 1.0e-5F;

template<>
  inline constexpr long double
  prec_fixed<long double> = 1.0e-14L;

template<typename _Tp>
  struct fixed_test
  {_Tp a, b, r, alpha, beta, gamma, delta;};

/**
 * 
 */
template<typename _Tp>
  struct quadrature_test
  {
    unsigned int num_tests = 0;
    unsigned int num_passed = 0;
    unsigned int num_failed = 0;
    bool verbose = false;

    enum STATUS
    {
      SUCCESS,
      FAILURE
    };

    quadrature_test()
    : num_tests{},
      num_passed{},
      num_failed{},
      verbose{true}
    { }

    ~quadrature_test()
    {
      tot_num_tests += this->num_tests;
      tot_num_passed += this->num_passed;
      tot_num_failed += this->num_failed;
    }

    void test_update(int status);
    void test_update(int status, const char* test_desc);
    void test_relative(_Tp result, _Tp expected, _Tp rel_error, const char* test_desc);
    void test_absolute(_Tp result, _Tp expected, _Tp abs_error, const char* test_desc);
    void test_factor(_Tp result, _Tp expected, _Tp factor, const char* test_desc);
    void test_integer(int result, int expected, const char* test_desc);

    static unsigned int tot_num_tests;
    static unsigned int tot_num_passed;
    static unsigned int tot_num_failed;
    static bool tot_verbose;

    static int test_summary();
  };

template<typename _Tp>
  unsigned int
  quadrature_test<_Tp>::tot_num_tests = 0;

template<typename _Tp>
  unsigned int
  quadrature_test<_Tp>::tot_num_passed = 0;

template<typename _Tp>
  unsigned int
  quadrature_test<_Tp>::tot_num_failed = 0;

template<typename _Tp>
  bool
  quadrature_test<_Tp>::tot_verbose = true;


template<typename _Tp>
  void
  quadrature_test<_Tp>::test_update(int status)
  {
    ++this->num_tests;

    if (status == 0)
      ++this->num_passed;
    else
      ++this->num_failed;
  }

template<typename _Tp>
  void
  quadrature_test<_Tp>::test_update(int status, const char* test_desc)
  {
    this->test_update(status);

    if (status || this->verbose)
      {
	std::cout << (status ? "FAIL: " : "PASS: ");

	std::cout << test_desc;

	if (status && !this->verbose)
	  std::cout << " [" << this->num_tests << ']';

	std::cout << '\n';
	std::cout << std::flush;
      }
  }

template<typename _Tp>
  void
  quadrature_test<_Tp>::test_relative(_Tp result, _Tp expected, _Tp rel_error,
				      const char* test_desc)
  {
    // Check for NaN vs inf vs number.
    int status = 0;
    if (std::isnan(result) || std::isnan(expected))
      status = std::isnan(result) || std::isnan(expected);
    else if (std::isinf(result) || std::isinf(expected))
      status = std::isinf(result) ^ std::isinf(expected);
    else if ((expected > _Tp{0} && expected < std::numeric_limits<_Tp>::min())
	  || (expected < _Tp{0} && expected > -(std::numeric_limits<_Tp>::min())))
      status = -1;
    else if (expected != _Tp{0})
      status = std::abs(result - expected) > rel_error * std::abs(expected);
    else
      status = std::abs(result) > rel_error;

    this->test_update(status);

    if (status || this->verbose)
      {
	std::cout << (status ? "FAIL: " : "PASS: ");

	std::cout << test_desc;

	if (status == 0)
	  std::cout << " (" << result << " observed vs " << expected << " expected)";
	else
	  {
	    const auto prec_old = std::cout.precision();
	    const auto prec_new = std::numeric_limits<_Tp>::max_digits10;
	    std::cout << std::setprecision(prec_new);
	    std::cout << " (" << result << " observed vs " << expected << " expected)";
	    std::cout << std::setprecision(prec_old);
	  }

	if (status == -1)
	  std::cout << " [test uses subnormal value]";

	if (status && !this->verbose)
	  std::cout << " [" << this->num_tests << ']';

	std::cout << '\n';
	std::cout << std::flush;
      }
  }

template<typename _Tp>
  void
  quadrature_test<_Tp>::test_absolute(_Tp result, _Tp expected, _Tp abs_error,
				      const char* test_desc)
  {
    // Check for NaN vs inf vs number.
    int status = 0;
    if (std::isnan(result) || std::isnan(expected))
      status = std::isnan(result) != std::isnan(expected);
    else if (std::isinf(result) || std::isinf(expected))
      status = std::isinf(result) != std::isinf(expected);
    else if ((expected > 0 && expected < std::numeric_limits<_Tp>::min())
	  || (expected < 0 && expected > -(std::numeric_limits<_Tp>::min())))
      status = -1;
    else
      status = std::abs(result - expected) > abs_error;

    this->test_update(status);

    if (status || this->verbose)
      {
	std::cout << (status ? "FAIL: " : "PASS: ");

	std::cout << test_desc;

	if (status == 0)
	  std::cout << " (" << result << " observed vs " << expected << " expected)";
	else
	  {
	    const auto prec_old = std::cout.precision();
	    const auto prec_new = std::numeric_limits<_Tp>::max_digits10;
	    std::cout << std::setprecision(prec_new);
	    std::cout << " (" << result << " observed vs " << expected << " expected)";
	    std::cout << std::setprecision(prec_old);
	  }

	if (status == -1)
	  std::cout << " [test uses subnormal value]";

	if (status && !this->verbose)
	  std::cout << " [" << this->num_tests << ']';

	std::cout << '\n';
	std::cout << std::flush;
      }
  }


template<typename _Tp>
  void
  quadrature_test<_Tp>::test_factor(_Tp result, _Tp expected, _Tp factor,
				    const char* test_desc)
  {
    int status = 0;
    if ((expected > 0 && expected < std::numeric_limits<_Tp>::min())
	|| (expected < 0 && expected > -(std::numeric_limits<_Tp>::min())))
      status = -1;
    else if (result == expected)
      status = 0;
    else if (expected == _Tp{0})
      status = (result > expected || result < expected);
    else
      {
	_Tp u = result / expected;
	status = (u > factor || u < _Tp{1} / factor);
      }

    this->test_update (status);

    if (status || this->verbose)
      {
	std::cout << (status ? "FAIL: " : "PASS: ");

	std::cout << test_desc;

	if (status == 0)
	  std::cout << " (" << result << " observed vs " << expected << " expected)";
	else
	  {
	    const auto prec_old = std::cout.precision();
	    const auto prec_new = std::numeric_limits<_Tp>::max_digits10;
	    std::cout << std::setprecision(prec_new);
	    std::cout << " (" << result << " observed vs " << expected << " expected)";
	    std::cout << std::setprecision(prec_old);
	  }

	if (status == -1)
	  std::cout << " [test uses subnormal value]";

	if (status && !this->verbose)
	  std::cout << " [" << this->num_tests << ']';

	std::cout << '\n';
	std::cout << std::flush;
      }
  }

template<typename _Tp>
  void
  quadrature_test<_Tp>::test_integer(int result, int expected,
				     const char* test_desc)
  {
    int status = (result != expected);

    this->test_update(status);

    if (status || this->verbose)
      {
	std::cout << (status ? "FAIL: " : "PASS: ");

	std::cout << test_desc;

	std::cout << " (" << result << " observed vs " << expected << " expected)";

	if (status && !this->verbose)
	  std::cout << " [" << this->num_tests << ']';

	std::cout << '\n';
	std::cout << std::flush;
      }
  }

template<typename _Tp>
  int
  quadrature_test<_Tp>::test_summary()
  {
    if (tot_verbose)
      std::cout << tot_num_tests << " tests, passed "
		<< tot_num_passed << ", failed "
		<< tot_num_failed << ".\n";

    if (tot_num_failed != 0)
      return FAILURE;

    if (tot_num_tests !=  tot_num_passed +  tot_num_failed)
      {
	if (tot_verbose)
	  std::cout << "TEST RESULTS DO NOT ADD UP " << tot_num_tests
		    << " != " << tot_num_passed << " + " << tot_num_failed << "\n";
	return FAILURE;
      }

    if (tot_num_passed ==  tot_num_tests)
      {
	if (! tot_verbose)
	  std::cout << "Completed [" << tot_num_passed << '/' << tot_num_tests << "]\n";

	return SUCCESS;
      }

    return FAILURE;
  }


/**
 * 
 * You want to hand in a rule or build one.
 * Have the params be variadic?
 * Expand the pack in the rule ctor.
 * Get sizeof...(Params) and loop for output?
 */
template<typename _Tp, typename _RuleTp, typename _FuncTp, typename... _Params>
  int
  test_quadrature_rule(_FuncTp f, _Tp a, _Tp b,
			_Tp tol, _Tp exact, const char* desc,
			_RuleTp /*rule*/, size_t n, _Params... params)
  {
    int status = 0;
    _RuleTp quad_rule(n, params...);

    std::ostringstream buf;
    buf << desc;
    buf << " a=" << a << " b=" << b;
    buf << " n=" << n;
    if constexpr (sizeof...(_Params) == 0)
      {
	// I think (params , ...) should expand to an empty list if sizeof...(_Params) == 0?
	// In any case, skip this problem.
      }
    else if constexpr (sizeof...(_Params) < 5)
      {
	_Tp param[sizeof...(_Params)]{params...};
	if constexpr (sizeof...(_Params) > 0)
	  buf << " alpha=" << param[0];
	if constexpr (sizeof...(_Params) > 1)
	  buf << " beta=" << param[1];
	if constexpr (sizeof...(_Params) > 2)
	  buf << " gamma=" << param[2];
	if constexpr (sizeof...(_Params) > 3)
	  buf << " delta=" << param[3];
      }
    else
      {
	_Tp param[sizeof...(_Params)]{params...};
	buf << " params=";
	const auto n = sizeof...(_Params);
	for (int p = 0; p < n; ++p)
	  buf << param[0] << (p < n-1 ? ", " : "");
      }

    auto result = quad_rule(f, a, b);
    quadrature_test<_Tp> qtest;
    qtest.test_relative(result, exact, tol, buf.str().c_str());

    return status;
  }


/**
 *
 */
template<typename _Tp>
  void
  belch(const __gnu_cxx::__integration_error<_Tp>& iex)
  {
    std::cout << "ERROR: " << iex.what()
	      << "       status = " << iex.error_code()
	      << "       result = " << iex.result()
	      << "       abserr = " << iex.abserr()
	      << std::endl;
    std::cerr << "ERROR: " << iex.what()
	      << "       status = " << iex.error_code()
	      << "       result = " << iex.result()
	      << "       abserr = " << iex.abserr()
	      << '\n';
  }


/**
 *
 */
void
belch(const std::exception& ex)
{
  std::cout << "ERROR: " << ex.what() << std::endl;
  std::cerr << "ERROR: " << ex.what() << '\n';
}


/**
 *
 */
template<typename _Tp>
  struct
  test_ival
  {
    _Tp a, b, r, e;
  };


/**
 * Main test function.
 */
template<typename _Tp>
int
test_quadrature()
{
  //feenableexcept(FE_OVERFLOW | FE_UNDERFLOW);
  //feenableexcept(FE_OVERFLOW);

  constexpr bool is_double = std::is_same_v<std::decay_t<_Tp>, double>;

  const auto _S_pi_2 = _S_pi<_Tp> / _Tp{2};
  const auto fpeps = _Tp{10} * std::numeric_limits<_Tp>::epsilon();//_Tp{1.0e-15};

  // Test the basic Gauss-Kronrod rules with a smooth positive function.
  const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 15 with a smooth positive function..." << std::endl;

      const auto exp_result = _Tp{7.716049357767090777e-02L};
      const auto exp_abserr = _Tp{2.990224871000550874e-06L};
      const auto exp_resabs = _Tp{7.716049357767090777e-02L};
      const auto exp_resasc = _Tp{4.434273814139995384e-02L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{2.6L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_cxx::Kronrod_15);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk15(f1) smooth result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk15(f1) smooth abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk15(f1) smooth resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk15(f1) smooth resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_cxx::Kronrod_15);

	qtest.test_relative(out.__result, -exp_result, fpeps, "qk15(f1) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk15(f1) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk15(f1) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk15(f1) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 21 with a smooth positive function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 21 with a smooth positive function..." << std::endl;

      const auto exp_result = _Tp{7.716049379303084599e-02L};
      const auto exp_abserr = _Tp{9.424302194248481445e-08L};
      const auto exp_resabs = _Tp{7.716049379303084599e-02L};
      const auto exp_resasc = _Tp{4.434311425038358484e-02L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{2.6L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_cxx::Kronrod_21);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk21(f1) smooth result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk21(f1) smooth abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk21(f1) smooth resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk21(f1) smooth resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_cxx::Kronrod_21);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk21(f1) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk21(f1) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk21(f1) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk21(f1) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 31 with a smooth positive function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 31 with a smooth positive function..." << std::endl;

      const auto exp_result = _Tp{7.716049382494900855e-02L};
      const auto exp_abserr = _Tp{1.713503193600029893e-09L};
      const auto exp_resabs = _Tp{7.716049382494900855e-02L};
      const auto exp_resasc = _Tp{4.427995051868838933e-02L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{2.6L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_cxx::Kronrod_31);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk31(f1) smooth result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk31(f1) smooth abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk31(f1) smooth resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk31(f1) smooth resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_cxx::Kronrod_31);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk31(f1) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk31(f1) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk31(f1) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk31(f1) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 41 with a smooth positive function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 41 with a smooth positive function..." << std::endl;

      const auto exp_result = _Tp{7.716049382681375302e-02L};
      const auto exp_abserr = _Tp{9.576386660975511224e-11L};
      const auto exp_resabs = _Tp{7.716049382681375302e-02L};
      const auto exp_resasc = _Tp{4.421521169637691873e-02L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{2.6L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_cxx::Kronrod_41);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk41(f1) smooth result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk41(f1) smooth abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk41(f1) smooth resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk41(f1) smooth resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_cxx::Kronrod_41);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk41(f1) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk41(f1) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk41(f1) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk41(f1) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 51 with a smooth positive function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 51 with a smooth positive function..." << std::endl;

      const auto exp_result = _Tp{7.716049382708510540e-02L};
      const auto exp_abserr = _Tp{1.002079980317363772e-11L};
      const auto exp_resabs = _Tp{7.716049382708510540e-02L};
      const auto exp_resasc = _Tp{4.416474291216854892e-02L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{2.6L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_cxx::Kronrod_51);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk51(f1) smooth result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qk51(f1) smooth abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk51(f1) smooth resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk51(f1) smooth resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_cxx::Kronrod_51);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk51(f1) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qk51(f1) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk51(f1) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk51(f1) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 61 with a smooth positive function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 61 with a smooth positive function..." << std::endl;

      const auto exp_result = _Tp{7.716049382713800753e-02L};
      const auto exp_abserr = _Tp{1.566060362296155616e-12L};
      const auto exp_resabs = _Tp{7.716049382713800753e-02L};
      const auto exp_resasc = _Tp{4.419287685934316506e-02L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{2.6L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_cxx::Kronrod_61);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk61(f1) smooth result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qk61(f1) smooth abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk61(f1) smooth resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk61(f1) smooth resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_cxx::Kronrod_61);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk61(f1) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qk61(f1) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk61(f1) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk61(f1) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Now test the basic rules with a positive function that has a
  // singularity. This should give large values of abserr which would
  // find discrepancies in the abserr calculation.

  // Test Gauss-Kronrod 15 with a singular positive function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 15 with a singular positive function..." << std::endl;

      const auto exp_result = _Tp{1.555688196612745777e+01L};
      const auto exp_abserr = _Tp{2.350164577239293706e+01L};
      const auto exp_resabs = _Tp{1.555688196612745777e+01L};
      const auto exp_resasc = _Tp{2.350164577239293706e+01L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{-0.9L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_cxx::Kronrod_15);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk15(f1) singular result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk15(f1) singular abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk15(f1) singular resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk15(f1) singular resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_cxx::Kronrod_15);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk15(f1) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk15(f1) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk15(f1) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk15(f1) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 21 with a singular positive function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 21 with a singular positive function..." << std::endl;

      const auto exp_result = _Tp{1.799045317938126232e+01L};
      const auto exp_abserr = _Tp{2.782360287710622515e+01L};
      const auto exp_resabs = _Tp{1.799045317938126232e+01L};
      const auto exp_resasc = _Tp{2.782360287710622515e+01L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{-0.9L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_cxx::Kronrod_21);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk21(f1) singular result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk21(f1) singular abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk21(f1) singular resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk21(f1) singular resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_cxx::Kronrod_21);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk21(f1) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk21(f1) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk21(f1) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk21(f1) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 31 with a singular positive function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 31 with a singular positive function..." << std::endl;

      const auto exp_result = _Tp{2.081873305159121657e+01L};
      const auto exp_abserr = _Tp{3.296500137482590276e+01L};
      const auto exp_resabs = _Tp{2.081873305159121301e+01L};
      const auto exp_resasc = _Tp{3.296500137482590276e+01L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{-0.9L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_cxx::Kronrod_31);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk31(f1) singular result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk31(f1) singular abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk31(f1) singular resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk31(f1) singular resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_cxx::Kronrod_31);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk31(f1) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk31(f1) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk31(f1) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk31(f1) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 41 with a singular positive function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 41 with a singular positive function..." << std::endl;

      const auto exp_result = _Tp{2.288677623903126701e+01L};
      const auto exp_abserr = _Tp{3.671538820274916048e+01L};
      const auto exp_resabs = _Tp{2.288677623903126701e+01L};
      const auto exp_resasc = _Tp{3.671538820274916048e+01L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{-0.9L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_cxx::Kronrod_41);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk41(f1) singular result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk41(f1) singular abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk41(f1) singular resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk41(f1) singular resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_cxx::Kronrod_41);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk41(f1) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk41(f1) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk41(f1) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk41(f1) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 51 with a singular positive function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 51 with a singular positive function..." << std::endl;

      const auto exp_result = _Tp{2.449953612016972215e+01L};
      const auto exp_abserr = _Tp{3.967771249391228849e+01L};
      const auto exp_resabs = _Tp{2.449953612016972215e+01L};
      const auto exp_resasc = _Tp{3.967771249391228849e+01L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{-0.9L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_cxx::Kronrod_51);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk51(f1) singular result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk51(f1) singular abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk51(f1) singular resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk51(f1) singular resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_cxx::Kronrod_51);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk51(f1) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk51(f1) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk51(f1) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk51(f1) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 61 with a singular positive function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 61 with a singular positive function..." << std::endl;

      const auto exp_result = _Tp{2.583030240976628988e+01L};
      const auto exp_abserr = _Tp{4.213750493076978643e+01L};
      const auto exp_resabs = _Tp{2.583030240976628988e+01L};
      const auto exp_resasc = _Tp{4.213750493076978643e+01L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{-0.9L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0}, _Tp{1}, __gnu_cxx::Kronrod_61);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk61(f1) singular result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk61(f1) singular abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk61(f1) singular resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk61(f1) singular resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{1}, _Tp{0}, __gnu_cxx::Kronrod_61);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk61(f1) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk61(f1) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk61(f1) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk61(f1) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test the basic Gauss-Kronrod rules with a smooth oscillating
  // function, over an unsymmetric range. This should find any
  // discrepancies in the abscissae.

  // Test Gauss-Kronrod 15 with a smooth oscillating function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 15 with a smooth oscillating function..." << std::endl;

      const auto exp_result = _Tp{-7.238969575483799046e-01L};
      const auto exp_abserr = _Tp{8.760080200939757174e-06L};
      const auto exp_resabs = _Tp{1.165564172429140788e+00L};
      const auto exp_resasc = _Tp{9.334560307787327371e-01L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{1.3L};
      auto f = make_function<_Tp, decltype(f1<_Tp>)>(f3<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0.3L}, _Tp{2.71L}, __gnu_cxx::Kronrod_15);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk15(f3) oscill result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk15(f3) oscill abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk15(f3) oscill resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk15(f3) oscill resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{2.71L}, _Tp{0.3L}, __gnu_cxx::Kronrod_15);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk15(f3) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk15(f3) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk15(f3) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk15(f3) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 21 with a smooth oscillating function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 21 with a smooth oscillating function..." << std::endl;

      const auto exp_result = _Tp{-7.238969575482959717e-01L};
      const auto exp_abserr = _Tp{7.999213141433641888e-11L};
      const auto exp_resabs = _Tp{1.150829032708484023e+00L};
      const auto exp_resasc = _Tp{9.297591249133687619e-01L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{1.3L};
      auto f = make_function<_Tp>(f3<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0.3L}, _Tp{2.71L}, __gnu_cxx::Kronrod_21);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk21(f3) oscill result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qk21(f3) oscill abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk21(f3) oscill resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk21(f3) oscill resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{2.71L}, _Tp{0.3L}, __gnu_cxx::Kronrod_21);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk21(f3) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qk21(f3) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk21(f3) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk21(f3) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 31 with a smooth oscillating function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 31 with a smooth oscillating function..." << std::endl;

      const auto exp_result = _Tp{-7.238969575482959717e-01L};
      const auto exp_abserr = _Tp{1.285805464427459261e-14L};
      const auto exp_resabs = _Tp{1.158150602093290571e+00L};
      const auto exp_resasc = _Tp{9.277828092501518853e-01L};
      quadrature_test<_Tp> qtest;

      const auto alpha =_Tp{1.3L};
      auto f = make_function<_Tp>(f3<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0.3L}, _Tp{2.71L}, __gnu_cxx::Kronrod_31);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk31(f3) oscill result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk31(f3) oscill abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk31(f3) oscill resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk31(f3) oscill resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{2.71L}, _Tp{0.3L}, __gnu_cxx::Kronrod_31);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk31(f3) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk31(f3) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk31(f3) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk31(f3) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 41 with a smooth oscillating function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 41 with a smooth oscillating function..." << std::endl;

      const auto exp_result = _Tp{-7.238969575482959717e-01L};
      const auto exp_abserr = _Tp{1.286535726271015626e-14L};
      const auto exp_resabs = _Tp{1.158808363486595328e+00L};
      const auto exp_resasc = _Tp{9.264382258645686985e-01L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{1.3L};
      auto f = make_function<_Tp>(f3<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0.3L}, _Tp{2.71L}, __gnu_cxx::Kronrod_41);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk41(f3) oscill result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk41(f3) oscill abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk41(f3) oscill resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk41(f3) oscill resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{2.71L}, _Tp{0.3L}, __gnu_cxx::Kronrod_41);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk41(f3) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk41(f3) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk41(f3) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk41(f3) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 51 with a smooth oscillating function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 51 with a smooth oscillating function..." << std::endl;

      const auto exp_result = _Tp{-7.238969575482961938e-01L};
      const auto exp_abserr = _Tp{1.285290995039385778e-14L};
      const auto exp_resabs = _Tp{1.157687209264406381e+00L};
      const auto exp_resasc = _Tp{9.264666884071264263e-01L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{1.3L};
      auto f = make_function<_Tp>(f3<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0.3L}, _Tp{2.71L}, __gnu_cxx::Kronrod_51);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk51(f3) oscill result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk51(f3) oscill abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk51(f3) oscill resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk51(f3) oscill resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{2.71L}, _Tp{0.3L}, __gnu_cxx::Kronrod_51);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk51(f3) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk51(f3) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk51(f3) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk51(f3) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Gauss-Kronrod 61 with a smooth oscillating function.
  try
    {
      //std::cout << ">>>> Test Gauss-Kronrod 61 with a smooth oscillating function..." << std::endl;

      const auto exp_result = _Tp{-7.238969575482959717e-01L};
      const auto exp_abserr = _Tp{1.286438572027470736e-14L};
      const auto exp_resabs = _Tp{1.158720854723590099e+00L};
      const auto exp_resasc = _Tp{9.270469641771273972e-01L};
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{1.3L};
      auto f = make_function<_Tp>(f3<_Tp>, alpha);
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{0.3L}, _Tp{2.71L}, __gnu_cxx::Kronrod_61);
	qtest.test_relative(out.__result, exp_result, fpeps, "qk61(f3) oscill result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk61(f3) oscill abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk61(f3) oscill resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk61(f3) oscill resasc");
      }
      {
	auto out = __gnu_cxx::qk_integrate(f, _Tp{2.71L}, _Tp{0.3L}, __gnu_cxx::Kronrod_61);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qk61(f3) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qk61(f3) reverse abserr");
	qtest.test_relative(out.__resabs, exp_resabs, fpeps, "qk61(f3) reverse resabs");
	qtest.test_relative(out.__resasc, exp_resasc, fpeps, "qk61(f3) reverse resasc");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test the non-adaptive Gaussian integrator.
  try
    {
      //std::cout << ">>>> Test non-adaptive Gaussian integrator..." << std::endl;

      int status = __gnu_cxx::NO_ERROR;
      const auto exp_result = _Tp{7.716049379303083211e-02L};
      const auto exp_abserr = _Tp{9.424302199601294244e-08L};
      const int exp_neval  = 21;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{2.6L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);
      {
	auto out = __gnu_cxx::qng_integrate(fc, _Tp{0}, _Tp{1}, _Tp{1.0e-1}, _Tp{0});
	qtest.test_relative(out.__result, exp_result, fpeps, "qng(f1) smooth result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qng(f1) smooth abserr");
	qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) smooth neval");
	qtest.test_integer(status, exp_ier, "qng(f1) smooth status");
      }
      {
	fc.num_evals(0);
	auto out = __gnu_cxx::qng_integrate(fc, _Tp{1}, _Tp{0}, _Tp{1.0e-1}, _Tp{0});
	qtest.test_relative(out.__result, -exp_result, fpeps, "qng(f1) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qng(f1) reverse abserr");
	qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) reverse neval");
	qtest.test_integer(status, exp_ier, "qng(f1) reverse status");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test the non-adaptive Gaussian integrator.
  try
    {
      //std::cout << ">>>> Test non-adaptive Gaussian integrator..." << std::endl;

      quadrature_test<_Tp> qtest;

      int status = 0;
      const auto exp_result = _Tp{7.716049382706505200e-02L};
      const auto exp_abserr = _Tp{2.666893044866214501e-12L};
      const int exp_neval  = 43;
      const int exp_ier    = __gnu_cxx::NO_ERROR;

      const auto alpha = _Tp{2.6L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);
      {
	auto out = __gnu_cxx::qng_integrate(fc, _Tp{0}, _Tp{1}, _Tp{0}, _Tp{1.0e-9L * eps_ratio<_Tp>});
	qtest.test_relative(out.__result, exp_result, fpeps, "qng(f1) smooth 43pt result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qng(f1) smooth 43pt abserr");
	qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) smooth 43pt neval");
	qtest.test_integer(status, exp_ier, "qng(f1) smooth 43pt status");
      }
      {
	fc.num_evals(0);
	auto out = __gnu_cxx::qng_integrate(fc, _Tp{1}, _Tp{0}, _Tp{0}, _Tp{1.0e-9L * eps_ratio<_Tp>});
	qtest.test_relative(out.__result, -exp_result, fpeps, "qng(f1) reverse 43pt result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qng(f1) reverse 43pt abserr");
	qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) reverse 43pt neval");
	qtest.test_integer(status, exp_ier, "qng(f1) reverse 43pt status");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test the non-adaptive Gaussian integrator.
  try
    {
      //std::cout << ">>>> Test non-adaptive Gaussian integrator..." << std::endl;

      int status = 0;
      const auto exp_result = _Tp{-7.238969575482961938e-01L};
      const auto exp_abserr = _Tp{1.277676889520056369e-14L};
      const int exp_neval  = 43;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      quadrature_test<_Tp> qtest;

      const auto alpha = _Tp{1.3L};
      auto f = make_function<_Tp>(f3<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);
      {
	auto out = __gnu_cxx::qng_integrate(fc, _Tp{0.3L}, _Tp{2.71L}, _Tp{0}, prec_fixed<_Tp>);
	qtest.test_relative(out.__result, exp_result, fpeps, "qnq(f3) oscill result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qng(f3) oscill abserr");
	qtest.test_integer(fc.num_evals(), exp_neval, "qng(f3) oscill neval");
	qtest.test_integer(status, exp_ier, "qng(f3) oscill status");
      }
      {
	fc.num_evals(0);
	auto out = __gnu_cxx::qng_integrate(fc, _Tp{2.71L}, _Tp{0.3L}, _Tp{0}, prec_fixed<_Tp>);
	qtest.test_relative(out.__result, -exp_result, fpeps, "qnq(f3) reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qng(f3) reverse abserr");
	qtest.test_integer(fc.num_evals(), exp_neval, "qng(f3) reverse neval");
	qtest.test_integer(status, exp_ier, "qng(f3) reverse status");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test the non-adaptive Gaussian integrator.
  try
    {
      //std::cout << ">>>> Test non-adaptive Gaussian integrator..." << std::endl;

      quadrature_test<_Tp> qtest;

      int status = 0;
      const auto exp_result = _Tp{7.716049382716029525e-02L};
      const auto exp_abserr = _Tp{8.566535680046930668e-16L};
      const int exp_neval  = 87;
      const int exp_ier    = __gnu_cxx::NO_ERROR;

      const auto alpha = _Tp{2.6L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);
      {
	auto out = __gnu_cxx::qng_integrate(fc, _Tp{0}, _Tp{1}, _Tp{0}, _Tp{1.0e-13L * eps_ratio<_Tp>});
	qtest.test_relative(out.__result, exp_result, fpeps, "qng(f1) 87pt smooth result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qng(f1) 87pt smooth abserr");
	qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) 87pt smooth neval");
	qtest.test_integer(status, exp_ier, "qng(f1) 87pt smooth status");
      }
      {
	fc.num_evals(0);
	auto out = __gnu_cxx::qng_integrate(fc, _Tp{1}, _Tp{0}, _Tp{0}, _Tp{1.0e-13L * eps_ratio<_Tp>});
	qtest.test_relative(out.__result, -exp_result, fpeps, "qng(f1) 87pt reverse result");
	if (is_double)
	  qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qng(f1) 87pt reverse abserr");
	qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) 87pt reverse neval");
	qtest.test_integer(status, exp_ier, "qng(f1) 87pt reverse status");
      }
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test the non-adaptive Gaussian integrator.
  try
    {
      //std::cout << ">>>> Test non-adaptive Gaussian integrator..." << std::endl;

      quadrature_test<_Tp> qtest;

      int status = 0;
      const auto exp_result = _Tp{3.222948711817264211e+01L};
      const auto exp_abserr = _Tp{2.782360287710622870e+01L};
      const int exp_neval  = 87;
      const int exp_ier    = __gnu_cxx::TOLERANCE_ERROR;

      const auto alpha = _Tp{-0.9L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::gauss_kronrod_integral_t<_Tp, decltype(fc(_Tp{}))> out;
      try
	{
	  out = __gnu_cxx::qng_integrate(fc, _Tp{0}, _Tp{1}, _Tp{0}, _Tp{1.0e-3L});
	}
      catch (__gnu_cxx::__integration_error<_Tp>& iex)
	{
	  status = iex.error_code();
	  out.__result = iex.result();
	  out.__abserr = iex.abserr();
	}
      qtest.test_relative(out.__result, exp_result, fpeps, "qng(f1) sing beyond 87pt result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qng(f1) sing beyond 87pt abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) sing beyond 87pt neval");
      qtest.test_integer(status, exp_ier, "qng(f1) sing beyond 87pt status");

      fc.num_evals(0);
      try
	{
	  out = __gnu_cxx::qng_integrate(fc, _Tp{1}, _Tp{0}, _Tp{0}, _Tp{1.0e-3L});
	}
      catch (__gnu_cxx::__integration_error<_Tp>& iex)
	{
	  status = iex.error_code();
	  out.__result = iex.result();
	  out.__abserr = iex.abserr();
	}
      qtest.test_relative(out.__result, -exp_result, fpeps, "qng(f1) reverse beyond 87pt result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-7L}, "qng(f1) rev beyond 87pt abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qng(f1) rev beyond 87pt neval");
      qtest.test_integer(status, exp_ier, "qng(f1) rev beyond 87pt status");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test the adaptive Gaussian integrator.
  try
    {
      //std::cout << ">>>> Test adaptive Gaussian integrator..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{7.716049382715854665e-02L};
      const auto exp_abserr = _Tp{6.679384885865053037e-12L};
      const int exp_neval  = 165;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 6;

      constexpr std::size_t num_test = 6;
      test_ival<_Tp> test[num_test]
      {
	{0.0L,     0.03125L, 3.966769831709074375e-06L, 6.678528276336181873e-12L},
	{0.5L,     1.0L,     5.491842501998222409e-02L, 6.097169993333454062e-16L},
	{0.125L,   0.25L,    2.776531175604360531e-03L, 3.082568839745514608e-17L},
	{0.0625L,  0.125L,   3.280661030752063693e-04L, 3.642265412331439511e-18L},
	{0.25L,    0.5L,     1.909827770934243926e-02L, 2.120334764359736934e-16L},
	{0.03125L, 0.0625L,  3.522704932261797744e-05L, 3.910988124757650942e-19L},
      };

      const auto alpha = _Tp{2.6L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      auto out = __gnu_cxx::qag_integrate(w, fc, _Tp{0}, _Tp{1},
					  _Tp{0}, _Tp{1.0e-10L},
					  __gnu_cxx::gauss_kronrod_integral<_Tp>(__gnu_cxx::Kronrod_15));

      qtest.test_relative(out.__result, exp_result, fpeps, "qag(f1) smooth result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-6L}, "qag(f1) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f1) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qag(f1) smooth last");
      qtest.test_integer(status, exp_ier, "qag(f1) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, fpeps, "qag(f1) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, fpeps, "qag(f1) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, fpeps, "qag(f1) smooth integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-6L}, "qag(f1) smooth abs error");

      fc.num_evals(0);
      out = __gnu_cxx::qag_integrate(w, fc, _Tp{1}, _Tp{0},
				     _Tp{0}, _Tp{1.0e-10L},
				     __gnu_cxx::gauss_kronrod_integral<_Tp>(__gnu_cxx::Kronrod_15));

      qtest.test_relative(out.__result, -exp_result, fpeps, "qag(f1) reverse result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-6L}, "qag(f1) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f1) reverse neval");
      qtest.test_integer(w.size(), exp_last, "qag(f1) reverse last");
      qtest.test_integer(status, exp_ier, "qag(f1) reverse status");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  try
    {
      //std::cout << ">>>> Test adaptive Gaussian integrator with absolute error bound..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{7.716049382716050342e-02L};
      const auto exp_abserr = _Tp{2.227969521869139532e-15L};
      const int exp_neval  = 315;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 8;

      constexpr std::size_t num_test = 8;
      test_ival<_Tp> test[num_test]
      {
	{0.0L,       0.0078125L, 3.696942726831556522e-08L, 1.371316364034059572e-15L},
	{0.25L,      0.5L,       1.909827770934243579e-02L, 2.120334764359736441e-16L},
	{0.5L,       1.0L,       5.491842501998223103e-02L, 6.097169993333454062e-16L},
	{0.0625L,    0.125L,     3.280661030752062609e-04L, 3.642265412331439511e-18L},
	{0.03125L,   0.0625L,    3.522704932261797744e-05L, 3.910988124757650460e-19L},
	{0.015625L,  0.03125L,   3.579060884684503576e-06L, 3.973555800712018091e-20L},
	{0.125L,     0.25L,      2.776531175604360097e-03L, 3.082568839745514608e-17L},
	{0.0078125L, 0.015625L,  3.507395216921808047e-07L, 3.893990926286736620e-21L},
      };

      const auto alpha = _Tp{2.6L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      auto out = __gnu_cxx::qag_integrate(w, fc, _Tp{0}, _Tp{1},
					  _Tp{1.0e-14L}, _Tp{0},
					  __gnu_cxx::gauss_kronrod_integral<_Tp>(__gnu_cxx::Kronrod_21));

      qtest.test_relative(out.__result, exp_result, fpeps, "qag(f1, 21pt) smooth result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-6L}, "qag(f1, 21pt) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f1, 21pt) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qag(f1, 21pt) smooth last");
      qtest.test_integer(status, exp_ier, "qag(f1, 21pt) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, fpeps, "qag(f1, 21pt) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, fpeps, "qag(f1, 21pt) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, fpeps, "qag(f1, 21pt) smooth integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-6L}, "qag(f1, 21pt) smooth abs error");

      fc.num_evals(0);
      out = __gnu_cxx::qag_integrate(w, fc, _Tp{1}, _Tp{0},
				     _Tp{1.0e-14L}, _Tp{0},
				     __gnu_cxx::gauss_kronrod_integral<_Tp>(__gnu_cxx::Kronrod_21));

      qtest.test_relative(out.__result, -exp_result, fpeps, "qag(f1, 21pt) reverse result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-6L}, "qag(f1, 21pt) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f1, 21pt) reverse neval");
      qtest.test_integer(w.size(), exp_last, "qag(f1, 21pt) reverse last");
      qtest.test_integer(status, exp_ier, "qag(f1, 21pt) reverse status");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Adaptive integration of an oscillatory function which terminates because
  // of roundoff error, uses the 31-pt rule.
  try
    {
      //std::cout << ">>>> Test adaptive integration of an oscillatory function\n"
	//	   ">>>> which terminates because of roundoff error, uses the 31-pt rule..." << std::endl;

      int status = __gnu_cxx::NO_ERROR;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{-7.238969575482959717e-01L};
      const auto exp_abserr = _Tp{1.285805464427459261e-14L};
      const int exp_neval   = 31;
      const int exp_ier     = __gnu_cxx::ROUNDOFF_ERROR;
      const int exp_last    = 1;

      const auto alpha = _Tp{1.3L};
      auto f = make_function<_Tp>(f3<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      __gnu_cxx::adaptive_integral_t<_Tp, decltype(fc(_Tp{}))> out;
      try
	{
	  out = __gnu_cxx::qag_integrate(w, fc, _Tp{0.3L}, _Tp{2.71L},
					 _Tp{1.0e-14L}, _Tp{0},
					 __gnu_cxx::gauss_kronrod_integral<_Tp>(__gnu_cxx::Kronrod_31));
	}
      catch (__gnu_cxx::__integration_error<_Tp>& iex)
	{
	  out.__result = iex.result();
	  out.__abserr = iex.abserr();
	  status = iex.error_code();
	}

      qtest.test_relative(out.__result, exp_result, fpeps, "qag(f3, 31pt) oscill result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-6L}, "qag(f3, 31pt) oscill abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f3, 31pt) oscill neval");
      qtest.test_integer(w.size(), exp_last, "qag(f3, 31pt) oscill last");
      qtest.test_integer(status, exp_ier, "qag(f3, 31pt) oscill status");

      fc.num_evals(0);
      try
	{
	  out = __gnu_cxx::qag_integrate(w, fc, _Tp{2.71L}, _Tp{0.3L},
					 _Tp{1.0e-14L}, _Tp{0},
					 __gnu_cxx::gauss_kronrod_integral<_Tp>(__gnu_cxx::Kronrod_31));
	}
      catch (__gnu_cxx::__integration_error<_Tp>& iex)
	{
	  out.__result = iex.result();
	  out.__abserr = iex.abserr();
	  status = iex.error_code();
	}

      qtest.test_relative(out.__result, -exp_result, fpeps, "qag(f3, 31pt) reverse result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-6L}, "qag(f3, 31pt) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f3, 31pt) reverse neval");
      qtest.test_integer(w.size(), exp_last, "qag(f3, 31pt) reverse last");
      qtest.test_integer(status, exp_ier, "qag(f3, 31pt) reverse status");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Check the singularity detection (singularity at x=-0.1 in this example).
  try
    {
      //std::cout << ">>>> Test singularity detection (singularity at x=-0.1 in this example)..." << std::endl;

      int status = __gnu_cxx::NO_ERROR;
      quadrature_test<_Tp> qtest;

      const int exp_neval  = 5151;
      const int exp_ier    = __gnu_cxx::SINGULAR_ERROR;
      const int exp_last   = 51;

      const auto alpha = _Tp{2};
      auto f = make_function<_Tp>(f16<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      __gnu_cxx::adaptive_integral_t<_Tp, decltype(fc(_Tp{}))> out;
      try
	{
	  out = __gnu_cxx::qag_integrate(w, fc, _Tp{-1}, _Tp{1},
					 _Tp{1.0e-14L}, _Tp{0},
					 __gnu_cxx::gauss_kronrod_integral<_Tp>(__gnu_cxx::Kronrod_51));
	}
      catch (__gnu_cxx::__integration_error<_Tp>& iex)
	{
	  out.__result = iex.result();
	  out.__abserr = iex.abserr();
	  status = iex.error_code();
	}

      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f16, 51pt) sing neval");
      qtest.test_integer(w.size(), exp_last, "qag(f16, 51pt) sing last");
      qtest.test_integer(status, exp_ier, "qag(f16, 51pt) sing status");

      fc.num_evals(0);
      try
	{
	  out = __gnu_cxx::qag_integrate(w, fc, _Tp{1}, _Tp{-1},
					 _Tp{1.0e-14L}, _Tp{0},
					 __gnu_cxx::gauss_kronrod_integral<_Tp>(__gnu_cxx::Kronrod_51));
	}
      catch (__gnu_cxx::__integration_error<_Tp>& iex)
	{
	  out.__result = iex.result();
	  out.__abserr = iex.abserr();
	  status = iex.error_code();
	}

      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f16, 51pt) rev neval");
      qtest.test_integer(w.size(), exp_last, "qag(f16, 51pt) rev last");
      qtest.test_integer(status, exp_ier, "qag(f16, 51pt) rev status");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Check for hitting the iteration limit.
  try
    {
      //std::cout << ">>>> Test hitting the iteration limit..." << std::endl;

      int status = __gnu_cxx::NO_ERROR;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{9.565151449233894709L};
      const auto exp_abserr = _Tp{1.570369823891028460e+01L};
      const int exp_neval  = 305;
      const int exp_ier    = __gnu_cxx::MAX_ITER_ERROR;
      const int exp_last   = 3;

      constexpr std::size_t num_test = 3;
      test_ival<_Tp> test[num_test]
      {
	{-5.000000000000000000e-01L,  0.000000000000000000L,     9.460353469435913709L,     1.570369823891028460e+01L},
	{-1.000000000000000000L,     -5.000000000000000000e-01L, 1.388888888888888812e-02L, 1.541976423090495140e-16L},
	{ 0.000000000000000000L,      1.000000000000000000L,     9.090909090909091161e-02L, 1.009293658750142399e-15L},
      };

      const auto alpha = _Tp{1};
      auto f = make_function<_Tp>(f16<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(3);

      __gnu_cxx::adaptive_integral_t<_Tp, decltype(fc(_Tp{}))> out;
      try
	{
	  out = __gnu_cxx::qag_integrate(w, fc, _Tp{-1}, _Tp{1},
					 _Tp{1.0e-14L}, _Tp{0},
					 __gnu_cxx::gauss_kronrod_integral<_Tp>(__gnu_cxx::Kronrod_61));
	}
      catch (__gnu_cxx::__integration_error<_Tp>& iex)
	{
	  out.__result = iex.result();
	  out.__abserr = iex.abserr();
	  status = iex.error_code();
	}

      qtest.test_relative(out.__result, exp_result, fpeps, "qag(f16, 61pt) limit result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-6L}, "qag(f16, 61pt) limit abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f16, 61pt) limit neval");
      qtest.test_integer(w.size(), exp_last, "qag(f16, 61pt) limit last");
      qtest.test_integer(status, exp_ier, "qag(f16, 61pt) limit status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, fpeps, "qag(f16, 61pt) limit lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, fpeps, "qag(f16, 61pt) limit upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, fpeps, "qag(f16, 61pt) limit integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-6L}, "qag(f16, 61pt) limit abs error");

      fc.num_evals(0);
      try
	{
	  out = __gnu_cxx::qag_integrate(w, fc, _Tp{1}, _Tp{-1},
					 _Tp{1.0e-14L}, _Tp{0},
					 __gnu_cxx::gauss_kronrod_integral<_Tp>(__gnu_cxx::Kronrod_61));
	}
      catch (__gnu_cxx::__integration_error<_Tp>& iex)
	{
	  out.__result = iex.result();
	  out.__abserr = iex.abserr();
	  status = iex.error_code();
	}

      qtest.test_relative(out.__result, -exp_result, fpeps, "qag(f16, 61pt) reverse result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-6L}, "qag(f16, 61pt) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qag(f16, 61pt) reverse neval");
      qtest.test_integer(w.size(), exp_last, "qag(f16, 61pt) reverse last");
      qtest.test_integer(status, exp_ier, "qag(f16, 61pt) reverse status");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test the adaptive integrator with extrapolation.
  try
    {
      //std::cout << ">>>> Test the adaptive integrator with extrapolation..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{7.716049382715789440e-02L};
      const auto exp_abserr = _Tp{2.216394961010438404e-12L};
      const int exp_neval  = 189;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 5;

      constexpr std::size_t num_test = 5;
      test_ival<_Tp> test[num_test]
      {
	0.0L,    0.0625L, 3.919381915366914693e-05L, 2.215538742580964735e-12L,
	0.5L,    1.0L,    5.491842501998223103e-02L, 6.097169993333454062e-16L,
	0.25L,   0.5L,    1.909827770934243579e-02L, 2.120334764359736441e-16L,
	0.0625L, 0.125L,  3.280661030752062609e-04L, 3.642265412331439511e-18L,
	0.125L,  0.25L,   2.776531175604360097e-03L, 3.082568839745514608e-17L,
      };

      const auto alpha = _Tp{2.6L};
      auto f = make_function<_Tp>(f1<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-10};
      auto out = __gnu_cxx::qags_integrate(w, fc, _Tp{0}, _Tp{1}, epsabs, epsrel);

      qtest.test_relative(out.__result, exp_result, fpeps, "qags(f1) smooth result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-6L}, "qags(f1) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qags(f1) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qags(f1) smooth last");
      qtest.test_integer(status, exp_ier, "qags(f1) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, fpeps, "qags(f1) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, fpeps, "qags(f1) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, fpeps, "qags(f1) smooth integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-6L}, "qags(f1) smooth abs error");

      fc.num_evals(0);
      out = __gnu_cxx::qags_integrate(w, fc, _Tp{1}, _Tp{0}, epsabs, epsrel);

      qtest.test_relative(out.__result, -exp_result, fpeps, "qags(f1) reverse result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-6L}, "qags(f1) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qags(f1) reverse neval");
      qtest.test_integer(w.size(), exp_last, "qags(f1) reverse last");
      qtest.test_integer(status, exp_ier, "qags(f1) reverse status");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test f11 using an absolute error bound.
  try
    {
      //std::cout << ">>>> Test f11 using an absolute error bound..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{-5.908755278982136588e+03L};
      const auto exp_abserr = _Tp{1.299646281053874554e-10L};
      const int exp_neval  = 357;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 9;

      constexpr std::size_t num_test = 9;
      test_ival<_Tp> test[num_test]
      {
	{1.000000000000000000e+00L, 4.902343750000000000e+00L, -3.890977835520834649e+00L, 6.448276035006137169e-11L},
	{5.005000000000000000e+02L, 1.000000000000000000e+03L, -3.297343675805121620e+03L, 3.660786868980994028e-11L},
	{1.258750000000000000e+02L, 2.507500000000000000e+02L, -6.517404019686431411e+02L, 7.235772003440423011e-12L},
	{2.507500000000000000e+02L, 5.005000000000000000e+02L, -1.475904154146372775e+03L, 1.638582774073219226e-11L},
	{3.221875000000000000e+01L, 6.343750000000000000e+01L, -1.201692001973227519e+02L, 1.334146129098576244e-12L},
	{6.343750000000000000e+01L, 1.258750000000000000e+02L, -2.829354222635842007e+02L, 3.141214202790722909e-12L},
	{1.660937500000000000e+01L, 3.221875000000000000e+01L, -4.959999906099650246e+01L, 5.506706097890446534e-13L},
	{4.902343750000000000e+00L, 8.804687500000000000e+00L, -7.457032710459004399e+00L, 8.278969410534525339e-14L},
	{8.804687500000000000e+00L, 1.660937500000000000e+01L, -1.971441499411640308e+01L, 2.188739744348345039e-13L},
      };

      auto alpha = _Tp{2};
      auto f = make_function<_Tp>(f11<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      const auto epsabs = _Tp{1.0e-7};
      const auto epsrel = _Tp{0};
      auto out = __gnu_cxx::qags_integrate(w, fc, _Tp{1}, _Tp{1000}, epsabs, epsrel);

      qtest.test_relative(out.__result, exp_result, fpeps, "qags(f11) smooth result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-3L}, "qags(f11) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qags(f11) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qags(f11) smooth last");
      qtest.test_integer(status, exp_ier, "qags(f11) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, fpeps, "qags(f11) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, fpeps, "qags(f11) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, fpeps, "qags(f11) smooth integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-5L}, "qags(f11) smooth abs error");

      fc.num_evals(0);
      out = __gnu_cxx::qags_integrate(w, fc, _Tp{1000}, _Tp{1}, epsabs, epsrel);

      qtest.test_relative(out.__result, -exp_result, fpeps, "qags(f11) reverse result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-3L}, "qags(f11) reverse abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qags(f11) reverse neval");
      qtest.test_integer(w.size(), exp_last, "qags(f11) reverse last");
      qtest.test_integer(status, exp_ier, "qags(f11) reverse status");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test infinite range integral f455 using a relative error bound.
  try
    {
      //std::cout << ">>>> Test infinite range integral f455 using a relative error bound..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{-3.616892186127022568e-01L};
      const auto exp_abserr = _Tp{3.016716913328831851e-06L};
      const int exp_neval  = 285;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 10;

      constexpr std::size_t num_test = 10;
      test_ival<_Tp> test[num_test]
      {
	{9.687500000000000000e-01L, 1.000000000000000000e+00L,  1.429785306003466313e-03L, 2.161214992172538524e-04L},
	{2.500000000000000000e-01L, 5.000000000000000000e-01L, -1.229943369113085765e-02L, 5.720644840858777846e-14L},
	{6.250000000000000000e-02L, 1.250000000000000000e-01L, -8.653752279614615461e-02L, 9.607595030230581153e-16L},
	{1.250000000000000000e-01L, 2.500000000000000000e-01L, -4.980050133751051655e-02L, 3.147380432198176412e-14L},
	{5.000000000000000000e-01L, 7.500000000000000000e-01L,  2.995321156568048898e-03L, 3.325474514168701167e-17L},
	{8.750000000000000000e-01L, 9.375000000000000000e-01L,  1.736218164975512294e-03L, 1.927589382528252344e-17L},
	{7.500000000000000000e-01L, 8.750000000000000000e-01L,  2.785385934678596704e-03L, 3.092399597147240624e-17L},
	{9.375000000000000000e-01L, 9.687500000000000000e-01L,  1.041689192004495576e-03L, 1.156507325466566521e-17L},
	{3.125000000000000000e-02L, 6.250000000000000000e-02L, -8.398745675010892142e-02L, 9.324480826368044019e-16L},
	{0.000000000000000000e+00L, 3.125000000000000000e-02L, -1.390003415539725340e-01L, 2.395037249893453013e-02L},
      };

      auto f = make_function<_Tp>(f455<_Tp>);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-3};
      auto out = __gnu_cxx::qagiu_integrate(w, fc, _Tp{0}, epsabs, epsrel);

      qtest.test_relative(out.__result, exp_result, /*1.0e-14*/epsrel, "qagiu(f455) smooth result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qagiu(f455) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qagiu(f455) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qagiu(f455) smooth last");
      qtest.test_integer(status, exp_ier, "qagiu(f455) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, fpeps, "qagiu(f455) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, fpeps, "qagiu(f455) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, /*fpeps*/epsrel, "qagiu(f455) smooth integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-4L}, "qagiu(f455) smooth abs error");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test infinite range integral f15 using a relative error bound.
  try
    {
      //std::cout << ">>>> Test infinite range integral f15 using a relative error bound..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{6.553600000000024738e+04L};
      const auto exp_abserr = _Tp{7.121667111456009280e-04L};
      const int exp_neval  = 285;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 10;

      constexpr std::size_t num_test = 10;
      test_ival<_Tp> test[num_test]
      {
	{0.000000000000000000e+00L, 1.953125000000000000e-03L, 1.099297665754340292e+00L, 7.101865971621337814e-04L},
	{6.250000000000000000e-02L, 1.250000000000000000e-01L, 6.977679035845269482e+02L, 6.973493131275552509e-07L},
	{3.906250000000000000e-03L, 7.812500000000000000e-03L, 1.498314766425578091e+04L, 3.408933028357320364e-07L},
	{2.500000000000000000e-01L, 5.000000000000000000e-01L, 8.064694554185326325e+00L, 9.167763417119923333e-08L},
	{1.250000000000000000e-01L, 2.500000000000000000e-01L, 8.873128656118993263e+01L, 3.769501719163865578e-07L},
	{3.125000000000000000e-02L, 6.250000000000000000e-02L, 4.096981198511257389e+03L, 1.205653952340679711e-07L},
	{1.562500000000000000e-02L, 3.125000000000000000e-02L, 1.574317583220441520e+04L, 1.380003928453846583e-07L},
	{5.000000000000000000e-01L, 1.000000000000000000e+00L, 3.256176475185617591e-01L, 1.912660677170175771e-08L},
	{1.953125000000000000e-03L, 3.906250000000000000e-03L, 9.225251570832365360e+02L, 2.132473175465897029e-09L},
	{7.812500000000000000e-03L, 1.562500000000000000e-02L, 2.899418134793237914e+04L, 1.934652413547325474e-07L},
      };

      auto alpha = _Tp{5};
      auto f = make_function<_Tp>(f15<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-7};
      auto out = __gnu_cxx::qagiu_integrate(w, fc, _Tp{0}, epsabs, epsrel);

      qtest.test_relative(out.__result, exp_result, epsrel, "qagiu(f15) smooth result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qagiu(f15) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qagiu(f15) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qagiu(f15) smooth last");
      qtest.test_integer(status, exp_ier, "qagiu(f15) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), _Tp{1} - test[i].b, fpeps, "qagiu(f15) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), _Tp{1} - test[i].a, fpeps, "qagiu(f15) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, epsrel, "qagiu(f15) smooth integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-4L}, "qagiu(f15) smooth abs error");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test infinite range integral f16 using an absolute error bound.
  try
    {
      //std::cout << ">>>> Test infinite range integral f16 using an absolute error bound..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{1.000000000006713292e-04L};
      const auto exp_abserr = _Tp{3.084062020905636316e-09L};
      const int exp_neval  = 165;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 6;

      constexpr std::size_t num_test = 6;
      test_ival<_Tp> test[num_test]
      {
	{0.000000000000000000e+00L, 3.125000000000000000e-02L, 7.633587786326674618e-05L, 3.084061858351569051e-09L},
	{5.000000000000000000e-01L, 1.000000000000000000e+00L, 9.900990099009899620e-07L, 3.112064814755089674e-17L},
	{1.250000000000000000e-01L, 2.500000000000000000e-01L, 3.629434715543053753e-06L, 4.908618166361344548e-17L},
	{6.250000000000000000e-02L, 1.250000000000000000e-01L, 6.501422186103209199e-06L, 3.014338672269481784e-17L},
	{3.125000000000000000e-02L, 6.250000000000000000e-02L, 1.062064387653501389e-05L, 6.795996738013555461e-18L},
	{2.500000000000000000e-01L, 5.000000000000000000e-01L, 1.922522349322310737e-06L, 4.543453652226561245e-17L},
      };

      const auto alpha = _Tp{1};
      auto f = make_function<_Tp>(f16<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      const auto epsabs = _Tp{1.0e-7};
      const auto epsrel = _Tp{0};
      auto out = __gnu_cxx::qagiu_integrate(w, fc, _Tp{99.9L}, epsabs, epsrel);

      qtest.test_relative(out.__result, exp_result, _Tp{1.0e-14L}, "qagiu(f16) smooth result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qagiu(f16) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qagiu(f16) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qagiu(f16) smooth last");
      qtest.test_integer(status, exp_ier, "qagiu(f16) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), _Tp{1} - test[i].b, fpeps, "qagiu(f16) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), _Tp{1} - test[i].a, fpeps, "qagiu(f16) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, epsabs, "qagiu(f16) smooth integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-4L}, "qagiu(f16) smooth abs error");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test infinite range integral myfn1 using an absolute error bound.
  try
    {
      //std::cout << ">>>> Test infinite range integral myfn1 using an absolute error bound..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{2.275875794468747770e+00L};
      const auto exp_abserr = _Tp{7.436490118267390744e-09L};
      const int exp_neval  = 270;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 5;

      constexpr std::size_t num_test = 5;
      test_ival<_Tp> test[num_test]
      {
	{5.000000000000000000e-01L, 1.000000000000000000e+00L, 1.691664195356748834e+00L, 4.265988974874425043e-09L},
	{1.250000000000000000e-01L, 2.500000000000000000e-01L, 4.639317228058405717e-04L, 3.169263960393051137e-09L},
	{0.000000000000000000e+00L, 1.250000000000000000e-01L, 4.379392477350953574e-20L, 8.360902986775307673e-20L},
	{2.500000000000000000e-01L, 3.750000000000000000e-01L, 1.146307471900291086e-01L, 1.231954072964969637e-12L},
	{3.750000000000000000e-01L, 5.000000000000000000e-01L, 4.691169201991640669e-01L, 5.208244060463541433e-15L},
      };

      auto f = make_function<_Tp>(myfn1<_Tp>);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      const auto epsabs = _Tp{1.0e-7};
      const auto epsrel = _Tp{0};
      //auto out = __gnu_cxx::qagi_integrate(w, fc, epsabs, epsrel);
      auto out = __gnu_cxx::qagis_integrate(w, fc, epsabs, epsrel);

      qtest.test_relative(out.__result, exp_result, _Tp{1.0e-14L}, "qagi(myfn1) smooth result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qagi(myfn1) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qagi(myfn1) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qagi(myfn1) smooth last");
      qtest.test_integer(status, exp_ier, "qagi(myfn1) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, fpeps, "qagi(myfn1) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, fpeps, "qagi(myfn1) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, _Tp{1.0e-14L}, "qagi(myfn1) smooth integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-4L}, "qagi(myfn1) smooth abs error");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test infinite range integral myfn2 using an absolute error bound.
  try
    {
      //std::cout << ">>>> Test infinite range integral myfn2 using an absolute error bound..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{2.718281828459044647e+00L};
      const auto exp_abserr = _Tp{1.588185109253204805e-10L};
      const int exp_neval  = 135;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 5;

      constexpr std::size_t num_test = 5;
      test_ival<_Tp> test[num_test]
      {
	{0.000000000000000000e+00L, 6.250000000000000000e-02L, 8.315287189746029816e-07L, 1.533437090413525935e-10L},
	{5.000000000000000000e-01L, 1.000000000000000000e+00L, 1.718281828459045091e+00L, 4.117868247943567505e-12L},
	{2.500000000000000000e-01L, 5.000000000000000000e-01L, 8.646647167633871867e-01L, 7.802455785301941044e-13L},
	{1.250000000000000000e-01L, 2.500000000000000000e-01L, 1.328565310599463256e-01L, 5.395586026138397182e-13L},
	{6.250000000000000000e-02L, 1.250000000000000000e-01L, 2.477920647947255521e-03L, 3.713312434866150125e-14L},
      };

      const auto alpha = _Tp{1};
      auto f = make_function<_Tp>(myfn2<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      const auto epsabs = _Tp{1.0e-7};
      const auto epsrel = _Tp{0};
      auto out = __gnu_cxx::qagil_integrate(w, fc, _Tp{1}, epsabs, epsrel);

      qtest.test_relative(out.__result, exp_result, _Tp{1.0e-14L}, "qagil(myfn2) smooth result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qagil(myfn2) smooth abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qagil(myfn2) smooth neval");
      qtest.test_integer(w.size(), exp_last, "qagil(myfn2) smooth last");
      qtest.test_integer(status, exp_ier, "qagil(myfn2) smooth status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, fpeps, "qagil(myfn2) smooth lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, fpeps, "qagil(myfn2) smooth upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, epsabs, "qagil(myfn2) smooth integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-4L}, "qagil(myfn2) smooth abs error");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test integral f454 with integrable singular points.
  try
    {
      //std::cout << ">>>> Test integral f454 with integrable singular points..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{5.274080611672716401e+01L};
      const auto exp_abserr = _Tp{1.755703848687062418e-04L};
      const int exp_neval  = 777;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 20;

      constexpr std::size_t num_test = 20;
      test_ival<_Tp> test[num_test]
      {
	{1.000000000000000000e+00L, 1.051776695296636976e+00L, -1.830392049835374568e-01L, 3.252808038935910834e-02L},
	{1.401269388548935790e+00L, 1.414213562373095145e+00L, -1.565132123531515207e-01L, 2.730454695485963826e-02L},
	{2.207106781186547462e+00L, 3.000000000000000000e+00L,  4.873920540843067783e+01L, 5.411138804637469780e-13L},
	{9.687500000000000000e-01L, 1.000000000000000000e+00L, -1.125078814079027711e-01L, 2.506431410088378817e-02L},
	{1.612436867076458391e+00L, 1.810660171779821415e+00L,  5.911661670635662835e-01L, 6.573952690524728748e-15L},
	{1.207106781186547462e+00L, 1.310660171779821415e+00L, -2.991531901645863023e-01L, 3.321267596107916554e-15L},
	{1.362436867076458391e+00L, 1.388325214724776657e+00L, -1.589345454585119055e-01L, 1.764527917763735212e-15L},
	{1.463769388548935790e+00L, 1.513325214724776657e+00L, -2.208730668000830344e-01L, 2.452183642810224359e-15L},
	{1.810660171779821415e+00L, 2.207106781186547462e+00L,  6.032891565603589079e+00L, 6.697855121200013106e-14L},
	{0.000000000000000000e+00L, 5.000000000000000000e-01L,  6.575875041899758092e-03L, 7.300687878575027348e-17L},
	{7.500000000000000000e-01L, 8.750000000000000000e-01L, -5.647871991778510847e-02L, 6.270397525408045936e-16L},
	{1.103553390593273731e+00L, 1.207106781186547462e+00L, -2.431894410706912923e-01L, 2.699945168224041491e-15L},
	{1.051776695296636976e+00L, 1.103553390593273731e+00L, -1.305470403178642658e-01L, 1.449363299575615261e-15L},
{0,0,0,0},
	{8.750000000000000000e-01L, 9.375000000000000000e-01L, -7.406626263352669715e-02L, 8.223007012367522077e-16L},
	{1.513325214724776657e+00L, 1.612436867076458391e+00L, -1.721363984401322045e-01L, 1.911097929242846383e-15L},
	{1.388325214724776657e+00L, 1.401269388548935790e+00L, -1.048692749517999567e-01L, 1.164282836272345215e-15L},
	{9.375000000000000000e-01L, 9.687500000000000000e-01L, -6.302287584527696551e-02L, 6.996944784151910810e-16L},
	{1.310660171779821415e+00L, 1.362436867076458391e+00L, -2.236786562536174916e-01L, 2.483331942899818875e-15L},
	{1.414213562373095145e+00L, 1.463769388548935790e+00L, -4.225328513207429193e-01L, 1.017446081816190118e-01L},
      };

      auto f = make_function<_Tp>(f454<_Tp>);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      std::vector<_Tp> pts{_Tp{0}, _Tp{1}, std::sqrt(_Tp{2}), _Tp{3}};

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-3};
      auto out = __gnu_cxx::qagp_integrate(w, fc, pts, epsabs, epsrel);

      qtest.test_relative(out.__result, exp_result, _Tp{1.0e-14L}, "qagp(f454) singular result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "qagp(f454) singular abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qagp(f454) singular neval");
      qtest.test_integer(w.size(), exp_last, "qagp(f454) singular last");
      qtest.test_integer(status, exp_ier, "qagp(f454) singular status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, fpeps, "qagp(f454) singular lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, fpeps, "qagp(f454) singular upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, epsrel, "qagp(f454) singular integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-4L}, "qagp(f454) singular abs error");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }


  // Test Cauchy integration using a relative error bound.
  try
    {
      //std::cout << ">>>> Test cauchy integration using a relative error bound..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{-8.994400695837000137e-02L};
      const auto exp_abserr = _Tp{1.185290176227023727e-06L};
      const int exp_neval  = 215;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 6;

      constexpr std::size_t num_test = 6;
      test_ival<_Tp> test[num_test]
      {
	{-1.000000000000000000e+00L, -7.500000000000000000e-01L, -1.234231128040012976e-01L, 1.172832717970022565e-06L},
	{-5.000000000000000000e-01L,  6.250000000000000000e-01L,  2.079093855884046535e-02L, 1.245463873006391609e-08L},
	{ 6.250000000000000000e-01L,  1.250000000000000000e+00L,  7.214232992127905808e-02L, 1.006998195150956048e-13L},
	{ 2.500000000000000000e+00L,  5.000000000000000000e+00L,  3.579970394639702888e-03L, 9.018232896137375412e-13L},
	{ 1.250000000000000000e+00L,  2.500000000000000000e+00L,  2.249831615049339983e-02L, 1.815172652101790755e-12L},
	{-7.500000000000000000e-01L, -5.000000000000000000e-01L, -8.553244917962132821e-02L, 1.833082948207153514e-15L},
      };

      auto f = make_function<_Tp>(f459<_Tp>);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-3};
      auto out = __gnu_cxx::qawc_integrate(w, fc, _Tp{-1}, _Tp{5}, _Tp{0}, epsabs, epsrel);

      qtest.test_relative(out.__result, exp_result, _Tp{1.0e-14L}, "qawc(f459) result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-6L}, "qawc(f459) abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qawc(f459) neval");
      qtest.test_integer(w.size(), exp_last, "qawc(f459) last");
      qtest.test_integer(status, exp_ier, "qawc(f459) status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, fpeps, "qawc(f459) lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, fpeps, "qawc(f459) upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, /*1.0e-14*/epsrel, "qawc(f459) integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-4L}, "qawc(f459) abs error");

      fc.num_evals(0);
      out = __gnu_cxx::qawc_integrate(w, fc, _Tp{5}, _Tp{-1}, _Tp{0}, epsabs, epsrel);

      qtest.test_relative(out.__result, -exp_result, _Tp{1.0e-14L}, "qawc(f459) rev result");
      qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-6L}, "qawc(f459) rev abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qawc(f459) rev neval");
      qtest.test_integer(w.size(), exp_last, "qawc(f459) rev last");
      qtest.test_integer(status, exp_ier, "qawc(f459) rev status");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test singular integration using a relative error bound.
  try
    {
      //std::cout << ">>>> Test adaptive singular integration using a relative error bound..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{-1.892751853489401670e-01L};
      const auto exp_abserr = _Tp{1.129133712015747658e-08L};
      const int exp_neval  = 280;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 8;

      constexpr std::size_t num_test = 8;
      test_ival<_Tp> test[num_test]
      {
	{0.000000000000000000e+00L, 7.812500000000000000e-03L, -4.126317299834445824e-05L, 1.129099387465713953e-08L},
	{2.500000000000000000e-01L, 5.000000000000000000e-01L, -6.240573216173390947e-02L, 6.928428071454762659e-16L},
	{5.000000000000000000e-01L, 1.000000000000000000e+00L, -1.076283950172247789e-01L, 3.423394967694403596e-13L},
	{6.250000000000000000e-02L, 1.250000000000000000e-01L, -3.408925115926728436e-03L, 3.784667152924835070e-17L},
	{3.125000000000000000e-02L, 6.250000000000000000e-02L, -8.914083918175634211e-04L, 9.896621209399419425e-18L},
	{1.562500000000000000e-02L, 3.125000000000000000e-02L, -2.574191402137795482e-04L, 2.857926564445496100e-18L},
	{1.250000000000000000e-01L, 2.500000000000000000e-01L, -1.456169844189576269e-02L, 1.616673288784094320e-16L},
	{7.812500000000000000e-03L, 1.562500000000000000e-02L, -8.034390712936630608e-05L, 8.919965558336773736e-19L},
      };

      auto f = make_function<_Tp>(f458<_Tp>);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::qaws_integration_table<_Tp> tb(_Tp{0}, _Tp{0}, 1, 0);
      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-7};
      auto out = __gnu_cxx::qaws_integrate(w, tb, fc, _Tp{0}, _Tp{1}, epsabs, epsrel);

      qtest.test_relative(out.__result, exp_result, _Tp{1.0e-14L}, "qaws(f458) ln(x-a) result");
      qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-6L}, "qaws(f458) ln(x-a) abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qaws(f458) ln(x-a) neval");
      qtest.test_integer(w.size(), exp_last, "qaws(f458) ln(x-a) last");
      qtest.test_integer(status, exp_ier, "qaws(f458) ln(x-a) status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, fpeps, "qaws(f458) ln(x-a) lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, fpeps, "qaws(f458) ln(x-a) upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, /*1.0e-14*/epsrel, "qaws(f458) ln(x-a) integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-4L}, "qaws(f458) ln(x-a) abs error");

      // Test without logs
      tb.set(_Tp{-0.5L}, _Tp{-0.3L}, 0, 0);
      out
	= __gnu_cxx::qaws_integrate(w, tb, fc, _Tp{0}, _Tp{1}, epsabs, epsrel);

      const auto exp_nolog_result = _Tp{9.896686656601706433e-01L};
      const auto exp_nolog_abserr = _Tp{5.888032513201251628e-08L};

      qtest.test_relative(out.__result, exp_nolog_result, _Tp{1.0e-14L}, "qaws(f458) AB result");
      qtest.test_relative(out.__abserr, exp_nolog_abserr, _Tp{1.0e-6L}, "qaws(f458) AB abserr");

      // Test with ln(x - a)
      tb.set(_Tp{-0.5L}, _Tp{-0.3L}, 1, 0);
      out = __gnu_cxx::qaws_integrate(w, tb, fc, _Tp{0}, _Tp{1}, epsabs, epsrel);

      const auto exp_logxma_result = _Tp{-3.636679470586539620e-01L};
      const auto exp_logxma_abserr = _Tp{2.851348775257054093e-08L};

      qtest.test_relative(out.__result, exp_logxma_result, _Tp{1.0e-14L}, "qaws(f458) AB ln(x-a) result");
      qtest.test_relative(out.__abserr, exp_logxma_abserr, _Tp{1.0e-6L}, "qaws(f458) AB ln(x-a) abserr");

      // Test with ln(b - x)
      tb.set(_Tp{-0.5L}, _Tp{-0.3L}, 0, 1);
      out
	= __gnu_cxx::qaws_integrate(w, tb, fc, _Tp{0}, _Tp{1}, epsabs, epsrel);

      const auto exp_logbmx_result = _Tp{-1.911489253363409802e+00L};
      const auto exp_logbmx_abserr = _Tp{9.854016753016499034e-09L};

      qtest.test_relative(out.__result, exp_logbmx_result, _Tp{1.0e-14L}, "qaws(f458) AB ln(b-x) result");
      qtest.test_relative(out.__abserr, exp_logbmx_abserr, _Tp{1.0e-6L}, "qaws(f458) AB ln(b-x) abserr");

      // Test with ln(x - a) ln(b - x)
      tb.set(_Tp{-0.5L}, _Tp{-0.3L}, 1, 1);
      out
	= __gnu_cxx::qaws_integrate(w, tb, fc, _Tp{0}, _Tp{1}, epsabs, epsrel);

      const auto exp_wlogs_result = _Tp{3.159922862811048172e-01L};
      const auto exp_wlogs_abserr = _Tp{2.336183482198144595e-08L};

      qtest.test_relative(out.__result, exp_wlogs_result, _Tp{1.0e-14L}, "qaws(f458) AB ln(x-a)ln(b-x) result");
      qtest.test_relative(out.__abserr, exp_wlogs_abserr, _Tp{1.0e-6L}, "qaws(f458) AB ln(x-a)ln(b-x) abserr");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test oscillatory integration using a relative error bound.
  try
    {
      //std::cout << ">>>> Test oscillatory integration using a relative error bound..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{-1.281368483991674190e-01L};
      const auto exp_abserr = _Tp{6.875028324415666248e-12L};
      const int exp_neval  = 305;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 9;

      constexpr std::size_t num_test = 9;
      test_ival<_Tp> test[num_test]
      {
	{5.000000000000000000e-01L, 1.000000000000000000e+00L,  2.190541162282139478e-02L, 1.302638552580516100e-13L},
	{2.500000000000000000e-01L, 5.000000000000000000e-01L, -2.587726479625663753e-02L, 7.259224351945759794e-15L},
	{1.250000000000000000e-01L, 2.500000000000000000e-01L,  5.483209176363500886e-02L, 1.249770395036711102e-14L},
	{1.562500000000000000e-02L, 3.125000000000000000e-02L, -3.886716016498160953e-02L, 4.315121611695628020e-16L},
	{3.125000000000000000e-02L, 6.250000000000000000e-02L, -9.178321994387816929e-02L, 1.018998440559284116e-15L},
	{6.250000000000000000e-02L, 1.250000000000000000e-01L, -3.081695575172510582e-02L, 7.832180081562836579e-16L},
	{7.812500000000000000e-03L, 1.562500000000000000e-02L, -1.242306301902117854e-02L, 1.379237060008662177e-16L},
	{3.906250000000000000e-03L, 7.812500000000000000e-03L, -3.659495117871544145e-03L, 4.062855738364339357e-17L},
	{0.000000000000000000e+00L, 3.906250000000000000e-03L, -1.447193692377651136e-03L, 8.326506625798146465e-07L},
      };

      auto f = make_function<_Tp>(f456<_Tp>);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);
      __gnu_cxx::oscillatory_integration_table<_Tp> wo(_Tp{10} * _S_pi<_Tp>, _Tp{1},
				       __gnu_cxx::oscillatory_integration_table<_Tp>::INTEG_SINE, 1000);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-7};
      auto out = __gnu_cxx::qawo_integrate(w, wo, fc, _Tp{0}, epsabs, epsrel);

      qtest.test_relative(out.__result, exp_result, _Tp{1.0e-14L}, "qawo(f456) result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-3L}, "qawo(f456) abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qawo(f456) neval");
      qtest.test_integer(w.size(), exp_last, "qawo(f456) last");
      qtest.test_integer(status, exp_ier, "qawo(f456) status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.lower_lim(i), test[i].a, fpeps, "qawo(f456) lower lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.upper_lim(i), test[i].b, fpeps, "qawo(f456) upper lim");

      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, _Tp{1.0e-14L}, "qawo(f456) integral");

      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{1.0e-2L}, "qawo(f456) abs error");

      // In reverse, flip limit and sign of length

      wo.set_length(_Tp{-1});
      fc.num_evals(0);
      out = qawo_integrate(w, wo, fc, _Tp{1}, epsabs, epsrel);

      qtest.test_relative(out.__result, -exp_result, _Tp{1.0e-14L}, "qawo(f456) rev result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-3L}, "qawo(f456) rev abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qawo(f456) rev neval");
      qtest.test_integer(w.size(), exp_last, "qawo(f456) rev last");
      qtest.test_integer(status, exp_ier, "qawo(f456) rev status");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test Fourier integration using an absolute error bound.
  try
    {
      //std::cout << ">>>> Test Fourier integration using an absolute error bound..." << std::endl;

      int status = 0;
      quadrature_test<_Tp> qtest;

      const auto exp_result = _Tp{9.999999999279802765e-01L};
      const auto exp_abserr = _Tp{1.556289974669056164e-08L};
      const int exp_neval  = 590;
      const int exp_ier    = __gnu_cxx::NO_ERROR;
      const int exp_last   = 12;

      constexpr std::size_t num_test = 12;
      test_ival<_Tp> test[num_test]
      {
	{0.0L, 0.0L,  1.013283128125232802e+00L, 1.224798040766472695e-12L},
	{0.0L, 0.0L, -1.810857954748607349e-02L, 1.396565155187268456e-13L},
	{0.0L, 0.0L,  7.466754034900931897e-03L, 1.053844511655910310e-16L},
	{0.0L, 0.0L, -1.352797860944863345e-03L, 5.854691731421723944e-18L},
	{0.0L, 0.0L,  8.136341270731781887e-04L, 2.439454888092388058e-17L},
	{0.0L, 0.0L, -7.093931338504278145e-04L, 2.130457268934021451e-17L},
	{0.0L, 0.0L,  1.680910783140869081e-03L, 9.757819552369539906e-18L},
	{0.0L, 0.0L, -4.360312526786496237e-03L, 6.505213034913026604e-19L},
	{0.0L, 0.0L,  1.119354921991485901e-03L, 4.553649124439220312e-18L},
	{0.0L, 0.0L,  2.950184068216192904e-03L, 7.155734338404329264e-18L},
	{0.0L, 0.0L, -9.462367583691360827e-04L, 7.643625316022806260e-18L},
	{0.0L, 0.0L, -2.168238443073697373e-03L, 1.105886215935214523e-17L},
      };

      auto f = make_function<_Tp>(f457<_Tp>);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);
      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> wc(1000);
      __gnu_cxx::oscillatory_integration_table<_Tp>
	wo(_S_pi<_Tp> / _Tp{2}, _Tp{1},
	  __gnu_cxx::oscillatory_integration_table<_Tp>::INTEG_COSINE, 1000);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-7};
      auto out = __gnu_cxx::qawf_integrate(w, wc, wo, fc, epsabs, epsrel);

      qtest.test_relative(out.__result, exp_result, _Tp{1.0e-14L}, "qawf(f457) result");
      if (is_double)
	qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-3L}, "qawf(f457) abserr");
      qtest.test_integer(fc.num_evals(), exp_neval, "qawf(f457) neval");
      qtest.test_integer(w.size(), exp_last, "qawf(f457) last");
      qtest.test_integer(status, exp_ier, "qawf(f457) status");

      const auto m = std::min(num_test, w.size());
      for (std::size_t i = 0; i < m; ++i)
	qtest.test_relative(w.result(i), test[i].r, _Tp{1.0e-12L}, "qawf(f457) integral");

      // We can only get within two orders of magnitude on the error
      // here, which is very sensitive to the floating point precision
      if (is_double)
	for (std::size_t i = 0; i < m; ++i)
	  qtest.test_relative(w.abs_error(i), test[i].e, _Tp{50}, "qawf(f457) abs error");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Sanity check monomial test function for fixed Gauss-Legendre rules.
  try
    {
      //std::cout << ">>>> Sanity check monomial test function for fixed Gauss-Legendre rules..." << std::endl;

      quadrature_test<_Tp> qtest;
      using dmon_t = monomial<_Tp>;

      qtest.test_absolute(dmon_t(2, _Tp{1})(_Tp{2}), _Tp{4}, 8*_S_eps, "monomial sanity check 1");

      qtest.test_absolute(dmon_t(1, _Tp{2})(_Tp{2}), _Tp{4}, 8*_S_eps, "monomial sanity check 2");

      qtest.test_absolute(integrate(dmon_t(2, _Tp{2}), _Tp{1}, _Tp{2}),
	  (_Tp{2}/_Tp{3})*(_Tp{2}*_Tp{2}*_Tp{2} - _Tp{1}*_Tp{1}*_Tp{1}), 8*_S_eps,
	  "integrate(monomial) sanity check");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test the fixed-order Gauss-Legendre rules with a monomial.
  try
    {
      //std::cout << ">>>> Test the fixed-order Gauss-Legendre rules with a monomial..." << std::endl;

      const auto a = _Tp{0}, b = _Tp{1.2};

      for (int n = 1; n < 1025; ++n)
	{
	  quadrature_test<_Tp> qtest;
	  __gnu_cxx::gauss_legendre_table<_Tp> tbl(n);

	  monomial<_Tp> mon(2*n-1, _Tp{1}); // n point rule exact for 2n-1 degree poly
	  auto expected = integrate(mon, a, b);
	  if (std::isinf(expected))
	    break;
	  auto result   = __gnu_cxx::glfixed_integrate(tbl, mon, a, b);

	  _Tp rel_tol;
	  if (tbl.precomputed)
	    rel_tol = prec_fixed<_Tp>;
	  else
	    rel_tol = _Tp{1.0e-7L};
	  std::ostringstream str;
	  str << "glfixed " << n << "-point: Integrating ("
	      << mon.constant << "*x^" << mon.degree
	      << ") over [" << a << "," << b << "]";
	  qtest.test_relative(result, expected, rel_tol, str.str().c_str());
	}
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Sanity check sin(x) test function for fixed Gauss-Legendre rules.
  try
    {
      //std::cout << ">>>> Sanity check sin(x) test function for fixed Gauss-Legendre rules..." << std::endl;

      quadrature_test<_Tp> qtest;

      qtest.test_absolute(f_sin(_Tp{2}), std::sin(_Tp{2}), _Tp{0}, "f_sin sanity check 1");
      qtest.test_absolute(f_sin(_Tp{7}), std::sin(_Tp{7}), _Tp{0}, "f_sin sanity check 2");
      qtest.test_absolute(integ_f_sin(_Tp{0}, _S_pi<_Tp>), _Tp{2}, _S_eps,
	  "integ_f_sin sanity check");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test the fixed-order Gauss-Legendre rules against sin(x) on [0, pi].
  try
    {
      //std::cout << ">>>> Test fixed-order Gauss-Legendre rules against sin(x) on [0, pi]..." << std::endl;

      const int n_max = 1024;
      const _Tp a = _Tp{0}, b = _S_pi<_Tp>;
      const auto expected = integ_f_sin(a, b);
      auto prev_abserr = _Tp{0};
      quadrature_test<_Tp> qtest;

      for (int n = 1; n <= n_max; ++n)
	{
	  __gnu_cxx::gauss_legendre_table<_Tp> tbl(n);

	  auto result = __gnu_cxx::glfixed_integrate(tbl, f_sin<_Tp>, a, b);
	  auto abserr = std::abs(expected - result);

	  if (n == 1)
	    {
	      std::ostringstream str;
	      str << "glfixed " << n << "-point: behavior for n == 1";
	      qtest.test_absolute(result, (b - a) * f_sin<_Tp>((b + a) / _Tp{2}), _Tp{0}, str.str().c_str());
	    }
	  else if (n < 9)
	    {
	      std::ostringstream str;
	      str << "glfixed " << n << "-point: observed drop in absolute error versus " << n-1 << "-points";
	      qtest.test_update(! (abserr < prev_abserr), str.str().c_str());
	    }
	  else if (tbl.precomputed)
	    {
	      std::ostringstream str;
	      str << "glfixed " << n << "-point: very low absolute error for high precision coefficients";
	      qtest.test_absolute(result, expected, _Tp{2} * n * _S_eps, str.str().c_str());
	    }
	  else
	    {
	      std::ostringstream str;
	      str << "glfixed " << n << "-point: acceptable absolute error for on-the-fly coefficients";
	      qtest.test_absolute(result, expected, _Tp{1.0e+6} * _S_eps, str.str().c_str());
	    }

	  prev_abserr = abserr;
	}
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test some fixed-order Gauss-Legendre rule points and weights on [-1, 1].
  // This verifies the (point, weight) retrieval API behaves sanely.
  try
    {
      //std::cout << ">>>> Test fixed-order Gauss-Legendre rule points and weights on [-1, 1]..." << std::endl;

      std::size_t n;
      quadrature_test<_Tp> qtest;

      // Analytical results for points and weights on [-1, 1]
      // Pulled from http://en.wikipedia.org/wiki/Gaussian_quadrature
      // Sorted in increasing order of Gauss points

      const _Tp
      e1[1][2]
      {
	{_Tp{0}, _Tp{2}}
      };

      const _Tp
      e2[2][2]
      {
	{_Tp{-1} / std::sqrt(_Tp{3}), _Tp{1}},
	{ _Tp{1} / std::sqrt(_Tp{3}), _Tp{1}}
      };

      const _Tp
      e3[3][2]
      {
	{-std::sqrt(_Tp{15}) / _Tp{5}, _Tp{5} / _Tp{9}},
	{		      _Tp{0}, _Tp{8} / _Tp{9}},
	{ std::sqrt(_Tp{15}) / _Tp{5}, _Tp{5} / _Tp{9}}
      };

      const _Tp e4c1 = _Tp{2} * std::sqrt(_Tp{6} / _Tp{5});
      const _Tp e4c2 = std::sqrt(_Tp{30});
      const _Tp
      e4[4][2]
      {
	{-std::sqrt((_Tp{3} + e4c1) / _Tp{7}), (_Tp{18} - e4c2) / _Tp{36}},
	{-std::sqrt((_Tp{3} - e4c1) / _Tp{7}), (_Tp{18} + e4c2) / _Tp{36}},
	{ std::sqrt((_Tp{3} - e4c1) / _Tp{7}), (_Tp{18} + e4c2) / _Tp{36}},
	{ std::sqrt((_Tp{3} + e4c1) / _Tp{7}), (_Tp{18} - e4c2) / _Tp{36}}
      };

      const _Tp e5c1 = std::sqrt(_Tp{10} / _Tp{7});
      const _Tp e5c2 = _Tp{13} * std::sqrt(_Tp{70});
      const _Tp
      e5[5][2]
      {
	{-std::sqrt((_Tp{5} + _Tp{2} * e5c1)) / _Tp{3}, (_Tp{322} - e5c2) / _Tp{900}},
	{-std::sqrt((_Tp{5} - _Tp{2} * e5c1)) / _Tp{3}, (_Tp{322} + e5c2) / _Tp{900}},
	{				       _Tp{0},	  _Tp{128} / _Tp{225}},
	{ std::sqrt((_Tp{5} - _Tp{2} * e5c1)) / _Tp{3}, (_Tp{322} + e5c2) / _Tp{900}},
	{ std::sqrt((_Tp{5} + _Tp{2} * e5c1)) / _Tp{3}, (_Tp{322} - e5c2) / _Tp{900}}
      };

      n = 1;
      __gnu_cxx::gauss_legendre_table<_Tp> tbl1(n);
      for (auto i = 0u; i < n; ++i)
	{
	  auto [xi, wi] = tbl1.get_point(_Tp{-1}, _Tp{1}, i);
	  std::ostringstream msg1, msg2;
	  msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
	  qtest.test_absolute(xi, e1[i][0], _S_eps, msg1.str().c_str());
	  msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
	  qtest.test_absolute(wi, e1[i][1], _S_eps, msg2.str().c_str());
	}

      n = 2;
      __gnu_cxx::gauss_legendre_table<_Tp> tbl2(n);
      for (auto i = 0u; i < n; ++i)
	{
	  auto [xi, wi] = tbl2.get_point(_Tp{-1}, _Tp{1}, i);
	  std::ostringstream msg1, msg2;
	  msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
	  qtest.test_absolute(xi, e2[i][0], _S_eps, msg1.str().c_str());
	  msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
	  qtest.test_absolute(wi, e2[i][1], _S_eps, msg2.str().c_str());
	}

      n = 3;
      __gnu_cxx::gauss_legendre_table<_Tp> tbl3(n);
      for (auto i = 0u; i < n; ++i)
	{
	  auto [xi, wi] = tbl3.get_point(_Tp{-1}, _Tp{1}, i);
	  std::ostringstream msg1, msg2;
	  msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
	  qtest.test_absolute(xi, e3[i][0], _S_eps, msg1.str().c_str());
	  msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
	  qtest.test_absolute(wi, e3[i][1], _S_eps, msg2.str().c_str());
	}

      n = 4;
      __gnu_cxx::gauss_legendre_table<_Tp> tbl4(n);
      for (auto i = 0u; i < n; ++i)
	{
	  auto [xi, wi] = tbl4.get_point(_Tp{-1}, _Tp{1}, i);
	  std::ostringstream msg1, msg2;
	  msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
	  qtest.test_absolute(xi, e4[i][0], _S_eps, msg1.str().c_str());
	  msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
	  qtest.test_absolute(wi, e4[i][1], _S_eps, msg2.str().c_str());
	}

      n = 5;
      __gnu_cxx::gauss_legendre_table<_Tp> tbl5(n);
      for (auto i = 0u; i < n; ++i)
	{
	  auto [xi, wi] = tbl5.get_point(_Tp{-1}, _Tp{1}, i);
	  std::ostringstream msg1, msg2;
	  msg1 << "glfixed " << n << "-point lookup: x(" << i << ')';
	  qtest.test_absolute(xi, e5[i][0], _S_eps, msg1.str().c_str());
	  msg2 << "glfixed " << n << "-point lookup: w(" << i << ')';
	  qtest.test_absolute(wi, e5[i][1], _S_eps, msg2.str().c_str());
	}
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test some fixed-order Gauss-Legendre rule points and weights on [-2, 3].
  // This verifies the (point, weight) retrieval API is okay on non-[-1,1].
  try
    {
      //std::cout << ">>>> Test some fixed-order Gauss-Legendre rule points and weights on [-2, 3]..." << std::endl;

      quadrature_test<_Tp> qtest;
      std::size_t n = 0;
      _Tp result;

      // Odd n = 3, f(x) = x**5 + x**4 + x**3 + x**2 + x**1 + 1
      n = 3;
      result = 0;
      __gnu_cxx::gauss_legendre_table<_Tp> tbl1(n);
      for (auto i = 0u; i < n; ++i)
	{
	  auto [x, w] = tbl1.get_point(_Tp{-2}, _Tp{3}, i);
	  result += w * (1 + x * (1 + x * (1 + x * (1 + x * (1 + x)))));
	}
      std::ostringstream msg1;
      msg1 << "glfixed " << n << "-point xi,wi eval";
      qtest.test_relative(result, _Tp{805} / _Tp{4}, _Tp{1.0e-8L}, msg1.str().c_str());

      // Even n = 4, f(x) = x**7 + x**6 + x**5 + x**4 + x**3 + x**2 + x**1 + 1
      n = 4;
      result = 0;
      __gnu_cxx::gauss_legendre_table<_Tp> tbl2(n);
      for (auto i = 0u; i < n; ++i)
	{
	  auto [x, w] = tbl2.get_point(_Tp{-2}, _Tp{3}, i);
	  result += w * (1 + x * (1 + x * (1 + x * (1 + x * (1 + x * (1 + x * (1 + x)))))));
	}
      std::ostringstream msg2;
      msg2 << "glfixed " << n << "-point xi,wi eval";
      qtest.test_relative(result, _Tp{73925} / _Tp{56}, _Tp{1.0e-8L}, msg2.str().c_str());
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test this newfangled cquad.
  try
    {
      //std::cout << ">>>> Test this newfangled cquad..." << std::endl;
      quadrature_test<_Tp> qtest;

      // Loop over the functions...
      for (int fid = 0; fid < 25; ++fid)
	{
	  auto f = make_function<_Tp>(func_tests<_Tp>[fid].fun);
	  auto a = func_tests<_Tp>[fid].a;
	  auto b = func_tests<_Tp>[fid].b;
	  auto exact = func_tests<_Tp>[fid].exact;
	  int status = 0;
	  //const auto rat = std::numeric_limits<_Tp>::epsilon() / std::numeric_limits<double>::epsilon();
	  auto rel_error = _Tp{1.0e-12L};

	  __gnu_cxx::cquad_workspace<_Tp, decltype(f(_Tp{}))> ws(200);

	  // Call our quadrature routine.
	  auto out = __gnu_cxx::cquad_integrate(ws, f, a, b, _Tp{0}, rel_error);

	  std::ostringstream rstr;
	  rstr << "cquad f" << (fid + 1);
	  qtest.test_relative(out.__result, exact, rel_error, rstr.str().c_str());

	  std::ostringstream upstr;
	  upstr << "cquad f" << (fid + 1) << " error("
			     << std::abs(out.__result - exact) << " actual vs "
			     << out.__abserr << " estimated)";
	  qtest.test_update(std::abs(out.__result - exact) > _Tp{5} * out.__abserr, upstr.str().c_str());

	  qtest.test_integer(status, 0, "cquad return code");
	  std::cout << std::flush;
	}
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test sinh-tanh.
  try
    {
      //std::cout << ">>>> Test this newfangled tanh-sinh..." << std::endl;
      quadrature_test<_Tp> qtest;

      // Loop over the functions...
      for (int fid = 0; fid < 25; ++fid)
	{
	  auto f = make_function<_Tp>(func_tests<_Tp>[fid].fun);
	  auto a = func_tests<_Tp>[fid].a;
	  auto b = func_tests<_Tp>[fid].b;
	  auto exact = func_tests<_Tp>[fid].exact;
	  int status = 0;
	  //const auto rat = std::numeric_limits<_Tp>::epsilon() / std::numeric_limits<double>::epsilon();
	  auto rel_error = _Tp{1.0e-12L};

	  // Call our quadrature routine.
	  auto out = __gnu_cxx::integrate_tanh_sinh(f, a, b, _Tp{0}, rel_error, 6);

	  std::ostringstream rstr;
	  rstr << "tanh_sinh f" << (fid + 1);
	  qtest.test_relative(out.__result, exact, rel_error, rstr.str().c_str());

	  std::ostringstream upstr;
	  upstr << "tanh_sinh f" << (fid + 1) << " error("
			     << std::abs(out.__result - exact) << " actual vs "
			     << out.__abserr << " estimated)";
	  qtest.test_update(std::abs(out.__result - exact) > _Tp{5} * out.__abserr, upstr.str().c_str());

	  qtest.test_integer(status, 0, "tanh_sinh return code");
	  std::cout << std::flush;
	}
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test exp_sinh integral f455 using a relative error bound.
  try
    {
      //std::cout << ">>>> Test exp_sinh integral f455 using a relative error bound..." << std::endl;

      const auto exp_result = _Tp{-3.616892186127022568e-01L};
      const auto exp_abserr = _Tp{3.016716913328831851e-06L};

      auto f = make_function<_Tp>(f455<_Tp>);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-3};
      auto out = __gnu_cxx::qagiu_integrate(w, fc, _Tp{0}, epsabs, epsrel);

      quadrature_test<_Tp> qtest;
      qtest.test_relative(out.__result, exp_result, /*1.0e-14*/epsrel, "exp_sinh(f455) smooth result");
      qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "exp_sinh(f455) smooth abserr");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test exp_sinh integral f15 using a relative error bound.
  try
    {
      //std::cout << ">>>> Test exp_sinh integral f15 using a relative error bound..." << std::endl;

      const auto exp_result = _Tp{6.553600000000024738e+04L};
      const auto exp_abserr = _Tp{7.121667111456009280e-04L};

      auto alpha = _Tp{5};
      auto f = make_function<_Tp>(f15<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      const auto epsabs = _Tp{0};
      const auto epsrel = _Tp{1.0e-7};
      auto out = __gnu_cxx::integrate_exp_sinh(fc, _Tp{0}, epsabs, epsrel);

      quadrature_test<_Tp> qtest;
      qtest.test_relative(out.__result, exp_result, epsrel, "exp_sinh(f15) smooth result");
      qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "exp_sinh(f15) smooth abserr");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test exp_sinh integral f16 using an absolute error bound.
  try
    {
      //std::cout << ">>>> Test exp_sinh integral f16 using an absolute error bound..." << std::endl;

      const auto exp_result = _Tp{1.000000000006713292e-04L};
      const auto exp_abserr = _Tp{3.084062020905636316e-09L};

      const auto alpha = _Tp{1};
      auto f = make_function<_Tp>(f16<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      const auto epsabs = _Tp{1.0e-7};
      const auto epsrel = _Tp{0};
      auto out = __gnu_cxx::integrate_exp_sinh(fc, _Tp{99.9L}, epsabs, epsrel);

      quadrature_test<_Tp> qtest;
      qtest.test_relative(out.__result, exp_result, _Tp{1.0e-14L}, "exp_sinh(f16) smooth result");
      qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "exp_sinh(f16) smooth abserr");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test sinh_sinh integral myfn1 using an absolute error bound.
  try
    {
      //std::cout << ">>>> Test exp_sinh integral myfn1 using an absolute error bound..." << std::endl;

      const auto exp_result = _Tp{2.275875794468747770e+00L};
      const auto exp_abserr = _Tp{7.436490118267390744e-09L};

      auto f = make_function<_Tp>(myfn1<_Tp>);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      const auto epsabs = _Tp{1.0e-7};
      const auto epsrel = _Tp{0};
      auto out = __gnu_cxx::integrate_sinh_sinh(fc, epsabs, epsrel);

      quadrature_test<_Tp> qtest;
      qtest.test_relative(out.__result, exp_result, _Tp{1.0e-14L}, "sinh_sinh(myfn1) smooth result");
      qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "sinh_sinh(myfn1) smooth abserr");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  // Test exp_sinh integral myfn2 using an absolute error bound.
  try
    {
      //std::cout << ">>>> Test exp_sinh integral myfn2 using an absolute error bound..." << std::endl;

      const auto exp_result = _Tp{2.718281828459044647e+00L};
      const auto exp_abserr = _Tp{1.588185109253204805e-10L};

      const auto alpha = _Tp{1};
      auto f = make_function<_Tp>(myfn2<_Tp>, alpha);
      auto fc = counted_function<_Tp, decltype(f)>(f);

      __gnu_cxx::integration_workspace<_Tp, decltype(fc(_Tp{}))> w(1000);

      const auto epsabs = _Tp{1.0e-7};
      const auto epsrel = _Tp{0};
      auto out = __gnu_cxx::integrate_exp_sinh(fc, _Tp{1}, epsabs, epsrel);

      quadrature_test<_Tp> qtest;
      qtest.test_relative(out.__result, exp_result, _Tp{1.0e-14L}, "exp_sinh(myfn2) smooth result");
      qtest.test_relative(out.__abserr, exp_abserr, _Tp{1.0e-5L}, "exp_sinh(myfn2) smooth abserr");
    }
  catch (__gnu_cxx::__integration_error<_Tp>& iex)
    {
      belch<_Tp>(iex);
    }
  catch (std::exception& ex)
    {
      belch(ex);
    }

  {
    using dmon_t = monomial<_Tp>;
    const dmon_t mon(5, _Tp{1});

    std::size_t n = 15;
    for (auto b = _Tp{1.1L}; b <= _Tp{4}; b += _Tp{0.1L})
      {
	const auto deg = mon.degree;
	const auto dterm = (deg % 2) == 0 ? _Tp{1} : _Tp{-1};

	// Test with a < b.
	// Then test with a > b.
	for (int k = -1; k != 3; k += 2)
	  {
	    auto a = b + _Tp(k);
	    auto bpa = b + a;
	    auto bma = b - a;

	    // Test Legendre quadrature.
	    auto exact = integrate(mon, a, b);
	    test_quadrature_rule(mon, a, b,
				 prec_fixed<_Tp>, exact, "legendre monomial",
				 __gnu_cxx::fixed_gauss_legendre_integral<_Tp>(n), n);

	    // Test Chebyshev T (first kind) quadrature.
	    exact = std::copysign(_Tp{1}, bma)
		  * _S_pi<_Tp> * std::pow(_Tp{0.5L} * bpa, _Tp(deg))
		  * __gnu_cxx::hyperg(_Tp(0.5L * (1 - deg)), _Tp(-0.5L * deg),
				      _Tp{1}, bma * bma / (bpa * bpa));
	    test_quadrature_rule(mon, a, b,
				 prec_fixed<_Tp>, exact, "chebyshev_t monomial",
				 __gnu_cxx::fixed_gauss_chebyshev_t_integral<_Tp>(n), n);

	    // Test Chebyshev U (second kind) quadrature.
	    exact = std::copysign(_Tp{1}, bma)
		  * _S_pi_2 * std::pow(_Tp{0.5L} * bpa, _Tp(deg))
		  * __gnu_cxx::hyperg(_Tp(0.5L * (1 - deg)), _Tp(-0.5L * deg),
				      _Tp{2}, bma * bma / (bpa * bpa))
		  * _Tp{0.25L} * bma * bma;
	    test_quadrature_rule(mon, a, b,
				 prec_fixed<_Tp>, exact, "chebyshev_u monomial",
				 __gnu_cxx::fixed_gauss_chebyshev_u_integral<_Tp>(n), n);

	    // Test Laguerre quadrature.
	    exact = std::pow(b, _Tp(-1 - deg))
		  * std::exp(a * b)
		  * __gnu_cxx::tgamma(_Tp(1 + deg), a * b);
	    test_quadrature_rule(mon, a, b,
				 prec_fixed<_Tp>, exact, "laguerre monomial",
				 __gnu_cxx::fixed_gauss_laguerre_integral<_Tp>(n, _Tp{0}),
				 n, _Tp{0});

	    // Test Hermite quadrature.
	    exact = _Tp{0.5L} * std::pow(b, _Tp(-0.5L * deg))
		  * (_Tp((1 - dterm) * deg) * a * std::tgamma(_Tp(0.5L * deg))
		      * __gnu_cxx::conf_hyperg(_Tp(0.5L * (1 - deg)), _Tp{1.5L}, -a * a * b)
		   + _Tp(1 + dterm) * std::tgamma(_Tp(0.5L * (1 + deg)))
		      * __gnu_cxx::conf_hyperg(_Tp(-0.5L * deg), _Tp{0.5L}, -a * a * b) / std::sqrt(b));
	    test_quadrature_rule(mon, a, b,
				 prec_fixed<_Tp>, exact, "hermite monomial",
				 __gnu_cxx::fixed_gauss_hermite_integral<_Tp>(n, _Tp{0}),
				 n, _Tp{0});
	  }
      }
  }

  {
    // Now test Gaussian rules on myfn1.
    int n = 200;

    test_quadrature_rule(myfn1<_Tp>, _Tp{1.2L}, _Tp{1.6L},
			 prec_fixed<_Tp>, _Tp{0.01505500344456001L}, "legendre myfn1",
			 __gnu_cxx::fixed_gauss_legendre_integral<_Tp>(n), n);
    test_quadrature_rule(myfn1<_Tp>, _Tp{1.2L}, _Tp{2.6L}, 
			 prec_fixed<_Tp>, _Tp{0.0582346516219999L}, "chebyshev_t myfn1",
			 __gnu_cxx::fixed_gauss_chebyshev_t_integral<_Tp>(n), n);
    test_quadrature_rule(myfn1<_Tp>, _Tp{1.2L}, _Tp{1.6L},
			 prec_fixed<_Tp>, _Tp{1.2279468957162412661311711271e-5L}, "gegenbauer myfn1",
			 __gnu_cxx::fixed_gauss_gegenbauer_integral<_Tp>(n, _Tp{2}),
			 n, _Tp{2});
    test_quadrature_rule(myfn1<_Tp>, _Tp{1.2L}, _Tp{1.6L},
			 prec_fixed<_Tp>, _Tp{1.228256086101808986e-1L}, "gegenbauer myfn1",
			 __gnu_cxx::fixed_gauss_gegenbauer_integral<_Tp>(n, _Tp{-0.5L}),
			 n, _Tp{-0.5L});
    test_quadrature_rule(myfn1<_Tp>, _Tp{1.2L}, _Tp{1.6L},
			 prec_fixed<_Tp>, _Tp{3.173064776410033e-5L}, "jacobi myfn1",
			 __gnu_cxx::fixed_gauss_jacobi_integral<_Tp>(n, _Tp{2}, _Tp{1.5L}),
			 n, _Tp{2}, _Tp{1.5L});
    test_quadrature_rule(myfn1<_Tp>, _Tp{1.2L}, _Tp{1.6L},
			 prec_fixed<_Tp>, _Tp{1.228256086101808986e-1L}, "jacobi myfn1",
			 __gnu_cxx::fixed_gauss_jacobi_integral<_Tp>(n, _Tp{-0.5L}, _Tp{-0.5L}),
			 n, _Tp{-0.5L}, _Tp{-0.5L});
    test_quadrature_rule(myfn1<_Tp>, _Tp{1.2L}, _Tp{0.6L},
			 prec_fixed<_Tp>, _Tp{0.006604180366378123L}, "laguerre myfn1",
			 __gnu_cxx::fixed_gauss_laguerre_integral<_Tp>(n, _Tp{0.5L}),
			 n, _Tp{0.5L});
    test_quadrature_rule(myfn1<_Tp>, _Tp{1.2L}, _Tp{0.6L},
			 prec_fixed<_Tp>, _Tp{0.6542819629825344L}, "hermite myfn1",
			 __gnu_cxx::fixed_gauss_hermite_integral<_Tp>(n, _Tp{1}), n, _Tp{1});
    test_quadrature_rule(myfn1<_Tp>, _Tp{1.2L}, _Tp{1.6L},
			 prec_fixed<_Tp>, _Tp{2.1315535492168832898083633e-4L}, "exponential myfn1",
			 __gnu_cxx::fixed_gauss_exponential_integral<_Tp>(n, _Tp{2}),
			 n, _Tp{2});
    test_quadrature_rule(myfn1<_Tp>, _Tp{1.2L}, _Tp{1.6L},
			 _Tp{1.0e-9L}, _Tp{4.8457468060064844e-20L}, "rational myfn1",
			 __gnu_cxx::fixed_gauss_rational_integral<_Tp>(15, _Tp{2}, _Tp{-33.4}),
			 15, _Tp{2}, _Tp{-33.4});
    test_quadrature_rule(myfn1<_Tp>, _Tp{1.2L}, _Tp{2.6L},
			 prec_fixed<_Tp>, _Tp{0.0081704088896491L}, "chebyshev_u myfn1",
			 __gnu_cxx::fixed_gauss_chebyshev_u_integral<_Tp>(n), n);
  }

  // Test Gegenbauer quadrature
  {
    using dmon_t = monomial<_Tp>;
    const dmon_t mon(5, _Tp{1});

    constexpr std::size_t num_test = 5;
    fixed_test<_Tp> test[num_test]
    {
      { 0.123L,  0.456L, 4.15933612154155020161400717857e-7L,  2.0L, 0.0L, 0.0L, 0.0L},
      { 7.747L, 12.0L,   7.44697808572324010134504819452e+5L,  0.5L, 0.0L, 0.0L, 0.0L},
      { 1.47L,   2.0L,   5.52024994284578980512106835228e+1L, -0.5L, 0.0L, 0.0L, 0.0L},
      {-1.47L,   2.0L,   7.95574829722734114107142857143L,     1.0L, 0.0L, 0.0L, 0.0L},
      { 0.0L,    0.47L,  1.79653588816666666666666666667e-3L,  0.0L, 0.0L, 0.0L, 0.0L}
    };

/* This is simpler - hyperg turns into a pochhammer.
    auto gegenbauer_moment
      = [n = mon.degree, alpha](auto a, auto b) -> auto
	{
	  using _Tp = decltype(alpha);
	  return n % 2 == 1
		 ? _Tp{0}
		 : _Tp{2} * std::tgamma(alpha + _Tp{1}) * std::tgamma(n + _Tp{1})
		   / std::tgamma(alpha + _Tp(n + 2))
		   * __gnu_cxx::hyperg(-alpha, _Tp(n + 1), alpha + _Tp(n + 2), _Tp{-1});
	};
*/

    std::size_t n = 50;
    for (std::size_t k = 0; k < num_test; ++k)
      {
	//exact = gegenbauer_moment(test[k].a, test[k].b);
	test_quadrature_rule(mon, test[k].a, test[k].b,
			     prec_fixed<_Tp>, test[k].r, "gegenbauer monomial",
			     __gnu_cxx::fixed_gauss_gegenbauer_integral<_Tp>(n, test[k].alpha),
			     n, test[k].alpha);
      }
  }

  // Test Jacobi quadrature
  {
    using dmon_t = monomial<_Tp>;
    const dmon_t mon(5, _Tp{1});

    std::size_t n = 50;
    const auto alpha = _Tp{2};
    const auto beta = _Tp{1.5L};

    constexpr std::size_t num_test = 5;
    fixed_test<_Tp> test[num_test]
    {
      { 0.123L,  0.456L, 9.052430592016123480501898e-7L,     alpha, beta, 0.0L, 0.0L},
      { 7.747L, 12.0L,   3.131716150347619771233591755e+6L,  alpha, beta, 0.0L, 0.0L},
      { 1.47L,   2.0L,   0.04435866422797298224404592896L,   alpha, beta, 0.0L, 0.0L},
      {-1.47L,   2.0L,   5.287059602300844442782407L,	alpha, beta, 0.0L, 0.0L},
      { 0.0L,    0.47L,  2.5337038518475893688512749675e-6L, alpha, beta, 0.0L, 0.0L}
    };

/*
    auto jacobi_moment
      = [n = mon.degree, alpha, beta](auto a, auto b) -> auto
	{
	  using _Tp = decltype(alpha + beta);
	  return std::tgamma(alpha + _Tp{1}) * std::tgamma(n + _Tp{1})
	       / std::tgamma(alpha + _Tp(n + 2))
	       * __gnu_cxx::hyperg(-beta, _Tp(n + 1), alpha + _Tp(n + 2), _Tp{-1})
	       + _Tp(n % 2 == 0 ? +1 : -1)
	       * std::tgamma(beta + _Tp{1}) * std::tgamma(n + _Tp{1})
	       / std::tgamma(beta + _Tp(n + 2))
	       * __gnu_cxx::hyperg(-alpha, _Tp(n + 1), beta + _Tp(n + 2), _Tp{-1});
	};
*/

    for (std::size_t k = 0; k < num_test; ++k)
      {
	//exact = jacobi_moment(test[k].a, test[k].b);
	test_quadrature_rule(mon, test[k].a, test[k].b,
			     prec_fixed<_Tp>, test[k].r, "jacobi monomial",
			     __gnu_cxx::fixed_gauss_jacobi_integral<_Tp>(n, test[k].alpha, test[k].beta),
			     n, test[k].alpha, test[k].beta);
      }
  }

  // Test Exponential quadrature
  {
    using dmon_t = monomial<_Tp>;
    const dmon_t mon(5, _Tp{1});

    constexpr std::size_t num_test = 5;
    fixed_test<_Tp> test[num_test]
    {
      { 0.123L,  0.456L, 1.598864206823942764921875e-4L,      1.0L, 0.0L, 0.0L, 0.0L},
      { 7.747L, 12.0L,   6.2461581848571833291063083819e+5L,  1.5L, 0.0L, 0.0L, 0.0L},
      { 1.47L,   2.0L,   2.22578063871903188095238095238e-1L, 2.0L, 0.0L, 0.0L, 0.0L},
      {-1.47L,   2.0L,   2.88968950008739567709168294271e+1L, 3.0L, 0.0L, 0.0L, 0.0L},
      { 0.0L,    0.47L,  4.62725113500425479890950520833e-7L, 5.0L, 0.0L, 0.0L, 0.0L}
    };

    std::size_t n = 50;
    for (std::size_t k = 0; k < num_test; ++k)
      {
	//exact = integrate(mon, test[k].a, test[k].b);
	test_quadrature_rule(mon, test[k].a, test[k].b,
			     prec_fixed<_Tp>, test[k].r, "exponential monomial",
			     __gnu_cxx::fixed_gauss_exponential_integral<_Tp>(n, test[k].alpha),
			     n, test[k].alpha);
      }
  }

  // Test Rational quadrature
  {
    using dmon_t = monomial<_Tp>;
    const dmon_t mon(5, _Tp{1});

    constexpr std::size_t num_test = 6;
    fixed_test<_Tp> test[num_test]
    {
      { 0.0L,    2.0L,    1.312245361412108703130374957e-10,    0.0L, -21.0L, 0.0L, 0.0L},
      { 0.123L,  0.456L,  1.70362044485924082779613124672e-2L,  1.0L, -12.0L, 0.0L, 0.0L},
      { 7.747L, 12.0L,    8.93065131938394658578136414201e-11L, 1.5L, -13.0L, 0.0L, 0.0L},
      { 1.47L,   2.0L,    7.17990217357447544326794457270e-13L, 2.0L, -22.0L, 0.0L, 0.0L},
      {-1.47L,   2.0L,   -1.10760676986664098133970869634e+1L,  3.0L, -21.0L, 0.0L, 0.0L},
      { 0.0L,    0.47L,   2.90392485414197833688178206557e-3L,  5.0L, -16.0L, 0.0L, 0.0L}
    };

    std::size_t n = 5;
    for (std::size_t k = 0; k < num_test; ++k)
      {
	//exact = integrate(mon, test[k].a, test[k].b);
	test_quadrature_rule(mon, test[k].a, test[k].b,
			     prec_fixed<_Tp>, test[k].r, "rational monomial",
			     __gnu_cxx::fixed_gauss_rational_integral<_Tp>(n, test[k].alpha, test[k].beta),
			     n, test[k].alpha, test[k].beta);
      }
  }

  return quadrature_test<_Tp>::test_summary();
}

int
main()
{
  std::cout << "\n\nTest double\n";
  std::cout << "=====================================\n";
  std::cerr << "\n\nTest double\n";
  std::cerr << "=====================================\n";
  test_quadrature<double>();

  std::cout << "\n\nTest long double\n";
  std::cout << "=====================================\n";
  std::cerr << "\n\nTest long double\n";
  std::cerr << "=====================================\n";
  test_quadrature<long double>();

  std::cout << "\n\nTest float\n";
  std::cout << "=====================================\n";
  std::cerr << "\n\nTest float\n";
  std::cerr << "=====================================\n";
  test_quadrature<float>();
}
