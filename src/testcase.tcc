/* cxx_integration/testcase.tcc
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * Copyright (C) 2016-2018 Free Software Foundation, Inc.
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

#ifndef QUADRATURE_TESTCASE_TCC
#define QUADRATURE_TESTCASE_TCC 1

#include <cmath>

/* These are the test functions from table 4.1 of the QUADPACK book */

/**
 * @f[
 *    f_1(x) = x^\alpha * \log(1/x)
 * @f]
 * @f[
 *    \int_{0}^{1} dx f_1(x) = 1/(\alpha + 1)^2
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f1(_Tp x, _Tp alpha)
  {
    return std::pow(x,alpha) * std::log(1/x);
  }

/**
 * @f[
 *    f_2(x) = 4^-alpha / ((x-pi/4)^2 + 16^-alpha)
 * @f]
 * @f[
 *    \int_{0}^{1} dx f_2(x) = \arctan((4-\pi)4^(\alpha-1))
 *                           + \arctan(\pi 4^(\alpha-1))
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f2(_Tp x, _Tp alpha)
  {
    const auto _S_pi_4 = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L} / _Tp{4};
    return std::pow(_Tp{4}, -alpha) / (std::pow((x - _S_pi_4), _Tp{2}) + std::pow(_Tp{16}, -alpha));
  }

/**
 * @f[
 *    f_3(x) = cos(2^alpha * sin(x))
 * @f]
 * @f[
 *    \int_{0}^{\pi} dx f_3(x) = \pi J_0(2^\alpha)
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f3(_Tp x, _Tp alpha)
  {
    return std::cos(std::pow(_Tp{2},alpha) * std::sin(x));
  }

/* Functions 4, 5 and 6 are duplicates of functions  1, 2 and 3 */
/* ....                                                         */

/**
 * @f[
 *    f_7(x) = |x - 1/3|^\alpha
 * @f]
 * @f[
 *    \int_{0}{1} dx f_7(x) = ((2/3)^(alpha+1) + (1/3)^(alpha+1))/(alpha + 1)
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f7(_Tp x, _Tp alpha)
  {
    return std::pow(std::abs(x - (_Tp{1}/_Tp{3})), alpha);
  }

/**
 * @f[
 *    f_8(x) = |x - \pi/4|^\alpha
 * @f]
 * @f[
 *    \int_{0}{1} dx f_8(x) = ((1 - \pi/4)^(\alpha + 1) + (\pi/4)^(\alpha + 1))
 *                          / (\alpha + 1)
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f8(_Tp x, _Tp alpha)
  {
    const auto _S_pi_4 = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L} / _Tp{4};
    return std::pow(std::abs(x - _S_pi_4), alpha);
  }

/**
 * @f[
 *    f_9(x) = \sqrt(1 - x^2) / (x + 1 + 2^{-\alpha})
 * @f]
 * @f[
 *    \int_{-1}{+1} dx f_9(x) = \pi/\sqrt((1+2^{-\alpha})^2-1)
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f9(_Tp x, _Tp alpha)
  {
    return 1 / ((x + 1 + std::pow(_Tp{2}, -alpha)) * std::sqrt(1 - x * x));
  }

/**
 * @f[
 *    f_{10}(x) = std::sin(x)^(alpha - 1)
 * @f]
 * @f[
 *    \int_{0}{pi/2} dx f_{10}(x) = 2^(\alpha-2) ((\Gamma(\alpha/2))^2)/\Gamma(\alpha)
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f10(_Tp x, _Tp alpha)
  {
    return std::pow(std::sin(x), alpha-1);
  }

/**
 * @f[
 *    f_{11}(x) = \log(1/x)^(\alpha - 1)
 * @f]
 * @f[
 *    \int_{0}{1} dx f_{11}(x) = \Gamma(\alpha)
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f11(_Tp x, _Tp alpha)
  {
    return std::pow(std::log(1/x), alpha-1);
  }

/**
 * @f[
 *    f_{12}(x) = \exp(20(x-1)) \sin(2^\alpha x)
 * @f]
 * @f[
 *    \int_{0}{1} dx f_{12}(x) =
 *      (20 \sin(2^\alpha) - 2^\alpha \cos(2^\alpha) + 2^\alpha \exp(-20))
 *       /(400 + 4^\alpha)
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f12(_Tp x, _Tp alpha)
  {
    return std::exp(20 * (x - 1)) * std::sin(std::pow(_Tp{2}, alpha) * x);
  }

/**
 * @f[
 *    f_{13}(x) = \cos(2^\alpha x)/\sqrt(x(1 - x))
 * @f]
 * @f[
 *    \int_{0}{1} dx f_{13}(x) = \pi \cos(2^(\alpha-1)) J_0(2^(\alpha-1))
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f13(_Tp x, _Tp alpha)
  {
    return std::cos(std::pow(_Tp{2}, alpha) * x) / std::sqrt(x * (1 - x));
  }

/**
 * @f[
 *    f_{14}(x) = \exp(-2^{-\alpha} x)\cos(x)/\sqrt(x)
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f14(_Tp x, _Tp alpha)
  {
    return std::exp(-std::pow(_Tp{2}, -alpha) * x) * std::cos(x) / std::sqrt(x);
  }

/**
 * @f[
 *    f_{15}(x) = x^2 \exp(-2^{-\alpha} x)
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f15(_Tp x, _Tp alpha)
  {
    return x * x * std::exp(-std::pow(_Tp{2}, -alpha) * x);
  }

/**
 * @f[
 *    f_{16}(x) = x^{\alpha - 1} / (1 + 10x)^2
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f16(_Tp x, _Tp alpha)
  {
    if (x == 0 && alpha == 1)
      return 1;  /* make the function continuous in x */
    if (x == 0 && alpha > 1)
      return 0;   /* avoid problems with pow(0,1) */
    return std::pow(x, alpha - 1) / std::pow((1 + 10 * x), _Tp{2});
  }

/**
 * @f[
 *    f_{17}(x) = 2^{-\alpha} / (((x - 1)^2 + 4^{-\alpha})(x-2))
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f17(_Tp x, _Tp alpha)
  {
    return std::pow(_Tp{2}, -alpha)
  	  / (((x - 1) * (x - 1) + std::pow(_Tp{4}, -alpha)) * (x - 2));
  }

/**
 * @f[
 *    f_{454}(x) = x^3 \log|(x^2-1)(x^2-2)|
 * @f]
 * @f[
 *    \int_{0}{\infty} dx f454(x) = 61 \log(2) + (77/4) \log(7) - 27
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f454(_Tp x)
  {
    _Tp x2 = x * x;
    _Tp x3 = x * x2;
    return x3 * std::log(std::abs((x2 - _Tp{1}) * (x2 - _Tp{2})));
  }

/**
 * @f[
 *    f_{455}(x) = log(x)/(1+100*x^2)
 * @f]
 * @f[
 *    \int_{0}{\infty} dx f455 = -\log(10)/20
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f455(_Tp x)
  {
    return std::log(x) / (_Tp{1} + _Tp{100} * x * x);
  }

/**
 * @f[
 *    f_{456}(x) = \log(x)
 * @f]
 * @f[
 *    \int_{0}{1} dx f456(x)\sin(10\pi x) = -(\gamma + \log(10pi) - Ci(10\pi))
 *                                        / (10 \pi)
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f456(_Tp x)
  {
    if (x == _Tp{0})
      return _Tp{0};
    return std::log(x);
  }

/**
 * @f[
 *    f_{457}(x) = 1/\sqrt(x)
 * @f]
 * @f[
 *    \int_{0}{+\infty} dx f457(x)\cos(\pi x/2) = 1
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f457(_Tp x)
  {
    if (x == _Tp{0})
      return _Tp{0};
    return 1 / std::sqrt(x);
  }

/**
 * @f[
 *    f_{458}(x) = 1/(1 + \log(x)^2)^2
 * @f]
 * @f[
 *    \int_{0}{1} dx log(x) f458(x) = (Ci(1) \sin(1) + (\pi/2 - Si(1)) \cos(1))
 *                                  / \pi
 *                             = -0.1892752
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f458(_Tp x)
  {
    if (x == _Tp{0})
      return _Tp{0};
    else
      {
	_Tp u = std::log(x);
	_Tp v = 1 + u * u;
	return _Tp{1} / (v * v);
      }
  }

/**
 * @f[
 *    f_{459}(x) = 1/(5 x^3 + 6)
 * @f]
 * @f[
 *    \int_{-1}{5} dx f459(x)/(x-0) = \log(125/631)/18
 * @f]
 */
template<typename _Tp>
  inline _Tp
  f459(_Tp x)
  {
    return _Tp{1} / (_Tp{5} * x * x * x + _Tp{6});
  }

/**
 * @f[
 *    myfn1(x) = \exp(-x - x^2)
 * @f]
 * @f[
 *    \int_{-\infty}{+\infty} dx myfn1(x) = \sqrt(\pi) \exp(-1/4)
 * @f]
 */
template<typename _Tp>
  inline _Tp
  myfn1(_Tp x)
  {
    return std::exp(-x - x * x);
  }

/**
 * @f[
 *    myfn2(x) = \exp(\alpha x)
 * @f]
 * @f[
 *    \int_{-\infty}{b} dx myfn2(x) = \exp(\alpha b)/\alpha
 * @f]
 */
template<typename _Tp>
  inline _Tp
  myfn2(_Tp x, _Tp alpha)
  {
    return std::exp(alpha * x);
  }

/* f_monomial = constant * x^degree

template<typename _Tp>
  inline _Tp
  f_monomial(_Tp x, _Tp alpha)
  {
    monomial<_Tp>* p = (monomial<_Tp> *)params;

    return p->constant * std::pow(x, p->degree);
  }
*/
template<typename _Tp>
  inline _Tp
  integ_f_monomial(_Tp a, _Tp b, monomial<_Tp>& p)
  {
    const int degreep1 = p.degree + 1;
    const _Tp bnp1 = std::pow(b, degreep1);
    const _Tp anp1 = std::pow(a, degreep1);
    return (p.constant / degreep1) * (bnp1 - anp1);
  }

template<typename _Tp>
  inline _Tp
  f_sin(_Tp x)
  {
    return std::sin(x);
  }

template<typename _Tp>
  inline _Tp
  integ_f_sin(_Tp a, _Tp b)
  {
    return -std::cos(b) + std::cos(a);
  }

/*
 * Another set of set functions..
 * A new quadrature routine for improper and oscillatory integrals
 * Applied Mathematics and Computation 189 (2007) pp. 452-461
 * E. Sermutlu, H.T. Eyyubo~glu
 */

template<typename _Tp>
  const auto _S_pi = _Tp{3.1415'92653'58979'32384'62643'38327'95028'84195L};

template<typename _Tp>
  inline _Tp
  cqf1(_Tp x)
  {
    return std::exp(x);
  }

template<typename _Tp>
  inline _Tp
  cqf2(_Tp x)
  {
    return x >= 0.3;
  }

template<typename _Tp>
  inline _Tp
  cqf3(_Tp x)
  {
    return std::sqrt(x);
  }


template<typename _Tp>
  inline _Tp
  cqf4(_Tp x)
  {
    return (_Tp{23} / _Tp{25}) * std::cosh(x) - std::cos(x);
  }

template<typename _Tp>
  inline _Tp
  cqf5(_Tp x)
  {
    _Tp x2 = x*x;
    return _Tp{1} / (x2 * (x2 + 1) + _Tp{0.9});
  }

template<typename _Tp>
  inline _Tp
  cqf6(_Tp x)
  {
    return x * std::sqrt(x);
  }

// Singular at 0.
template<typename _Tp>
  inline _Tp
  cqf7(_Tp x)
  {
    return _Tp{1} / std::sqrt(x);
  }

template<typename _Tp>
  inline _Tp
  cqf8(_Tp x)
  {
    _Tp x2 = x * x;
    return _Tp{1} / (_Tp{1} + x2*x2);
  }

template<typename _Tp>
  inline _Tp
  cqf9(_Tp x)
  {
    return _Tp{2} / (_Tp{2} + std::sin(_Tp{10} * _S_pi<_Tp> * x));
  }

template<typename _Tp>
  inline _Tp
  cqf10(_Tp x)
  {
    return _Tp{1} / (_Tp{1} + x);
  }

template<typename _Tp>
  inline _Tp
  cqf11(_Tp x)
  {
    return _Tp{1} / (_Tp{1} + std::exp(x));
  }

template<typename _Tp>
  inline _Tp
  cqf12(_Tp x)
  {
    if (x == _Tp{0})
      return _Tp{1};
    else
      return x / (std::exp(x) - _Tp{1});
  }

template<typename _Tp>
  inline _Tp
  cqf13(_Tp x)
  {
    return std::sin(_Tp{100} * _S_pi<_Tp> * x) / (_S_pi<_Tp> * x);
  }

template<typename _Tp>
  inline _Tp
  cqf14(_Tp x)
  {
    return std::sqrt(_Tp{5}) * std::exp(-_Tp{10} * _S_pi<_Tp> * x * x);
  }

template<typename _Tp>
  inline _Tp
  cqf15(_Tp x)
  {
    return _Tp{25} * std::exp(-_Tp{25} * x);
  }

template<typename _Tp>
  inline _Tp
  cqf16(_Tp x)
  {
    return _Tp{50} / (_S_pi<_Tp> * (_Tp{2500} * x * x + _Tp{1}));
  }

template<typename _Tp>
  inline _Tp
  cqf17(_Tp x)
  {
    _Tp t1 = _Tp{50} * _S_pi<_Tp> * x, t2;
    t2 = std::sin(t1) / t1;
    return _Tp{50} * t2 * t2;
  }

template<typename _Tp>
  inline _Tp
  cqf18(_Tp x)
  {
    return std::cos(std::cos(x)
		  + _Tp{3} * std::sin(x)
		  + _Tp{2} * std::cos(_Tp{2} * x)
		  + _Tp{3} * std::cos(_Tp{3} * x));
  }

// The paper regularizes this to 0 for x < 1e-15.
template<typename _Tp>
  inline _Tp
  cqf19(_Tp x)
  {
    return std::log(x);
  }

template<typename _Tp>
  inline _Tp
  cqf20(_Tp x)
  {
    return _Tp{1} / (x * x + _Tp{1.005});
  }

template<typename _Tp>
  inline _Tp
  cqf21(_Tp x)
  {
    return _Tp{1} / std::cosh(_Tp{20} * x - _Tp{4})
         + _Tp{1} / std::cosh(_Tp{20} * x - _Tp{8})
         + _Tp{1} / std::cosh(_Tp{20} * x - _Tp{12});
  }

template<typename _Tp>
  inline _Tp
  cqf22(_Tp x)
  {
    return _Tp{4} * _S_pi<_Tp> * _S_pi<_Tp> * x
	  * std::sin(_Tp{20} * _S_pi<_Tp> * x)
	  * std::cos(_Tp{2} * _S_pi<_Tp> * x);
  }

template<typename _Tp>
  inline _Tp
  cqf23(_Tp x)
  {
    _Tp t = _Tp{230} * x - _Tp{30};
    return _Tp{1} / (_Tp{1} + t * t);
  }

template<typename _Tp>
  inline _Tp
  cqf24(_Tp x)
  {
    return std::floor(std::exp(x));
  }

template<typename _Tp>
  inline _Tp
  cqf25(_Tp x)
  {
    return (x < _Tp{1}) * (x + _Tp{1})
         + (_Tp{1} <= x && x <= _Tp{3}) * (_Tp{3} - x)
         + (x > _Tp{3}) * _Tp{2};
  }

// The limits of integration and the expected results fot the 
template<typename _Tp>
  const static func_test<_Tp>
  func_tests[25]
  {
    &cqf1,  _Tp{0},  _Tp{1},	 _Tp{1.718281828459045235360287471352662497759L},
    &cqf2,  _Tp{0},  _Tp{1},	 _Tp{0.7L},
    &cqf3,  _Tp{0},  _Tp{1},	 _Tp{2.0L} / _Tp{3.0L},
    &cqf4,  _Tp{-1}, _Tp{1},	 _Tp{0.4794282266888016673585779618353075006421L},
    &cqf5,  _Tp{-1}, _Tp{1},	 _Tp{1.582232963729672933117468949026169067924L},
    &cqf6,  _Tp{0},  _Tp{1},	 _Tp{0.4L},
    &cqf7,  _Tp{0},  _Tp{1},	 _Tp{2.0L},
    &cqf8,  _Tp{0},  _Tp{1},	 _Tp{0.8669729873399110375739951638828707136522L},
    &cqf9,  _Tp{0},  _Tp{1},	 _Tp{1.154700538379251529018297561003914911295L},
    &cqf10, _Tp{0},  _Tp{1},	 _Tp{0.6931471805599453094172321214581765680755L},
    &cqf11, _Tp{0},  _Tp{1},	 _Tp{0.3798854930417224753682366264903209261602L},
    &cqf12, _Tp{0},  _Tp{1},	 _Tp{0.7775046341122482764175865454257105071925L},
    &cqf13, _Tp{0.1},_Tp{1},	 _Tp{0.9098637539166842915557830641141434835684e-2L},
    &cqf14, _Tp{0},  _Tp{10},	 _Tp{0.3535533905932737622004221810524245196424L},
    &cqf15, _Tp{0},  _Tp{10},	 _Tp{1.0L},
    &cqf16, _Tp{0},  _Tp{10},	 _Tp{0.4993633810764567446362485183117640508837L},
    &cqf17, _Tp{0},  _Tp{1},	 _Tp{0.4989868086930455024989853136560590525530L},
    &cqf18, _Tp{0},  _S_pi<_Tp>, _Tp{0.2910187828600526985238845968602020712077L},
    &cqf19, _Tp{0},  _Tp{1},	 _Tp{-1.0L},
    &cqf20, _Tp{-1}, _Tp{1},	 _Tp{1.564396444069049773091493015808472813088L},
    &cqf21, _Tp{0},  _Tp{1},	 _Tp{0.4693392062964089464236343930140556287142L},
    &cqf22, _Tp{0},  _Tp{1},	 _Tp{-0.6346651825433925734267966430867682594338L},
    &cqf23, _Tp{0},  _Tp{1},	 _Tp{0.1349248564946777269188547624864782167827e-1L},
    &cqf24, _Tp{0},  _Tp{3},	 _Tp{17.66438353924651497034012402929007814264L},
    &cqf25, _Tp{0},  _Tp{5},	 _Tp{7.5L}
  };

#endif // QUADRATURE_TESTCASE_TCC
