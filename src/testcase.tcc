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
template<typename Tp>
  inline Tp
  f1(Tp x, Tp alpha)
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
template<typename Tp>
  inline Tp
  f2(Tp x, Tp alpha)
  {
    const auto s_pi_4 = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L} / Tp{4};
    return std::pow(Tp{4}, -alpha) / (std::pow((x - s_pi_4), Tp{2}) + std::pow(Tp{16}, -alpha));
  }

/**
 * @f[
 *    f_3(x) = cos(2^alpha * sin(x))
 * @f]
 * @f[
 *    \int_{0}^{\pi} dx f_3(x) = \pi J_0(2^\alpha)
 * @f]
 */
template<typename Tp>
  inline Tp
  f3(Tp x, Tp alpha)
  {
    return std::cos(std::pow(Tp{2},alpha) * std::sin(x));
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
template<typename Tp>
  inline Tp
  f7(Tp x, Tp alpha)
  {
    return std::pow(std::abs(x - (Tp{1}/Tp{3})), alpha);
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
template<typename Tp>
  inline Tp
  f8(Tp x, Tp alpha)
  {
    const auto s_pi_4 = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L} / Tp{4};
    return std::pow(std::abs(x - s_pi_4), alpha);
  }

/**
 * @f[
 *    f_9(x) = \sqrt(1 - x^2) / (x + 1 + 2^{-\alpha})
 * @f]
 * @f[
 *    \int_{-1}{+1} dx f_9(x) = \pi/\sqrt((1+2^{-\alpha})^2-1)
 * @f]
 */
template<typename Tp>
  inline Tp
  f9(Tp x, Tp alpha)
  {
    return 1 / ((x + 1 + std::pow(Tp{2}, -alpha)) * std::sqrt(1 - x * x));
  }

/**
 * @f[
 *    f_{10}(x) = std::sin(x)^(alpha - 1)
 * @f]
 * @f[
 *    \int_{0}{pi/2} dx f_{10}(x) = 2^(\alpha-2) ((\Gamma(\alpha/2))^2)/\Gamma(\alpha)
 * @f]
 */
template<typename Tp>
  inline Tp
  f10(Tp x, Tp alpha)
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
template<typename Tp>
  inline Tp
  f11(Tp x, Tp alpha)
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
template<typename Tp>
  inline Tp
  f12(Tp x, Tp alpha)
  {
    return std::exp(20 * (x - 1)) * std::sin(std::pow(Tp{2}, alpha) * x);
  }

/**
 * @f[
 *    f_{13}(x) = \cos(2^\alpha x)/\sqrt(x(1 - x))
 * @f]
 * @f[
 *    \int_{0}{1} dx f_{13}(x) = \pi \cos(2^(\alpha-1)) J_0(2^(\alpha-1))
 * @f]
 */
template<typename Tp>
  inline Tp
  f13(Tp x, Tp alpha)
  {
    return std::cos(std::pow(Tp{2}, alpha) * x) / std::sqrt(x * (1 - x));
  }

/**
 * @f[
 *    f_{14}(x) = \exp(-2^{-\alpha} x)\cos(x)/\sqrt(x)
 * @f]
 */
template<typename Tp>
  inline Tp
  f14(Tp x, Tp alpha)
  {
    return std::exp(-std::pow(Tp{2}, -alpha) * x) * std::cos(x) / std::sqrt(x);
  }

/**
 * @f[
 *    f_{15}(x) = x^2 \exp(-2^{-\alpha} x)
 * @f]
 */
template<typename Tp>
  inline Tp
  f15(Tp x, Tp alpha)
  {
    return x * x * std::exp(-std::pow(Tp{2}, -alpha) * x);
  }

/**
 * @f[
 *    f_{16}(x) = x^{\alpha - 1} / (1 + 10x)^2
 * @f]
 */
template<typename Tp>
  inline Tp
  f16(Tp x, Tp alpha)
  {
    if (x == 0 && alpha == 1)
      return 1;  /* make the function continuous in x */
    if (x == 0 && alpha > 1)
      return 0;   /* avoid problems with pow(0,1) */
    return std::pow(x, alpha - 1) / std::pow((1 + 10 * x), Tp{2});
  }

/**
 * @f[
 *    f_{17}(x) = 2^{-\alpha} / (((x - 1)^2 + 4^{-\alpha})(x-2))
 * @f]
 */
template<typename Tp>
  inline Tp
  f17(Tp x, Tp alpha)
  {
    return std::pow(Tp{2}, -alpha)
  	  / (((x - 1) * (x - 1) + std::pow(Tp{4}, -alpha)) * (x - 2));
  }

/**
 * @f[
 *    f_{454}(x) = x^3 \log|(x^2-1)(x^2-2)|
 * @f]
 * @f[
 *    \int_{0}{\infty} dx f454(x) = 61 \log(2) + (77/4) \log(7) - 27
 * @f]
 */
template<typename Tp>
  inline Tp
  f454(Tp x)
  {
    Tp x2 = x * x;
    Tp x3 = x * x2;
    return x3 * std::log(std::abs((x2 - Tp{1}) * (x2 - Tp{2})));
  }

/**
 * @f[
 *    f_{455}(x) = log(x)/(1+100*x^2)
 * @f]
 * @f[
 *    \int_{0}{\infty} dx f455 = -\log(10)/20
 * @f]
 */
template<typename Tp>
  inline Tp
  f455(Tp x)
  {
    return std::log(x) / (Tp{1} + Tp{100} * x * x);
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
template<typename Tp>
  inline Tp
  f456(Tp x)
  {
    if (x == Tp{0})
      return Tp{0};
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
template<typename Tp>
  inline Tp
  f457(Tp x)
  {
    if (x == Tp{0})
      return Tp{0};
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
template<typename Tp>
  inline Tp
  f458(Tp x)
  {
    if (x == Tp{0})
      return Tp{0};
    else
      {
	Tp u = std::log(x);
	Tp v = 1 + u * u;
	return Tp{1} / (v * v);
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
template<typename Tp>
  inline Tp
  f459(Tp x)
  {
    return Tp{1} / (Tp{5} * x * x * x + Tp{6});
  }

/**
 * @f[
 *    myfn1(x) = \exp(-x - x^2)
 * @f]
 * @f[
 *    \int_{-\infty}{+\infty} dx myfn1(x) = \sqrt(\pi) \exp(-1/4)
 * @f]
 */
template<typename Tp>
  inline Tp
  myfn1(Tp x)
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
template<typename Tp>
  inline Tp
  myfn2(Tp x, Tp alpha)
  {
    return std::exp(alpha * x);
  }

/* f_monomial = constant * x^degree

template<typename Tp>
  inline Tp
  f_monomial(Tp x, Tp alpha)
  {
    monomial<Tp>* p = (monomial<Tp> *)params;

    return p->constant * std::pow(x, p->degree);
  }
*/
template<typename Tp>
  inline Tp
  integ_f_monomial(Tp a, Tp b, monomial<Tp>& p)
  {
    const int degreep1 = p.degree + 1;
    const Tp bnp1 = std::pow(b, degreep1);
    const Tp anp1 = std::pow(a, degreep1);
    return (p.constant / degreep1) * (bnp1 - anp1);
  }

template<typename Tp>
  inline Tp
  f_sin(Tp x)
  {
    return std::sin(x);
  }

template<typename Tp>
  inline Tp
  integ_f_sin(Tp a, Tp b)
  {
    return -std::cos(b) + std::cos(a);
  }

/*
 * Another set of set functions..
 * A new quadrature routine for improper and oscillatory integrals
 * Applied Mathematics and Computation 189 (2007) pp. 452-461
 * E. Sermutlu, H.T. Eyyubo~glu
 */

template<typename Tp>
  const auto s_pi = Tp{3.1415'92653'58979'32384'62643'38327'95028'84195L};

template<typename Tp>
  inline Tp
  cqf1(Tp x)
  {
    return std::exp(x);
  }

template<typename Tp>
  inline Tp
  cqf2(Tp x)
  {
    return x >= 0.3;
  }

template<typename Tp>
  inline Tp
  cqf3(Tp x)
  {
    return std::sqrt(x);
  }


template<typename Tp>
  inline Tp
  cqf4(Tp x)
  {
    return (Tp{23} / Tp{25}) * std::cosh(x) - std::cos(x);
  }

template<typename Tp>
  inline Tp
  cqf5(Tp x)
  {
    Tp x2 = x*x;
    return Tp{1} / (x2 * (x2 + 1) + Tp{0.9});
  }

template<typename Tp>
  inline Tp
  cqf6(Tp x)
  {
    return x * std::sqrt(x);
  }

// Singular at 0.
template<typename Tp>
  inline Tp
  cqf7(Tp x)
  {
    return Tp{1} / std::sqrt(x);
  }

template<typename Tp>
  inline Tp
  cqf8(Tp x)
  {
    Tp x2 = x * x;
    return Tp{1} / (Tp{1} + x2*x2);
  }

template<typename Tp>
  inline Tp
  cqf9(Tp x)
  {
    return Tp{2} / (Tp{2} + std::sin(Tp{10} * s_pi<Tp> * x));
  }

template<typename Tp>
  inline Tp
  cqf10(Tp x)
  {
    return Tp{1} / (Tp{1} + x);
  }

template<typename Tp>
  inline Tp
  cqf11(Tp x)
  {
    return Tp{1} / (Tp{1} + std::exp(x));
  }

template<typename Tp>
  inline Tp
  cqf12(Tp x)
  {
    if (x == Tp{0})
      return Tp{1};
    else
      return x / (std::exp(x) - Tp{1});
  }

template<typename Tp>
  inline Tp
  cqf13(Tp x)
  {
    return std::sin(Tp{100} * s_pi<Tp> * x) / (s_pi<Tp> * x);
  }

template<typename Tp>
  inline Tp
  cqf14(Tp x)
  {
    return std::sqrt(Tp{5}) * std::exp(-Tp{10} * s_pi<Tp> * x * x);
  }

template<typename Tp>
  inline Tp
  cqf15(Tp x)
  {
    return Tp{25} * std::exp(-Tp{25} * x);
  }

template<typename Tp>
  inline Tp
  cqf16(Tp x)
  {
    return Tp{50} / (s_pi<Tp> * (Tp{2500} * x * x + Tp{1}));
  }

template<typename Tp>
  inline Tp
  cqf17(Tp x)
  {
    Tp t1 = Tp{50} * s_pi<Tp> * x, t2;
    t2 = std::sin(t1) / t1;
    return Tp{50} * t2 * t2;
  }

template<typename Tp>
  inline Tp
  cqf18(Tp x)
  {
    return std::cos(std::cos(x)
		  + Tp{3} * std::sin(x)
		  + Tp{2} * std::cos(Tp{2} * x)
		  + Tp{3} * std::cos(Tp{3} * x));
  }

// The paper regularizes this to 0 for x < 1e-15.
template<typename Tp>
  inline Tp
  cqf19(Tp x)
  {
    return std::log(x);
  }

template<typename Tp>
  inline Tp
  cqf20(Tp x)
  {
    return Tp{1} / (x * x + Tp{1.005});
  }

template<typename Tp>
  inline Tp
  cqf21(Tp x)
  {
    return Tp{1} / std::cosh(Tp{20} * x - Tp{4})
         + Tp{1} / std::cosh(Tp{20} * x - Tp{8})
         + Tp{1} / std::cosh(Tp{20} * x - Tp{12});
  }

template<typename Tp>
  inline Tp
  cqf22(Tp x)
  {
    return Tp{4} * s_pi<Tp> * s_pi<Tp> * x
	  * std::sin(Tp{20} * s_pi<Tp> * x)
	  * std::cos(Tp{2} * s_pi<Tp> * x);
  }

template<typename Tp>
  inline Tp
  cqf23(Tp x)
  {
    Tp t = Tp{230} * x - Tp{30};
    return Tp{1} / (Tp{1} + t * t);
  }

template<typename Tp>
  inline Tp
  cqf24(Tp x)
  {
    return std::floor(std::exp(x));
  }

template<typename Tp>
  inline Tp
  cqf25(Tp x)
  {
    return (x < Tp{1}) * (x + Tp{1})
         + (Tp{1} <= x && x <= Tp{3}) * (Tp{3} - x)
         + (x > Tp{3}) * Tp{2};
  }

// The limits of integration and the expected results fot the 
template<typename Tp>
  const static func_test<Tp>
  func_tests[25]
  {
    &cqf1,  Tp{0},  Tp{1},	 Tp{1.718281828459045235360287471352662497759L},
    &cqf2,  Tp{0},  Tp{1},	 Tp{0.7L},
    &cqf3,  Tp{0},  Tp{1},	 Tp{2.0L} / Tp{3.0L},
    &cqf4,  Tp{-1}, Tp{1},	 Tp{0.4794282266888016673585779618353075006421L},
    &cqf5,  Tp{-1}, Tp{1},	 Tp{1.582232963729672933117468949026169067924L},
    &cqf6,  Tp{0},  Tp{1},	 Tp{0.4L},
    &cqf7,  Tp{0},  Tp{1},	 Tp{2.0L},
    &cqf8,  Tp{0},  Tp{1},	 Tp{0.8669729873399110375739951638828707136522L},
    &cqf9,  Tp{0},  Tp{1},	 Tp{1.154700538379251529018297561003914911295L},
    &cqf10, Tp{0},  Tp{1},	 Tp{0.6931471805599453094172321214581765680755L},
    &cqf11, Tp{0},  Tp{1},	 Tp{0.3798854930417224753682366264903209261602L},
    &cqf12, Tp{0},  Tp{1},	 Tp{0.7775046341122482764175865454257105071925L},
    &cqf13, Tp{0.1},Tp{1},	 Tp{0.9098637539166842915557830641141434835684e-2L},
    &cqf14, Tp{0},  Tp{10},	 Tp{0.3535533905932737622004221810524245196424L},
    &cqf15, Tp{0},  Tp{10},	 Tp{1.0L},
    &cqf16, Tp{0},  Tp{10},	 Tp{0.4993633810764567446362485183117640508837L},
    &cqf17, Tp{0},  Tp{1},	 Tp{0.4989868086930455024989853136560590525530L},
    &cqf18, Tp{0},  s_pi<Tp>, Tp{0.2910187828600526985238845968602020712077L},
    &cqf19, Tp{0},  Tp{1},	 Tp{-1.0L},
    &cqf20, Tp{-1}, Tp{1},	 Tp{1.564396444069049773091493015808472813088L},
    &cqf21, Tp{0},  Tp{1},	 Tp{0.4693392062964089464236343930140556287142L},
    &cqf22, Tp{0},  Tp{1},	 Tp{-0.6346651825433925734267966430867682594338L},
    &cqf23, Tp{0},  Tp{1},	 Tp{0.1349248564946777269188547624864782167827e-1L},
    &cqf24, Tp{0},  Tp{3},	 Tp{17.66438353924651497034012402929007814264L},
    &cqf25, Tp{0},  Tp{5},	 Tp{7.5L}
  };

#endif // QUADRATURE_TESTCASE_TCC
