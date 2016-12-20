/* integration/qng_integrate.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * Copyright (C) 2016 Edward Smith-Rowland
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

/**
   Gauss-Kronrod-Patterson quadrature coefficients for use in
   quadpack routine qng. These coefficients were calculated with
   101 decimal digit arithmetic by L. W. Fullerton, Bell Labs, Nov
   1981.
 */

namespace __gnu_test
{

/* x1, abscissae common to the 10-, 21-, 43- and 87-point rule */
static const long double
x1[5]
{
  0.973906528517171720077964012084452L,
  0.865063366688984510732096688423493L,
  0.679409568299024406234327365114874L,
  0.433395394129247190799265943165784L,
  0.148874338981631210884826001129720L,
};

/* w10, weights of the 10-point formula */
static const long double
w10[5]
{
  0.066671344308688137593568809893332L,
  0.149451349150580593145776339657697L,
  0.219086362515982043995534934228163L,
  0.269266719309996355091226921569469L,
  0.295524224714752870173892994651338L,
};

/* x2, abscissae common to the 21-, 43- and 87-point rule */
static const long double
x2[5]
{
  0.995657163025808080735527280689003L,
  0.930157491355708226001207180059508L,
  0.780817726586416897063717578345042L,
  0.562757134668604683339000099272694L,
  0.294392862701460198131126603103866L,
};

/* w21a, weights of the 21-point formula for abscissae x1 */
static const long double
w21a[5]
{
  0.032558162307964727478818972459390L,
  0.075039674810919952767043140916190L,
  0.109387158802297641899210590325805L,
  0.134709217311473325928054001771707L,
  0.147739104901338491374841515972068L,
};

/* w21b, weights of the 21-point formula for abscissae x2 */
static const long double
w21b[6]
{
  0.011694638867371874278064396062192L,
  0.054755896574351996031381300244580L,
  0.093125454583697605535065465083366L,
  0.123491976262065851077958109831074L,
  0.142775938577060080797094273138717L,
  0.149445554002916905664936468389821L,
};

/* x3, abscissae common to the 43- and 87-point rule */
static const long double
x3[11]
{
  0.999333360901932081394099323919911L,
  0.987433402908088869795961478381209L,
  0.954807934814266299257919200290473L,
  0.900148695748328293625099494069092L,
  0.825198314983114150847066732588520L,
  0.732148388989304982612354848755461L,
  0.622847970537725238641159120344323L,
  0.499479574071056499952214885499755L,
  0.364901661346580768043989548502644L,
  0.222254919776601296498260928066212L,
  0.074650617461383322043914435796506L,
};

/* w43a, weights of the 43-point formula for abscissae x1, x3 */
static const long double
w43a[10]
{
  0.016296734289666564924281974617663L,
  0.037522876120869501461613795898115L,
  0.054694902058255442147212685465005L,
  0.067355414609478086075553166302174L,
  0.073870199632393953432140695251367L,
  0.005768556059769796184184327908655L,
  0.027371890593248842081276069289151L,
  0.046560826910428830743339154433824L,
  0.061744995201442564496240336030883L,
  0.071387267268693397768559114425516L,
};

/* w43b, weights of the 43-point formula for abscissae x3 */
static const long double
w43b[12]
{
  0.001844477640212414100389106552965L,
  0.010798689585891651740465406741293L,
  0.021895363867795428102523123075149L,
  0.032597463975345689443882222526137L,
  0.042163137935191811847627924327955L,
  0.050741939600184577780189020092084L,
  0.058379395542619248375475369330206L,
  0.064746404951445885544689259517511L,
  0.069566197912356484528633315038405L,
  0.072824441471833208150939535192842L,
  0.074507751014175118273571813842889L,
  0.074722147517403005594425168280423L,
};

/* x4, abscissae of the 87-point rule */
static const long double
x4[22]
{
  0.999902977262729234490529830591582L,
  0.997989895986678745427496322365960L,
  0.992175497860687222808523352251425L,
  0.981358163572712773571916941623894L,
  0.965057623858384619128284110607926L,
  0.943167613133670596816416634507426L,
  0.915806414685507209591826430720050L,
  0.883221657771316501372117548744163L,
  0.845710748462415666605902011504855L,
  0.803557658035230982788739474980964L,
  0.757005730685495558328942793432020L,
  0.706273209787321819824094274740840L,
  0.651589466501177922534422205016736L,
  0.593223374057961088875273770349144L,
  0.531493605970831932285268948562671L,
  0.466763623042022844871966781659270L,
  0.399424847859218804732101665817923L,
  0.329874877106188288265053371824597L,
  0.258503559202161551802280975429025L,
  0.185695396568346652015917141167606L,
  0.111842213179907468172398359241362L,
  0.037352123394619870814998165437704L,
};

/* w87a, weights of the 87-point formula for abscissae x1, x2, x3 */
static const long double
w87a[21]
{
  0.008148377384149172900002878448190L,
  0.018761438201562822243935059003794L,
  0.027347451050052286161582829741283L,
  0.033677707311637930046581056957588L,
  0.036935099820427907614589586742499L,
  0.002884872430211530501334156248695L,
  0.013685946022712701888950035273128L,
  0.023280413502888311123409291030404L,
  0.030872497611713358675466394126442L,
  0.035693633639418770719351355457044L,
  0.000915283345202241360843392549948L,
  0.005399280219300471367738743391053L,
  0.010947679601118931134327826856808L,
  0.016298731696787335262665703223280L,
  0.021081568889203835112433060188190L,
  0.025370969769253827243467999831710L,
  0.029189697756475752501446154084920L,
  0.032373202467202789685788194889595L,
  0.034783098950365142750781997949596L,
  0.036412220731351787562801163687577L,
  0.037253875503047708539592001191226L,
};

/* w87b, weights of the 87-point formula for abscissae x4 */
static const long double
w87b[23]
{
  0.000274145563762072350016527092881L,
  0.001807124155057942948341311753254L,
  0.004096869282759164864458070683480L,
  0.006758290051847378699816577897424L,
  0.009549957672201646536053581325377L,
  0.012329447652244853694626639963780L,
  0.015010447346388952376697286041943L,
  0.017548967986243191099665352925900L,
  0.019938037786440888202278192730714L,
  0.022194935961012286796332102959499L,
  0.024339147126000805470360647041454L,
  0.026374505414839207241503786552615L,
  0.028286910788771200659968002987960L,
  0.030052581128092695322521110347341L,
  0.031646751371439929404586051078883L,
  0.033050413419978503290785944862689L,
  0.034255099704226061787082821046821L,
  0.035262412660156681033782717998428L,
  0.036076989622888701185500318003895L,
  0.036698604498456094498018047441094L,
  0.037120549269832576114119958413599L,
  0.037334228751935040321235449094698L,
  0.037361073762679023410321241766599L,
};

template<typename _Tp>
  _Tp
  rescale_error(_Tp err, const _Tp result_abs, const _Tp result_asc)
  {
    err = std::abs(err);

    if (result_asc != 0 && err != 0)
      {
	auto scale = std::pow((200 * err / result_asc), 1.5);

	if (scale < 1)
          err = result_asc * scale;
	else 
          err = result_asc;
      }
    if (result_abs > std::numeric_limits<_Tp>::min() / (50 * std::numeric_limits<_Tp>::epsilon()))
      {
	auto min_err = 50 * std::numeric_limits<_Tp>::epsilon() * result_abs;

	if (min_err > err) 
          err = min_err;
      }

    return err;
  }

template<typename _FuncTp, typename _Tp>
  int
  qng_integrate(const _FuncTp& __func,
        	_Tp a, _Tp b,
        	_Tp epsabs, _Tp epsrel,
        	_Tp& result, _Tp& abserr, size_t& neval)
  {
    _Tp fv1[5], fv2[5], fv3[5], fv4[5];
    _Tp savfun[21];  /* array of function values which have been computed */
    _Tp res10, res21, res43, res87; /* 10, 21, 43 and 87 point results */
    _Tp result_kronrod, err; 
    _Tp resabs; /* approximation to the integral of abs(f) */
    _Tp resasc; /* approximation to the integral of abs(f-i/(b-a)) */

    const _Tp half_length = 0.5 * (b - a);
    const _Tp abs_half_length = std::abs(half_length);
    const _Tp center = 0.5 * (b + a);
    const _Tp f_center = __func(center);

    if (epsabs <= 0 && (epsrel < 50 * std::numeric_limits<_Tp>::epsilon() || epsrel < 0.5e-28))
      {
	result = 0;
	abserr = 0;
	neval = 0;
	std::__throw_runtime_error("tolerance cannot be achieved with given epsabs and epsrel");
      };

    /* Compute the integral using the 10- and 21-point formula. */

    res10 = 0;
    res21 = _Tp(w21b[5]) * f_center;
    resabs = _Tp(w21b[5]) * std::abs(f_center);

    for (int k = 0; k < 5; ++k)
      {
	const _Tp abscissa = half_length * _Tp(x1[k]);
	const _Tp fval1 = __func(center + abscissa);
	const _Tp fval2 = __func(center - abscissa);
	const _Tp fval = fval1 + fval2;
	res10 += _Tp(w10[k]) * fval;
	res21 += _Tp(w21a[k]) * fval;
	resabs += _Tp(w21a[k]) * (std::abs(fval1) + std::abs(fval2));
	savfun[k] = fval;
	fv1[k] = fval1;
	fv2[k] = fval2;
      }

    for (int k = 0; k < 5; ++k)
      {
	const _Tp abscissa = half_length * _Tp(x2[k]);
	const _Tp fval1 = __func(center + abscissa);
	const _Tp fval2 = __func(center - abscissa);
	const _Tp fval = fval1 + fval2;
	res21 += _Tp(w21b[k]) * fval;
	resabs += _Tp(w21b[k]) * (std::abs(fval1) + std::abs(fval2));
	savfun[k + 5] = fval;
	fv3[k] = fval1;
	fv4[k] = fval2;
      }

    resabs *= abs_half_length;

    { 
      const _Tp mean = 0.5 * res21;

      resasc = _Tp(w21b[5]) * std::abs(f_center - mean);

      for (int k = 0; k < 5; ++k)
	{
          resasc +=
            (_Tp(w21a[k]) * (std::abs(fv1[k] - mean) + std::abs(fv2[k] - mean))
            + _Tp(w21b[k]) * (std::abs(fv3[k] - mean) + std::abs(fv4[k] - mean)));
	}
      resasc *= abs_half_length;
    }

    result_kronrod = res21 * half_length;

    err = rescale_error((res21 - res10) * half_length, resabs, resasc);

    /*   Test for convergence. */

    if (err < epsabs || err < epsrel * std::abs(result_kronrod))
      {
	result = result_kronrod;
	abserr = err;
	neval = 21;
	return 0;
      }

    /* Compute the integral using the 43-point formula. */

    res43 = _Tp(w43b[11]) * f_center;

    for (int k = 0; k < 10; ++k)
      res43 += savfun[k] * _Tp(w43a[k]);

    for (int k = 0; k < 11; ++k)
      {
	const _Tp abscissa = half_length * _Tp(x3[k]);
	const _Tp fval = __func(center + abscissa) 
                       + __func(center - abscissa);
	res43 += fval * _Tp(w43b[k]);
	savfun[k + 10] = fval;
      }

    /* Test for convergence. */

    result_kronrod = res43 * half_length;
    err = rescale_error ((res43 - res21) * half_length, resabs, resasc);

    if (err < epsabs || err < epsrel * std::abs(result_kronrod))
      {
	result = result_kronrod;
	abserr = err;
	neval = 43;
	return 0;
      }

    /* Compute the integral using the 87-point formula. */

    res87 = _Tp(w87b[22]) * f_center;

    for (int k = 0; k < 21; ++k)
      res87 += savfun[k] * _Tp(w87a[k]);

    for (int k = 0; k < 22; ++k)
      {
	const _Tp abscissa = half_length * _Tp(x4[k]);
	res87 += _Tp(w87b[k]) * (__func(center + abscissa)
			       + __func(center - abscissa));
      }

    /* Test for convergence */

    result_kronrod = res87 * half_length;

    err = rescale_error ((res87 - res43) * half_length, resabs, resasc);

    if (err < epsabs || err < epsrel * std::abs(result_kronrod))
      {
	result = result_kronrod;
	abserr = err;
	neval = 87;
	return 0;
      }

    /* Failed to converge */

    result = result_kronrod;
    abserr = err;
    neval = 87;

    std::__throw_runtime_error("failed to reach tolerance with highest-order rule");
  }

} // namespace __gnu_test