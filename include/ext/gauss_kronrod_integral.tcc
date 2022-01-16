// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2011-2020 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.
//
// Ported from GSL by Edward Smith-Rowland
// Originally written by Brian Gaugh
//
// Implements Gauss-Kronrod integration
// Based on gsl/integration/qk.c

#ifndef GAUSS_KRONROD_INTERGAL_TCC
#define GAUSS_KRONROD_INTERGAL_TCC 1

#include <type_traits>
#include <vector>
#include <cmath>
#include <array>
#include <stdexcept>

#include <ext/integration_error.h>
#include <ext/gauss_kronrod_rule.tcc>

namespace __gnu_cxx
{

  template<typename Tp>
    gauss_kronrod_integral<Tp>::gauss_kronrod_integral(unsigned gk_rule)
    : m_rule{gk_rule},
      m_x_kronrod{}, m_w_gauss{}, m_w_kronrod{}
    {
      const int n = (this->m_rule - 1) / 2;
      const auto eps = 4 * std::numeric_limits<Tp>::epsilon();
      build_gauss_kronrod(n, eps,
			    this->m_x_kronrod, this->m_w_gauss,
			    this->m_w_kronrod);
    }

  template<typename AreaTp, typename AbsAreaTp>
    inline bool
    test_positivity(AreaTp result, AbsAreaTp resabs)
    {
      return(std::abs(result)
	  >= (1 - 50 * std::numeric_limits<AbsAreaTp>::epsilon()) * resabs);
    }

  // Class template for internal implementation of
  // each individual integration rule.
  template<typename Tp, typename FuncTp, Kronrod_Rule _GK_rule>
    class qk_integrator;

  /**
   * Gauss-Kronrod 15-point rule implementation.
   */
  template<typename Tp, typename FuncTp>
    struct qk_integrator<Tp, FuncTp, Kronrod_15>
    {
      // Abscissae of the 15-point Kronrod rule
      static constexpr std::array<Tp, 8>
      s_x_kronrod =
      {
	Tp{0.991455371120812639206854697526329L},
	Tp{0.949107912342758524526189684047851L},
	Tp{0.864864423359769072789712788640926L},
	Tp{0.741531185599394439863864773280788L},
	Tp{0.586087235467691130294144838258730L},
	Tp{0.405845151377397166906606412076961L},
	Tp{0.207784955007898467600689403773245L},
	Tp{0.000000000000000000000000000000000L}
      };
      // Weights of the 7-point Gauss rule
      static constexpr std::array<Tp, 4>
      s_w_gauss =
      {
	Tp{0.129484966168869693270611432679082L},
	Tp{0.279705391489276667901467771423780L},
	Tp{0.381830050505118944950369775488975L},
	Tp{0.417959183673469387755102040816327L}
      };
      // Weights of the 15-point Kronrod rule
      static constexpr std::array<Tp, 8>
      s_w_kronrod =
      {
	Tp{0.022935322010529224963732008058970L},
	Tp{0.063092092629978553290700663189204L},
	Tp{0.104790010322250183839876322541518L},
	Tp{0.140653259715525918745189590510238L},
	Tp{0.169004726639267902826583426598550L},
	Tp{0.190350578064785409913256402421014L},
	Tp{0.204432940075298892414161999234649L},
	Tp{0.209482141084727828012999174891714L}
      };
    };

  /**
   * Gauss-Kronrod 21-point rule implementation.
   */
  template<typename Tp, typename FuncTp>
    struct qk_integrator<Tp, FuncTp, Kronrod_21>
    {
      // Abscissae of the 21-point Kronrod rule
      static constexpr std::array<Tp, 11>
      s_x_kronrod =
      {
	Tp{0.995657163025808080735527280689003L},
	Tp{0.973906528517171720077964012084452L},
	Tp{0.930157491355708226001207180059508L},
	Tp{0.865063366688984510732096688423493L},
	Tp{0.780817726586416897063717578345042L},
	Tp{0.679409568299024406234327365114874L},
	Tp{0.562757134668604683339000099272694L},
	Tp{0.433395394129247190799265943165784L},
	Tp{0.294392862701460198131126603103866L},
	Tp{0.148874338981631210884826001129720L},
	Tp{0.000000000000000000000000000000000L}
      };
      // Weights of the 10-point Gauss rule
      static constexpr std::array<Tp, 5>
      s_w_gauss =
      {
	Tp{0.066671344308688137593568809893332L},
	Tp{0.149451349150580593145776339657697L},
	Tp{0.219086362515982043995534934228163L},
	Tp{0.269266719309996355091226921569469L},
	Tp{0.295524224714752870173892994651338L}
      };
      // Weights of the 21-point Kronrod rule
      static constexpr std::array<Tp, 11>
      s_w_kronrod =
      {
	Tp{0.011694638867371874278064396062192L},
	Tp{0.032558162307964727478818972459390L},
	Tp{0.054755896574351996031381300244580L},
	Tp{0.075039674810919952767043140916190L},
	Tp{0.093125454583697605535065465083366L},
	Tp{0.109387158802297641899210590325805L},
	Tp{0.123491976262065851077958109831074L},
	Tp{0.134709217311473325928054001771707L},
	Tp{0.142775938577060080797094273138717L},
	Tp{0.147739104901338491374841515972068L},
	Tp{0.149445554002916905664936468389821L}
      };
    };

  /**
   * Gauss-Kronrod 31-point rule implementation.
   */
  template<typename Tp, typename FuncTp>
    struct qk_integrator<Tp, FuncTp, Kronrod_31>
    {
      // Abscissae of the 31-point Kronrod rule
      static constexpr std::array<Tp, 16>
      s_x_kronrod =
      {
	Tp{0.998002298693397060285172840152271L},
	Tp{0.987992518020485428489565718586613L},
	Tp{0.967739075679139134257347978784337L},
	Tp{0.937273392400705904307758947710209L},
	Tp{0.897264532344081900882509656454496L},
	Tp{0.848206583410427216200648320774217L},
	Tp{0.790418501442465932967649294817947L},
	Tp{0.724417731360170047416186054613938L},
	Tp{0.650996741297416970533735895313275L},
	Tp{0.570972172608538847537226737253911L},
	Tp{0.485081863640239680693655740232351L},
	Tp{0.394151347077563369897207370981045L},
	Tp{0.299180007153168812166780024266389L},
	Tp{0.201194093997434522300628303394596L},
	Tp{0.101142066918717499027074231447392L},
	Tp{0.000000000000000000000000000000000L}
      };
      // Weights of the 15-point Gauss rule
      static constexpr std::array<Tp, 8>
      s_w_gauss =
      {
	Tp{0.030753241996117268354628393577204L},
	Tp{0.070366047488108124709267416450667L},
	Tp{0.107159220467171935011869546685869L},
	Tp{0.139570677926154314447804794511028L},
	Tp{0.166269205816993933553200860481209L},
	Tp{0.186161000015562211026800561866423L},
	Tp{0.198431485327111576456118326443839L},
	Tp{0.202578241925561272880620199967519L}
      };
      // Weights of the 31-point Kronrod rule
      static constexpr std::array<Tp, 16>
      s_w_kronrod =
      {
	Tp{0.005377479872923348987792051430128L},
	Tp{0.015007947329316122538374763075807L},
	Tp{0.025460847326715320186874001019653L},
	Tp{0.035346360791375846222037948478360L},
	Tp{0.044589751324764876608227299373280L},
	Tp{0.053481524690928087265343147239430L},
	Tp{0.062009567800670640285139230960803L},
	Tp{0.069854121318728258709520077099147L},
	Tp{0.076849680757720378894432777482659L},
	Tp{0.083080502823133021038289247286104L},
	Tp{0.088564443056211770647275443693774L},
	Tp{0.093126598170825321225486872747346L},
	Tp{0.096642726983623678505179907627589L},
	Tp{0.099173598721791959332393173484603L},
	Tp{0.100769845523875595044946662617570L},
	Tp{0.101330007014791549017374792767493L}
      };
    };

  /**
   * Gauss-Kronrod 41-point rule implementation.
   */
  template<typename Tp, typename FuncTp>
    struct qk_integrator<Tp, FuncTp, Kronrod_41>
    {
      // Abscissae of the 41-point Kronrod rule
      static constexpr std::array<Tp, 21>
      s_x_kronrod =
      {
	Tp{0.998859031588277663838315576545863L},
	Tp{0.993128599185094924786122388471320L},
	Tp{0.981507877450250259193342994720217L},
	Tp{0.963971927277913791267666131197277L},
	Tp{0.940822633831754753519982722212443L},
	Tp{0.912234428251325905867752441203298L},
	Tp{0.878276811252281976077442995113078L},
	Tp{0.839116971822218823394529061701521L},
	Tp{0.795041428837551198350638833272788L},
	Tp{0.746331906460150792614305070355642L},
	Tp{0.693237656334751384805490711845932L},
	Tp{0.636053680726515025452836696226286L},
	Tp{0.575140446819710315342946036586425L},
	Tp{0.510867001950827098004364050955251L},
	Tp{0.443593175238725103199992213492640L},
	Tp{0.373706088715419560672548177024927L},
	Tp{0.301627868114913004320555356858592L},
	Tp{0.227785851141645078080496195368575L},
	Tp{0.152605465240922675505220241022678L},
	Tp{0.076526521133497333754640409398838L},
	Tp{0.000000000000000000000000000000000L}
      };
      // Weights of the 20-point Gauss rule
      static constexpr std::array<Tp, 10>
      s_w_gauss =
      {
	Tp{0.017614007139152118311861962351853L},
	Tp{0.040601429800386941331039952274932L},
	Tp{0.062672048334109063569506535187042L},
	Tp{0.083276741576704748724758143222046L},
	Tp{0.101930119817240435036750135480350L},
	Tp{0.118194531961518417312377377711382L},
	Tp{0.131688638449176626898494499748163L},
	Tp{0.142096109318382051329298325067165L},
	Tp{0.149172986472603746787828737001969L},
	Tp{0.152753387130725850698084331955098L}
      };
      // Weights of the 41-point Kronrod rule
      static constexpr std::array<Tp, 21>
      s_w_kronrod =
      {
	Tp{0.003073583718520531501218293246031L},
	Tp{0.008600269855642942198661787950102L},
	Tp{0.014626169256971252983787960308868L},
	Tp{0.020388373461266523598010231432755L},
	Tp{0.025882133604951158834505067096153L},
	Tp{0.031287306777032798958543119323801L},
	Tp{0.036600169758200798030557240707211L},
	Tp{0.041668873327973686263788305936895L},
	Tp{0.046434821867497674720231880926108L},
	Tp{0.050944573923728691932707670050345L},
	Tp{0.055195105348285994744832372419777L},
	Tp{0.059111400880639572374967220648594L},
	Tp{0.062653237554781168025870122174255L},
	Tp{0.065834597133618422111563556969398L},
	Tp{0.068648672928521619345623411885368L},
	Tp{0.071054423553444068305790361723210L},
	Tp{0.073030690332786667495189417658913L},
	Tp{0.074582875400499188986581418362488L},
	Tp{0.075704497684556674659542775376617L},
	Tp{0.076377867672080736705502835038061L},
	Tp{0.076600711917999656445049901530102L}
      };
    };

  /**
   * Gauss-Kronrod 51-point rule implementation.
   */
  template<typename Tp, typename FuncTp>
    struct qk_integrator<Tp, FuncTp, Kronrod_51>
    {
      // Abscissae of the 51-point Kronrod rule
      static constexpr std::array<Tp, 26>
      s_x_kronrod =
      {
	Tp{0.999262104992609834193457486540341L},
	Tp{0.995556969790498097908784946893902L},
	Tp{0.988035794534077247637331014577406L},
	Tp{0.976663921459517511498315386479594L},
	Tp{0.961614986425842512418130033660167L},
	Tp{0.942974571228974339414011169658471L},
	Tp{0.920747115281701561746346084546331L},
	Tp{0.894991997878275368851042006782805L},
	Tp{0.865847065293275595448996969588340L},
	Tp{0.833442628760834001421021108693570L},
	Tp{0.797873797998500059410410904994307L},
	Tp{0.759259263037357630577282865204361L},
	Tp{0.717766406813084388186654079773298L},
	Tp{0.673566368473468364485120633247622L},
	Tp{0.626810099010317412788122681624518L},
	Tp{0.577662930241222967723689841612654L},
	Tp{0.526325284334719182599623778158010L},
	Tp{0.473002731445714960522182115009192L},
	Tp{0.417885382193037748851814394594572L},
	Tp{0.361172305809387837735821730127641L},
	Tp{0.303089538931107830167478909980339L},
	Tp{0.243866883720988432045190362797452L},
	Tp{0.183718939421048892015969888759528L},
	Tp{0.122864692610710396387359818808037L},
	Tp{0.061544483005685078886546392366797L},
	Tp{0.000000000000000000000000000000000L}
      };
      // Weights of the 25-point Gauss rule
      static constexpr std::array<Tp, 13>
      s_w_gauss =
      {
	Tp{0.011393798501026287947902964113235L},
	Tp{0.026354986615032137261901815295299L},
	Tp{0.040939156701306312655623487711646L},
	Tp{0.054904695975835191925936891540473L},
	Tp{0.068038333812356917207187185656708L},
	Tp{0.080140700335001018013234959669111L},
	Tp{0.091028261982963649811497220702892L},
	Tp{0.100535949067050644202206890392686L},
	Tp{0.108519624474263653116093957050117L},
	Tp{0.114858259145711648339325545869556L},
	Tp{0.119455763535784772228178126512901L},
	Tp{0.122242442990310041688959518945852L},
	Tp{0.123176053726715451203902873079050L}
      };
      // Weights of the 51-point Kronrod rule
      static constexpr std::array<Tp, 26>
      s_w_kronrod =
      {
	Tp{0.001987383892330315926507851882843L},
	Tp{0.005561932135356713758040236901066L},
	Tp{0.009473973386174151607207710523655L},
	Tp{0.013236229195571674813656405846976L},
	Tp{0.016847817709128298231516667536336L},
	Tp{0.020435371145882835456568292235939L},
	Tp{0.024009945606953216220092489164881L},
	Tp{0.027475317587851737802948455517811L},
	Tp{0.030792300167387488891109020215229L},
	Tp{0.034002130274329337836748795229551L},
	Tp{0.037116271483415543560330625367620L},
	Tp{0.040083825504032382074839284467076L},
	Tp{0.042872845020170049476895792439495L},
	Tp{0.045502913049921788909870584752660L},
	Tp{0.047982537138836713906392255756915L},
	Tp{0.050277679080715671963325259433440L},
	Tp{0.052362885806407475864366712137873L},
	Tp{0.054251129888545490144543370459876L},
	Tp{0.055950811220412317308240686382747L},
	Tp{0.057437116361567832853582693939506L},
	Tp{0.058689680022394207961974175856788L},
	Tp{0.059720340324174059979099291932562L},
	Tp{0.060539455376045862945360267517565L},
	Tp{0.061128509717053048305859030416293L},
	Tp{0.061471189871425316661544131965264L},
	Tp{0.061580818067832935078759824240066L}
      };
    };

  /**
   * Gauss-Kronrod 51-point rule implementation.
   */
  template<typename Tp, typename FuncTp>
    struct qk_integrator<Tp, FuncTp, Kronrod_61>
    {
      // Abscissae of the 61-point Kronrod rule
      static constexpr std::array<Tp, 31>
      s_x_kronrod =
      {
	Tp{0.999484410050490637571325895705811L},
	Tp{0.996893484074649540271630050918695L},
	Tp{0.991630996870404594858628366109486L},
	Tp{0.983668123279747209970032581605663L},
	Tp{0.973116322501126268374693868423707L},
	Tp{0.960021864968307512216871025581798L},
	Tp{0.944374444748559979415831324037439L},
	Tp{0.926200047429274325879324277080474L},
	Tp{0.905573307699907798546522558925958L},
	Tp{0.882560535792052681543116462530226L},
	Tp{0.857205233546061098958658510658944L},
	Tp{0.829565762382768397442898119732502L},
	Tp{0.799727835821839083013668942322683L},
	Tp{0.767777432104826194917977340974503L},
	Tp{0.733790062453226804726171131369528L},
	Tp{0.697850494793315796932292388026640L},
	Tp{0.660061064126626961370053668149271L},
	Tp{0.620526182989242861140477556431189L},
	Tp{0.579345235826361691756024932172540L},
	Tp{0.536624148142019899264169793311073L},
	Tp{0.492480467861778574993693061207709L},
	Tp{0.447033769538089176780609900322854L},
	Tp{0.400401254830394392535476211542661L},
	Tp{0.352704725530878113471037207089374L},
	Tp{0.304073202273625077372677107199257L},
	Tp{0.254636926167889846439805129817805L},
	Tp{0.204525116682309891438957671002025L},
	Tp{0.153869913608583546963794672743256L},
	Tp{0.102806937966737030147096751318001L},
	Tp{0.051471842555317695833025213166723L},
	Tp{0.000000000000000000000000000000000L}
      };
      // Weights of the 30-point Gauss rule
      static constexpr std::array<Tp, 15>
      s_w_gauss =
      {
	Tp{0.007968192496166605615465883474674L},
	Tp{0.018466468311090959142302131912047L},
	Tp{0.028784707883323369349719179611292L},
	Tp{0.038799192569627049596801936446348L},
	Tp{0.048402672830594052902938140422808L},
	Tp{0.057493156217619066481721689402056L},
	Tp{0.065974229882180495128128515115962L},
	Tp{0.073755974737705206268243850022191L},
	Tp{0.080755895229420215354694938460530L},
	Tp{0.086899787201082979802387530715126L},
	Tp{0.092122522237786128717632707087619L},
	Tp{0.096368737174644259639468626351810L},
	Tp{0.099593420586795267062780282103569L},
	Tp{0.101762389748405504596428952168554L},
	Tp{0.102852652893558840341285636705415L}
      };
      // Weights of the 61-point Kronrod rule
      static constexpr std::array<Tp, 31>
      s_w_kronrod =
      {
	Tp{0.001389013698677007624551591226760L},
	Tp{0.003890461127099884051267201844516L},
	Tp{0.006630703915931292173319826369750L},
	Tp{0.009273279659517763428441146892024L},
	Tp{0.011823015253496341742232898853251L},
	Tp{0.014369729507045804812451432443580L},
	Tp{0.016920889189053272627572289420322L},
	Tp{0.019414141193942381173408951050128L},
	Tp{0.021828035821609192297167485738339L},
	Tp{0.024191162078080601365686370725232L},
	Tp{0.026509954882333101610601709335075L},
	Tp{0.028754048765041292843978785354334L},
	Tp{0.030907257562387762472884252943092L},
	Tp{0.032981447057483726031814191016854L},
	Tp{0.034979338028060024137499670731468L},
	Tp{0.036882364651821229223911065617136L},
	Tp{0.038678945624727592950348651532281L},
	Tp{0.040374538951535959111995279752468L},
	Tp{0.041969810215164246147147541285970L},
	Tp{0.043452539701356069316831728117073L},
	Tp{0.044814800133162663192355551616723L},
	Tp{0.046059238271006988116271735559374L},
	Tp{0.047185546569299153945261478181099L},
	Tp{0.048185861757087129140779492298305L},
	Tp{0.049055434555029778887528165367238L},
	Tp{0.049795683427074206357811569379942L},
	Tp{0.050405921402782346840893085653585L},
	Tp{0.050881795898749606492297473049805L},
	Tp{0.051221547849258772170656282604944L},
	Tp{0.051426128537459025933862879215781L},
	Tp{0.051494729429451567558340433647099L}
      };
    };

  // Integrates func from a to b using integration rule qkintrule
  // returns a tuple with the results of a single Gauss-Kronrod integration
  // Based on GSL function gsl_integration_qk()
  // Return tuple slots are as follows:
  // 0: result - result of integration using Kronrod scheme
  // 1: abserr - Estimated error as difference between Gauss and Kronrod
  // 2: resabs - Integral of absolute value of function
  // 3: resasc - Integral of absolute value of difference between function
  //             and weighted mean function value
  template<typename Tp, typename FuncTp>
    auto
    qk_integrate(FuncTp func, Tp lower, Tp upper,
		 Kronrod_Rule qkintrule)
    -> gauss_kronrod_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
    {
      gauss_kronrod_integral<Tp> gk_integ(qkintrule);
      return gk_integ.integrate(func, lower, upper);
    }

  template<typename Tp>
    template<typename FuncTp,
	     typename KronrodIter, typename GaussIter>
      auto
      gauss_kronrod_integral<Tp>::
      s_integrate(const KronrodIter& x_kronrod,
		   const GaussIter& w_gauss,
		   const KronrodIter& w_kronrod,
		   FuncTp func, Tp lower, Tp upper)
      -> gauss_kronrod_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
      {
	using RetTp = std::invoke_result_t<FuncTp, Tp>;
	using AreaTp = decltype(RetTp{} * Tp{});

	const auto KronrodSz = std::size(x_kronrod);
	//const auto GaussSz = std::size(w_gauss);
	//static_assert(KronrodSz == 2 * GaussSz + (KronrodSz & 1));
	//std::array<AreaTp, KronrodSz> fv1;
	//std::array<AreaTp, KronrodSz> fv2;
	std::vector<AreaTp> fv1(KronrodSz);
	std::vector<AreaTp> fv2(KronrodSz);

	const auto center = (lower + upper) / Tp{2};
	const auto half_length = (upper - lower) / Tp{2};
	const auto abs_half_length = std::abs(half_length);
	const auto f_center = func(center);

	auto result_gauss = AreaTp{0};
	auto result_kronrod = f_center * w_kronrod[KronrodSz - 1];
	auto result_abs = std::abs(result_kronrod);

	if (KronrodSz % 2 == 0)
	  result_gauss = f_center * w_gauss[KronrodSz / 2 - 1];

	for (std::size_t jj = 0; jj < (KronrodSz - 1) / 2; ++jj)
	  {
	    const std::size_t jtw = jj * 2 + 1;
	    const auto abscissa = half_length * x_kronrod[jtw];
	    const auto fval1 = func(center - abscissa);
	    const auto fval2 = func(center + abscissa);
	    const auto fsum = fval1 + fval2;
	    fv1[jtw] = fval1;
	    fv2[jtw] = fval2;

	    result_gauss += w_gauss[jj] * fsum;
	    result_kronrod += w_kronrod[jtw] * fsum;
	    result_abs += w_kronrod[jtw]
			  * (std::abs(fval1) + std::abs(fval2));
	  }

	for (std::size_t jj = 0; jj < KronrodSz / 2; ++jj)
	  {
	    std::size_t jtwm1 = jj * 2;
	    const auto abscissa = half_length * x_kronrod[jtwm1];
	    const auto fval1 = func(center - abscissa);
	    const auto fval2 = func(center + abscissa);
	    fv1[jtwm1] = fval1;
	    fv2[jtwm1] = fval2;

	    result_kronrod += w_kronrod[jtwm1] * (fval1 + fval2);
	    result_abs += w_kronrod[jtwm1]
			  * (std::abs(fval1) + std::abs(fval2));
	  }

	auto mean = result_kronrod / Tp{2};
	auto result_asc = w_kronrod[KronrodSz - 1]
			  * std::abs(f_center - mean);

	for (std::size_t jj = 0; jj < KronrodSz - 1; ++jj)
	  result_asc += w_kronrod[jj]
			* (std::abs(fv1[jj] - mean)
			 + std::abs(fv2[jj] - mean));

	auto err = (result_kronrod - result_gauss) * half_length;

	result_kronrod *= half_length;
	result_abs *= abs_half_length;
	result_asc *= abs_half_length;

	return {result_kronrod,
		  rescale_error(err, result_abs, result_asc),
		  result_abs, result_asc};
      }

  template<typename Tp>
    template<typename FuncTp>
      auto
      gauss_kronrod_integral<Tp>::
      integrate(FuncTp func, Tp lower, Tp upper) const
      -> gauss_kronrod_integral_t<Tp, std::invoke_result_t<FuncTp, Tp>>
      {
	switch(this->m_rule)
	  {
	  case Kronrod_15:
	    {
	      using _GK = qk_integrator<Tp, FuncTp, Kronrod_15>;
	      return s_integrate(_GK::s_x_kronrod, _GK::s_w_gauss,
				  _GK::s_w_kronrod, func, lower, upper);
	    }
	  case Kronrod_21:
	    {
	      using _GK = qk_integrator<Tp, FuncTp, Kronrod_21>;
	      return s_integrate(_GK::s_x_kronrod, _GK::s_w_gauss,
				  _GK::s_w_kronrod, func, lower, upper);
	    }
	  case Kronrod_31:
	    {
	      using _GK = qk_integrator<Tp, FuncTp, Kronrod_31>;
	      return s_integrate(_GK::s_x_kronrod, _GK::s_w_gauss,
				  _GK::s_w_kronrod, func, lower, upper);
	    }
	  case Kronrod_41:
	    {
	      using _GK = qk_integrator<Tp, FuncTp, Kronrod_41>;
	      return s_integrate(_GK::s_x_kronrod, _GK::s_w_gauss,
				  _GK::s_w_kronrod, func, lower, upper);
	    }
	  case Kronrod_51:
	    {
	      using _GK = qk_integrator<Tp, FuncTp, Kronrod_51>;
	      return s_integrate(_GK::s_x_kronrod, _GK::s_w_gauss,
				  _GK::s_w_kronrod, func, lower, upper);
	    }
	  case Kronrod_61:
	    {
	      using _GK = qk_integrator<Tp, FuncTp, Kronrod_61>;
	      return s_integrate(_GK::s_x_kronrod, _GK::s_w_gauss,
				  _GK::s_w_kronrod, func, lower, upper);
	    }
	  default:
	    {
	      return s_integrate(this->m_x_kronrod, this->m_w_gauss,
				  this->m_w_kronrod, func, lower, upper);
	    }
	  }
      }

} // namespace __gnu_cxx

#endif // GAUSS_KRONROD_INTERGAL_TCC
