/*
$HOME/bin/bin/g++ -std=c++17 help.cpp

clang-5.0 -std=c++17 -stdlib=libc++ help.cpp

clang-5.0 -std=c++17 -stdlib=libstdc++ help.cpp
*/

#include <tuple>
#include <functional>
#include <memory>
#include <type_traits>

template<typename _Tp, typename _FuncTp, typename... _Parms,
	 typename _Ret = std::invoke_result_t<_FuncTp, _Tp, _Parms...>>
  std::function<_Ret(_Tp)>
  make_function(_FuncTp f, _Parms... p)
  { return [f, p...](_Tp x)->_Ret{ return f(x, p...); }; }

template<typename _Tp, typename _FuncTp,
	 typename _Ret = std::invoke_result_t<_FuncTp, _Tp>>
  struct counted_function
  {
    counted_function(_FuncTp f)
    : m_func(f), m_neval(new int{0})
    { }

    _Ret
    operator()(_Tp x) const
    {
      ++(*this->m_neval);
      return this->m_func(x);
    }

    _FuncTp m_func;
    mutable std::shared_ptr<int> m_neval;
  };

enum Kronrod_Rule
{
  QK_15 = 15,
  QK_21 = 21,
  QK_31 = 31,
  QK_41 = 41,
  QK_51 = 51,
  QK_61 = 61
};

template<typename _FuncTp, typename _Tp, typename _Integrator>
  std::tuple<_Tp, _Tp, bool>
  qc25c(_FuncTp func, _Tp lower, _Tp upper, _Tp center,
	_Integrator quad)
  {
    using quad_ret = std::tuple<_Tp&, _Tp&, _Tp&, _Tp&>;

    auto func_cauchy = [func, center](_Tp x)
		       -> _Tp
		       { return func(x) / (x - center); };
    _Tp result, abserr, resabs, resasc;
    quad_ret{result, abserr, resabs, resasc}
      = quad(func_cauchy, lower, upper);

    return std::make_tuple(_Tp{}, _Tp{}, true);
  }

template<typename _Tp, typename _FuncTp>
  std::tuple<_Tp, _Tp, _Tp, _Tp>
  qk_integrate(_FuncTp __func, _Tp __lower, _Tp __upper,
	       Kronrod_Rule __qkintrule)
  {
    return std::make_tuple(_Tp{}, _Tp{}, _Tp{}, _Tp{});
  }

template<typename _FuncTp, typename _Tp, typename _Integrator>
  std::tuple<_Tp, _Tp>
  qawc_integrate(_FuncTp __func,
		 _Tp __lower, _Tp __upper, _Tp __center,
		 _Tp __max_abs_err, _Tp __max_rel_err,
		 _Integrator __quad)
  {
    _Tp __result0, __abserr0;
    bool __err_reliable;
    std::tie(__result0, __abserr0, __err_reliable)
      = qc25c(__func, __lower, __upper, __center, __quad);
    return std::make_tuple(_Tp{}, _Tp{});
  }


template<typename _FuncTp, typename _Tp>
  std::tuple<_Tp, _Tp>
  qawc_integrate(_FuncTp __func,
		 _Tp __lower, _Tp __upper, _Tp __center,
		 _Tp __max_abs_err, _Tp __max_rel_err,
		 Kronrod_Rule __qk_rule = QK_15)
  {
    auto __quad
      = [__qk_rule]
	(_FuncTp __func, _Tp __lower, _Tp __upper)
	-> std::tuple<_Tp, _Tp, _Tp, _Tp>
	{ return qk_integrate(__func, __lower, __upper, __qk_rule); };

    return qawc_integrate(__func, __lower, __upper, __center,
			  __max_abs_err, __max_rel_err, __quad);
  }

template<typename _Tp>
  _Tp
  f459(_Tp x)
  { return x / (x - 0.5); }

int
main()
{
  auto f = make_function<double>(f459<double>);
  auto fc = counted_function<double, decltype(f)>(f);
  qawc_integrate(f459<double>,
                 0.0, 1.0, 0.5, 0.00001, 0.00001);
// This fails too.
//  qawc_integrate([](double x)->double{return x;},
//                 0.0, 1.0, 0.5, 0.00001, 0.00001);
}

/*
help.cpp: In instantiation of 'std::tuple<_Tp, _Tp, bool> qc25c(_FuncTp, _Tp, _Tp, _Tp, _Integrator) [with _FuncTp = main()::<lambda(double)>; _Tp = double; _Integrator = qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]::<lambda(main()::<lambda(double)>, double, double)>]':
help.cpp:84:14:   required from 'std::tuple<_Tp, _Tp> qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, _Integrator) [with _FuncTp = main()::<lambda(double)>; _Tp = double; _Integrator = qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]::<lambda(main()::<lambda(double)>, double, double)>]'
help.cpp:102:26:   required from 'std::tuple<_Tp, _Tp> qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]'
help.cpp:117:49:   required from here
help.cpp:61:13: error: no match for call to '(qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]::<lambda(main()::<lambda(double)>, double, double)>) (qc25c(_FuncTp, _Tp, _Tp, _Tp, _Integrator) [with _FuncTp = main()::<lambda(double)>; _Tp = double; _Integrator = qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]::<lambda(main()::<lambda(double)>, double, double)>]::<lambda(double)>&, double&, double&)'
       = quad(func_cauchy, lower, upper);
         ~~~~^~~~~~~~~~~~~~~~~~~~~~~~~~~
help.cpp:99:10: note: candidate: 'qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]::<lambda(main()::<lambda(double)>, double, double)>'
  -> std::tuple<_Tp, _Tp, _Tp, _Tp>
          ^~~~~~~~~~~~~~~~~~~~~~~~~
help.cpp:99:10: note:   no known conversion for argument 1
 from 'qc25c(_FuncTp, _Tp, _Tp, _Tp, _Integrator) [with _FuncTp = main()::<lambda(double)>; _Tp = double; _Integrator = qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]::<lambda(main()::<lambda(double)>, double, double)>]::<lambda(double)>'
   to 'main()::<lambda(double)>'
*/
