/*
$HOME/bin/bin/g++ -std=c++17 help.cpp

clang-5.0 -std=c++17 -stdlib=libc++ help.cpp

clang-5.0 -std=c++17 -stdlib=libstdc++ -I$HOME/bin/include/c++/9.0.0 -I/home/ed/bin/include/c++/9.0.0/x86_64-pc-linux-gnu help.cpp
*/

#include <tuple>
#include <functional>
#include <memory>
#include <type_traits>

enum Kronrod_Rule
{
  QK_15 = 15,
  QK_21 = 21,
  QK_31 = 31,
  QK_41 = 41,
  QK_51 = 51,
  QK_61 = 61
};

template<typename _Ret, typename _Tp, typename _Integrator>
  std::tuple<_Tp, _Tp, bool>
  qc25c(std::function<_Ret(_Tp)> func, _Tp lower, _Tp upper, _Tp center,
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

template<typename _Ret, typename _Tp>
  std::tuple<_Tp, _Tp, _Tp, _Tp>
  qk_integrate(std::function<_Ret(_Tp)> __func, _Tp __lower, _Tp __upper,
	       Kronrod_Rule __qkintrule)
  {
    return std::make_tuple(_Tp{}, _Tp{}, _Tp{}, _Tp{});
  }

template<typename _Ret, typename _Tp, typename _Integrator>
  std::tuple<_Tp, _Tp>
  qawc_integrate(std::function<_Ret(_Tp)> __func,
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

template<typename _Ret, typename _Tp>
  std::tuple<_Tp, _Tp>
  qawc_integrate(std::function<_Ret(_Tp)> __func,
		 _Tp __lower, _Tp __upper, _Tp __center,
		 _Tp __max_abs_err, _Tp __max_rel_err,
		 Kronrod_Rule __qk_rule = QK_15)
  {
    auto __quad
      = [__qk_rule]
	(std::function<_Ret(_Tp)> __func, _Tp __lower, _Tp __upper)
	-> std::tuple<_Tp, _Tp, _Tp, _Tp>
	{ return qk_integrate(__func, __lower, __upper, __qk_rule); };

    return qawc_integrate(__func, __lower, __upper, __center,
			  __max_abs_err, __max_rel_err, __quad);
  }

int
main()
{
  qawc_integrate([](double x)->double{return x;},
                 0.0, 1.0, 0.5, 0.00001, 0.00001);
}

/*
help.cpp: In instantiation of 'std::tuple<_Tp, _Tp, bool> qc25c(_FuncTp, _Tp, _Tp, _Tp, _Integrator) [with _FuncTp = main()::<lambda(double)>; _Tp = double; _Integrator = qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]::<lambda(main()::<lambda(double)>, double, double)>]':
help.cpp:84:14:   required from 'std::tuple<_Tp, _Tp> qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, _Integrator) [with _FuncTp = main()::<lambda(double)>; _Tp = double; _Integrator = qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]::<lambda(main()::<lambda(double)>, double, double)>]'
help.cpp:101:26:   required from 'std::tuple<_Tp, _Tp> qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]'
help.cpp:109:49:   required from here
help.cpp:61:13: error: no match for call to '(qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]::<lambda(main()::<lambda(double)>, double, double)>) (qc25c(_FuncTp, _Tp, _Tp, _Tp, _Integrator) [with _FuncTp = main()::<lambda(double)>; _Tp = double; _Integrator = qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]::<lambda(main()::<lambda(double)>, double, double)>]::<lambda(double)>&, double&, double&)'
       = quad(func_cauchy, lower, upper);
         ~~~~^~~~~~~~~~~~~~~~~~~~~~~~~~~
help.cpp:98:10: note: candidate: 'qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]::<lambda(main()::<lambda(double)>, double, double)>'
  -> std::tuple<_Tp, _Tp, _Tp, _Tp>
          ^~~~~~~~~~~~~~~~~~~~~~~~~
help.cpp:98:10: note:   no known conversion for argument 1
 from 'qc25c(_FuncTp, _Tp, _Tp, _Tp, _Integrator) [with _FuncTp = main()::<lambda(double)>; _Tp = double; _Integrator = qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule) [with _FuncTp = main()::<lambda(double)>; _Tp = double]::<lambda(main()::<lambda(double)>, double, double)>]::<lambda(double)>'
   to 'main()::<lambda(double)>'

How does _Integrator get to be qawc_integrate(_FuncTp, _Tp, _Tp, _Tp, _Tp, _Tp, Kronrod_Rule)?
Why isn't it qk_integrate?

Clang-5.0:
----------
help.cpp:61:9: error: no matching function for call to object of type '(lambda at help.cpp:96:9)'
      = quad(func_cauchy, lower, upper);
        ^~~~
help.cpp:84:9: note: in instantiation of function template specialization 'qc25c<(lambda at help.cpp:108:18), double, (lambda at help.cpp:96:9)>' requested here
      = qc25c(__func, __lower, __upper, __center, __quad);
        ^
help.cpp:101:12: note: in instantiation of function template specialization 'qawc_integrate<(lambda at help.cpp:108:18), double, (lambda at help.cpp:96:9)>' requested here
    return qawc_integrate(__func, __lower, __upper, __center,
           ^
help.cpp:108:3: note: in instantiation of function template specialization 'qawc_integrate<(lambda at help.cpp:108:18), double>' requested here
  qawc_integrate([](double x)->double{return x;},
  ^
help.cpp:96:9: note: candidate function not viable: no known conversion from '(lambda at help.cpp:56:24)' to '(lambda at help.cpp:108:18)' for 1st argument
      = [__qk_rule]
        ^

I need type erasure? Sigh.
*/
