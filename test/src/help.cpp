
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

template<typename Ret, typename Tp, typename Integrator>
  std::tuple<_Tp, _Tp, bool>
  qc25c(std::function<Ret(Tp)> func, Tp lower, Tp upper, Tp center,
	Integrator quad)
  {
    using quad_ret = std::tuple<Tp&, Tp&, Tp&, Tp&>;

    auto func_cauchy = [func, center](Tp x)
		       -> Tp
		       { return func(x) / (x - center); };
    Tp result, abserr, resabs, resasc;
    quad_ret{result, abserr, resabs, resasc}
      = quad(func_cauchy, lower, upper);

    return std::make_tuple(Tp{}, Tp{}, true);
  }

template<typename Ret, typename Tp>
  std::tuple<Tp, Tp, Tp, Tp>
  qk_integrate(std::function<Ret(Tp)> func, Tp lower, Tp upper,
	       Kronrod_Rule qkintrule)
  {
    return std::make_tuple(Tp{}, Tp{}, Tp{}, Tp{});
  }

template<typename Ret, typename Tp, typename Integrator>
  std::tuple<Tp, Tp>
  qawc_integrate(std::function<Ret(Tp)> func,
		 Tp lower, Tp upper, Tp center,
		 Tp max_abs_err, Tp max_rel_err,
		 Integrator quad)
  {
    Tp result0, abserr0;
    bool err_reliable;
    std::tie(result0, abserr0, err_reliable)
      = qc25c(func, lower, upper, center, quad);
    //return std::make_tuple(Tp{}, Tp{}); // ???
    return std::make_tuple(result0, abserr0);
  }

template<typename Ret, typename Tp>
  std::tuple<Tp, Tp>
  qawc_integrate(std::function<Ret(Tp)> func,
		 Tp lower, Tp upper, Tp center,
		 Tp max_abs_err, Tp max_rel_err,
		 Kronrod_Rule qk_rule = QK_15)
  {
    auto quad
      = [qk_rule]
	(std::function<Ret(Tp)> func, Tp lower, Tp upper)
	-> std::tuple<Tp, Tp, Tp, Tp>
	{ return qk_integrate(func, lower, upper, qk_rule); };

    return qawc_integrate(func, lower, upper, center,
			  max_abs_err, max_rel_err, quad);
  }

int
main()
{
  qawc_integrate([](double x)->double{return x;},
                 0.0, 1.0, 0.5, 0.00001, 0.00001);
}
