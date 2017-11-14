
// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_chebyshev_v(int n1, int n2, _Tp x)
  {
    const auto _S_eps = __gnu_cxx::__epsilon(x);
    const auto _S_inf = __gnu_cxx::__infinity(x);
    if (std::abs(x - _Tp{1}) < _S_eps)
      return (n1 + n2) & 1 ? -_S_inf : _S_inf;
    else
      return __gnu_cxx::chebyshev_v(n1, x)
	   * __gnu_cxx::chebyshev_v(n2, x)
	   * std::sqrt((_Tp{1} + x) / (_Tp{1} - x))
	   / _Tp{2};
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_chebyshev_v()
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto integ_precision = _Tp{100} * eps;
    const auto comp_precision = _Tp{10} * integ_precision;

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [n1, n2](_Tp x)->_Tp{return normalized_chebyshev_v(n1, n2, x);};

	    auto [result, error]
		= integrate_singular_endpoints(func,
				 _Tp{-1}, _Tp{1},
				 _Tp{0.5}, _Tp{-0.5}, 0, 0,
				 integ_precision, _Tp{0});

	    assert(std::abs(delta<_Tp>(n1, n2) - result) < abs_precision);
	  }
      }
  }

int
main()
{
  test_chebyshev_v<float>();

  test_chebyshev_v<double>();

  test_chebyshev_v<long double>();
}
