
// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  normalized_chebyshev_u(int n1, int n2, _Tp x)
  {
    const auto _S_pi_2 = __gnu_cxx::__const_pi_half(x);
    return __gnu_cxx::chebyshev_u(n2, x)
	 * __gnu_cxx::chebyshev_u(n1, x)
	 * std::sqrt(_Tp{1} - x * x)
	 / _S_pi_2;
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_chebyshev_u()
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto integ_precision = _Tp{100} * eps;
    const auto comp_precision = _Tp{10} * integ_precision;

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [n1, n2](_Tp x)->_Tp{return normalized_chebyshev_u(n1, n2, x);};

	    auto [result, error]
		= integrate(func, _Tp{-1}, _Tp{1}, integ_precision, _Tp{0});

	    assert(std::abs(delta<_Tp>(n1, n2) - result) < abs_precision);
	  }
      }
  }

int
main()
{
  test_chebyshev_u<float>();

  test_chebyshev_u<double>();

  test_chebyshev_u<long double>();
}
