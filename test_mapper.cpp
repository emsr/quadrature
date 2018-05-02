/*
$HOME/bin/bin/g++ -std=c++2a -o test_mapper test_mapper.cpp
*/

#include <iostream>
#include "integration_transform.h"

template<typename _Tp>
  void
  test_map_minf_pinf()
  {
    auto fun = [](_Tp x)->_Tp{ return x; };
    using func_t = decltype(fun);

    for (int i = 0; i <= 100; ++i)
      {
	auto x = _Tp{0.01} * i;
        auto y = __gnu_cxx::map_minf_pinf<func_t, _Tp>(fun)(x);
	std::cout << ' ' << x << ' ' << y << '\n';
      }
  }

int
main()
{
  test_map_minf_pinf<double>();
}
