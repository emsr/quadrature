cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  cxx_integration
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -Wno-psabi)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_executable(test_phase_iterator test_phase_iterator.cpp)
target_include_directories(test_phase_iterator PRIVATE
  include
  ../cxx_math_const/include)
target_link_libraries(test_phase_iterator quadmath)

add_executable(test_quadrature test_quadrature.cpp)
target_include_directories(test_quadrature PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(test_quadrature quadmath)

add_executable(test_trapezoid_integral test_trapezoid_integral.cpp)
target_include_directories(test_trapezoid_integral PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(test_trapezoid_integral quadmath)

add_executable(test_midpoint_integral test_midpoint_integral.cpp)
target_include_directories(test_midpoint_integral PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(test_midpoint_integral quadmath)

add_executable(test_simpson_integral test_simpson_integral.cpp)
target_include_directories(test_simpson_integral PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(test_simpson_integral quadmath)

add_executable(test_double_exp_integrate test_double_exp_integrate.cpp)
target_include_directories(test_double_exp_integrate PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(test_double_exp_integrate quadmath)

add_executable(test_gauss_hermite test_gauss_hermite.cpp)
target_include_directories(test_gauss_hermite PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(test_gauss_hermite quadmath)

add_executable(test_gauss_laguerre test_gauss_laguerre.cpp)
target_include_directories(test_gauss_laguerre PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(test_gauss_laguerre quadmath)

add_executable(test_mapper test_mapper.cpp)
target_include_directories(test_mapper PRIVATE
  include
  ../include
  ../cxx_fp_utils/include)
target_link_libraries(test_mapper quadmath)

add_executable(test_composite_trapezoid_integral test_composite_trapezoid_integral.cpp)
target_include_directories(test_composite_trapezoid_integral PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(test_composite_trapezoid_integral quadmath)

add_executable(test_gauss_kronrod_rule test_gauss_kronrod_rule.cpp)
target_include_directories(test_gauss_kronrod_rule PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(test_gauss_kronrod_rule quadmath)

# Orthogonal polynomials (C++17) tests.

add_executable(assoc_laguerre_test assoc_laguerre_test.cpp)
target_include_directories(assoc_laguerre_test PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(assoc_laguerre_test quadmath)

add_executable(assoc_legendre_test assoc_legendre_test.cpp)
target_include_directories(assoc_legendre_test PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(assoc_legendre_test quadmath)

add_executable(sph_legendre_test sph_legendre_test.cpp)
target_include_directories(sph_legendre_test PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(sph_legendre_test quadmath)

add_executable(hermite_test hermite_test.cpp)
target_include_directories(hermite_test PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(hermite_test quadmath)

add_executable(laguerre_test laguerre_test.cpp)
target_include_directories(laguerre_test PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(laguerre_test quadmath)

add_executable(legendre_test legendre_test.cpp)
target_include_directories(legendre_test PRIVATE
  include
  ../cxx_math_const/include
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include)
target_link_libraries(legendre_test quadmath)

# Orthogonal polynomials (tr29124_test project extensions) tests.

if (FOUND_TR29124)
  add_executable(gegenbauer_test gegenbauer_test.cpp)
  target_include_directories(gegenbauer_test PRIVATE
    include
    ../cxx_math_const/include
    ../include
    ../cxx_fp_utils/include
    ../polynomial/include
    ../cxx_summation/include)
  target_link_libraries(gegenbauer_test quadmath)

  add_executable(jacobi_test jacobi_test.cpp)
  target_include_directories(jacobi_test PRIVATE
    include
    ../cxx_math_const/include
    ../include
    ../cxx_fp_utils/include
    ../polynomial/include
    ../cxx_summation/include)
  target_link_libraries(jacobi_test quadmath)

  add_executable(chebyshev_t_test chebyshev_t_test.cpp)
  target_include_directories(chebyshev_t_test PRIVATE
    include
    ../cxx_math_const/include
    ../include
    ../cxx_fp_utils/include
    ../polynomial/include
    ../cxx_summation/include)
  target_link_libraries(chebyshev_t_test quadmath)

  add_executable(chebyshev_u_test chebyshev_u_test.cpp)
  target_include_directories(chebyshev_u_test PRIVATE
    include
    ../cxx_math_const/include
    ../include
    ../cxx_fp_utils/include
    ../polynomial/include
    ../cxx_summation/include)
  target_link_libraries(chebyshev_u_test quadmath)

  add_executable(chebyshev_v_test chebyshev_v_test.cpp)
  target_include_directories(chebyshev_v_test PRIVATE
    include
    ../cxx_math_const/include
    ../include
    ../cxx_fp_utils/include
    ../polynomial/include
    ../cxx_summation/include)
  target_link_libraries(chebyshev_v_test quadmath)

  add_executable(chebyshev_w_test chebyshev_w_test.cpp)
  target_include_directories(chebyshev_w_test PRIVATE
    include
    ../cxx_math_const/include
    ../include
    ../cxx_fp_utils/include
    ../polynomial/include
    ../cxx_summation/include)
  target_link_libraries(chebyshev_w_test quadmath)

  add_executable(radpoly_test radpoly_test.cpp)
  target_include_directories(radpoly_test PRIVATE
    include
    ../cxx_math_const/include
    ../include
    ../cxx_fp_utils/include
    ../polynomial/include
    ../cxx_summation/include)
  target_link_libraries(radpoly_test quadmath)

  add_executable(zernike_test zernike_test.cpp)
  target_include_directories(zernike_test PRIVATE
    include
    ../cxx_math_const/include
    ../include
    ../cxx_fp_utils/include
    ../polynomial/include
    ../cxx_summation/include)
  target_link_libraries(zernike_test quadmath)
endif (FOUND_TR29124)

