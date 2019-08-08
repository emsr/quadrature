cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  cxx_integration
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

set(TR29124_LIBRARIES
  cxx_math_constants
  cxx_fp_utils
  cxx_complex_utils
  cxx_complex_math
  cxx_polynomial
  cxx_summation
  cxx_continued_fractions
  cxx_integration
  cxx_special_functions)

find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_custom_target(docs_integration
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation for cxx_integration" VERBATIM
  )
endif (DOXYGEN_FOUND)

# Header-only library.
add_library(cxx_integration INTERFACE)
target_include_directories(cxx_integration INTERFACE include)

add_executable(test_phase_iterator test_phase_iterator.cpp)
target_link_libraries(test_phase_iterator
  cxx_math_constants
  cxx_fp_utils
  cxx_complex_utils
  cxx_complex_math
  cxx_polynomial
  cxx_summation
  cxx_continued_fractions
  cxx_integration
  quadmath)

add_executable(test_quadrature test_quadrature.cpp)
target_link_libraries(test_quadrature
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(test_trapezoid_integral test_trapezoid_integral.cpp)
target_link_libraries(test_trapezoid_integral
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(test_midpoint_integral test_midpoint_integral.cpp)
target_link_libraries(test_midpoint_integral
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(test_simpson_integral test_simpson_integral.cpp)
target_link_libraries(test_simpson_integral
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(test_double_exp_integrate test_double_exp_integrate.cpp)
target_link_libraries(test_double_exp_integrate
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(test_gauss_hermite test_gauss_hermite.cpp)
target_link_libraries(test_gauss_hermite
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(test_gauss_laguerre test_gauss_laguerre.cpp)
target_link_libraries(test_gauss_laguerre
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(test_mapper test_mapper.cpp)
target_link_libraries(test_mapper
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(test_composite_trapezoid_integral test_composite_trapezoid_integral.cpp)
target_link_libraries(test_composite_trapezoid_integral
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(test_gauss_kronrod_rule test_gauss_kronrod_rule.cpp)
target_link_libraries(test_gauss_kronrod_rule
  ${TR29124_LIBRARIES}
  quadmath)

# Orthogonal polynomials (C++17) tests.

add_executable(assoc_laguerre_test assoc_laguerre_test.cpp)
target_link_libraries(assoc_laguerre_test
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(assoc_legendre_test assoc_legendre_test.cpp)
target_link_libraries(assoc_legendre_test
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(sph_legendre_test sph_legendre_test.cpp)
target_link_libraries(sph_legendre_test
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(hermite_test hermite_test.cpp)
target_link_libraries(hermite_test
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(laguerre_test laguerre_test.cpp)
target_link_libraries(laguerre_test
  ${TR29124_LIBRARIES}
  quadmath)

add_executable(legendre_test legendre_test.cpp)
target_link_libraries(legendre_test
  ${TR29124_LIBRARIES}
  quadmath)

# Orthogonal polynomials (tr29124_test project extensions) tests.

if (FOUND_TR29124)
  add_executable(gegenbauer_test gegenbauer_test.cpp)
  target_link_libraries(gegenbauer_test
    ${TR29124_LIBRARIES}
    quadmath)

  add_executable(jacobi_test jacobi_test.cpp)
  target_link_libraries(jacobi_test
    ${TR29124_LIBRARIES}
    quadmath)

  add_executable(chebyshev_t_test chebyshev_t_test.cpp)
  target_link_libraries(chebyshev_t_test
    ${TR29124_LIBRARIES}
    quadmath)

  add_executable(chebyshev_u_test chebyshev_u_test.cpp)
  target_link_libraries(chebyshev_u_test
    ${TR29124_LIBRARIES}
    quadmath)

  add_executable(chebyshev_v_test chebyshev_v_test.cpp)
  target_link_libraries(chebyshev_v_test
    ${TR29124_LIBRARIES}
    quadmath)

  add_executable(chebyshev_w_test chebyshev_w_test.cpp)
  target_link_libraries(chebyshev_w_test
    ${TR29124_LIBRARIES}
    quadmath)

  add_executable(radpoly_test radpoly_test.cpp)
  target_link_libraries(radpoly_test
    ${TR29124_LIBRARIES}
    quadmath)

  add_executable(zernike_test zernike_test.cpp)
  target_link_libraries(zernike_test
    ${TR29124_LIBRARIES}
    quadmath)
endif (FOUND_TR29124)

