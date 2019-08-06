cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  cxx_integration
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_custom_target(make_cxx_integration_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

set(TR29124_LIBRARIES
  cxx_math_constants
  cxx_fp_utils
  cxx_complex_utils
  cxx_complex_math
  cxx_polynomial
  cxx_special_functions
  cxx_continued_fractions
  cxx_summation
  cxx_integration)

find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_custom_target(docs_cxx_integration
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
  cxx_integration
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
