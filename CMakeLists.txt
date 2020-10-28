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

find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_custom_target(docs_cxx_integration
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation for cxx_integration" VERBATIM
  )
endif (DOXYGEN_FOUND)

add_library(cxx_integration_polynomial INTERFACE)
target_include_directories(cxx_integration_polynomial INTERFACE 3rdparty/cxx_polynomial/include)

add_library(cxx_integration_complex_utils INTERFACE)
target_include_directories(cxx_integration_complex_utils INTERFACE 3rdparty/cxx_complex_utils/include)

add_library(cxx_integration_fp_utils INTERFACE)
target_include_directories(cxx_integration_fp_utils INTERFACE 3rdparty/cxx_fp_utils/include)

add_library(cxx_integration_numeric_limits INTERFACE)
target_include_directories(cxx_integration_numeric_limits INTERFACE 3rdparty/cxx_numeric_limits/include)

add_library(cxx_integration_special_functions INTERFACE)
target_include_directories(cxx_integration_special_functions INTERFACE 3rdparty/cxx_special_functions/include)
target_link_libraries(cxx_integration_special_functions
  INTERFACE cxx_integration_numeric_limits)

add_library(cxx_integration INTERFACE)
target_include_directories(cxx_integration INTERFACE include)

add_executable(test_phase_iterator src/test_phase_iterator.cpp)
target_link_libraries(test_phase_iterator
  cxx_integration)

add_executable(test_quadrature src/test_quadrature.cpp)
target_link_libraries(test_quadrature
  cxx_integration
  cxx_integration_fp_utils
  cxx_integration_special_functions)

add_executable(test_trapezoid_integral src/test_trapezoid_integral.cpp)
target_link_libraries(test_trapezoid_integral
  cxx_integration
  cxx_integration_polynomial
  cxx_integration_complex_utils
  cxx_integration_fp_utils)

add_executable(test_midpoint_integral src/test_midpoint_integral.cpp)
target_link_libraries(test_midpoint_integral
  cxx_integration
  cxx_integration_polynomial)

add_executable(test_simpson_integral src/test_simpson_integral.cpp)
target_link_libraries(test_simpson_integral
  cxx_integration
  cxx_integration_polynomial
  cxx_integration_complex_utils
  cxx_integration_fp_utils)

add_executable(test_double_exp_integrate src/test_double_exp_integrate.cpp)
target_link_libraries(test_double_exp_integrate
  cxx_integration
  cxx_integration_polynomial)

add_executable(test_gauss_hermite src/test_gauss_hermite.cpp)
target_link_libraries(test_gauss_hermite
  cxx_integration)

add_executable(test_gauss_laguerre src/test_gauss_laguerre.cpp)
target_link_libraries(test_gauss_laguerre
  cxx_integration)

add_executable(test_mapper src/test_mapper.cpp)
target_link_libraries(test_mapper
  cxx_integration)

add_executable(test_composite_trapezoid_integral src/test_composite_trapezoid_integral.cpp)
target_link_libraries(test_composite_trapezoid_integral
  cxx_integration
  cxx_integration_polynomial)

add_executable(test_gauss_kronrod_rule src/test_gauss_kronrod_rule.cpp)
target_link_libraries(test_gauss_kronrod_rule
  cxx_integration)

add_executable(test_factorial_integration src/test_factorial.cpp)
target_link_libraries(test_factorial_integration
  cxx_integration_special_functions
  cxx_integration_fp_utils)
