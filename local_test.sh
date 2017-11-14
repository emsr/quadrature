#!  /bin/bash

tool="cp -f"
makedir="mkdir -p"

#rm -rf "testsuite"

#utildir="testsuite/util"
test_dir="testsuite/special_functions"
text_dir="testsuite/ext/special_functions"
quad_dir="testsuite/ext/integration"


${makedir} ${quad_dir}/gauss_kronrod
${makedir} ${quad_dir}/midpoint
${makedir} ${quad_dir}/multi_singular
${makedir} ${quad_dir}/singular
${makedir} ${quad_dir}/trapezoid

#${tool} check_gauss_kronrod.cc    ${quad_dir}/gauss_kronrod/check.cc
#${tool} check_midpoint.cc         ${quad_dir}/midpoint/check.cc
#${tool} check_multi_singular.cc   ${quad_dir}/multi_singular/check.cc
#${tool} check_singular.cc         ${quad_dir}/singular/check.cc
#${tool} check_trapezoid.cc        ${quad_dir}/trapezoid/check.cc
##${tool} check_.cc   ${quad_dir}//check.cc


${makedir} ${test_dir}/01_assoc_laguerre
${makedir} ${test_dir}/02_assoc_legendre
${makedir} ${test_dir}/15_hermite
${makedir} ${test_dir}/16_laguerre
${makedir} ${test_dir}/17_legendre
${makedir} ${test_dir}/20_sph_legendre

${tool} orthonorm_assoc_laguerre.cc    ${test_dir}/01_assoc_laguerre/orthonorm.cc
${tool} orthonorm_assoc_legendre.cc    ${test_dir}/02_assoc_legendre/orthonorm.cc
${tool} orthonorm_hermite.cc           ${test_dir}/15_hermite/orthonorm.cc
${tool} orthonorm_laguerre.cc          ${test_dir}/16_laguerre/orthonorm.cc
${tool} orthonorm_legendre.cc          ${test_dir}/17_legendre/orthonorm.cc
${tool} orthonorm_sph_legendre.cc      ${test_dir}/20_sph_legendre/orthonorm.cc


${makedir} ${text_dir}/chebyshev_t
${makedir} ${text_dir}/chebyshev_u
${makedir} ${text_dir}/chebyshev_v
${makedir} ${text_dir}/chebyshev_w
${makedir} ${text_dir}/gegenbauer
${makedir} ${text_dir}/jacobi
${makedir} ${text_dir}/radpoly
${makedir} ${text_dir}/zernike

${tool} orthonorm_chebyshev_t.cc    ${text_dir}/chebyshev_t/orthonorm.cc
${tool} orthonorm_chebyshev_u.cc    ${text_dir}/chebyshev_u/orthonorm.cc
${tool} orthonorm_chebyshev_v.cc    ${text_dir}/chebyshev_v/orthonorm.cc
${tool} orthonorm_chebyshev_w.cc    ${text_dir}/chebyshev_w/orthonorm.cc
${tool} orthonorm_gegenbauer.cc     ${text_dir}/gegenbauer/orthonorm.cc
${tool} orthonorm_jacobi.cc         ${text_dir}/jacobi/orthonorm.cc
${tool} orthonorm_radpoly.cc        ${text_dir}/radpoly/orthonorm.cc
${tool} orthonorm_zernike.cc        ${text_dir}/zernike/orthonorm.cc
#${tool} orthonorm_.cc    ${text_dir}//orthonorm.cc
