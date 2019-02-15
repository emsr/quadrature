
CXX_VER = -std=gnu++17
CXX_INST_DIR = $(HOME)/bin$(SUFFIX)
ifeq ("$(wildcard $(CXX_INST_DIR)/bin/g++)","")
  CXX_VER = -std=gnu++2a
  CXX_INST_DIR = $(HOME)/bin
  ifeq ("$(wildcard $(CXX_INST_DIR)/bin/g++)","")
    CXX_VER = -std=gnu++17
    ifeq ($(wildcard "/mingw64"),"")
      CXX_INST_DIR = /mingw64
    else
      CXX_INST_DIR = /usr
    endif
  endif
endif

#OPT = -O3
OPT = -g
#OPT = -g -fsanitize=signed-integer-overflow -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fsanitize=alignment
GCC = $(CXX_INST_DIR)/bin/gcc $(OPT) -Wall -Wextra
CXX17 = $(CXX_INST_DIR)/bin/g++ -std=gnu++17 -fconcepts $(OPT) -Wall -Wextra -Wno-psabi
CXX20 = $(CXX_INST_DIR)/bin/g++ -std=gnu++2a $(OPT) -Wall -Wextra -Wno-psabi
#CXXMAX = $(CXX20)
CXXMAX = $(CXX_INST_DIR)/bin/g++ $(CXX_VER) $(OPT) -Wall -Wextra -Wno-psabi
CXX_INC_DIR = $(CXX_INST_DIR)/include/c++/8.0.0/bits
CXX_LIB_DIR = $(CXX_INST_DIR)/lib64

#WRAPPER_LIB_DIR = ../wrappers/release
WRAPPER_LIB_DIR = ../wrappers/debug
WRAPPER_LIBS = -L$(WRAPPER_LIB_DIR) -lwrap_burkhardt -lgfortran

INC_DIR = include/ext
INCLUDES =  -I../include -Iinclude -I../polynomial/include

TEST_OUT_DIR = test_output

BIN_DIR = bin

INCS = \
  $(INC_DIR)/cquad_const.tcc \
  $(INC_DIR)/cquad_integrate.tcc \
  $(INC_DIR)/cquad_workspace.h \
  $(INC_DIR)/double_exp_integrate.tcc \
  $(INC_DIR)/extrapolation_table.h \
  $(INC_DIR)/extrapolation_table.tcc \
  $(INC_DIR)/fourier_transform.h \
  $(INC_DIR)/fourier_transform.tcc \
  $(INC_DIR)/gauss_hermite_integrate.h \
  $(INC_DIR)/gauss_jacobi_integrate.tcc \
  $(INC_DIR)/gauss_jacobi_interface.tcc \
  $(INC_DIR)/gauss_kronrod_integral.h \
  $(INC_DIR)/gauss_kronrod_integral.tcc \
  $(INC_DIR)/gauss_kronrod_rule.tcc \
  $(INC_DIR)/gauss_laguerre_integrate.h \
  $(INC_DIR)/gauss_legendre_table.h \
  $(INC_DIR)/gauss_legendre_table.tcc \
  $(INC_DIR)/gauss_quadrature.h \
  $(INC_DIR)/gauss_quadrature.tcc \
  $(INC_DIR)/glfixed_integrate.tcc \
  $(INC_DIR)/integration_error.h \
  $(INC_DIR)/integration.h \
  $(INC_DIR)/integration.tcc \
  $(INC_DIR)/integration_transform.h \
  $(INC_DIR)/integration_workspace.h \
  $(INC_DIR)/integration_workspace.tcc \
  $(INC_DIR)/jacobi.h \
  $(INC_DIR)/matrix.h \
  $(INC_DIR)/matrix.tcc \
  $(INC_DIR)/midpoint_integral.h \
  $(INC_DIR)/midpoint_integral.tcc \
  $(INC_DIR)/oscillatory_integration_table.h \
  $(INC_DIR)/oscillatory_integration_table.tcc \
  $(INC_DIR)/qag_integrate.tcc \
  $(INC_DIR)/qagp_integrate.tcc \
  $(INC_DIR)/qags_integrate.tcc \
  $(INC_DIR)/qawc_integrate.tcc \
  $(INC_DIR)/qawf_integrate.tcc \
  $(INC_DIR)/qawo_integrate.tcc \
  $(INC_DIR)/qaws_integrate.tcc \
  $(INC_DIR)/qaws_integration_table.h \
  $(INC_DIR)/qaws_integration_table.tcc \
  $(INC_DIR)/qcheb_integrate.tcc \
  $(INC_DIR)/qng_integrate.tcc \
  $(INC_DIR)/simpson_integral.h \
  $(INC_DIR)/simpson_integral.tcc \
  $(INC_DIR)/test_integral.tcc \
  $(INC_DIR)/trapezoid_integral.h \
  $(INC_DIR)/trapezoid_integral.tcc \
  $(INC_DIR)/triangle_rules.h

BINS = \
  $(BIN_DIR)/test_phase_iterator \
  $(BIN_DIR)/test_quadrature \
  $(BIN_DIR)/test_trapezoid_integral \
  $(BIN_DIR)/test_midpoint_integral \
  $(BIN_DIR)/test_simpson_integral \
  $(BIN_DIR)/test_double_exp_integrate \
  $(BIN_DIR)/test_gauss_hermite \
  $(BIN_DIR)/test_gauss_laguerre \
  $(BIN_DIR)/test_mapper \
  $(BIN_DIR)/test_composite_trapezoid_integral \
  $(BIN_DIR)/assoc_laguerre_test \
  $(BIN_DIR)/assoc_legendre_test \
  $(BIN_DIR)/sph_legendre_test \
  $(BIN_DIR)/hermite_test \
  $(BIN_DIR)/laguerre_test \
  $(BIN_DIR)/legendre_test \
  $(BIN_DIR)/gegenbauer_test \
  $(BIN_DIR)/jacobi_test \
  $(BIN_DIR)/chebyshev_t_test \
  $(BIN_DIR)/chebyshev_u_test \
  $(BIN_DIR)/chebyshev_v_test \
  $(BIN_DIR)/chebyshev_w_test \
  $(BIN_DIR)/radpoly_test \
  $(BIN_DIR)/zernike_test \
  $(BIN_DIR)/test_gauss_kronrod_rule

all: $(BIN_DIR) $(BINS)


# These require tr29124_test wrapper project libs.
BUILDERS = \
  $(BIN_DIR)/build_clenshaw_curtis \
  $(BIN_DIR)/build_double_exp_rules

builders: $(BIN_DIR) $(BUILDERS)


ortho_test: $(TEST_OUT_DIR)
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/legendre_test > $(TEST_OUT_DIR)/legendre_test.txt 2> $(TEST_OUT_DIR)/legendre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/chebyshev_t_test > $(TEST_OUT_DIR)/chebyshev_t_test.txt 2> $(TEST_OUT_DIR)/chebyshev_t_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/chebyshev_u_test > $(TEST_OUT_DIR)/chebyshev_u_test.txt 2> $(TEST_OUT_DIR)/chebyshev_u_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/chebyshev_v_test > $(TEST_OUT_DIR)/chebyshev_v_test.txt 2> $(TEST_OUT_DIR)/chebyshev_v_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/chebyshev_w_test > $(TEST_OUT_DIR)/chebyshev_w_test.txt 2> $(TEST_OUT_DIR)/chebyshev_w_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/gegenbauer_test > $(TEST_OUT_DIR)/gegenbauer_test.txt 2> $(TEST_OUT_DIR)/gegenbauer_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/jacobi_test > $(TEST_OUT_DIR)/jacobi_test.txt 2> $(TEST_OUT_DIR)/jacobi_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/assoc_laguerre_test > $(TEST_OUT_DIR)/assoc_laguerre_test.txt 2> $(TEST_OUT_DIR)/assoc_laguerre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/assoc_legendre_test > $(TEST_OUT_DIR)/assoc_legendre_test.txt 2> $(TEST_OUT_DIR)/assoc_legendre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/hermite_test > $(TEST_OUT_DIR)/hermite_test.txt 2> $(TEST_OUT_DIR)/hermite_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/sph_legendre_test > $(TEST_OUT_DIR)/sph_legendre_test.txt 2> $(TEST_OUT_DIR)/sph_legendre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/laguerre_test > $(TEST_OUT_DIR)/laguerre_test.txt 2> $(TEST_OUT_DIR)/laguerre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/radpoly_test > $(TEST_OUT_DIR)/radpoly_test.txt 2> $(TEST_OUT_DIR)/radpoly_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/zernike_test > $(TEST_OUT_DIR)/zernike_test.txt 2> $(TEST_OUT_DIR)/zernike_test.err


test: $(TEST_OUT_DIR) $(BINS)
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_gauss_kronrod_rule > $(TEST_OUT_DIR)/test_gauss_kronrod_rule.txt 2> $(TEST_OUT_DIR)/test_gauss_kronrod_rule.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_double_exp_integrate > $(TEST_OUT_DIR)/test_double_exp_integrate.txt 2> $(TEST_OUT_DIR)/test_double_exp_integrate.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_trapezoid_integral > $(TEST_OUT_DIR)/test_trapezoid_integral.txt 2> $(TEST_OUT_DIR)/test_trapezoid_integral.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_midpoint_integral > $(TEST_OUT_DIR)/test_midpoint_integral.txt 2> $(TEST_OUT_DIR)/test_midpoint_integral.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_simpson_integral > $(TEST_OUT_DIR)/test_simpson_integral.txt 2> $(TEST_OUT_DIR)/test_simpson_integral.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_phase_iterator > $(TEST_OUT_DIR)/test_phase_iterator.txt 2> $(TEST_OUT_DIR)/test_phase_iterator.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_mapper > $(TEST_OUT_DIR)/test_mapper.txt 2> $(TEST_OUT_DIR)/test_mapper.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_quadrature > $(TEST_OUT_DIR)/test_quadrature.txt 2> $(TEST_OUT_DIR)/test_quadrature.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_gauss_hermite > $(TEST_OUT_DIR)/test_gauss_hermite.txt 2> $(TEST_OUT_DIR)/test_gauss_hermite.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_gauss_laguerre > $(TEST_OUT_DIR)/test_gauss_laguerre.txt 2> $(TEST_OUT_DIR)/test_gauss_laguerre.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_composite_trapezoid_integral > $(TEST_OUT_DIR)/test_composite_trapezoid_integral.txt 2> $(TEST_OUT_DIR)/test_composite_trapezoid_integral.err

builds: $(TEST_OUT_DIR) $(BUILDERS)
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/build_double_exp_rules > $(TEST_OUT_DIR)/build_double_exp_rules.txt 2> $(TEST_OUT_DIR)/build_double_exp_rules.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/build_clenshaw_curtis > $(TEST_OUT_DIR)/build_clenshaw_curtis.txt 2> $(TEST_OUT_DIR)/build_clenshaw_curtis.err


docs:
	#rm -rf html/*
	#rm -rf latex/*
	doxygen
	cd docs/latex && make


# Binaries...

$(BIN_DIR)/build_clenshaw_curtis: $(INCS) build_clenshaw_curtis.cpp
	$(CXXMAX) $(INCLUDES) -I../wrappers -o $(BIN_DIR)/build_clenshaw_curtis build_clenshaw_curtis.cpp -lquadmath $(WRAPPER_LIBS)

$(BIN_DIR)/build_double_exp_rules: $(INCS) build_double_exp_rules.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/build_double_exp_rules build_double_exp_rules.cpp -lquadmath $(WRAPPER_LIBS)

$(BIN_DIR)/test_gauss_kronrod_rule: $(INCS) test_gauss_kronrod_rule.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/test_gauss_kronrod_rule test_gauss_kronrod_rule.cpp -lquadmath

$(BIN_DIR)/assoc_laguerre_test: $(INCS) assoc_laguerre_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/assoc_laguerre_test assoc_laguerre_test.cpp -lquadmath

$(BIN_DIR)/assoc_legendre_test: $(INCS) assoc_legendre_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/assoc_legendre_test assoc_legendre_test.cpp -lquadmath

$(BIN_DIR)/sph_legendre_test: $(INCS) sph_legendre_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/sph_legendre_test sph_legendre_test.cpp -lquadmath

$(BIN_DIR)/test_phase_iterator: test_phase_iterator.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/test_phase_iterator test_phase_iterator.cpp -lquadmath

$(BIN_DIR)/test_quadrature: test_quadrature.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/test_quadrature test_quadrature.cpp -lquadmath -lubsan

$(BIN_DIR)/test_trapezoid_integral: test_trapezoid_integral.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -I../polynomial -o $(BIN_DIR)/test_trapezoid_integral test_trapezoid_integral.cpp -lquadmath

$(BIN_DIR)/test_midpoint_integral: test_midpoint_integral.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -I../polynomial -o $(BIN_DIR)/test_midpoint_integral test_midpoint_integral.cpp -lquadmath

$(BIN_DIR)/test_simpson_integral: test_simpson_integral.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -I../polynomial -o $(BIN_DIR)/test_simpson_integral test_simpson_integral.cpp -lquadmath

$(BIN_DIR)/test_double_exp_integrate: test_double_exp_integrate.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/test_double_exp_integrate test_double_exp_integrate.cpp -lquadmath

$(BIN_DIR)/test_gauss_hermite: test_gauss_hermite.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/test_gauss_hermite test_gauss_hermite.cpp -lquadmath

$(BIN_DIR)/test_gauss_laguerre: test_gauss_laguerre.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/test_gauss_laguerre test_gauss_laguerre.cpp -lquadmath

$(BIN_DIR)/test_mapper: test_mapper.cpp include/ext/integration_transform.h
	$(CXXMAX) -Iinclude -o $(BIN_DIR)/test_mapper test_mapper.cpp -lquadmath

$(BIN_DIR)/test_composite_trapezoid_integral: test_composite_trapezoid_integral.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -I../polynomial -o $(BIN_DIR)/test_composite_trapezoid_integral test_composite_trapezoid_integral.cpp -lquadmath

$(BIN_DIR)/hermite_test: $(INCS) hermite_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/hermite_test hermite_test.cpp -lquadmath

$(BIN_DIR)/laguerre_test: $(INCS) laguerre_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/laguerre_test laguerre_test.cpp -lquadmath

$(BIN_DIR)/legendre_test: $(INCS) legendre_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/legendre_test legendre_test.cpp -lquadmath

$(BIN_DIR)/gegenbauer_test: $(INCS) gegenbauer_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/gegenbauer_test gegenbauer_test.cpp -lquadmath

$(BIN_DIR)/jacobi_test: $(INCS) jacobi_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/jacobi_test jacobi_test.cpp -lquadmath

$(BIN_DIR)/chebyshev_t_test: $(INCS) chebyshev_t_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/chebyshev_t_test chebyshev_t_test.cpp -lquadmath

$(BIN_DIR)/chebyshev_u_test: $(INCS) chebyshev_u_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/chebyshev_u_test chebyshev_u_test.cpp -lquadmath

$(BIN_DIR)/chebyshev_v_test: $(INCS) chebyshev_v_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/chebyshev_v_test chebyshev_v_test.cpp -lquadmath

$(BIN_DIR)/chebyshev_w_test: $(INCS) chebyshev_w_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/chebyshev_w_test chebyshev_w_test.cpp -lquadmath

$(BIN_DIR)/radpoly_test: $(INCS) radpoly_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/radpoly_test radpoly_test.cpp -lquadmath

$(BIN_DIR)/zernike_test: $(INCS) zernike_test.cpp
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/zernike_test zernike_test.cpp -lquadmath



$(TEST_OUT_DIR): $(TEST_OUT_DIR)
	if test ! -d $(TEST_OUT_DIR); then \
	  mkdir $(TEST_OUT_DIR); \
	fi


$(BIN_DIR):
	if test ! -d $(BIN_DIR); then \
	  mkdir $(BIN_DIR); \
	fi

clean:
	rm -f a.out
	rm -f *.stackdump
	rm -f $(BINS) $(BUILDERS)

