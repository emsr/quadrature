
#SUFFIX = _tr29124
#SUFFIX = _specfun
CXX_INST_DIR = $(HOME)/bin$(SUFFIX)
ifeq ("$(wildcard $(CXX_INST_DIR))","")
  SUFFIX = 
  CXX_INST_DIR = $(HOME)/bin
  ifeq ("$(wildcard $(CXX_INST_DIR))","")
    ifneq ($(wildcard "/mingw64"),"")
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
CXXMAX = $(CXX20)
CXX_INC_DIR = $(CXX_INST_DIR)/include/c++/8.0.0/bits
CXX_LIB_DIR = $(CXX_INST_DIR)/lib64

#WRAPPER_LIB_DIR = ../wrappers/release
WRAPPER_LIB_DIR = ../wrappers/debug
WRAPPER_LIBS = -L$(WRAPPER_LIB_DIR) -lwrap_burkhardt -lgfortran

INC_DIR = include/ext
INCLUDES =  -I../include -Iinclude -I../polynomial/include

OUTPUT_DIR = output

OBJ_DIR = obj
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
  $(BIN_DIR)/build_clenshaw_curtis \
  $(BIN_DIR)/build_double_exp_rules \
  $(BIN_DIR)/test_gauss_kronrod_rule


all: $(OBJ_DIR) $(BINS)


ortho_test: $(OUTPUT_DIR)
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/legendre_test > $(OUTPUT_DIR)/legendre_test.txt 2> $(OUTPUT_DIR)/legendre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/chebyshev_t_test > $(OUTPUT_DIR)/chebyshev_t_test.txt 2> $(OUTPUT_DIR)/chebyshev_t_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/chebyshev_u_test > $(OUTPUT_DIR)/chebyshev_u_test.txt 2> $(OUTPUT_DIR)/chebyshev_u_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/chebyshev_v_test > $(OUTPUT_DIR)/chebyshev_v_test.txt 2> $(OUTPUT_DIR)/chebyshev_v_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/chebyshev_w_test > $(OUTPUT_DIR)/chebyshev_w_test.txt 2> $(OUTPUT_DIR)/chebyshev_w_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/gegenbauer_test > $(OUTPUT_DIR)/gegenbauer_test.txt 2> $(OUTPUT_DIR)/gegenbauer_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/jacobi_test > $(OUTPUT_DIR)/jacobi_test.txt 2> $(OUTPUT_DIR)/jacobi_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/assoc_laguerre_test > $(OUTPUT_DIR)/assoc_laguerre_test.txt 2> $(OUTPUT_DIR)/assoc_laguerre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/assoc_legendre_test > $(OUTPUT_DIR)/assoc_legendre_test.txt 2> $(OUTPUT_DIR)/assoc_legendre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/hermite_test > $(OUTPUT_DIR)/hermite_test.txt 2> $(OUTPUT_DIR)/hermite_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/sph_legendre_test > $(OUTPUT_DIR)/sph_legendre_test.txt 2> $(OUTPUT_DIR)/sph_legendre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/laguerre_test > $(OUTPUT_DIR)/laguerre_test.txt 2> $(OUTPUT_DIR)/laguerre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/radpoly_test > $(OUTPUT_DIR)/radpoly_test.txt 2> $(OUTPUT_DIR)/radpoly_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/zernike_test > $(OUTPUT_DIR)/zernike_test.txt 2> $(OUTPUT_DIR)/zernike_test.err


test: $(OUTPUT_DIR) $(BIN_DIR)/test_quadrature
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/build_double_exp_rules > $(OUTPUT_DIR)/build_double_exp_rules.txt 2> $(OUTPUT_DIR)/build_double_exp_rules.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/build_clenshaw_curtis > $(OUTPUT_DIR)/build_clenshaw_curtis.txt 2> $(OUTPUT_DIR)/build_clenshaw_curtis.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_gauss_kronrod_rule > $(OUTPUT_DIR)/test_gauss_kronrod_rule.txt 2> $(OUTPUT_DIR)/test_gauss_kronrod_rule.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_double_exp_integrate > $(OUTPUT_DIR)/test_double_exp_integrate.txt 2> $(OUTPUT_DIR)/test_double_exp_integrate.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_trapezoid_integral > $(OUTPUT_DIR)/test_trapezoid_integral.txt 2> $(OUTPUT_DIR)/test_trapezoid_integral.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_midpoint_integral > $(OUTPUT_DIR)/test_midpoint_integral.txt 2> $(OUTPUT_DIR)/test_midpoint_integral.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_simpson_integral > $(OUTPUT_DIR)/test_simpson_integral.txt 2> $(OUTPUT_DIR)/test_simpson_integral.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_phase_iterator > $(OUTPUT_DIR)/test_phase_iterator.txt 2> $(OUTPUT_DIR)/test_phase_iterator.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_mapper > $(OUTPUT_DIR)/test_mapper.txt 2> $(OUTPUT_DIR)/test_mapper.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_quadrature > $(OUTPUT_DIR)/test_quadrature.txt 2> $(OUTPUT_DIR)/test_quadrature.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_gauss_hermite > $(OUTPUT_DIR)/test_gauss_hermite.txt 2> $(OUTPUT_DIR)/test_gauss_hermite.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_gauss_laguerre > $(OUTPUT_DIR)/test_gauss_laguerre.txt 2> $(OUTPUT_DIR)/test_gauss_laguerre.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAPPER_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_composite_trapezoid_integral > $(OUTPUT_DIR)/test_composite_trapezoid_integral.txt 2> $(OUTPUT_DIR)/test_composite_trapezoid_integral.err


docs:
	#rm -rf html/*
	#rm -rf latex/*
	doxygen
	cd docs/latex && make


# Binaries...

$(BIN_DIR)/build_clenshaw_curtis: $(BIN_DIR) $(OBJ_DIR)/build_clenshaw_curtis.o
	$(CXXMAX) -o $(BIN_DIR)/build_clenshaw_curtis $(OBJ_DIR)/build_clenshaw_curtis.o -lquadmath $(WRAPPER_LIBS)

$(BIN_DIR)/build_double_exp_rules: $(BIN_DIR) $(OBJ_DIR)/build_double_exp_rules.o
	$(CXXMAX) -o $(BIN_DIR)/build_double_exp_rules $(OBJ_DIR)/build_double_exp_rules.o -lquadmath $(WRAPPER_LIBS)

$(BIN_DIR)/test_gauss_kronrod_rule: $(BIN_DIR) $(OBJ_DIR)/test_gauss_kronrod_rule.o
	$(CXXMAX) -o $(BIN_DIR)/test_gauss_kronrod_rule $(OBJ_DIR)/test_gauss_kronrod_rule.o -lquadmath

$(BIN_DIR)/assoc_laguerre_test: $(BIN_DIR) $(OBJ_DIR)/assoc_laguerre_test.o
	$(CXXMAX) -o $(BIN_DIR)/assoc_laguerre_test $(OBJ_DIR)/assoc_laguerre_test.o -lquadmath

$(BIN_DIR)/assoc_legendre_test: $(BIN_DIR) $(OBJ_DIR)/assoc_legendre_test.o
	$(CXXMAX) -o $(BIN_DIR)/assoc_legendre_test $(OBJ_DIR)/assoc_legendre_test.o -lquadmath

$(BIN_DIR)/sph_legendre_test: $(BIN_DIR) $(OBJ_DIR)/sph_legendre_test.o
	$(CXXMAX) -o $(BIN_DIR)/sph_legendre_test $(OBJ_DIR)/sph_legendre_test.o -lquadmath

$(BIN_DIR)/test_phase_iterator: $(BIN_DIR) test_phase_iterator.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/test_phase_iterator test_phase_iterator.cpp -lquadmath

$(BIN_DIR)/test_quadrature: $(BIN_DIR) test_quadrature.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/test_quadrature test_quadrature.cpp -lquadmath -lubsan

$(BIN_DIR)/test_trapezoid_integral: $(BIN_DIR) test_trapezoid_integral.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -I../polynomial -o $(BIN_DIR)/test_trapezoid_integral test_trapezoid_integral.cpp -lquadmath

$(BIN_DIR)/test_midpoint_integral: $(BIN_DIR) test_midpoint_integral.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -I../polynomial -o $(BIN_DIR)/test_midpoint_integral test_midpoint_integral.cpp -lquadmath

$(BIN_DIR)/test_simpson_integral: $(BIN_DIR) test_simpson_integral.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -I../polynomial -o $(BIN_DIR)/test_simpson_integral test_simpson_integral.cpp -lquadmath

$(BIN_DIR)/test_double_exp_integrate: $(BIN_DIR) test_double_exp_integrate.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/test_double_exp_integrate test_double_exp_integrate.cpp -lquadmath

$(BIN_DIR)/test_gauss_hermite: $(BIN_DIR) test_gauss_hermite.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/test_gauss_hermite test_gauss_hermite.cpp -lquadmath

$(BIN_DIR)/test_gauss_laguerre: $(BIN_DIR) test_gauss_laguerre.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -o $(BIN_DIR)/test_gauss_laguerre test_gauss_laguerre.cpp -lquadmath

$(BIN_DIR)/test_mapper: $(BIN_DIR) test_mapper.cpp include/ext/integration_transform.h
	$(CXXMAX) -Iinclude -o $(BIN_DIR)/test_mapper test_mapper.cpp -lquadmath

$(BIN_DIR)/test_composite_trapezoid_integral: $(BIN_DIR) test_composite_trapezoid_integral.cpp $(INCS)
	$(CXXMAX) $(INCLUDES) -I../polynomial -o $(BIN_DIR)/test_composite_trapezoid_integral test_composite_trapezoid_integral.cpp -lquadmath

$(BIN_DIR)/hermite_test: $(BIN_DIR) $(OBJ_DIR)/hermite_test.o
	$(CXXMAX) -o $(BIN_DIR)/hermite_test $(OBJ_DIR)/hermite_test.o -lquadmath

$(BIN_DIR)/laguerre_test: $(BIN_DIR) $(OBJ_DIR)/laguerre_test.o
	$(CXXMAX) -o $(BIN_DIR)/laguerre_test $(OBJ_DIR)/laguerre_test.o -lquadmath

$(BIN_DIR)/legendre_test: $(BIN_DIR) $(OBJ_DIR)/legendre_test.o
	$(CXXMAX) -o $(BIN_DIR)/legendre_test $(OBJ_DIR)/legendre_test.o -lquadmath

$(BIN_DIR)/gegenbauer_test: $(BIN_DIR) $(OBJ_DIR)/gegenbauer_test.o
	$(CXXMAX) -o $(BIN_DIR)/gegenbauer_test $(OBJ_DIR)/gegenbauer_test.o -lquadmath

$(BIN_DIR)/jacobi_test: $(BIN_DIR) $(OBJ_DIR)/jacobi_test.o
	$(CXXMAX) -o $(BIN_DIR)/jacobi_test $(OBJ_DIR)/jacobi_test.o -lquadmath

$(BIN_DIR)/chebyshev_t_test: $(BIN_DIR) $(OBJ_DIR)/chebyshev_t_test.o
	$(CXXMAX) -o $(BIN_DIR)/chebyshev_t_test $(OBJ_DIR)/chebyshev_t_test.o -lquadmath

$(BIN_DIR)/chebyshev_u_test: $(BIN_DIR) $(OBJ_DIR)/chebyshev_u_test.o
	$(CXXMAX) -o $(BIN_DIR)/chebyshev_u_test $(OBJ_DIR)/chebyshev_u_test.o -lquadmath

$(BIN_DIR)/chebyshev_v_test: $(BIN_DIR) $(OBJ_DIR)/chebyshev_v_test.o
	$(CXXMAX) -o $(BIN_DIR)/chebyshev_v_test $(OBJ_DIR)/chebyshev_v_test.o -lquadmath

$(BIN_DIR)/chebyshev_w_test: $(BIN_DIR) $(OBJ_DIR)/chebyshev_w_test.o
	$(CXXMAX) -o $(BIN_DIR)/chebyshev_w_test $(OBJ_DIR)/chebyshev_w_test.o -lquadmath

$(BIN_DIR)/radpoly_test: $(BIN_DIR) $(OBJ_DIR)/radpoly_test.o
	$(CXXMAX) -o $(BIN_DIR)/radpoly_test $(OBJ_DIR)/radpoly_test.o -lquadmath

$(BIN_DIR)/zernike_test: $(BIN_DIR) $(OBJ_DIR)/zernike_test.o
	$(CXXMAX) -o $(BIN_DIR)/zernike_test $(OBJ_DIR)/zernike_test.o -lquadmath

# Objects...

$(OBJ_DIR)/build_clenshaw_curtis.o: $(OBJ_DIR) $(INCS) build_clenshaw_curtis.cpp
	$(CXXMAX) -c $(INCLUDES) -I../wrappers -o $(OBJ_DIR)/build_clenshaw_curtis.o build_clenshaw_curtis.cpp

$(OBJ_DIR)/build_double_exp_rules.o: $(OBJ_DIR) $(INCS) build_double_exp_rules.cpp
	$(CXXMAX) -c $(INCLUDES) -I../wrappers -o $(OBJ_DIR)/build_double_exp_rules.o build_double_exp_rules.cpp

$(OBJ_DIR)/test_gauss_kronrod_rule.o: $(OBJ_DIR) $(INCS) test_gauss_kronrod_rule.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/test_gauss_kronrod_rule.o test_gauss_kronrod_rule.cpp

$(OBJ_DIR)/assoc_laguerre_test.o: $(OBJ_DIR) $(INCS) assoc_laguerre_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/assoc_laguerre_test.o assoc_laguerre_test.cpp

$(OBJ_DIR)/assoc_legendre_test.o: $(OBJ_DIR) $(INCS) assoc_legendre_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/assoc_legendre_test.o assoc_legendre_test.cpp

$(OBJ_DIR)/sph_legendre_test.o: $(OBJ_DIR) $(INCS) sph_legendre_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/sph_legendre_test.o sph_legendre_test.cpp

$(OBJ_DIR)/hermite_test.o: $(OBJ_DIR) $(INCS) hermite_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/hermite_test.o hermite_test.cpp

$(OBJ_DIR)/laguerre_test.o: $(OBJ_DIR) $(INCS) laguerre_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/laguerre_test.o laguerre_test.cpp

$(OBJ_DIR)/legendre_test.o: $(OBJ_DIR) $(INCS) legendre_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/legendre_test.o legendre_test.cpp

$(OBJ_DIR)/gegenbauer_test.o: $(OBJ_DIR) $(INCS) gegenbauer_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/gegenbauer_test.o gegenbauer_test.cpp

$(OBJ_DIR)/jacobi_test.o: $(OBJ_DIR) $(INCS) jacobi_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/jacobi_test.o jacobi_test.cpp

$(OBJ_DIR)/chebyshev_t_test.o: $(OBJ_DIR) $(INCS) chebyshev_t_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/chebyshev_t_test.o chebyshev_t_test.cpp

$(OBJ_DIR)/chebyshev_u_test.o: $(OBJ_DIR) $(INCS) chebyshev_u_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/chebyshev_u_test.o chebyshev_u_test.cpp

$(OBJ_DIR)/chebyshev_v_test.o: $(OBJ_DIR) $(INCS) chebyshev_v_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/chebyshev_v_test.o chebyshev_v_test.cpp

$(OBJ_DIR)/chebyshev_w_test.o: $(OBJ_DIR) $(INCS) chebyshev_w_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/chebyshev_w_test.o chebyshev_w_test.cpp

$(OBJ_DIR)/radpoly_test.o: $(OBJ_DIR) $(INCS) radpoly_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/radpoly_test.o radpoly_test.cpp

$(OBJ_DIR)/zernike_test.o: $(OBJ_DIR) $(INCS) zernike_test.cpp
	$(CXXMAX) -c $(INCLUDES) -o $(OBJ_DIR)/zernike_test.o zernike_test.cpp


$(OUTPUT_DIR): $(OUTPUT_DIR)
	if test ! -d $(OUTPUT_DIR); then \
	  mkdir $(OUTPUT_DIR); \
	fi


$(OBJ_DIR):
	if test ! -d $(OBJ_DIR); then \
	  mkdir $(OBJ_DIR); \
	fi


$(BIN_DIR):
	if test ! -d $(BIN_DIR); then \
	  mkdir $(BIN_DIR); \
	fi

clean:
	rm -f a.out
	rm -f *.stackdump
	rm -rf $(OBJ_DIR)/*
	rm -f $(BINS)

