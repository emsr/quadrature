
#SUFFIX = _tr29124
SUFFIX = _specfun
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
OPT = -g -fsanitize=signed-integer-overflow -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fsanitize=alignment
GCC = $(CXX_INST_DIR)/bin/gcc $(OPT) -Wall -Wextra
CXX = $(CXX_INST_DIR)/bin/g++ -std=gnu++14 $(OPT) -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -Wall -Wextra -Wno-psabi -I..
CXX17 = $(CXX_INST_DIR)/bin/g++ -std=gnu++17 -fconcepts $(OPT) -Wall -Wextra -Wno-psabi -I..
CXX_INC_DIR = $(CXX_INST_DIR)/include/c++/8.0.0/bits
CXX_LIB_DIR = $(CXX_INST_DIR)/lib64

OBJ_DIR = obj
BIN_DIR = .

BINS = \
  $(BIN_DIR)/test_phase_iterator \
  $(BIN_DIR)/test_quadrature \
  $(BIN_DIR)/test_trapezoid_integral \
  $(BIN_DIR)/test_midpoint_integral \
  $(BIN_DIR)/test_double_exp_integrate \
  $(BIN_DIR)/test_gauss_hermite \
  $(BIN_DIR)/test_mapper \
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
  $(BIN_DIR)/zernike_test


all: $(OBJ_DIR) $(BINS)


check:
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/assoc_laguerre_test > assoc_laguerre_test.txt 2> assoc_laguerre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/assoc_legendre_test > assoc_legendre_test.txt 2> assoc_legendre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/sph_legendre_test > sph_legendre_test.txt 2> sph_legendre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/chebyshev_t_test > chebyshev_t_test.txt 2> chebyshev_t_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/chebyshev_u_test > chebyshev_u_test.txt 2> chebyshev_u_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/chebyshev_v_test > chebyshev_v_test.txt 2> chebyshev_v_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/chebyshev_w_test > chebyshev_w_test.txt 2> chebyshev_w_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/gegenbauer_test > gegenbauer_test.txt 2> gegenbauer_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/hermite_test > hermite_test.txt 2> hermite_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/jacobi_test > jacobi_test.txt 2> jacobi_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/laguerre_test > laguerre_test.txt 2> laguerre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/legendre_test > legendre_test.txt 2> legendre_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/radpoly_test > radpoly_test.txt 2> radpoly_test.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/zernike_test > zernike_test.txt 2> zernike_test.err


test: $(BIN_DIR)/test_quadrature
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_phase_iterator > test_phase_iterator.txt 2> test_phase_iterator.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_mapper > test_mapper.txt 2> test_mapper.err
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_quadrature > test_quadrature.txt 2> test_quadrature.err


# Binaries...

$(BIN_DIR)/assoc_laguerre_test: $(OBJ_DIR)/assoc_laguerre_test.o
	$(CXX17) -o $(BIN_DIR)/assoc_laguerre_test $(OBJ_DIR)/assoc_laguerre_test.o -lquadmath

$(BIN_DIR)/assoc_legendre_test: $(OBJ_DIR)/assoc_legendre_test.o
	$(CXX17) -o $(BIN_DIR)/assoc_legendre_test $(OBJ_DIR)/assoc_legendre_test.o -lquadmath

$(BIN_DIR)/sph_legendre_test: $(OBJ_DIR)/sph_legendre_test.o
	$(CXX17) -o $(BIN_DIR)/sph_legendre_test $(OBJ_DIR)/sph_legendre_test.o -lquadmath

$(BIN_DIR)/test_phase_iterator: test_phase_iterator.cpp *.h *.tcc
	$(CXX17) -o $(BIN_DIR)/test_phase_iterator test_phase_iterator.cpp -lquadmath

$(BIN_DIR)/test_quadrature: test_quadrature.cpp *.h *.tcc
	$(CXX17) -o $(BIN_DIR)/test_quadrature test_quadrature.cpp -lquadmath -lubsan

$(BIN_DIR)/test_trapezoid_integral: test_trapezoid_integral.cpp trapezoid_integral.h trapezoid_integral.tcc
	$(CXX17) -o $(BIN_DIR)/test_trapezoid_integral test_trapezoid_integral.cpp -lquadmath

$(BIN_DIR)/test_simpson_integral: test_simpson_integral.cpp simpson_integral.h simpson_integral.tcc
	$(CXX17) -o $(BIN_DIR)/test_simpson_integral test_simpson_integral.cpp -lquadmath

$(BIN_DIR)/test_midpoint_integral: test_midpoint_integral.cpp midpoint_integral.h midpoint_integral.tcc
	$(CXX17) -o $(BIN_DIR)/test_midpoint_integral test_midpoint_integral.cpp -lquadmath

$(BIN_DIR)/test_double_exp_integrate: test_double_exp_integrate.cpp double_exp_integrate.tcc
	$(CXX17) -o $(BIN_DIR)/test_double_exp_integrate test_double_exp_integrate.cpp -lquadmath

$(BIN_DIR)/test_gauss_hermite: test_gauss_hermite.cpp gauss_hermite_integrate.h
	$(CXX17) -o $(BIN_DIR)/test_gauss_hermite test_gauss_hermite.cpp -lquadmath

$(BIN_DIR)/test_mapper: test_mapper.cpp integration_transform.h
	$(CXX17) -o $(BIN_DIR)/test_mapper test_mapper.cpp -lquadmath

$(BIN_DIR)/hermite_test: $(OBJ_DIR)/hermite_test.o
	$(CXX17) -o $(BIN_DIR)/hermite_test $(OBJ_DIR)/hermite_test.o -lquadmath

$(BIN_DIR)/laguerre_test: $(OBJ_DIR)/laguerre_test.o
	$(CXX17) -o $(BIN_DIR)/laguerre_test $(OBJ_DIR)/laguerre_test.o -lquadmath

$(BIN_DIR)/legendre_test: $(OBJ_DIR)/legendre_test.o
	$(CXX17) -o $(BIN_DIR)/legendre_test $(OBJ_DIR)/legendre_test.o -lquadmath

$(BIN_DIR)/gegenbauer_test: $(OBJ_DIR)/gegenbauer_test.o
	$(CXX17) -o $(BIN_DIR)/gegenbauer_test $(OBJ_DIR)/gegenbauer_test.o -lquadmath

$(BIN_DIR)/jacobi_test: $(OBJ_DIR)/jacobi_test.o
	$(CXX17) -o $(BIN_DIR)/jacobi_test $(OBJ_DIR)/jacobi_test.o -lquadmath

$(BIN_DIR)/chebyshev_t_test: $(OBJ_DIR)/chebyshev_t_test.o
	$(CXX17) -o $(BIN_DIR)/chebyshev_t_test $(OBJ_DIR)/chebyshev_t_test.o -lquadmath

$(BIN_DIR)/chebyshev_u_test: $(OBJ_DIR)/chebyshev_u_test.o
	$(CXX17) -o $(BIN_DIR)/chebyshev_u_test $(OBJ_DIR)/chebyshev_u_test.o -lquadmath

$(BIN_DIR)/chebyshev_v_test: $(OBJ_DIR)/chebyshev_v_test.o
	$(CXX17) -o $(BIN_DIR)/chebyshev_v_test $(OBJ_DIR)/chebyshev_v_test.o -lquadmath

$(BIN_DIR)/chebyshev_w_test: $(OBJ_DIR)/chebyshev_w_test.o
	$(CXX17) -o $(BIN_DIR)/chebyshev_w_test $(OBJ_DIR)/chebyshev_w_test.o -lquadmath

$(BIN_DIR)/radpoly_test: $(OBJ_DIR)/radpoly_test.o
	$(CXX17) -o $(BIN_DIR)/radpoly_test $(OBJ_DIR)/radpoly_test.o -lquadmath

$(BIN_DIR)/zernike_test: $(OBJ_DIR)/zernike_test.o
	$(CXX17) -o $(BIN_DIR)/zernike_test $(OBJ_DIR)/zernike_test.o -lquadmath

# Objects...

$(OBJ_DIR)/assoc_laguerre_test.o: *.h *.tcc assoc_laguerre_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/assoc_laguerre_test.o assoc_laguerre_test.cpp

$(OBJ_DIR)/assoc_legendre_test.o: *.h *.tcc assoc_legendre_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/assoc_legendre_test.o assoc_legendre_test.cpp

$(OBJ_DIR)/sph_legendre_test.o: *.h *.tcc sph_legendre_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/sph_legendre_test.o sph_legendre_test.cpp

$(OBJ_DIR)/hermite_test.o: *.h *.tcc hermite_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/hermite_test.o hermite_test.cpp

$(OBJ_DIR)/laguerre_test.o: *.h *.tcc laguerre_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/laguerre_test.o laguerre_test.cpp

$(OBJ_DIR)/legendre_test.o: *.h *.tcc legendre_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/legendre_test.o legendre_test.cpp

$(OBJ_DIR)/gegenbauer_test.o: *.h *.tcc gegenbauer_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/gegenbauer_test.o gegenbauer_test.cpp

$(OBJ_DIR)/jacobi_test.o: *.h *.tcc jacobi_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/jacobi_test.o jacobi_test.cpp

$(OBJ_DIR)/chebyshev_t_test.o: *.h *.tcc chebyshev_t_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/chebyshev_t_test.o chebyshev_t_test.cpp

$(OBJ_DIR)/chebyshev_u_test.o: *.h *.tcc chebyshev_u_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/chebyshev_u_test.o chebyshev_u_test.cpp

$(OBJ_DIR)/chebyshev_v_test.o: *.h *.tcc chebyshev_v_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/chebyshev_v_test.o chebyshev_v_test.cpp

$(OBJ_DIR)/chebyshev_w_test.o: *.h *.tcc chebyshev_w_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/chebyshev_w_test.o chebyshev_w_test.cpp

$(OBJ_DIR)/radpoly_test.o: *.h *.tcc radpoly_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/radpoly_test.o radpoly_test.cpp

$(OBJ_DIR)/zernike_test.o: *.h *.tcc zernike_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/zernike_test.o zernike_test.cpp


$(OBJ_DIR): $(OUT_DIR)
	if test ! -d $(OBJ_DIR); then \
	  mkdir $(OBJ_DIR); \
	fi

clean:
	rm -rf $(OBJ_DIR)/*
	rm -f $(BINS)

