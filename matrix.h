
#ifndef MATRIX_H
#define MATRIX_H 1

namespace __gnu_cxx
{

  template<typename _RandAccIter, typename _RandAccIterRHS>
    int
    _S_tridiag(std::size_t __n,
	       _RandAccIter __dsub, _RandAccIter __diag, _RandAccIter __dsup,
	       _RandAccIterRHS __b);

  template<typename _RandAccIter, typename _RandAccIterRHS>
    int
    _S_tridiag_symm(std::size_t __n,
		    _RandAccIter __diag, _RandAccIter __dsub,
		    _RandAccIterRHS __b);

} // namespace __gnu_cxx

#include "matrix.tcc"

#endif // MATRIX_H
