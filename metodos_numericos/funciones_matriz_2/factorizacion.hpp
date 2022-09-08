#ifndef FACTORIZACION_HPP
#define FACTORIZACION_HPP

#include "matrix_like/matrix_like.tcc"
#include "matrix_like/real_matrix.hpp"
#include <exception>
#include <cmath>

using mymtx::matrix_like;
using mymtx::array_like;
using mymtx::RealMatrix;
using mymtx::RealVector;

void crout(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_superior);
void crout_tridiagonal(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_superior);
void doolittle(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_superior);
void doolittle_tridiagonal(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_superior);
void LDU_factor(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_diagonal, RealMatrix &matriz_superior);
void LDU_factor_tridiagonal(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_diagonal, RealMatrix &matriz_superior);
void factor_cholesky(mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular);
void factor_cholesky_tridiag(mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular);

class cant_factor_exception : public std::exception{
    using std::exception::exception;
    public:
    virtual const char* 
        what() 
        const throw();
};

#endif /* FACTORIZACION_HPP */
