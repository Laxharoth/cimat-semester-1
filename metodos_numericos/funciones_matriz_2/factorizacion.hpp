#ifndef FACTORIZACION_HPP
#define FACTORIZACION_HPP

#include "matrix_like/real_matrix.hpp"
#include <exception>
#include <cmath>

using mymtx::RealMatrix;

void crout(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_superior);
void crout_as_band(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx,int height, int width);
void crout_tridiagonal(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_superior);
void doolittle(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_superior);
void doolittle_as_band(RealMatrix &A_mtx, RealMatrix &L_mtx, RealMatrix &U_mtx,int height, int width);
void doolittle_tridiagonal(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_superior);
void LDU_factor(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_diagonal, RealMatrix &matriz_superior);
void LDU_factor_tridiagonal(RealMatrix &matriz, RealMatrix &matriz_inferior, RealMatrix &matriz_diagonal, RealMatrix &matriz_superior);
void factor_cholesky(const mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular);
void factor_cholesky_as_band(const mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular, int heigh);
void factor_cholesky_tridiag(mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular);
void qr_decomposition(const mymtx::RealMatrix& A, mymtx::RealMatrix&Q, mymtx::RealMatrix&R);

class cant_factor_exception : public std::exception{
    using std::exception::exception;
    public:
    virtual const char* 
        what() 
        const throw();
};

#endif /* FACTORIZACION_HPP */
