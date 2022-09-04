#include "matrix_like/real_matrix.hpp"
#include "funcion_matriz.hpp"

void factor_cholesky(mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular);
void factor_cholesky_tridiag(mymtx::RealMatrix &matix, mymtx::RealMatrix &triangular);
void solve_jacobi(mymtx::RealMatrix &matix, mymtx::RealVector &variables, mymtx::RealVector &solutions, double *error);
void solve_gauss_seidel(mymtx::RealMatrix &matix, mymtx::RealVector &variables, mymtx::RealVector &solutions, double *error);
void solve_cholesky(mymtx::RealMatrix &cholesky_factored,mymtx::RealVector &variables, mymtx::RealVector &solutions);

