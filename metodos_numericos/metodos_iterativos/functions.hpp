#include "matrix_like/matrix.hpp"
#include "matrix_like/funcion_matriz.hpp"

void factor_cholesky(mymtx::matrix &matix, mymtx::matrix &triangular);
void factor_cholesky_tridiag(mymtx::matrix &matix, mymtx::matrix &triangular);
void solve_jacobi(mymtx::matrix &matix, mymtx::vector &variables, mymtx::vector &solutions, double *error);
void solve_gauss_seidel(mymtx::matrix &matix, mymtx::vector &variables, mymtx::vector &solutions, double *error);
void solve_cholesky(mymtx::matrix &cholesky_factored,mymtx::vector &variables, mymtx::vector &solutions);

