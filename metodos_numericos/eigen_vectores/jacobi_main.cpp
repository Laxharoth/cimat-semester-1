#include "matrix_like/real_matrix.hpp"
#include "funcion_matriz.hpp"
#include "macros.hpp"
#include <algorithm>

#define SIZE 1000

const double h = 2.0 / (SIZE + 1);

int main(int argc, char const *argv[]){
    auto mtx = mymtx::RealMatrix::tridiag(SIZE,1,-2,1);
    // mtx*=(1/(h*h));
    mymtx::RealMatrix eigvec(SIZE,SIZE);
    mymtx::RealVector eigval(SIZE);
    ANNOUNCE_TEST("jacobi eigen")
    measure_time(
        jacobi_eigen(mtx, eigval, eigvec,10000);
    )
    mymtx::RealVector::sort(eigval);
    for(int i = 0; i < SIZE; i++)
        strm_out("Eigen"<<i+1<<":" << eigval[i])
    return 0;
}
