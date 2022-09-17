#include "matrix_like/real_matrix.hpp"
#include "funcion_matriz.hpp"
#include "macros.hpp"
#include "matrix_like/print.cpp"
#include <algorithm>

#define SIZE 5

const double h = 2.0 / (SIZE + 1);

struct eigen_index{
    int index;
    double val;
};

int main(int argc, char const *argv[]){
    // auto f = std::ofstream("bathe_out.txt");
    // macros::out = &f;
    auto mtx = mymtx::RealMatrix::tridiag(SIZE,-1,2,-1);
{
    // mtx*=(1/(h*h));
    const unsigned amt = 2;
    mymtx::RealMatrix eigvec(SIZE,SIZE);
    mymtx::RealVector eigval(SIZE);
    ANNOUNCE_TEST("bathe eigen")
    mymtx::RealMatrix I(mtx.shape_y, amt);
    mymtx::RealVector eig(amt);
    measure_time(
        bathe_subspace(mtx, I ,eig);
    )
    strm_out("vectors");
    printm( mymtx::RealMatrix::traspose(I) );
    strm_out("values");
    printv(eig);
}

//     ANNOUNCE_TEST("Ralergh")
// {
//     mymtx::RealVector V0{0, 0, -0.49, 0.62, -0.60};
//     double val = 1.5;
//     printm(mtx);
//     ralergh_method(mtx, V0, val);
//     strm_out("val" << std::endl << val);
//     strm_out("vec");
//     printv(V0);
// }
    return 0;
}
