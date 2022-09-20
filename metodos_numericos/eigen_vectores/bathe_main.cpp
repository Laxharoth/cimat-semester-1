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
    ANNOUNCE_TEST("bathe eigen"){
    // mtx*=(1/(h*h));
    const unsigned amt = 2;
    mymtx::RealMatrix eigvec(SIZE,SIZE);
    mymtx::RealVector eigval(SIZE);
    
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
    /*
    ANNOUNCE_TEST("Rayleigh"){
    mymtx::RealVector V0{1, 0, -1, 0.0, 1};
    double val = 2.1;
    printm(mtx);
    rayleigh_method(mtx, V0, val);
    strm_out("val" << std::endl << val);
    strm_out("vec");
    printv(V0);
    strm_out("mtx");
    printm(mtx);
    }
    ANNOUNCE_TEST("QR"){
        auto mtx = mymtx::RealMatrix::tridiag(5,-1,2,-1);
        mymtx::RealMatrix Q(mtx.shape_y, mtx.shape_x);
        mymtx::RealMatrix R(mtx.shape_y, mtx.shape_x);
        qr_decomposition(mtx,Q,R);
        strm_out("Q");
        printm(Q);
        strm_out("R");
        printm(R);
        strm_out("QR");
        printm(Q*R);
    }
    ANNOUNCE_TEST("Gradient"){
        mymtx::RealMatrix mtx{
		{4,1,0},
		{1,3,1},
		{1,3,4}
        };
        mymtx::RealVector result{5.0,8.0,20.0};
        mymtx::RealVector expected{1.0, 1.0,4.0};
        mymtx::RealVector actual(result.size);
        conjugate_gradient(mtx,actual,result);
        strm_out("expected");
        printv(expected,*macros::out);
        strm_out("actual");
        printv(actual,*macros::out);
    }
    ANNOUNCE_TEST("Gradient Jacobi"){
        mymtx::RealMatrix mtx{
		{4,1,0},
		{1,3,1},
		{1,3,4}
        };
        mymtx::RealVector result{5.0,8.0,20.0};
        mymtx::RealVector expected{1.0, 1.0,4.0};
        mymtx::RealVector actual(result.size);
        conjugate_gradient_jacobi(mtx,actual,result);
        strm_out("expected");
        printv(expected,*macros::out);
        strm_out("actual");
        printv(actual,*macros::out);
    }
    return 0;
    //*/
}
