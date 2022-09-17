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
    // auto f = std::ofstream("jacobi_out.txt");
    // macros::out = &f;

    auto mtx = mymtx::RealMatrix::tridiag(SIZE,1,-2,1);
    mtx*=-1;
    mymtx::RealMatrix eigvec(SIZE,SIZE);
    mymtx::RealVector eigval(SIZE);
    ANNOUNCE_TEST("jacobi eigen")
    measure_time(
        jacobi_eigen(mtx, eigval, eigvec,10000);
    )
    eigen_index indexs[SIZE];
    for(int i = 0; i < SIZE; i++){
        indexs[i].index = i;
        indexs[i].val = eigval[i];
    }
    eigvec = mymtx::RealMatrix::traspose(eigvec);
    std::sort(indexs,indexs + SIZE, [](const eigen_index&a, const eigen_index&b){ return std::abs(a.val) < std::abs(b.val); } );
    strm_out( "min value:" << eigval[indexs[0].index] );
    strm_out( "max value:" << eigval[indexs[SIZE - 1].index] );
    mymtx::RealVector::fwrite("min.vec",eigvec[indexs[0].index]);
    mymtx::RealVector::fwrite("max.vec",eigvec[indexs[SIZE - 1].index]);
    if(SIZE < 10){
        printf("values\n");
        printv(eigval);
        printf("vectors\n");
        printm( mymtx::MatrixTraspose(eigvec) );
    }
    // strm_out("");
    // printm(mtx);
    return 0;
}
