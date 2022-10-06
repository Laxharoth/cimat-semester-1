#include "matrix/funcion_matriz.hpp"
#include "matrix/real_matrix.hpp"
#include <cmath>
#include <cstdio>
#include <cmath>
#include <ostream>
#include <iostream>

#include "print_time.hpp"

using namespace mymtx;

#define OUT_DIR "mtx.out/"

void printm(const mymtx::RealMatrix &m, std::ostream& out){
    for (size_t i = 0; i < m.shape_y; ++i)
    {
        for( auto j = m.begin(i); j < m.end(i); ++j){
            out << *j << " ";
        }
        out << std::endl;
    }
    
}
double norm(const RealMatrix &m);
const int size = 5000;
int main(int argc, const char** argv) {
    //generar matriz pendiagonal
    double coef[] = {-4,-8,40,-8,-4};
    auto A = RealMatrix::pendiag(size, coef);
    RealMatrix L(A.shape_y,A.shape_y);
    //generar vector solucion 
    RealVector a(size);
    RealVector b(size);
    b=size;
    b[0]=20;
    b[1]=50;
    b[b.size-2] = 50;
    b[b.size-1] = 20;
    measure_time(
        factor_cholesky_as_band(A, L,2);
        solve_cholesky(L, a, b);
    )
    RealVector::fwrite(OUT_DIR "sol_sistema_ec.vec", a);
    //Error
    printf("Error de la solucion %lf\n", (A*a - b).distance() );
    return 0;
}
double norm(const RealMatrix &m){
    double sum = 0;
    for(int i = 0; i < m.shape_y; ++i){
        for(auto j = m.begin(i); j < m.end(i); ++j){
            sum += *j**j;
        }
    }
    return std::sqrt(sum);
}