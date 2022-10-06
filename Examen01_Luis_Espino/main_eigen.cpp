#include "matrix/funcion_matriz.hpp"
#include "matrix/real_matrix.hpp"
#include <cmath>
#include <cstdio>
#include <cmath>
#include <ostream>

#include "print_time.hpp"

using namespace mymtx;
#include <iostream>

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
    
    //eigen vectores
    
    RealMatrix vectores(10,size);
    RealVector valores(10);
    RealVector V1(size); double val;
    // 10 mayores vectores
    std::cout << "10 mayores vectores" << std::endl;
    measure_time(
    power_iteration(A, V1, 1E-5, val, 10, &vectores, &valores,10000);
    )
    RealVector::fwrite( "mayores_valores.vec", valores);
    RealMatrix::fwrite( "mayores_vectores.mtx", vectores);

    // 10 menores vectores
    std::cout << "10 menores vectores" << std::endl;
    measure_time(
    inverse_power_iteration(A, V1, 1E-5, val, 10, &vectores, &valores,10000);
    )
    RealVector::fwrite( "menores_valores.vec", valores);
    RealMatrix::fwrite( "menores_vectores.mtx", vectores);
    
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