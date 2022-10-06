#include "matrix_like/matrix.hpp"
#include "matrix_like/funcion_matriz.hpp"
#include "macros.hpp"

#include <fstream>

int main(int argc, char **argv){
    auto file = std::ofstream("solutions/MediumMatrixTest.txt");
	macros::out = &file;
    auto f_matriz = std::ifstream("./A.txt");
    size_t rows;
    size_t cols;
    f_matriz >> rows;
    f_matriz >> cols;
    mymtx::matrix matriz(rows, cols);
    for(size_t i = 0; i < rows; ++i){
        for(size_t j = 0; j < cols; ++j){
            f_matriz >> matriz[i][j];
        }
    }
    f_matriz.close();

    auto f_vector = std::ifstream("b.txt");
    
    f_matriz >> rows;
    f_matriz >> cols;

    mymtx::vector vec(rows);
    for(size_t i = 0; i < rows; ++i){
        f_matriz >> vec[i];
    }
    f_vector.close();
    mymtx::vector variables(vec.size);
    {
        auto cpy = matriz;
        ANNOUNCE_TEST("Metodo: Gauss")
        measure_time( gauss( cpy, variables, vec ) );
        out_vector( variables, "solutions/long_gauss.txt")
    }
    {
        auto cpy = matriz;
        ANNOUNCE_TEST("Metodo: Factorizacion Crout")
        measure_time( crout( cpy,cpy,cpy ) );
        ANNOUNCE_TEST("Metodo: Solucion Crout")
        measure_time( solucion_crout( cpy,variables,vec ) );
        out_vector( variables, "solutions/long_crout.txt")
    }
    {
        auto cpy = matriz;
        ANNOUNCE_TEST("Metodo: Factorizacion Doolittle")
        measure_time( doolittle( cpy,cpy,cpy ) );
        ANNOUNCE_TEST("Metodo: Solucion Doolittle")
        measure_time( solucion_doolittle( cpy, variables, vec ) );
        out_vector( variables, "solutions/long_doolittle.txt")
    }
    {
        auto cpy = matriz;
        ANNOUNCE_TEST("Metodo: Factorizacion LDU")
        measure_time( doolittle( cpy,cpy,cpy ) );
        ANNOUNCE_TEST("Metodo: Solucion LDU")
        measure_time( solucion_LDU( cpy,variables,vec ) );
        out_vector( variables, "solutions/long_ldu.txt")
    }

    return 0;
}
