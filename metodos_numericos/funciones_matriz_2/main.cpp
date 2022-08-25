#include "matrix_like/matrix_like.tcc"
#include "funcion_matriz.hpp"
#include "factorizacion.hpp"
#include "macros.hpp"

#include <fstream>

int main(int argc, char **argv){
    auto file = std::ofstream("MediumMatrixTest.txt");
	macros::out = &file;
    auto f_matriz = std::ifstream("./A.txt");
    size_t rows;
    size_t cols;
    f_matriz >> rows;
    f_matriz >> cols;
    matrix<double> matriz(rows, cols);
    for(size_t i = 0; i < rows; ++i){
        for(size_t j = 0; j < cols; ++j){
            f_matriz >> matriz[i][j];
        }
    }
    f_matriz.close();

    auto f_vector = std::ifstream("b.txt");
    
    f_matriz >> rows;
    f_matriz >> cols;

    vector<double> vec(rows);
    for(size_t i = 0; i < rows; ++i){
        f_matriz >> vec[i];
    }
    f_vector.close();
    vector<double> variables(vec.get_size());
    {
        auto cpy = matriz;
        ANNOUNCE_TEST("Metodo: Gauss")
        measure_time( gauss( cpy, variables, vec, vec.get_size() ) );
    }
    {
        auto cpy = matriz;
        ANNOUNCE_TEST("Metodo: Factorizacion Crout")
        measure_time( metodo_de_crout( cpy,cpy,cpy,cpy.get_shape_y() ) );
        ANNOUNCE_TEST("Metodo: Solucion Crout")
        measure_time( solucion_crout( cpy, variables, vec, vec.get_size() ) );
    }
    {
        auto cpy = matriz;
        ANNOUNCE_TEST("Metodo: Factorizacion Doolittle")
        measure_time( metodo_de_doolittle( cpy,cpy,cpy,cpy.get_shape_y() ) );
        ANNOUNCE_TEST("Metodo: Solucion Doolittle")
        measure_time( solucion_doolittle( cpy, variables, vec, vec.get_size() ) );
    }

    return 0;
}