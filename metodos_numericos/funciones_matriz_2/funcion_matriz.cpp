#include "funcion_matriz.hpp"

void solucion_diagonal(matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result){
    const size_t size = matriz.get_shape_x();
    for(int i=0; i<size; ++i){
        incognitas[i] = result[i] / matriz[i][i];
    }
}
double determinante_diagonal(matrix_like<double> &matriz_diagonal){
    const size_t size = matriz_diagonal.get_shape_x();
    double result{1};
    for(int i=0; i<size; ++i) result*=matriz_diagonal[i][i];
    return result;
}
void inversa_diagonal(matrix_like<double> &matriz_diagonal,array_like<double> &inversa){
    const size_t size = matriz_diagonal.get_shape_x();
    for(int i=0; i<size; ++i){
        inversa[i] = 1 / matriz_diagonal[i][i];
    }
}
double determinante_triangular(matrix_like<double> &matriz_triangular){
    const size_t size = matriz_triangular.get_shape_x();
    double result{1};
    for(int i=0; i<size; ++i) result*=matriz_triangular[i][i];
    return result;
}
void solucion_triangular_inf( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result){
    const size_t size = matriz.get_shape_x();
    for(int i=0; i<size; ++i){
        incognitas[i] = result[i];
        for(int j=0; j<=i - 1; ++j){
            incognitas[i] -= matriz[i][j] * incognitas[j];
        }
        incognitas[i] /= matriz[i][i];
    }
}
void solucion_triangular_sup( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result){
    const size_t size = matriz.get_shape_x();
    for(int i = size - 1; i>=0; --i){
        incognitas[i] = result[i];
        for(int j=i+1; j< size; ++j){
            incognitas[i] -= matriz[i][j] * incognitas[j];
        }
        incognitas[i] /= matriz[i][i];
    }
}
void gauss( matrix_like<double> &matriz, array_like<double> &variables, array_like<double> &resultados, const int &size ){
    const size_t size = matriz.get_shape_x();
    for( int i = 0 ; i < size ; ++i) variables[i] = resultados[i];
    for(int i=0; i<size; ++i){
        const double divide_privote = matriz[i][i];
        for(int j=i; j<size; ++j){
            matriz[i][j] /= divide_privote;
        }
        resultados[i] /= divide_privote;
        for(int j=i+1; j < size; ++j){
            const double coeficiente_eliminar = matriz[j][i];
            matriz[j][i] = 0;
            for(int k=i+1; k < size; ++k){
                matriz[j][k] -= coeficiente_eliminar * matriz[i][k];
            }
            resultados[j] -= coeficiente_eliminar * resultados[i] / matriz[i][i];
        }
    }
    solucion_triangular_sup(matriz, variables, resultados);
}
void solucion_LDU( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result){
    auto matriz_doolittle = LDU_wrapper::from(&matriz, LDU_wrapper::DOOLITTLE);
    auto matriz_crout = LDU_wrapper::from(&matriz, LDU_wrapper::CROUT);
    mymtx::vector<double> aux1(incognitas.get_size());
    mymtx::vector<double> aux2(incognitas.get_size());
    solucion_triangular_inf( matriz_doolittle, aux1, result);
    solucion_diagonal( matriz, aux2, aux1);
    solucion_triangular_sup( matriz_crout, incognitas, aux2);
}
void solucion_crout( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result, const size_t &size){
    auto matriz_crout = LDU_wrapper::from(&matriz, LDU_wrapper::CROUT);
    mymtx::vector<double> aux1(incognitas.get_size());
    solucion_triangular_inf( matriz, aux1, result);
    solucion_triangular_sup( matriz_crout, incognitas, aux1);
}
void solucion_doolittle( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result, const size_t &size){
    auto matriz_doolittle = LDU_wrapper::from(&matriz, LDU_wrapper::DOOLITTLE);
    mymtx::vector<double> aux1(incognitas.get_size());
    solucion_triangular_inf( matriz_doolittle, aux1, result );
    solucion_triangular_sup( matriz, incognitas, aux1);
}
