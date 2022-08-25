#include "funcion_matriz.hpp"

void solucion_diagonal(matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result, size_t size){
    for(int i=0; i<size; ++i){
        incognitas[i] = result[i] / matriz[i][i];
    }
}
double determinante_diagonal(matrix_like<double> &matriz_diagonal, const size_t &size){
    double result{1};
    for(int i=0; i<size; ++i) result*=matriz_diagonal[i][i];
    return result;
}
void inversa_diagonal(matrix_like<double> &matriz_diagonal,array_like<double> &inversa, const size_t &size){
    for(int i=0; i<size; ++i){
        inversa[i] = 1 / matriz_diagonal[i][i];
    }
}
double determinante_triangular(matrix_like<double> &matriz_triangular, const size_t &size){
    double result{1};
    for(int i=0; i<size; ++i) result*=matriz_triangular[i][i];
    return result;
}
void solucion_triangular_inf( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result, const size_t &size){
    for(int i=0; i<size; ++i){
        incognitas[i] = result[i];
        for(int j=0; j<=i - 1; ++j){
            incognitas[i] -= matriz[i][j] * incognitas[j];
        }
        incognitas[i] /= matriz[i][i];
    }
}
void solucion_triangular_sup( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result, const size_t &size){
    for(int i = size - 1; i>=0; --i){
        incognitas[i] = result[i];
        for(int j=i+1; j<=i - 1; ++j){
            incognitas[i] -= matriz[i][j] * incognitas[j];
        }
        incognitas[i] /= matriz[i][i];
    }
void gauss( matrix_like<double> &matriz, array_like<double> &variables, array_like<double> &resultados, const int &size ){
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
    solucion_triangular_sup(matriz, variables, resultados, size);
}
void solucion_LDU( matrix_like<double> &matriz, array_like<double> &incognitas, array_like<double> &result, const size_t &size){
