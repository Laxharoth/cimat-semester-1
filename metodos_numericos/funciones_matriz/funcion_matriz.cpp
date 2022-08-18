#include "funcion_matriz.hpp"

void solucion_diagonal(double *matriz, double *incognitas, double *result, int size){
    for(int i=0; i<size; ++i){
        incognitas[i] = result[i] / matriz[i];
    }
}
double determinante_diagonal(double *matriz_diagonal, const int &size){
    double result{1};
    for(int i=0; i<size; ++i) result*=matriz_diagonal[i];
    return result;
}
void inversa_diagonal(double *matriz_diagonal,double *inversa, const int &size){
    for(int i=0; i<size; ++i){
        inversa[i] = 1 / matriz_diagonal[i];
    }
}
double determinante_triangular(double **matriz_triangular, const int &size){
    double result{1};
    for(int i=0; i<size; ++i) result*=matriz_triangular[i][i];
    return result;
}
void solucion_triangular_inf( double **matriz, double *incognitas, double *result, const int &size){
    for(int i=0; i<size; ++i){
        incognitas[i] = result[i];
        for(int j=0; j<=i - 1; ++j){
            incognitas[i] -= matriz[i][j] * incognitas[j];
        }
        incognitas[i] /= matriz[i][i];
    }
}
void solucion_triangular_sup( double **matriz, double *incognitas, double *result, const int &row){
    for(int i = row - 1; i>=0; --i){
        incognitas[i] = result[i];
        for(int j=i+1; j<=i - 1; ++j){
            incognitas[i] -= matriz[i][j] * incognitas[j];
        }
        incognitas[i] /= matriz[i][i];
    }
}