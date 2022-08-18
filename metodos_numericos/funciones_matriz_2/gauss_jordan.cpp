#include "gauss_jordan.hpp"

void gauss_jordan( double **matriz, double *variables, double *resultados, const int &size ){
    for(int i=0; i<size; ++i){
        const double divide_privote = matriz[i][i];
        for(int j=i; j<size; ++j){
            matriz[i][j] /= divide_privote;
        }
        resultados[i] /= divide_privote;
        for(int j=0; j < size; ++j){
            if(i == j) continue;
            const double coeficiente_eliminar = matriz[j][i];
            matriz[j][i] = 0;
            for(int k=i+1; k < size; ++k){
                matriz[j][k] -= coeficiente_eliminar * matriz[i][k];
            }
            resultados[j] -= coeficiente_eliminar * resultados[i] / matriz[i][i];
        }
    }
    memcpy(variables, resultados, sizeof(double) * size);
    solucion_triangular_sup(matriz, variables, resultados, size);
}
