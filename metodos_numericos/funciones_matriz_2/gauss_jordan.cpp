#include "gauss_jordan.hpp"

void gauss_jordan( matrix_like<double> &matriz, array_like<double> &variables, array_like<double> &resultados, const int &size ){
    for( int i = 0 ; i < size ; ++i) variables[i] = resultados[i];
    for(int i=0; i<size; ++i){
        const double divide_privote = matriz[i][i];
        for(int j=i; j<size; ++j){
            matriz[i][j] /= divide_privote;
        }
        variables[i] /= divide_privote;
        for(int j=0; j < size; ++j){
            if(i == j) continue;
            const double coeficiente_eliminar = matriz[j][i];
            matriz[j][i] = 0;
            for(int k=i+1; k < size; ++k){
                matriz[j][k] -= coeficiente_eliminar * matriz[i][k];
            }
            variables[j] -= coeficiente_eliminar * variables[i] / matriz[i][i];
        }
    }
    solucion_triangular_sup(matriz, variables, resultados, size);
}
