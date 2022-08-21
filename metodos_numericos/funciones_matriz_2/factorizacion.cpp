#include "matrix_like/matrix_like.tcc"
#include "factorizacion.hpp"

void metodo_de_crout(matrix_like<double> &matriz, matrix_like<double> &matriz_inferior, matrix_like<double> &matriz_superior, const int &size){
    factorizacion_LDU(matriz, matriz_inferior, matriz_inferior, matriz_superior, size);
}
void metodo_de_doolittle(matrix_like<double> &matriz, matrix_like<double> &matriz_inferior, matrix_like<double> &matriz_superior, const int &size){
    auto calcular_factor_superior = [&](const int &i, const int &j){
        matriz_superior[i][j] = matriz[i][j];
        for(int k = 0; k <=i-1; ++k)
            matriz_superior[i][j] -= matriz_inferior[i][k] * matriz_superior[k][j]; 
    };
    auto calcular_factor_inferior = [&](const int &i, const int &j){
        matriz_inferior[i][j] = matriz[i][j];
        for(int k = 0; k <=j-1; ++k)
            matriz_inferior[i][j] -= matriz_inferior[i][k] * matriz_superior[k][j]; 
        matriz_inferior[i][j] /= matriz_superior[j][j];
    };
    for(int j=0; j<size; ++j){
        for(int i=0; i<j; ++i){
            calcular_factor_superior(i,j);
            calcular_factor_inferior(j,i);
        }
        calcular_factor_superior(j,j);
        if(matriz_inferior[j][j] == 0) throw cant_factor_exception();
    }
}
void factorizacion_LDU(matrix_like<double> &matriz, matrix_like<double> &matriz_inferior, matrix_like<double> &matriz_diagonal, matrix_like<double> &matriz_superior, const int &size){
    auto calcular_factor_inferior = [&](const int &i, const int &j, matrix_like<double> &matriz_objetivo){
        matriz_objetivo[i][j] = matriz[i][j];
        for(int k = 0; k <= j-1; ++k)
            matriz_objetivo[i][j] -= matriz_inferior[i][k] * matriz_superior[k][j]; 
    };
    auto calcular_factor_superior = [&](const int &i, const int &j){
        matriz_superior[i][j] = matriz[i][j];
        for(int k = 0; k <= i-1; ++k)
            matriz_superior[i][j] -= matriz_inferior[i][k] * matriz_superior[k][j];         
        matriz_superior[i][j] /= matriz_inferior[i][i];
    };
    for(int i=0; i<size; ++i){
        for(int j=0; j<i; ++j){
            calcular_factor_inferior(i,j,matriz_inferior);
            calcular_factor_superior(j,i);
        }
        calcular_factor_inferior(i,i,matriz_diagonal);
        if(matriz_inferior[i][i] == 0) throw cant_factor_exception();
    }
}
const char* cant_factor_exception::what() const throw(){
    return "zero found in diagonal";
}

