#include "factorizacion.hpp"

void metodo_de_crout(double **matriz, double **matriz_inferior, double **matriz_superior, const int &size){
    auto obtener_suma_inferior = [&](const int &i, const int &j){
        double sum{}; 
        for(int k = 0; k <= j-1; ++k)
            sum += matriz_inferior[i][k] * matriz_superior[k][j]; 
        return sum;
    };
    auto obtener_suma_superior = [&](const int &i, const int &j){
        double sum{}; 
        for(int k = 0; k <= i-1; ++k)
            sum += matriz_inferior[i][k] * matriz_superior[k][j];         
        return sum;
    };
    for(int i=0; i<size; ++i){
        for(int j=0; j<i; ++j){
            matriz_inferior[i][j] = matriz[i][j] - obtener_suma_inferior(i,j);
            matriz_superior[j][i] = (matriz[j][i] - obtener_suma_superior(j,i)) / matriz_inferior[j][j];            
        }
        matriz_inferior[i][i] = matriz[i][i] - obtener_suma_inferior(i,i);
        if(matriz_inferior[i][i] == 0) throw cant_factor_exception();
    }
}
void metodo_de_doolittle(double **matriz, double **matriz_inferior, double **matriz_superior, const int &size){
    auto obtener_suma_superior = [&](const int &i, const int &j){
        double sum{}; 
        for(int k = 0; k <=i-1; ++k)
            sum += matriz_inferior[i][k] * matriz_superior[k][j]; 
        return sum;
    };
    auto obtener_suma_inferior = [&](const int &i, const int &j){
        double sum{}; 
        for(int k = 0; k <=j-1; ++k)
            sum += matriz_inferior[i][k] * matriz_superior[k][j]; 
        return sum;
    };
    for(int j=0; j<size; ++j){
        for(int i=0; i<j; ++i){
            matriz_superior[i][j] =  matriz[i][j] - obtener_suma_superior(i,j);
            matriz_inferior[j][i] = (matriz[j][i] - obtener_suma_inferior(j,i)) / matriz_superior[i][i];
        }
        matriz_superior[j][j] = matriz[j][j] - obtener_suma_superior(j,j);
        if(matriz_inferior[j][j] == 0) throw cant_factor_exception();
    }
}

const char* cant_factor_exception::what() const throw(){
    return "zero found in diagonal";
}

