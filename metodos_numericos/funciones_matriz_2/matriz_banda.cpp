#include "matriz_banda.hpp"

MatrizBanda::MatrizBanda(const int &left, const int &right, const int &size){
    matriz = new double *[size];
    for(int i=0; i<size; ++i){
        matriz[i] = new double[left + right +1];
    }
}
MatrizBanda::~MatrizBanda(){
    for(int i=0; i<size; ++i){
        delete[] matriz[i];
    }
    delete[] matriz;
}

MatrizBanda::row_wrapper::row_wrapper(double** &matriz, int &left, int &right, int row)
    :matriz(matriz),right(right),left(left),row(row){};

double MatrizBanda::row_wrapper::default_value = 0;

double &MatrizBanda::get(const int &row, const int &col){
    return row_wrapper(matriz, left, right, row).get(col);
}
MatrizBanda::row_wrapper MatrizBanda::operator[](const int &row){
    return row_wrapper(matriz, left, right, row);
}

double &MatrizBanda::row_wrapper::get(const int &col){
    if(row > left - col || col > right- row){ 
        default_value = 0;
        return default_value;
    }
    if(row == col){ return matriz[row][left]; }
    if(row < col){return matriz[col + 1][ left + col - row ]; }
    return matriz[col][ left - row  + col ];
}
double &MatrizBanda::row_wrapper::operator[](const int &col){
    return get(col);
}