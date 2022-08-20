#include "matriz_banda.hpp"

MatrizBanda::MatrizBanda(const int &left, const int &right, const int &size):
    left(left),right(right),size(size){
    matriz = new double *[size];
    for(int i=0; i<size; ++i){
        matriz[i] = new double[left + right +1];
    }
    this->wrapper = new row_wrapper(this->matriz, this->left, this->right, this->size, 0);
}
MatrizBanda::~MatrizBanda(){
    delete this->wrapper;
    for(int i=0; i<size; ++i){
        delete[] matriz[i];
    }
    delete[] matriz;
}

MatrizBanda::row_wrapper::row_wrapper(double** &matriz, int &left, int &right, int &size, int row)
    :matriz(matriz),right(right),left(left),size(size),row(row){};

double MatrizBanda::row_wrapper::default_value = 0;
MatrizBanda::row_wrapper &MatrizBanda::get(const int &row){
    this->wrapper->row = row;
    return *(this->wrapper);
}
MatrizBanda::row_wrapper &MatrizBanda::operator[](const int &row){
    return this->get(row);
}
double &MatrizBanda::row_wrapper::get(const int &col){
    if( (row > col && row - col > left ) || ( col > row && col - row > right ) ){ 
        default_value = 0;
        return default_value;
    }
    if(row == col){ return matriz[row][left]; }
    if(row < col){return matriz[col][left + col - row ]; }
    return matriz[col][ left - row  + col ];
}
double &MatrizBanda::row_wrapper::operator[](const int &col){
    return get(col);
}