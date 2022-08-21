#include "matriz_banda.hpp"

MatrizBanda::MatrizBanda(const size_t &left, const size_t &right, const size_t &size):
    left(left),right(right),size(size){ 
    this->shape_x  = this->shape_y = size;
    matriz = new matrix<double>(size, left + 1 + right);
    this->wrapper = new row_wrapper(*(this->matriz), this->left, this->right, this->size, 0);
}
MatrizBanda::~MatrizBanda(){
    delete this->wrapper;
    delete matriz;
}

MatrizBanda::row_wrapper::row_wrapper(matrix<double> &matriz, size_t &left, size_t &right, size_t &size, size_t row)
    :matriz(matriz),right(right),left(left),size(size),row(row){ this->size = size; };

double MatrizBanda::row_wrapper::default_value = 0;
MatrizBanda::row_wrapper &MatrizBanda::get(const size_t &row){
    this->wrapper->row = row;
    return *(this->wrapper);
}
MatrizBanda::row_wrapper &MatrizBanda::operator[](const size_t &row){
    return this->get(row);
}
double &MatrizBanda::row_wrapper::get(const size_t &col){
    if( (row > col && row - col > left ) || ( col > row && col - row > right ) ){ 
        default_value = 0;
        return default_value;
    }
    if(row == col){ return matriz[row][left]; }
    if(row < col){return matriz[col][left + col - row ]; }
    return matriz[col][ left - row  + col ];
}
double &MatrizBanda::row_wrapper::operator[](const size_t &col){
    return get(col);
}