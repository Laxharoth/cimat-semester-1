#include "matriz_banda.hpp"

namespace mymtx{

double MatrizBanda::row_wrapper::default_value = 0;
double MatrizDiagonal::row_wrapper::default_value = 0;
MatrizBanda::MatrizBanda(const size_t &left, const size_t &right, const size_t &size):
    matrix_like<double>(size,size),
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
    :array_like<double>(size),matriz(matriz),right(right),left(left),row(row){};

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
size_t MatrizBanda::row_wrapper::get_rbegin_n() const{
    if( row <= left ) return 0;
    return  row - left;
};
size_t MatrizBanda::row_wrapper::get_rend_n() const{
    if(row > right) return get_size();
    return get_size() - right + row;
};

double &MatrizBanda::row_wrapper::operator[](const size_t &col){
    return get(col);
}
MatrizDiagonal::MatrizDiagonal(size_t size):matrix_like<double>(size,size),size(size){
    this->shape_x = this->size;
    this->shape_y = this->size;
    this->data = new double[this->size];
    this->wrapper = new MatrizDiagonal::row_wrapper(this->data, size,0);
}
MatrizDiagonal::MatrizDiagonal(std::initializer_list<double> init):matrix_like<double>(init.size(),init.size()),size(init.size()){
    this->shape_x = this->size;
    this->shape_y = this->size;
    this->data = new double[this->size];
    this->wrapper = new MatrizDiagonal::row_wrapper(this->data, size,0);
}
MatrizDiagonal::row_wrapper & MatrizDiagonal::operator[](const size_t &col){
    this->wrapper->row = col;
    return *(this->wrapper);
}
MatrizDiagonal::~MatrizDiagonal(){
    delete this->wrapper;
    delete[] this->data;
}
MatrizDiagonal::row_wrapper::row_wrapper(double *data, size_t size, size_t row):array_like<double>(size),data(data),row(row){}
double &MatrizDiagonal::row_wrapper::operator[](const size_t &col){
    if( row == row ) return this->data[row];
    default_value = 0;
    return default_value;
}
size_t MatrizDiagonal::row_wrapper::get_row() const {return row;}
size_t MatrizBanda::row_wrapper::get_row() const {return row;}
size_t MatrizDiagonal::row_wrapper::get_rbegin_n() const{
    return row;
};
size_t MatrizDiagonal::row_wrapper::get_rend_n() const{
    return row + 1;
};
MatrizBanda::row_wrapper *MatrizBanda::row_wrapper::allocate_this_cpy(){
    return new row_wrapper(this->matriz, this->left, this->right, this->size, this->row);
}
MatrizDiagonal::row_wrapper *MatrizDiagonal::row_wrapper::allocate_this_cpy(){
    return new row_wrapper(this->data, this->size, this->row);
}
}