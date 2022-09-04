#include "real_matrix.hpp"

namespace mymtx{

RealMatrix::RealMatrix(size_t shape_y,size_t shape_x):shape_x(shape_x),shape_y(shape_y){
    this->data = new double[this->shape_x * this->shape_y];
}
RealMatrix::RealMatrix(const RealMatrix &other):shape_x(other.shape_x),shape_y(other.shape_y){
    this->data = new double[this->shape_x * this->shape_y];
    memcpy(this->data,other.data,this->shape_x*this->shape_y*sizeof(double));
}
RealMatrix::RealMatrix(matrix_like<double> &other):shape_x(other.get_shape_x()),shape_y(other.get_shape_y()){
    this->data = new double[this->shape_x * this->shape_y];
    for(int i = 0; i < this->shape_x; ++i)
        for(int j = 0; j < this->shape_y; ++j)
            (*this)[i][j] = other[i][j];
}
RealMatrix::RealMatrix(std::initializer_list<std::initializer_list<double>> initial):
    shape_y(initial.size()),shape_x((initial.begin())->size()){
    this->data = new double[this->shape_x * this->shape_y];
    size_t i{0},j{0};
    for(auto ii = initial.begin(); ii != initial.end(); ++i,++ii, j=0){
        for(auto ij = (*ii).begin(); ij != (*ii).end(); ++j,++ij)
            (*this)[i][j] = *ij;
    }
}
RealMatrix::~RealMatrix(){
    delete[] data;
}
RealMatrix RealMatrix::traspose(const RealMatrix &m){
    RealMatrix traposed( m.shape_x, m.shape_y );
    for( size_t i = 0; i < m.shape_y; ++i )
        for( size_t j = 0; j < m.shape_x; ++j )
            traposed[j][i] = m[i][j];
    return traposed;
}
RealVector RealMatrix::operator[](const size_t &row){ return RealVector(data + row * shape_x, shape_x); }
const RealVector RealMatrix::operator[](const size_t &row) const { return RealVector(data + row * shape_x, shape_x); }
vector_iterator RealMatrix::begin(const size_t row){
    return vector_iterator(data + row * shape_x, row, 0);
}
vector_iterator RealMatrix::end(const size_t row){
    return vector_iterator(data + row * shape_x, row, shape_x);
}
const_vector_iterator RealMatrix::begin(const size_t row) const{
    return const_vector_iterator(data + row * shape_x, row, 0);
}
const_vector_iterator RealMatrix::end(const size_t row) const{
    return const_vector_iterator(data + row * shape_x, row, shape_x);
}
RealVector &RealMatrix::operator*=(RealVector &vec) const{
    auto cpy = (*this)*vec;
    vec = cpy;
    return vec;
}
RealVector RealMatrix::operator*(const RealVector &vec) const{
    RealVector result(vec.size);
    for(size_t i = 0; i < shape_y; ++i){
        result[i] = this->operator[](i) * vec;
    }
    return result;
}
RealMatrix &RealMatrix::operator*=(const RealMatrix &other){
    auto cpy = (*this)*other;
    (*this) = cpy;
    return *this;
}
RealMatrix &RealMatrix::operator*=(double coef){
    for (size_t i = 0; i < this->shape_y; i++){
        (*this)[i]*=coef;
    }
    return *this;
}
RealMatrix RealMatrix::operator*(const RealMatrix &other) const{
    auto t = traspose(other);
    auto cpy = RealMatrix(this->shape_y,this->shape_x);
    for(size_t i = 0; i < this->shape_y; ++i){
        for(size_t j = 0; j < this->shape_x; ++j){
            cpy[i][j] = this->operator[](i) * t[j];
        }
    }
    return cpy;
}
RealMatrix &RealMatrix::operator=(const RealMatrix &other){
    memcpy(this->data, other.data, this->shape_x * this->shape_y *sizeof(double));
    return *this;
}

}
#include "real_vector.cpp"
#include "real_vector_iterator.cpp"