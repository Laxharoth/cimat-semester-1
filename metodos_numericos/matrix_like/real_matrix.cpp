#include "real_matrix.hpp"

namespace mymtx{

RealMatrix::RealMatrix(size_t shape_y,size_t shape_x):shape_x(shape_x),shape_y(shape_y){
    this->data = new double[this->shape_x * this->shape_y]{0};
}
RealMatrix::RealMatrix(const RealMatrix &other):shape_x(other.shape_x),shape_y(other.shape_y){
    this->data = new double[this->shape_x * this->shape_y];
    memcpy(this->data,other.data,this->shape_x*this->shape_y*sizeof(double));
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
RealMatrix RealMatrix::identity(const size_t n){
    RealMatrix identity(n,n);
    for( size_t i = 0; i < n; ++i )identity[i][i]=1;
    return identity;
}
RealMatrix RealMatrix::tridiag(const size_t n, double (*low)(const int i), double (*dig)(const int i), double (*up)(const int i)){
    RealMatrix trid(n,n);
    for( size_t i = 0; i < n; ++i ){
        trid[i][i] = dig(i);
        if(i>0)trid[i][i-1] = low(i);
        if(i<n-1)trid[i][i+1] = up(i);
    }
    return trid;
}
RealMatrix RealMatrix::tridiag(const size_t n, const double low,const double dig, const double up){
    RealMatrix trid(n,n);
    for( size_t i = 0; i < n; ++i ){
        trid[i][i] = dig;
        if(i>0)trid[i][i-1] = low;
        if(i<n-1)trid[i][i+1] = up;
    }
    return trid;
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
RealMatrix &RealMatrix::operator*=(const double coef){
    for (size_t i = 0; i < this->shape_y; i++){
        (*this)[i]*=coef;
    }
    return *this;
}
RealMatrix RealMatrix::operator*(const double coef) const{
    RealMatrix cpy = *this;
    for (size_t i = 0; i < this->shape_y; i++){
        cpy[i]*=coef;
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
RealMatrix &RealMatrix::operator-=(const RealMatrix &other){
    for(int i = 0; i < this->shape_y; i++){
        for(int j = 0; j < this->shape_x; j++){
            (*this)[i][j] -= other[i][j];
        }
    }
    return *this;
}
RealMatrix &RealMatrix::prod_as_band(const double coef, const size_t height, const size_t width){
    size_t begin,end;
    for (size_t i = 0; i < this->shape_y; i++){
        begin = i - height;
        end = 1 + i + width;
        if( i < height ) begin = 0;
        if( i > this->shape_y - height ) end = this->shape_x;
        for (size_t j = begin; j < end; j++){
            (*this)[i][j] *= coef;
        }
    }
    return *this;
}
RealVector RealMatrix::prod_as_band(RealVector &vec , const size_t height, const size_t width) const{
    RealVector result(vec.size);
    size_t begin,end;
    for (size_t i = 0; i < this->shape_y; i++){
        begin = i - height;
        end = 1 + i + width;
        if( i < height ) begin = 0;
        if( i > this->shape_y - height ) end = this->shape_x;
        result[i] = 0;
        for (size_t j = begin; j < end; j++){
            result[i] += (*this)[i][j] * vec[j];
        }
    }
    return result;
}
void RealMatrix::fwrite( const char* filename, const RealMatrix& matrix){
    std::ofstream f(filename, std::ofstream::out | std::ofstream::binary);
    f.write( (char*)(&(matrix.shape_y)), sizeof(size_t) );
    f.write( (char*)(&(matrix.shape_x)), sizeof(size_t) );
    f.write( (char*)(matrix.data), sizeof(double)*matrix.shape_y*matrix.shape_y );
    f.close();
}
void RealVector::fwrite( const char* filename, const RealVector& vector){
    std::ofstream f(filename, std::ofstream::out | std::ofstream::binary);
    f.write( (char*)(&(vector.size)), sizeof(size_t) );
    f.write( (char*)(vector.data), sizeof(double)*vector.size );
    f.close();
}
RealMatrix RealMatrix::fread( const char* filename){
    size_t shape_x, shape_y;
    std::ifstream f(filename, std::ifstream::in | std::ifstream::binary);
    f.read( (char*)(&(shape_y)), sizeof(size_t) );
    f.read( (char*)(&(shape_x)), sizeof(size_t) );
    RealMatrix matrix(shape_y, shape_x);
    f.read( (char*)(matrix.data), sizeof(double)*shape_y*shape_y );
    f.close();
    return matrix;
}
RealVector RealVector::fread( const char* filename){
    size_t size;
    std::ifstream f(filename, std::ifstream::in | std::ifstream::binary);
    f.read( (char*)(&(size)), sizeof(size_t) );
    RealVector vector(size);
    f.read( (char*)(vector.data), sizeof(double)*size );
    f.close();
    return vector;
}
void abs(RealMatrix &A){
    for(size_t k=0; k<A.shape_y; ++k){
        for(auto i = A[k].begin(); i<A[k].end(); ++i){
            *i = std::abs(*i);
        }
    }
}
void abs(RealVector &V){
    for(auto i = V.begin(); i<V.end(); ++i){
        *i = std::abs(*i);
    }
}
}
#include "real_vector.cpp"
#include "real_vector_iterator.cpp"