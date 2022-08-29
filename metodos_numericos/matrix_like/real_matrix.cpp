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
RealVector RealMatrix::operator[](const size_t &row){ return RealVector(data + row * shape_x, row); }
const RealVector RealMatrix::operator[](const size_t &row) const { return RealVector(data + row * shape_x, row); }
vector_iterator RealMatrix::begin(const size_t row){
    return vector_iterator(data + row * shape_x, row, 0);
}
vector_iterator RealMatrix::end(const size_t row){
    return vector_iterator(data + row * shape_x, row, shape_x);
}
RealVector::RealVector(const size_t size):size(size),allocated(true){
    this->data = new double[this->size];
}
RealVector::RealVector(const RealVector &other):size(other.size),allocated(true){
    this->data = new double[this->size];
    for(size_t i = 0; i < this->size; ++i){
        this->data[i] = other[i];
    }
}
RealVector::RealVector(array_like<double> &other):size(other.get_size()),allocated(true){
    this->data = new double[this->size];
    for(size_t i = 0; i < this->size; ++i){
        this->data[i] = other[i];
    }
}
RealVector::RealVector(std::initializer_list<double> initial):size(initial.size()),allocated(true){
    this->data = new double[this->size];
    size_t i{0};
    for(auto j= initial.begin(); j!= initial.end(); ++j,++i)
        data[i] = *j;
}
RealVector::RealVector(double *data, size_t size):size(size),allocated(false),data(data){}
RealVector::~RealVector(){
    if(allocated) delete[] data;
}
double &RealVector::operator[](const size_t col){return data[col];}
const double &RealVector::operator[](const size_t col) const {return data[col];}
vector_iterator RealVector::begin(){
    return vector_iterator(data,0,0);
}
vector_iterator RealVector::end(){
        return vector_iterator(data,0,size);
}

vector_iterator::vector_iterator(double *data, const size_t row, const size_t col):data(data),row(row),col(col){}
vector_iterator::vector_iterator(const vector_iterator &other):data(other.data),row(other.row),col(other.col){}
size_t vector_iterator::get_row() const { return row; }
size_t vector_iterator::get_col() const { return col; }
vector_iterator& vector_iterator::operator++()     {++col;return *this;}
vector_iterator vector_iterator::operator++(int)   {vector_iterator tmp(*this); operator++(); return tmp;}
vector_iterator& vector_iterator::operator--()     {--col;return *this;}
vector_iterator vector_iterator::operator--(int)   {vector_iterator tmp(*this); operator--(); return tmp;}
vector_iterator& vector_iterator::operator+=(int c){col+=c;return *this;}
vector_iterator vector_iterator::operator+(int c) const {vector_iterator tmp(*this); tmp+=c; return tmp;}
vector_iterator& vector_iterator::operator-=(int c){col-=c;return *this; }
vector_iterator vector_iterator::operator-(int c) const {vector_iterator tmp(*this); tmp-=c; return tmp;}
bool vector_iterator::operator==(const vector_iterator& rhs) const { return col==rhs.col; }
bool vector_iterator::operator!=(const vector_iterator& rhs) const { return col!=rhs.col; }
bool vector_iterator::operator>(const vector_iterator& rhs ) const { return col>rhs.col;  }
bool vector_iterator::operator<(const vector_iterator& rhs ) const { return col<rhs.col;  }
bool vector_iterator::operator>=(const vector_iterator& rhs) const { return col>=rhs.col; }
bool vector_iterator::operator<=(const vector_iterator& rhs) const { return col<=rhs.col; }
double& vector_iterator::operator*() {return data[col];}
}