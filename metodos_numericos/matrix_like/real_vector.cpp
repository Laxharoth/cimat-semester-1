#include "real_matrix.hpp"
namespace mymtx{
RealVector::RealVector(const size_t size):size(size),allocated(true){
    this->data = new double[this->size];
}
RealVector::RealVector(const RealVector &other):size(other.size),allocated(true){
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
const_vector_iterator RealVector::begin() const{
    return const_vector_iterator(data,0,0);
}
const_vector_iterator RealVector::end() const{
        return const_vector_iterator(data,0,size);
}
RealVector RealVector::operator*(const double coef){
    RealVector cpy = *this;
    cpy*=coef;
    return cpy;
}
RealVector &RealVector::operator*=(const double coef){
    for( auto i = this->begin(); i != this->end(); ++i )
        (*i)*=coef;
    return *this;
}
RealVector RealVector::operator/(const double coef){
    RealVector cpy = *this;
    cpy/=coef;
    return cpy;
}
RealVector &RealVector::operator/=(const double coef){
    for( auto i = this->begin(); i != this->end(); ++i )
        (*i)/=coef;
    return *this;
}
double RealVector::operator*(const RealVector &other) const {
    double sum{0};
    auto it = other.begin();
    auto ti = this->begin();
    while(it != other.end()){
        sum += (*it)*(*ti);
        ++it;++ti;
    }
    return sum;
}
RealVector &RealVector::operator+=(const RealVector &other){
    auto it = other.begin();
    auto ti = this->begin();
    while(it != other.end()){
        (*ti) += (*it);
        ++it;++ti;
    }
    return *this;
}
RealVector RealVector::operator+(const RealVector &other){
    auto cpy = *this;
    cpy+=other;
    return cpy;
}
RealVector &RealVector::operator-=(const RealVector &other){
    auto it = other.begin();
    auto ti = this->begin();
    while(it != other.end()){
        (*ti) -= (*it);
        ++it;++ti;
    }
    return *this;
}
RealVector RealVector::operator-(const RealVector &other){
    auto cpy = *this;
    cpy-=other;
    return cpy;
}
RealVector &RealVector::operator=(const RealVector &other){
    memcpy(this->data,other.data,sizeof(double)*this->size);
    return *this;
}
RealMatrix RealVector::cross_product(const RealVector &other) const {
    RealMatrix result(other.size,other.size);
    for (size_t i = 0; i < other.size; i++){
        for (size_t j = 0; j < other.size; j++){
            result[i][j] = other[i]*other[j];
        }
    }
    return result;
}
}