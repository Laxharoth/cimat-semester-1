#include "real_matrix.hpp"
namespace mymtx{
RealVector::RealVector(const size_t size):size(size),allocated(true){
    this->data = new double[this->size];
}
RealVector::RealVector(const RealVector &other):size(other.size),allocated(true){
    this->data = new double[this->size];
    memcpy(this->data,other.data,this->size*sizeof(double));
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
RealVector RealVector::normal(const size_t size){
    RealVector norm(size);
    const double val = 1 / std::sqrt(size);
    for( auto i = norm.begin(); i != norm.end(); ++i )
        *i = val;
    return norm;
}
void RealVector::sort(RealVector &v){
    std::sort(v.data, v.data+v.size);
}
double RealVector::distance() const{
    double sum = 0;
    for(auto i=this->begin(); i!=this->end(); ++i) sum += (*i) * (*i);
    return std::sqrt(sum);
}
RealMatrix RealVector::as_matrix() const {
    RealMatrix m(this->size,1);
    memcmp(m.data, this->data, sizeof(double)*this->size);
    return m;
}
RealMatrix::Column::Column(const RealMatrix::Column &other):size(other.size),increment(other.increment),data(other.data){}
RealMatrix::Column::Column(double *data, size_t size,const size_t increment):
    data(data),size(size),increment(increment){}
double &RealMatrix::Column::operator[](const size_t row){
    return *(data+row*increment);
}
const double &RealMatrix::Column::operator[](const size_t row) const{
    return *(data+row*increment);
}
double RealMatrix::Column::distance() const{
    double v;
    for (size_t i = 0; i < size; i++)
        v = (*this)[i] * (*this)[i];
    return v;
}
RealMatrix RealMatrix::Column::as_matrix() const{
    RealMatrix result(size,1);
    for (size_t i = 0; i < size; i++)
        result[i][0] = (*this)[i];
    return result;
}
RealVector RealMatrix::Column::as_vector() const{
    RealVector v(size);
    for(int i=0;i<this->size;i++) v[i] = (*this)[i];
    return v;
}
double RealMatrix::Column::operator*(const RealMatrix::Column&other) const{
    double r=0;
    for (size_t i = 0; i < size; i++){
        r+= (*this)[i]*other[i];
    }
    return r;
}
RealMatrix::Column &RealMatrix::Column::operator=(const RealVector &other){
    for (size_t i = 0; i < size; i++){
        (*this)[i] = other[i];
    }
    return *this;
}
RealMatrix::Column &RealMatrix::Column::operator-=(const RealVector &other){
    for (size_t i = 0; i < size; i++){
        (*this)[i] -= other[i];
    }
    return *this;
}
RealVector RealMatrix::Column::operator-(const RealMatrix::Column &other)const{
    RealVector v(this->size);
    for (size_t i = 0; i < this->size; i++){
        v[i] = (*this)[i]-other[i];
    }
    return v;
}
}
mymtx::RealVector operator*(const mymtx::RealVector &v, const double c){
    mymtx::RealVector cpy = v;
    cpy*=c;
    return cpy;
}
mymtx::RealVector operator*(const double c, const mymtx::RealVector &v){
    return v*c;
}
double operator*(const mymtx::RealMatrix::Column &c, const mymtx::RealVector &v){
    double r=0;
    for( auto i = 0; i < v.size; ++i )
        r += c[i] * v[i];
    return r;
}
double operator*(const mymtx::RealVector &v, const mymtx::RealMatrix::Column &c){
    return c * v;
}
mymtx::RealVector operator*(const mymtx::RealMatrix::Column &c, double &coef){
    mymtx::RealVector v(c.size);
    for (size_t i = 0; i < c.size; i++){
        v[i] = c[i] * coef;
    }
    return v;
}
mymtx::RealVector operator*(double &coef, const mymtx::RealMatrix::Column &c){
    return c * coef;
}
mymtx::RealVector operator+(const mymtx::RealMatrix::Column &c, const mymtx::RealVector &v){
    mymtx::RealVector vr(c.size);
    for (size_t i = 0; i < c.size; i++)
        vr[i] = c[i] + v[i];
    return vr;
}
mymtx::RealVector operator+(const mymtx::RealVector &v, const mymtx::RealMatrix::Column &c){
    return c+v;
}
mymtx::RealVector operator-(const mymtx::RealMatrix::Column &c, const mymtx::RealVector &v){
    mymtx::RealVector vr(c.size);
    for (size_t i = 0; i < c.size; i++)
        vr[i] = c[i] - v[i];
    return vr;
}
mymtx::RealVector operator-(const mymtx::RealVector &v, const mymtx::RealMatrix::Column &c){
    mymtx::RealVector vr(c.size);
    for (size_t i = 0; i < c.size; i++)
        vr[i] = v[i] + c[i];
    return vr;
}
