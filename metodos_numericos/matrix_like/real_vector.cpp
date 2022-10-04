#include "real_matrix.hpp"
#include <cmath>
#include <cstddef>
#include <future>
namespace mymtx{
RealVector::RealVector(const size_t size):size(size),allocated(true){
    this->data = new double[this->size];
    memset(this->data, 0, this->size * sizeof(double));
}
RealVector::RealVector(const RealVector &other):size(other.size),allocated(true){
    this->data = new double[this->size];
    memcpy(this->data,other.data,this->size*sizeof(double));
}
RealVector::RealVector(RealVector &&other):size(other.size),allocated(other.allocated){
    this->data = other.data;
    other.data = nullptr;
}
RealVector::RealVector(std::initializer_list<double> initial):size(initial.size()),allocated(true){
    this->data = new double[this->size];
    size_t i{0};
    for(auto j= initial.begin(); j!= initial.end(); ++j,++i)
        data[i] = *j;
}
RealVector::RealVector(double *data, size_t size):size(size),allocated(false),data(data){}
RealVector::~RealVector(){
    if(allocated && data!=nullptr) delete[] data;
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
struct vector_coef_data{size_t begin; size_t end; RealVector *v;double coef;};
RealVector &RealVector::operator*=(const double coef){
    #ifndef NO_ASYNC
    if(size > MIN_OPER_FOR_THREAD){
    std::vector<std::future<void>> futures;
    for(size_t i= 0; i<3; ++i){
        size_t begin = (size *i) / 3;
        size_t end = (size *(i+1)) / 3;
        futures.push_back(std::async(std::launch::async,[](vector_coef_data d){
            for (size_t i = d.begin; i < d.end; i++){
                (*d.v)[i] *= d.coef;
            }
        }, vector_coef_data{begin,end,this,coef}));
    }
    futures.clear();
    }else
    #endif
    {
    for( auto i = this->begin(); i != this->end(); ++i )
        (*i)*=coef;
    }

    return *this;
}
RealVector RealVector::operator/(const double coef){
    RealVector cpy = *this;
    cpy/=coef;
    return cpy;
}
RealVector &RealVector::operator/=(const double coef){
    (*this) *= 1/coef;
    return *this;
}
struct vector_vector_data{size_t begin; size_t end; const RealVector* const v1; const RealVector* const v2; double * const sum;};
double RealVector::operator*(const RealVector &other) const {
    double sum{0};
    #ifndef NO_ASYNC
    double sums[3];
    if(size > MIN_OPER_FOR_THREAD){
    std::vector<std::future<void>> futures;
    for(size_t i= 0; i<3; ++i){
        size_t begin = (size*i) / 3;
        size_t end = (size *(i+1)) / 3;
        futures.push_back(std::async(std::launch::async,[](vector_vector_data d){
            (*d.sum) = 0;
            for (size_t i = d.begin; i < d.end; i++){
                (*d.sum) += (*d.v1)[i] * (*d.v2)[i];
            }
        }, vector_vector_data{begin,end,this,&other,sums+i}));
    }
    futures.clear();
    sum = sums[0]+sums[1]+sums[2];
    }
    else
    #endif
    {
    auto it = other.begin();
    auto ti = this->begin();
    while(it != other.end()){
        sum += (*it)*(*ti);
        ++it;++ti;
    }
    }
    return sum;
}
struct vector_sum_vector_data{size_t begin; size_t end; RealVector* const v1; const RealVector* const v2;};
RealVector &RealVector::operator+=(const RealVector &other){
    #ifndef NO_ASYNC
    if(size > MIN_OPER_FOR_THREAD){
    std::vector<std::future<void>> futures;
    for(size_t i= 0; i<3; ++i){
        size_t begin = (size*i) / 3;
        size_t end = (size *(i+1)) / 3;
        futures.push_back(std::async(std::launch::async,[](vector_sum_vector_data d){
            for (size_t i = d.begin; i < d.end; i++){
                (*d.v1)[i] += (*d.v2)[i];
            }
        }, vector_sum_vector_data{begin,end,this,&other}));
    }}
    else
    #endif 
    {
    auto it = other.begin();
    auto ti = this->begin();
    while(it != other.end()){
        (*ti) += (*it);
        ++it;++ti;
    }
    }
    return *this;
}
RealVector RealVector::operator+(const RealVector &other){
    auto cpy = *this;
    cpy+=other;
    return cpy;
}
RealVector &RealVector::operator-=(const RealVector &other){
    #ifndef NO_ASYNC
    if(size > MIN_OPER_FOR_THREAD){
    std::vector<std::future<void>> futures;
    for(size_t i= 0; i<3; ++i){
        size_t begin = (size*i) / 3;
        size_t end = (size *(i+1)) / 3;
        futures.push_back(std::async(std::launch::async,[](vector_sum_vector_data d){
            for (size_t i = d.begin; i < d.end; i++){
                (*d.v1)[i] -= (*d.v2)[i];
            }
        }, vector_sum_vector_data{begin,end,this,&other}));
    }}
    else
    #endif 
    {
    auto it = other.begin();
    auto ti = this->begin();
    while(it != other.end()){
        (*ti) -= (*it);
        ++it;++ti;
    }
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
RealVector &RealVector::operator=(const double coef){
    memset(this->data,coef,sizeof(double)*this->size);
    return *this;
}
struct vector_cross_vector_data{size_t begin; size_t end; const RealVector* const v1; const RealVector* const v2; RealMatrix *res;};
RealMatrix RealVector::cross_product(const RealVector &other) const {
    RealMatrix result(other.size,other.size);
    #ifndef NO_ASYNC
    if(size > MIN_OPER_FOR_THREAD){
    std::vector<std::future<void>> futures;
    for(size_t i= 0; i<3; ++i){
        size_t begin = (size*i) / 3;
        size_t end = (size *(i+1)) / 3;
        futures.push_back(std::async(std::launch::async,[](vector_cross_vector_data d){
            for (size_t i = d.begin; i < d.end; i++){
                for (size_t j = 0; j < (*d.v2).size; j++){
                    (*d.res)(i,j) = (*d.v1)[i]*(*d.v2)[j];
                }
            }
        }, vector_cross_vector_data{begin,end,this,&other,&result}));
    }}
    else
    #endif 
    {
    for (size_t i = 0; i < other.size; i++){
        for (size_t j = 0; j < other.size; j++){
            result[i][j] = other[i]*other[j];
        }
    }}
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
    double sum = (*this)*(*this);
    return std::sqrt(sum);
}
RealMatrix RealVector::as_matrix() const {
    RealMatrix m(this->size,1);
    memcpy(m.data, this->data, sizeof(double)*this->size);
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
    double v=0;
    for (size_t i = 0; i < size; i++)
        v += (*this)[i] * (*this)[i];
    return std::sqrt(v);
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
struct column_sum_vector_data{size_t begin; size_t end; RealMatrix::Column* const v1; const RealVector* const v2;};
RealMatrix::Column &RealMatrix::Column::operator-=(const RealVector &other){
    #ifndef NO_ASYNC
    if(size > MIN_OPER_FOR_THREAD){
    std::vector<std::future<void>> futures;
    for(size_t i= 0; i<3; ++i){
        size_t begin = (size*i) / 3;
        size_t end = (size *(i+1)) / 3;
        futures.push_back(std::async(std::launch::async,[](column_sum_vector_data d){
            for (size_t i = d.begin; i < d.end; i++){
                (*d.v1)[i] -= (*d.v2)[i];
            }
        }, column_sum_vector_data{begin,end,this,&other}));
    }}
    else
    #endif 
    {
    for (size_t i = 0; i < size; i++){
        (*this)[i] -= other[i];
    }
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
RealVector map(const RealVector &v,double (*callback)(const double)){
    RealVector v_new(v.size);
    auto i = v.begin();
    for(auto j =v_new.begin(); j !=v_new.end();++i,++j)
        *j = callback(*i);
    return v_new;
}
double reduce(const RealVector &v,double (*callback)(double acc, const double cur),const double start){
    double acc = start;
    for(auto i=v.begin();i!=v.end();++i)
        acc = callback(acc,*i);
    return acc;
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
