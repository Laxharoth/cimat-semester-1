#ifndef MATRIX_LIKE_HPP
#define MATRIX_LIKE_HPP

#include "matrix_exception.hpp"
#include <initializer_list>
#include <cstdlib>
#include <string.h>

namespace mymtx{

template <class T>
class matrix_like;
template <class T>
class matrix;
template <class T>
class array_like;
}
#include "matrix_iterator.tcc"
namespace mymtx{
template <class T>
class vector;
template <class T>
class marray;

template <class T>
class matrix_like{
    matrix_like(){};
protected:
    size_t shape_y;
    size_t shape_x;
public:
    virtual array_like<T> &operator[](const size_t &row) = 0;
    size_t get_shape_y() const { return shape_y; };
    size_t get_shape_x() const { return shape_x; };
    vector<T> operator*(array_like<T>& vec){
        if( vec.get_size() != this->shape_y )
            throw wrong_shape_exception(this->get_shape_x(), this->get_shape_y(), vec.get_size(), 1);
        vector<T> result(vec.get_size());
        for(size_t i = 0; i < this->get_shape_y(); ++i) {
            result[i] = 0;
            for(size_t j = 0; j < this->get_shape_y(); ++j) {
                T a = (*this)[j][i];
                T b = vec[j];
                result[i] += a * b;
            }
        }
        return result;
    };
    matrix_like(size_t shape_y, size_t shape_x):shape_y(shape_y),shape_x(shape_x) {};
};

template <class T>
class matrix : public matrix_like<T>{
    T *data;
    marray<T> *array;
    public:
    matrix(const size_t &rows, const size_t &cols):matrix_like<T>(rows,cols){
        data = new T[rows*cols];
        array = new marray<T>(this->data, this->shape_y, this->shape_x);
    };
    matrix(std::initializer_list<std::initializer_list<T>> l)
        : matrix(l.size(), 0) 
    {
        for(auto i = l.begin(); i != l.end(); ++i){
            if( (*i).size() > this->shape_x ){
                this->shape_x = (*i).size();
            }
        }
        data = new T[this->shape_y*this->shape_x];
        this->array = new marray<T>(this->data, this->shape_y, this->shape_x);
        size_t i{0}, j{0};
        static size_t asdfds = 0;
        for(auto it_i = l.begin(); i < this->shape_y; ++i, ++it_i){
            for(auto it_j = (*it_i).begin(); j < this->shape_x; ++j, ++it_j){
                (*this)[i][j] = *it_j;
            }
            j = 0;
        }
    }
    matrix(const matrix<T> &other):matrix_like<T>(other.shape_y,other.shape_x){
        data = new T[this->shape_y*this->shape_x];
        array = new marray<T>(this->data, this->shape_y, this->shape_x);
        memcpy(this->data, other.data, sizeof(T) * this->shape_y * this->shape_x);
    }
    ~matrix(){ delete array; delete[] data; }
    marray<T> &operator[](const size_t &row) { array->row = row; return *array; }
};

template <class T>
class array_like{
    array_like(){};
protected:
    size_t size;
    size_t get_begin_n() const {return 0;};
    size_t get_end_n() const {return size;};
    
public:
    array_like(size_t size):size(size){};
    size_t get_size() const { return size; };
    virtual T &operator[](const size_t &row) = 0;
    array_like_iterator<T> begin(){
        return array_like_iterator<T>(get_begin_n(),this);
    }
    array_like_iterator<T> rbegin(){
        return array_like_iterator<T>(get_rbegin_n(),this);
    }
    array_like_iterator<T> end(){
        return array_like_iterator<T>(get_end_n(),this);
    }
    array_like_iterator<T> rend(){
        return array_like_iterator<T>(get_rend_n(),this);
    }
protected:
    virtual size_t get_row() const {return 0;}; 
    virtual size_t get_rbegin_n() const {return 0;}; 
    virtual size_t get_rend_n() const {return size;}; 
};

template <class T>
class vector : public array_like<T>{
    T *data;
    public:
    vector(size_t size):array_like<T>(size) { this->size = size ;data = new T[this->size]; }
    vector(std::initializer_list<T> initial):array_like<T>(initial.size()) { 
        data = new T[this->size];
        auto it = initial.begin();
        for(size_t i = 0; i < this->size; ++i) {
            data[i] = *it++;
        }
    }
    vector(vector<T> &other):array_like<T>(other.get_size()){
        data = new T[this->size];
        for(size_t i = 0; i < this->size; ++i) {
            data[i] = other[i];
        }
    }
    ~vector() { delete[] data; }
    T &operator[](const size_t &row) { return data[row]; }
};

template <class T>
class marray: public array_like<T>{
    public:
    size_t &rows,&cols;
    size_t row;
    T *data;
    marray(T *data, size_t &rows, size_t &cols)
        : array_like<T>(cols),data(data), rows(rows), cols(cols){};
    T &operator[](const size_t &col){
        return data[ row * cols + col];
    }
    size_t get_row() const {return row;}
};

}
#endif /* MATRIX_LIKE_HPP */
