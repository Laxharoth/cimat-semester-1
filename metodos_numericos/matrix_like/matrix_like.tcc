#ifndef MATRIX_LIKE_HPP
#define MATRIX_LIKE_HPP

#include "matrix_exception.hpp"
#include <initializer_list>
#include <cstdlib>
#include <string.h>

template <class T>
class array_like{
protected:
    size_t size;
public:
    virtual T &operator[](const size_t &row) = 0;
    size_t get_size() const { return size; };
};
template <class T>
class vector : public array_like<T>{
    T *data;
    public:
    vector(size_t size) { this->size = size ;data = new T[this->size]; }
    vector(std::initializer_list<T> initial) { 
        this->size = initial.size();
        data = new T[this->size];
        auto it = initial.begin();
        for(size_t i = 0; i < this->size; ++i) {
            data[i] = *it++;
        }
    }
    vector(vector<T> &other){
        this->size = other.get_size();
        data = new T[this->size];
        for(size_t i = 0; i < this->size; ++i) {
            data[i] = other[i];
        }
    }
    ~vector() { delete[] data; }
    T &operator[](const size_t &row) { return data[row]; }
};
template <class T>
class matrix_like{
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
};
template <class T>
class marray: public array_like<T>{
    public:
    size_t &rows,&cols;
    size_t row;
    T *data;
    marray(T *data, size_t &rows, size_t &cols)
        : data(data), rows(rows), cols(cols){};
    T &operator[](const size_t &col){
        return data[ row * cols + col];
    }
};
template <class T>
class matrix : public matrix_like<T>{
    T *data;
    marray<T> *array;
    public:
    matrix(const size_t &rows, const size_t &cols){
        this->shape_y = rows;
        this->shape_x = cols;
        data = new T[rows*cols];
        array = new marray<T>(this->data, this->shape_y, this->shape_x);
    };
    matrix(std::initializer_list<std::initializer_list<T>> l){
        this->shape_y = l.size();
        this->shape_x = 0;
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
    matrix(const matrix<T> &other){
        this->shape_y = other.get_shape_y();
        this->shape_x = other.get_shape_x();
        data = new T[this->shape_y*this->shape_x];
        array = new marray<T>(this->data, this->shape_y, this->shape_x);
        memcpy(this->data, other.data, sizeof(T) * this->shape_y * this->shape_x);
    }
    ~matrix(){ delete array; delete[] data; }
    marray<T> &operator[](const size_t &row) { array->row = row; return *array; }
};

#endif /* MATRIX_LIKE_HPP */
