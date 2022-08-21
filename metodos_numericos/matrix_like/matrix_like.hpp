#ifndef MATRIX_LIKE_HPP
#define MATRIX_LIKE_HPP

#include <initializer_list>
#include <cstdlib>

template <class T>
class array_like{
protected:
    size_t size;
public:
    virtual T &operator[](const size_t &row) = 0;
};
template <class T>
class matrix_like{
protected:
    size_t shape_y;
    size_t shape_x;
public:
    virtual array_like<T> &operator[](const size_t &row) = 0;
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
    ~matrix(){ delete array; delete[] data; }
    marray<T> &operator[](const size_t &row) { array->row = row; return *array; }
};

#endif /* MATRIX_LIKE_HPP */
