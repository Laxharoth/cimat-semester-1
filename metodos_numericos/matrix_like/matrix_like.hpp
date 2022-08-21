#ifndef MATRIX_LIKE_HPP
#define MATRIX_LIKE_HPP

#include <initializer_list>

template <class T>
class array_like{
public:
    virtual T &operator[](const int &row) = 0;
};
template <class T>
class matrix_like{
public:
    virtual array_like<T> &operator[](const int &row) = 0;
};
template <class T>
class marray: public array_like<T>{
    public:
    int &rows,&cols;
    int row;
    T *data;
    marray(T *data, int &rows, int &cols)
        : data(data), rows(rows), cols(cols){};
    T &operator[](const int &col){
        return data[ row * cols + col];
    }
};
template <class T>
class matrix : public matrix_like<T>{
    T *data;
    int rows,cols;
    marray<T> *array;
    public:
    matrix(const int &rows, const int &cols) : rows(rows), cols(cols){
        data = new T[rows*cols];
        array = new marray<T>(this->data, this->rows, this->cols);
    };
    matrix(std::initializer_list<std::initializer_list<T>> l){
        this->rows = l.size();
        this->cols = 0;
        for(auto i = l.begin(); i != l.end(); ++i){
            if( (*i).size() > this->cols ){
                this->cols = (*i).size();
            }
        }
        data = new T[rows*cols];
        this->array = new marray<T>(this->data, this->rows, this->cols);
        int i{0}, j{0};
        static int asdfds = 0;
        for(auto it_i = l.begin(); i < this->rows; ++i, ++it_i){
            for(auto it_j = (*it_i).begin(); j < this->cols; ++j, ++it_j){
                (*this)[i][j] = *it_j;
            }
            j = 0;
        }
    }
    ~matrix(){ delete array; delete[] data; }
    marray<T> &operator[](const int &row) { array->row = row; return *array; }
};

#endif /* MATRIX_LIKE_HPP */
