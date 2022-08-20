#ifndef MATRIX_LIKE_HPP
#define MATRIX_LIKE_HPP
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
    ~matrix(){ delete array; delete[] data; }
    marray<T> &operator[](const int &row) { array->row = row; return *array; }
};

#endif /* MATRIX_LIKE_HPP */
