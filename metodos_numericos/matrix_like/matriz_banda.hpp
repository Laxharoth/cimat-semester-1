#ifndef MATRIZ_BANDA_HPP
#define MATRIZ_BANDA_HPP

#include "matrix_like.tcc"
#include <initializer_list>
#include <cstdlib>

namespace mymtx{

class MatrizBanda: public matrix_like<double>{
    class row_wrapper;
    matrix<double> *matriz;
    row_wrapper *wrapper;
    size_t left{}, right{},size{};
class row_wrapper: public array_like<double>{
    static double default_value;
    matrix<double> &matriz;
    size_t &left, &right;
    size_t row;
public:
    row_wrapper(matrix<double> &matriz, size_t &left, size_t &right, size_t &size, size_t row);
    double &get(const size_t &col);
    double & operator[](const size_t &col);
    friend class MatrizBanda;
    size_t get_row() const;
    size_t get_rbegin_n() const;
    size_t get_rend_n() const;
};
public:
    row_wrapper &get(const size_t &row);
    row_wrapper &operator[](const size_t &row);
    MatrizBanda(const size_t &left, const size_t &right, const size_t &size);
    ~MatrizBanda();
};
class MatrizDiagonal : public matrix_like<double> {
    class row_wrapper;
    size_t size;
    double *data;
    row_wrapper *wrapper;
    class row_wrapper: public array_like<double>{
    static double default_value;
    double *data;
    size_t row;
    public:
        row_wrapper(double *data, size_t size);
        double & operator[](const size_t &col);
        friend class MatrizDiagonal;
        size_t get_row() const;
        size_t get_rbegin_n() const;
        size_t get_rend_n() const;
    };
    MatrizDiagonal(size_t size);
    MatrizDiagonal(std::initializer_list<double> init);
    ~MatrizDiagonal();
    row_wrapper & operator[](const size_t &col);

};
}
#endif /* MATRIZ_BANDA_HPP */
