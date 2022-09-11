#ifndef REAL_MATRIX_HPP
#define REAL_MATRIX_HPP

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <fstream>
#include <initializer_list>

namespace mymtx{
    class RealMatrix;
    class RealVector;
    class vector_iterator;
    class const_vector_iterator;

    class RealMatrix{
        RealMatrix();
        double *data;
        public:
        const size_t shape_x,shape_y;
        RealMatrix(size_t shape_y,size_t shape_x);
        RealMatrix(const RealMatrix &other);
        RealMatrix(std::initializer_list<std::initializer_list<double>> initial);
        ~RealMatrix();
        RealVector operator[](const size_t &row);
        const RealVector operator[](const size_t &row) const;
        vector_iterator begin(const size_t row);
        vector_iterator end(const size_t row);
        const_vector_iterator begin(const size_t row) const;
        const_vector_iterator end(const size_t row) const;
        RealMatrix &operator*=(const double coef);
        RealMatrix operator*(const double coef) const;
        RealVector &operator*=(RealVector &vec) const;
        RealVector operator*(const RealVector &vec) const;
        RealMatrix &operator*=(const RealMatrix &other);
        RealMatrix operator*(const RealMatrix &other) const;
        RealMatrix &operator=(const RealMatrix &other);
        RealMatrix &operator-=(const RealMatrix &other);
        static RealMatrix traspose(const RealMatrix &m);
        static RealMatrix identity(const size_t n);
        static RealMatrix tridiag(const size_t n, double (*low)(const int i), double (*dig)(const int i), double (*up)(const int i));
        static RealMatrix tridiag(const size_t n, const double low, const double dig, const double up);
        RealMatrix &prod_as_band(const double coef, const size_t height, const size_t width);
        RealVector prod_as_band(RealVector &vec , const size_t height, const size_t width) const;
        RealMatrix prod_as_band(const RealMatrix &other, const size_t height, const size_t width);
        static void fwrite( const char* filename, const RealMatrix& matrix);
        static RealMatrix fread( const char* filename);
    };
    class RealVector{
        RealVector();
        double *data;
        const bool allocated;
        public:
        const size_t size;
        RealVector(const size_t size);
        RealVector(const RealVector &other);
        RealVector(std::initializer_list<double> initial);
        RealVector(double *data, size_t size);
        ~RealVector();
        double &operator[](const size_t col);
        const double &operator[](const size_t col) const;
        vector_iterator begin();
        vector_iterator end();
        const_vector_iterator begin() const;
        const_vector_iterator end() const;
        RealVector operator*(const double coef);
        RealVector &operator*=(const double coef);
        RealVector operator/(const double coef);
        RealVector &operator/=(const double coef);
        double operator*(const RealVector &other) const;
        RealVector operator+(const RealVector &other);
        RealVector &operator+=(const RealVector &other);
        RealVector operator-(const RealVector &other);
        RealVector &operator-=(const RealVector &other);
        RealVector &operator=(const RealVector &other);
        RealMatrix cross_product(const RealVector &other) const;
        static RealVector normal(const size_t size);
        static void sort(RealVector &v);
        static void fwrite( const char* filename, const RealVector& vector);
        static RealVector fread( const char* filename);
    };
    class vector_iterator{
        vector_iterator(){}
protected:
    double *data;
    size_t row,col;
public:
    vector_iterator(double *data, const size_t row, const size_t col);
    vector_iterator(const vector_iterator &other);
    size_t get_row() const;
    size_t get_col() const;
    vector_iterator& operator++();
    vector_iterator operator++(int);
    vector_iterator& operator--();
    vector_iterator operator--(int);
    vector_iterator& operator+=(int c);
    vector_iterator operator+(int c) const;
    vector_iterator& operator-=(int c);
    vector_iterator operator-(int c) const ;
    bool operator==(const vector_iterator& rhs) const;
    bool operator!=(const vector_iterator& rhs) const;
    bool operator>(const vector_iterator& rhs)  const;
    bool operator<(const vector_iterator& rhs)  const;
    bool operator>=(const vector_iterator& rhs) const;
    bool operator<=(const vector_iterator& rhs) const;
    double& operator[](int);
    double& operator*();
    };
class const_vector_iterator{
        const_vector_iterator(){}
protected:
    const double *data;
    size_t row,col;
public:
    const_vector_iterator(const double *data, const size_t row, const size_t col);
    const_vector_iterator(const const_vector_iterator &other);
    size_t get_row() const;
    size_t get_col() const;
    const_vector_iterator& operator++();
    const_vector_iterator operator++(int);
    const_vector_iterator& operator--();
    const_vector_iterator operator--(int);
    const_vector_iterator& operator+=(int c);
    const_vector_iterator operator+(int c) const;
    const_vector_iterator& operator-=(int c);
    const_vector_iterator operator-(int c) const ;
    bool operator==(const const_vector_iterator& rhs) const;
    bool operator!=(const const_vector_iterator& rhs) const;
    bool operator>(const const_vector_iterator& rhs)  const;
    bool operator<(const const_vector_iterator& rhs)  const;
    bool operator>=(const const_vector_iterator& rhs) const;
    bool operator<=(const const_vector_iterator& rhs) const;
    const double& operator[](int) const;
    const double& operator*() const;
    };
void abs(RealMatrix &A);
}
#endif /* REAL_MATRIX_HPP */
