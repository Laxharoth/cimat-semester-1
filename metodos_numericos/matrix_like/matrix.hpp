#ifndef REAL_MATRIX_HPP
#define REAL_MATRIX_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <future>
#include <initializer_list>
#include <string.h>
#include <vector>

#ifndef MIN_OPER_FOR_THREAD
#define MIN_OPER_FOR_THREAD 500
#endif

namespace mymtx {
class matrix;
class vector;
class vector_iterator;
class const_vector_iterator;

class matrix {
  matrix();
  double *data;

public:
  const size_t shape_x, shape_y;
  matrix(size_t shape_y, size_t shape_x);
  matrix(const matrix &other);
  matrix(std::initializer_list<std::initializer_list<double>> initial);
  matrix(matrix &&other);
  ~matrix();
  vector operator[](const size_t &row);
  double &operator()(const size_t row, const size_t col);
  const double &operator()(const size_t row, const size_t col) const;
  const vector operator[](const size_t &row) const;
  vector_iterator begin(const size_t row);
  vector_iterator end(const size_t row);
  const_vector_iterator begin(const size_t row) const;
  const_vector_iterator end(const size_t row) const;
  matrix &operator*=(const double coef);
  vector &operator*=(vector &vec) const;
  vector operator*(const vector &vec) const;
  matrix &operator*=(const matrix &other);
  matrix operator*(const matrix &other) const;
  matrix &operator=(const matrix &other);
  matrix &operator-=(const matrix &other);
  matrix operator-(const matrix &other) const;
  static matrix traspose(const matrix &m);
  static matrix identity(const size_t n);
  static matrix tridiag(const size_t n, double (*low)(const int i),
                        double (*dig)(const int i), double (*up)(const int i));
  static matrix tridiag(const size_t n, const double low, const double dig,
                        const double up);
  static matrix pendiag(const size_t n, const double coef[5]);
  matrix &prod_as_band(const double coef, const size_t height,
                       const size_t width);
  vector prod_as_band(vector &vec, const size_t height,
                      const size_t width) const;
  matrix prod_as_band(const matrix &other, const size_t height,
                      const size_t width);
  static void fwrite(const char *filename, const matrix &matrix);
  static matrix fread(const char *filename);
  class Column {
    Column();
    double *data;

  public:
    const size_t size;
    const size_t increment;
    Column(const matrix::Column &other);
    Column(double *data, size_t size, const size_t increment);
    double &operator[](const size_t row);
    const double &operator[](const size_t row) const;
    double operator*(const Column &other) const;
    double distance() const;
    matrix as_matrix() const;
    vector as_vector() const;
    Column &operator=(const vector &other);
    vector operator-(const Column &other) const;
    matrix::Column &operator-=(const vector &other);
  };
  Column column(const size_t col);
  const Column column(const size_t col) const;
  vector operator*(const Column &other) const;
  friend class vector;
};
class vector {
  vector();
  double *data;
  const bool allocated;

public:
  const size_t size;
  vector(const size_t size);
  vector(const vector &other);
  vector(vector &&other);
  vector(std::initializer_list<double> initial);
  vector(double *data, size_t size);
  ~vector();
  double &operator[](const size_t col);
  const double &operator[](const size_t col) const;
  vector_iterator begin();
  vector_iterator end();
  const_vector_iterator begin() const;
  const_vector_iterator end() const;
  vector &operator*=(const double coef);
  vector operator/(const double coef);
  vector &operator/=(const double coef);
  double operator*(const vector &other) const;
  vector operator+(const vector &other);
  vector &operator+=(const vector &other);
  vector operator-(const vector &other);
  vector &operator-=(const vector &other);
  vector &operator=(const vector &other);
  vector &operator=(const double coef);
  vector subvector(size_t from, size_t size);
  matrix cross_product(const vector &other) const;
  double distance() const;
  matrix as_matrix() const;
  static vector normal(const size_t size);
  static void sort(vector &v);
  static void fwrite(const char *filename, const vector &vector);
  static vector fread(const char *filename);
  friend class matrix;
};
class vector_iterator {
  vector_iterator() {}

protected:
  double *data;
  size_t row, col;

public:
  vector_iterator(double *data, const size_t row, const size_t col);
  vector_iterator(const vector_iterator &other);
  size_t get_row() const;
  size_t get_col() const;
  vector_iterator &operator++();
  vector_iterator operator++(int);
  vector_iterator &operator--();
  vector_iterator operator--(int);
  vector_iterator &operator+=(int c);
  vector_iterator operator+(int c) const;
  vector_iterator &operator-=(int c);
  vector_iterator operator-(int c) const;
  bool operator==(const vector_iterator &rhs) const;
  bool operator!=(const vector_iterator &rhs) const;
  bool operator>(const vector_iterator &rhs) const;
  bool operator<(const vector_iterator &rhs) const;
  bool operator>=(const vector_iterator &rhs) const;
  bool operator<=(const vector_iterator &rhs) const;
  double &operator[](int);
  double &operator*();
};
class const_vector_iterator {
  const_vector_iterator() {}

protected:
  const double *data;
  size_t row, col;

public:
  const_vector_iterator(const double *data, const size_t row, const size_t col);
  const_vector_iterator(const const_vector_iterator &other);
  size_t get_row() const;
  size_t get_col() const;
  const_vector_iterator &operator++();
  const_vector_iterator operator++(int);
  const_vector_iterator &operator--();
  const_vector_iterator operator--(int);
  const_vector_iterator &operator+=(int c);
  const_vector_iterator operator+(int c) const;
  const_vector_iterator &operator-=(int c);
  const_vector_iterator operator-(int c) const;
  bool operator==(const const_vector_iterator &rhs) const;
  bool operator!=(const const_vector_iterator &rhs) const;
  bool operator>(const const_vector_iterator &rhs) const;
  bool operator<(const const_vector_iterator &rhs) const;
  bool operator>=(const const_vector_iterator &rhs) const;
  bool operator<=(const const_vector_iterator &rhs) const;
  const double &operator[](int) const;
  const double &operator*() const;
};
class MatrixTraspose {
  const matrix &t;

public:
  const size_t shape_x;
  const size_t shape_y;
  const double &operator()(size_t row, size_t col) const;
  vector operator*(const vector &vec) const;
  MatrixTraspose(const matrix &m);
  MatrixTraspose(const MatrixTraspose &other);
};

void abs(matrix &A);
vector map(const vector &v, double (*callback)(const double));
double reduce(const vector &v, double (*callback)(double acc, const double cur),
              const double start);
} // namespace mymtx
mymtx::matrix operator*(const mymtx::matrix &mtx, const double c);
mymtx::matrix operator*(const double c, const mymtx::matrix &mtx);
mymtx::matrix operator*(const mymtx::matrix &mtx,
                        const mymtx::MatrixTraspose &mtxt);
mymtx::matrix operator*(const mymtx::MatrixTraspose &mtxt,
                        const mymtx::matrix &mtx);
mymtx::vector operator*(const mymtx::vector &v, const double c);
mymtx::vector operator*(const double c, const mymtx::vector &v);
double operator*(const mymtx::matrix::Column &c, const mymtx::vector &v);
double operator*(const mymtx::vector &v, const mymtx::matrix::Column &c);
mymtx::vector operator*(const mymtx::matrix::Column &c, double &coef);
mymtx::vector operator*(double &coef, const mymtx::matrix::Column &c);
mymtx::vector operator+(const mymtx::matrix::Column &c, const mymtx::vector &v);
mymtx::vector operator+(const mymtx::vector &v, const mymtx::matrix::Column &c);
mymtx::vector operator-(const mymtx::matrix::Column &c, const mymtx::vector &v);
mymtx::vector operator-(const mymtx::vector &v, const mymtx::matrix::Column &c);
#endif /* REAL_MATRIX_HPP */
