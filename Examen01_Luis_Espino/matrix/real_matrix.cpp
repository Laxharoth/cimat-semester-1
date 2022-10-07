#include "matrix.hpp"

struct mrow_data {
  size_t begin;
  size_t end;
  const size_t row;
  const mymtx::matrix *mtx;
  const mymtx::matrix *other;
  mymtx::matrix *cpy;
};
struct mrowMT_data {
  size_t begin;
  size_t end;
  const size_t row;
  const mymtx::matrix *mtx;
  const mymtx::MatrixTraspose *other;
  mymtx::matrix *cpy;
};
struct mrowMV_data {
  const mymtx::vector *vec;
  mymtx::vector *result;
  const mymtx::matrix *mtx;
  const size_t row, begin, end;
};
struct mrowMC_data {
  const mymtx::matrix::Column *vec;
  mymtx::vector *result;
  const mymtx::matrix *mtx;
  const size_t row;
};
struct mrowMTV_data {
  const mymtx::vector *vec;
  mymtx::vector *result;
  const mymtx::MatrixTraspose *mtx;
  const size_t row;
};
void mult_row(mrow_data d) {
  for (size_t k = d.begin; k < d.end; ++k) {
    for (size_t j = 0; j < d.other->shape_x; ++j) {
      (*d.cpy)[d.row][j] += (*d.mtx)(d.row, k) * (*d.other)(k, j);
    }
  }
};
void mult_matrix_MT_row(mrowMT_data d) {
  for (size_t k = d.begin; k < d.end; ++k) {
    for (size_t j = 0; j < d.other->shape_x; ++j) {
      (*d.cpy)[d.row][j] += (*d.mtx)(d.row, k) * (*d.other)(k, j);
    }
  }
};
void mult_matrix_TM_row(mrowMT_data d) {
  for (size_t k = d.begin; k < d.end; ++k) {
    for (size_t j = 0; j < d.mtx->shape_x; ++j) {
      (*d.cpy)[d.row][j] += (*d.other)(d.row, k) * (*d.mtx)(k, j);
    }
  }
};
void mult_mtx_vec_Row(mrowMV_data d) {
  for (size_t j = d.begin; j < d.end; j++) {
    (*d.result)[d.row] += (*d.mtx)[d.row][j] * (*d.vec)[j];
  }
}
void mult_mtx_col_Row(mrowMC_data d) {
  (*d.result)[d.row] = d.mtx->operator[](d.row) * (*(d.vec));
}
void mult_mtxt_vec_Row(mrowMTV_data d) {
  (*d.result)[d.row] = 0;
  for (size_t i = 0; i < d.vec->size; i++) {
    (*d.result)[d.row] += d.mtx->operator()(d.row, i) * d.vec->operator[](i);
  }
}
namespace mymtx {

matrix::matrix(size_t shape_y, size_t shape_x)
    : shape_x(shape_x), shape_y(shape_y) {
  this->data = new double[this->shape_x * this->shape_y];
  memset(this->data, 0, this->shape_x * this->shape_y * sizeof(double));
}
matrix::matrix(const matrix &other)
    : shape_x(other.shape_x), shape_y(other.shape_y) {
  this->data = new double[this->shape_x * this->shape_y];
  memcpy(this->data, other.data,
         this->shape_x * this->shape_y * sizeof(double));
}
matrix::matrix(matrix &&other)
    : shape_x(other.shape_x), shape_y(other.shape_y) {
  this->data = other.data;
  other.data = nullptr;
}
matrix::matrix(std::initializer_list<std::initializer_list<double>> initial)
    : shape_y(initial.size()), shape_x((initial.begin())->size()) {
  this->data = new double[this->shape_x * this->shape_y];
  size_t i{0}, j{0};
  for (auto ii = initial.begin(); ii != initial.end(); ++i, ++ii, j = 0) {
    for (auto ij = (*ii).begin(); ij != (*ii).end(); ++j, ++ij)
      (*this)[i][j] = *ij;
  }
}
matrix::~matrix() {
  if (data != nullptr)
    delete[] data;
}
matrix matrix::traspose(const matrix &m) {
  matrix traposed(m.shape_x, m.shape_y);
  for (size_t i = 0; i < m.shape_y; ++i)
    for (size_t j = 0; j < m.shape_x; ++j)
      traposed[j][i] = m[i][j];
  return traposed;
}
matrix matrix::identity(const size_t n) {
  matrix identity(n, n);
  for (size_t i = 0; i < n; ++i)
    identity[i][i] = 1;
  return identity;
}
matrix matrix::tridiag(const size_t n, double (*low)(const int i),
                       double (*dig)(const int i), double (*up)(const int i)) {
  matrix trid(n, n);
  for (size_t i = 0; i < n; ++i) {
    if (i > 0)
      trid[i][i - 1] = low(i);
    trid[i][i] = dig(i);
    if (i < n - 1)
      trid[i][i + 1] = up(i);
  }
  return trid;
}
matrix matrix::tridiag(const size_t n, const double low, const double dig,
                       const double up) {
  matrix trid(n, n);
  for (size_t i = 0; i < n; ++i) {
    if (i > 0)
      trid[i][i - 1] = low;
    trid[i][i] = dig;
    if (i < n - 1)
      trid[i][i + 1] = up;
  }
  return trid;
}
matrix matrix::pendiag(const size_t n, const double coef[5]) {
  matrix trid(n, n);
  for (size_t i = 0; i < n; ++i) {
    if (i > 1)
      trid(i, i - 2) = coef[0];
    if (i > 0)
      trid(i, i - 1) = coef[1];
    trid(i, i) = coef[2];
    if (i + 1 <= n - 1)
      trid(i, i + 1) = coef[3];
    if (i + 2 <= n - 1)
      trid(i, i + 2) = coef[4];
  }
  return trid;
}
vector matrix::operator[](const size_t &row) {
  return vector(data + row * shape_x, shape_x);
}
double &matrix::operator()(const size_t row, const size_t col) {
  return this->data[row * shape_x + col];
}
const double &matrix::operator()(const size_t row, const size_t col) const {
  return this->data[row * shape_x + col];
}
const vector matrix::operator[](const size_t &row) const {
  return vector(data + row * shape_x, shape_x);
}
vector_iterator matrix::begin(const size_t row) {
  return vector_iterator(data + row * shape_x, row, 0);
}
vector_iterator matrix::end(const size_t row) {
  return vector_iterator(data + row * shape_x, row, shape_x);
}
const_vector_iterator matrix::begin(const size_t row) const {
  return const_vector_iterator(data + row * shape_x, row, 0);
}
const_vector_iterator matrix::end(const size_t row) const {
  return const_vector_iterator(data + row * shape_x, row, shape_x);
}
vector &matrix::operator*=(vector &vec) const {
  auto cpy = (*this) * vec;
  vec = cpy;
  return vec;
}
vector matrix::operator*(const vector &vec) const {
  vector result(vec.size);
#ifndef NO_ASYNC
  std::vector<std::future<void>> futures;
#endif
  for (size_t i = 0; i < shape_y; ++i) {
#ifndef NO_ASYNC
    if (shape_y > MIN_OPER_FOR_THREAD)
      futures.push_back(
          std::async(std::launch::async, mult_mtx_vec_Row,
                     mrowMV_data{&vec, &result, this, i, 0, shape_x}));
    else
#endif
      for (size_t j = 0; j < shape_x; j++) {
        result[i] += (*this)[i][j] * vec[j];
      }
  }
  return result;
}
vector matrix::operator*(const matrix::Column &vec) const {
  vector result(vec.size);
#ifndef NO_ASYNC
  std::vector<std::future<void>> futures;
#endif
  for (size_t i = 0; i < shape_y; ++i) {
#ifndef NO_ASYNC
    if (shape_y > MIN_OPER_FOR_THREAD)
      futures.push_back(std::async(std::launch::async, mult_mtx_col_Row,
                                   mrowMC_data{&vec, &result, this, i}));
    else
#endif
      for (size_t j = 0; j < vec.size; j++) {
        result[i] += (*this)(i, j) * vec[j];
      }
  }
  return result;
}
vector MatrixTraspose::operator*(const vector &vec) const {
  vector result(vec.size);
#ifndef NO_ASYNC
  std::vector<std::future<void>> futures;
#endif
  for (size_t i = 0; i < shape_y; ++i) {
#ifndef NO_ASYNC
    if (shape_y > MIN_OPER_FOR_THREAD)
      futures.push_back(std::async(std::launch::async, mult_mtxt_vec_Row,
                                   mrowMTV_data{&vec, &result, this, i}));
    else
#endif
      for (size_t j = 0; j < vec.size; j++) {
        result[i] += (*this)(i, j) * vec[j];
      }
  }
  return result;
}
matrix &matrix::operator*=(const matrix &other) {
  auto cpy = (*this) * other;
  (*this) = cpy;
  return *this;
}
matrix &matrix::operator*=(const double coef) {
  for (size_t i = 0; i < this->shape_y; i++) {
    (*this)[i] *= coef;
  }
  return *this;
}
// RealMatrix RealMatrix::operator*(const double coef) const{
//     RealMatrix cpy = *this;
//     for (size_t i = 0; i < this->shape_y; i++){
//         cpy[i]*=coef;
//     }
//     return *this;
// }
matrix matrix::operator*(const matrix &other) const {
  auto cpy = matrix(this->shape_y, other.shape_x);
  std::vector<std::future<void>> futures;
  for (size_t i = 0; i < this->shape_y; ++i) {
#ifndef NO_ASYNC
    if (shape_y > MIN_OPER_FOR_THREAD)
      futures.push_back(
          std::async(std::launch::async, mult_row,
                     mrow_data{0, other.shape_y, i, this, &other, &cpy}));
    else
#endif
      for (size_t k = 0; k < other.shape_y; ++k) {
        for (size_t j = 0; j < other.shape_x; ++j) {
          cpy(i, j) += (*this)(i, k) * other(k, j);
        }
      }
  }
  return cpy;
}
matrix matrix::prod_as_band(const matrix &other, const size_t height,
                            const size_t width) {
  auto cpy = matrix(this->shape_y, other.shape_x);
  size_t begin, end;
#ifndef NO_ASYNC
  std::vector<std::future<void>> futures;
#endif
  for (size_t i = 0; i < this->shape_y; ++i) {
    begin = i - height;
    end = 1 + i + width;
    if (i < height)
      begin = 0;
    if (i > this->shape_y - height)
      end = this->shape_x;
#ifndef NO_ASYNC
    if (end - begin > 500)
      futures.push_back(
          std::async(std::launch::async, mult_row,
                     mrow_data{begin, end, i, this, &other, &cpy}));
    else
#endif
      for (size_t k = begin; k < end; ++k) {
        for (size_t j = 0; j < other.shape_x; ++j) {
          cpy[i][j] += (*this)(i, k) * other(k, j);
        }
      }
  }
  return cpy;
}
matrix &matrix::operator=(const matrix &other) {
  memcpy(this->data, other.data,
         this->shape_x * this->shape_y * sizeof(double));
  return *this;
}
matrix &matrix::operator-=(const matrix &other) {
  for (int i = 0; i < this->shape_y; i++) {
    for (int j = 0; j < this->shape_x; j++) {
      (*this)[i][j] -= other[i][j];
    }
  }
  return *this;
}
matrix matrix::operator-(const matrix &other) const {
  auto cpy = *this;
  cpy -= other;
  return cpy;
}
matrix &matrix::prod_as_band(const double coef, const size_t height,
                             const size_t width) {
  size_t begin, end;
  for (size_t i = 0; i < this->shape_y; i++) {
    begin = i - height;
    end = 1 + i + width;
    if (i < height)
      begin = 0;
    if (i > this->shape_y - height)
      end = this->shape_x;
    for (size_t j = begin; j < end; j++) {
      (*this)[i][j] *= coef;
    }
  }
  return *this;
}
vector matrix::prod_as_band(vector &vec, const size_t height,
                            const size_t width) const {
  vector result(vec.size);
  size_t begin, end;
#ifndef NO_ASYNC
  std::vector<std::future<void>> futures;
#endif
  for (size_t i = 0; i < this->shape_y; i++) {
    begin = i - height;
    end = 1 + i + width;
    if (i < height)
      begin = 0;
    if (i > this->shape_y - height)
      end = this->shape_x;
    result[i] = 0;
#ifndef NO_ASYNC
    if (shape_y > MIN_OPER_FOR_THREAD)
      futures.push_back(
          std::async(std::launch::async, mult_mtx_vec_Row,
                     mrowMV_data{&vec, &result, this, i, begin, end}));
    else
#endif
      for (size_t j = begin; j < end; j++) {
        result[i] += (*this)[i][j] * vec[j];
      }
  }
  return result;
}
void matrix::fwrite(const char *filename, const matrix &matrix) {
  std::ofstream f(filename, std::ofstream::out | std::ofstream::binary);
  f.write((char *)(&(matrix.shape_y)), sizeof(size_t));
  f.write((char *)(&(matrix.shape_x)), sizeof(size_t));
  f.write((char *)(matrix.data),
          sizeof(double) * matrix.shape_y * matrix.shape_x);
  f.close();
}
void vector::fwrite(const char *filename, const vector &vector) {
  std::ofstream f(filename, std::ofstream::out | std::ofstream::binary);
  f.write((char *)(&(vector.size)), sizeof(size_t));
  f.write((char *)(vector.data), sizeof(double) * vector.size);
  f.close();
}
matrix matrix::fread(const char *filename) {
  size_t shape_x, shape_y;
  std::ifstream f(filename, std::ifstream::in | std::ifstream::binary);
  f.read((char *)(&(shape_y)), sizeof(size_t));
  f.read((char *)(&(shape_x)), sizeof(size_t));
  matrix matrix(shape_y, shape_x);
  f.read((char *)(matrix.data), sizeof(double) * shape_y * shape_y);
  f.close();
  return matrix;
}
vector vector::fread(const char *filename) {
  size_t size;
  std::ifstream f(filename, std::ifstream::in | std::ifstream::binary);
  f.read((char *)(&(size)), sizeof(size_t));
  vector vector(size);
  f.read((char *)(vector.data), sizeof(double) * size);
  f.close();
  return vector;
}
void abs(matrix &A) {
  for (size_t k = 0; k < A.shape_y; ++k) {
    for (auto i = A[k].begin(); i < A[k].end(); ++i) {
      *i = std::abs(*i);
    }
  }
}
void abs(vector &V) {
  for (auto i = V.begin(); i < V.end(); ++i) {
    *i = std::abs(*i);
  }
}

const double &MatrixTraspose::operator()(size_t row, size_t col) const {
  return t[col][row];
}
MatrixTraspose::MatrixTraspose(const matrix &m)
    : t(m), shape_y(m.shape_x), shape_x(m.shape_y) {}
MatrixTraspose::MatrixTraspose(const MatrixTraspose &other)
    : t(other.t), shape_y(other.shape_y), shape_x(other.shape_x) {}
matrix::Column matrix::column(const size_t col) {
  return matrix::Column(data + col, shape_y, shape_x);
}
const matrix::Column matrix::column(const size_t col) const {
  return matrix::Column(data + col, shape_y, shape_x);
}
} // namespace mymtx
mymtx::matrix operator*(const mymtx::matrix &mtx, const double c) {
  mymtx::matrix cpy = mtx;
  cpy *= c;
  return cpy;
}
mymtx::matrix operator*(const double c, const mymtx::matrix &mtx) {
  return mtx * c;
}
mymtx::matrix operator*(const mymtx::matrix &mtx,
                        const mymtx::MatrixTraspose &mtxt) {
  mymtx::matrix res(mtx.shape_y, mtxt.shape_x);
#ifndef NO_ASYNC
  std::vector<std::future<void>> futures;
#endif
  for (size_t i = 0; i < res.shape_y; ++i) {
#ifndef NO_ASYNC
    if (res.shape_y > MIN_OPER_FOR_THREAD)
      futures.push_back(
          std::async(std::launch::async, mult_matrix_MT_row,
                     mrowMT_data{0, mtxt.shape_x, i, &mtx, &mtxt, &res}));
    else
#endif
      for (size_t k = 0; k < mtxt.shape_y; ++k) {
        for (size_t j = 0; j < res.shape_x; ++j) {
          res[i][j] += mtx[i][k] * mtxt(k, j);
        }
      }
  }
  return res;
}
mymtx::matrix operator*(const mymtx::MatrixTraspose &mtxt,
                        const mymtx::matrix &mtx) {
  mymtx::matrix res(mtxt.shape_y, mtx.shape_x);
#ifndef NO_ASYNC
  std::vector<std::future<void>> futures;
#endif
  for (size_t i = 0; i < res.shape_y; ++i) {
#ifndef NO_ASYNC
    if (res.shape_y > MIN_OPER_FOR_THREAD)
      futures.push_back(
          std::async(std::launch::async, mult_matrix_TM_row,
                     mrowMT_data{0, mtxt.shape_x, i, &mtx, &mtxt, &res}));
    else
#endif
      for (size_t k = 0; k < mtx.shape_y; ++k) {
        for (size_t j = 0; j < res.shape_x; ++j) {
          res[i][j] += mtxt(i, k) * mtx(k, j);
        }
      }
  }
  return res;
}
#include "vector.cpp"
#include "vector_iterator.cpp"