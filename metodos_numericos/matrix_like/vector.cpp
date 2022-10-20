#include "matrix.hpp"
#include <cmath>
#include <cstddef>
#include <future>
namespace mymtx {
vector::vector(const size_t size) : size(size), allocated(true) {
  this->data = new double[this->size];
  memset(this->data, 0, this->size * sizeof(double));
}
vector::vector(const vector &other) : size(other.size), allocated(true) {
  this->data = new double[this->size];
  memcpy(this->data, other.data, this->size * sizeof(double));
}
vector::vector(vector &&other) : size(other.size), allocated(other.allocated) {
  this->data = other.data;
  other.data = nullptr;
}
vector::vector(std::initializer_list<double> initial)
    : size(initial.size()), allocated(true) {
  this->data = new double[this->size];
  size_t i{0};
  for (auto j = initial.begin(); j != initial.end(); ++j, ++i)
    data[i] = *j;
}
vector::vector(double *data, size_t size)
    : size(size), allocated(false), data(data) {}
vector::~vector() {
  if (allocated && data != nullptr)
    delete[] data;
}
double &vector::operator[](const size_t col) { return data[col]; }
const double &vector::operator[](const size_t col) const { return data[col]; }
vector_iterator vector::begin() { return vector_iterator(data, 0, 0); }
vector_iterator vector::end() { return vector_iterator(data, 0, size); }
const_vector_iterator vector::begin() const {
  return const_vector_iterator(data, 0, 0);
}
const_vector_iterator vector::end() const {
  return const_vector_iterator(data, 0, size);
}
struct vector_coef_data {
  size_t begin;
  size_t end;
  vector *v;
  double coef;
};
vector &vector::operator*=(const double coef) {
#ifndef NO_ASYNC
  if (size > MIN_OPER_FOR_THREAD) {
    std::vector<std::future<void>> futures;
    for (size_t i = 0; i < 3; ++i) {
      size_t begin = (size * i) / 3;
      size_t end = (size * (i + 1)) / 3;
      futures.push_back(std::async(
          std::launch::async,
          [](vector_coef_data d) {
            for (size_t i = d.begin; i < d.end; i++) {
              (*d.v)[i] *= d.coef;
            }
          },
          vector_coef_data{begin, end, this, coef}));
    }
    futures.clear();
  } else
#endif
  {
    for (auto i = this->begin(); i != this->end(); ++i)
      (*i) *= coef;
  }

  return *this;
}
vector vector::operator/(const double coef) {
  vector cpy = *this;
  cpy /= coef;
  return cpy;
}
vector &vector::operator/=(const double coef) {
  (*this) *= 1 / coef;
  return *this;
}
struct vector_vector_data {
  size_t begin;
  size_t end;
  const vector *const v1;
  const vector *const v2;
  double *const sum;
};
double vector::operator*(const vector &other) const {
  double sum{0};
#ifndef NO_ASYNC
  double sums[3];
  if (size > MIN_OPER_FOR_THREAD) {
    std::vector<std::future<void>> futures;
    for (size_t i = 0; i < 3; ++i) {
      size_t begin = (size * i) / 3;
      size_t end = (size * (i + 1)) / 3;
      futures.push_back(std::async(
          std::launch::async,
          [](vector_vector_data d) {
            (*d.sum) = 0;
            for (size_t i = d.begin; i < d.end; i++) {
              (*d.sum) += (*d.v1)[i] * (*d.v2)[i];
            }
          },
          vector_vector_data{begin, end, this, &other, sums + i}));
    }
    futures.clear();
    sum = sums[0] + sums[1] + sums[2];
  } else
#endif
  {
    auto it = other.begin();
    auto ti = this->begin();
    while (it != other.end()) {
      sum += (*it) * (*ti);
      ++it;
      ++ti;
    }
  }
  return sum;
}
struct vector_sum_vector_data {
  size_t begin;
  size_t end;
  vector *const v1;
  const vector *const v2;
};
vector &vector::operator+=(const vector &other) {
#ifndef NO_ASYNC
  if (size > MIN_OPER_FOR_THREAD) {
    std::vector<std::future<void>> futures;
    for (size_t i = 0; i < 3; ++i) {
      size_t begin = (size * i) / 3;
      size_t end = (size * (i + 1)) / 3;
      futures.push_back(std::async(
          std::launch::async,
          [](vector_sum_vector_data d) {
            for (size_t i = d.begin; i < d.end; i++) {
              (*d.v1)[i] += (*d.v2)[i];
            }
          },
          vector_sum_vector_data{begin, end, this, &other}));
    }
  } else
#endif
  {
    auto it = other.begin();
    auto ti = this->begin();
    while (it != other.end()) {
      (*ti) += (*it);
      ++it;
      ++ti;
    }
  }
  return *this;
}
vector vector::operator+(const vector &other) {
  auto cpy = *this;
  cpy += other;
  return cpy;
}
vector &vector::operator-=(const vector &other) {
#ifndef NO_ASYNC
  if (size > MIN_OPER_FOR_THREAD) {
    std::vector<std::future<void>> futures;
    for (size_t i = 0; i < 3; ++i) {
      size_t begin = (size * i) / 3;
      size_t end = (size * (i + 1)) / 3;
      futures.push_back(std::async(
          std::launch::async,
          [](vector_sum_vector_data d) {
            for (size_t i = d.begin; i < d.end; i++) {
              (*d.v1)[i] -= (*d.v2)[i];
            }
          },
          vector_sum_vector_data{begin, end, this, &other}));
    }
  } else
#endif
  {
    auto it = other.begin();
    auto ti = this->begin();
    while (it != other.end()) {
      (*ti) -= (*it);
      ++it;
      ++ti;
    }
  }
  return *this;
}
vector vector::operator-(const vector &other) const {
  auto cpy = *this;
  cpy -= other;
  return cpy;
}
vector &vector::operator=(const vector &other) {
  memcpy(this->data, other.data, sizeof(double) * this->size);
  return *this;
}
vector &vector::operator=(const double coef) {
  memset(this->data, coef, sizeof(double) * this->size);
  return *this;
}
struct vector_cross_vector_data {
  size_t begin;
  size_t end;
  const vector *const v1;
  const vector *const v2;
  matrix *res;
};
matrix vector::cross_product(const vector &other) const {
  matrix result(other.size, other.size);
#ifndef NO_ASYNC
  if (size > MIN_OPER_FOR_THREAD) {
    std::vector<std::future<void>> futures;
    for (size_t i = 0; i < 3; ++i) {
      size_t begin = (size * i) / 3;
      size_t end = (size * (i + 1)) / 3;
      futures.push_back(std::async(
          std::launch::async,
          [](vector_cross_vector_data d) {
            for (size_t i = d.begin; i < d.end; i++) {
              for (size_t j = 0; j < (*d.v2).size; j++) {
                (*d.res)(i, j) = (*d.v1)[i] * (*d.v2)[j];
              }
            }
          },
          vector_cross_vector_data{begin, end, this, &other, &result}));
    }
  } else
#endif
  {
    for (size_t i = 0; i < other.size; i++) {
      for (size_t j = 0; j < other.size; j++) {
        result[i][j] = other[i] * other[j];
      }
    }
  }
  return result;
}
vector vector::normal(const size_t size) {
  vector norm(size);
  const double val = 1 / std::sqrt(size);
  for (auto i = norm.begin(); i != norm.end(); ++i)
    *i = val;
  return norm;
}
void vector::sort(vector &v) { std::sort(v.data, v.data + v.size); }
double vector::distance() const {
  double sum = (*this) * (*this);
  return std::sqrt(sum);
}
matrix vector::as_matrix() const {
  matrix m(this->size, 1);
  memcpy(m.data, this->data, sizeof(double) * this->size);
  return m;
}
matrix::Column::Column(const matrix::Column &other)
    : size(other.size), increment(other.increment), data(other.data) {}
matrix::Column::Column(double *data, size_t size, const size_t increment)
    : data(data), size(size), increment(increment) {}
double &matrix::Column::operator[](const size_t row) {
  return *(data + row * increment);
}
const double &matrix::Column::operator[](const size_t row) const {
  return *(data + row * increment);
}
double matrix::Column::distance() const {
  double v = 0;
  for (size_t i = 0; i < size; i++)
    v += (*this)[i] * (*this)[i];
  return std::sqrt(v);
}
matrix matrix::Column::as_matrix() const {
  matrix result(size, 1);
  for (size_t i = 0; i < size; i++)
    result[i][0] = (*this)[i];
  return result;
}
vector matrix::Column::as_vector() const {
  vector v(size);
  for (int i = 0; i < this->size; i++)
    v[i] = (*this)[i];
  return v;
}
double matrix::Column::operator*(const matrix::Column &other) const {
  double r = 0;
  for (size_t i = 0; i < size; i++) {
    r += (*this)[i] * other[i];
  }
  return r;
}
matrix::Column &matrix::Column::operator=(const vector &other) {
  for (size_t i = 0; i < size; i++) {
    (*this)[i] = other[i];
  }
  return *this;
}
struct column_sum_vector_data {
  size_t begin;
  size_t end;
  matrix::Column *const v1;
  const vector *const v2;
};
matrix::Column &matrix::Column::operator-=(const vector &other) {
#ifndef NO_ASYNC
  if (size > MIN_OPER_FOR_THREAD) {
    std::vector<std::future<void>> futures;
    for (size_t i = 0; i < 3; ++i) {
      size_t begin = (size * i) / 3;
      size_t end = (size * (i + 1)) / 3;
      futures.push_back(std::async(
          std::launch::async,
          [](column_sum_vector_data d) {
            for (size_t i = d.begin; i < d.end; i++) {
              (*d.v1)[i] -= (*d.v2)[i];
            }
          },
          column_sum_vector_data{begin, end, this, &other}));
    }
  } else
#endif
  {
    for (size_t i = 0; i < size; i++) {
      (*this)[i] -= other[i];
    }
  }
  return *this;
}
vector matrix::Column::operator-(const matrix::Column &other) const {
  vector v(this->size);
  for (size_t i = 0; i < this->size; i++) {
    v[i] = (*this)[i] - other[i];
  }
  return v;
}
vector map(const vector &v, double (*callback)(const double)) {
  vector v_new(v.size);
  auto i = v.begin();
  for (auto j = v_new.begin(); j != v_new.end(); ++i, ++j)
    *j = callback(*i);
  return v_new;
}
double reduce(const vector &v, double (*callback)(double acc, const double cur),
              const double start) {
  double acc = start;
  for (auto i = v.begin(); i != v.end(); ++i)
    acc = callback(acc, *i);
  return acc;
}
vector vector::subvector(size_t from, size_t size) {
  return vector(data + from, size);
}
} // namespace mymtx
mymtx::vector operator*(const mymtx::vector &v, const double c) {
  mymtx::vector cpy = v;
  cpy *= c;
  return cpy;
}
mymtx::vector operator*(const double c, const mymtx::vector &v) {
  return v * c;
}
double operator*(const mymtx::matrix::Column &c, const mymtx::vector &v) {
  double r = 0;
  for (auto i = 0; i < v.size; ++i)
    r += c[i] * v[i];
  return r;
}
double operator*(const mymtx::vector &v, const mymtx::matrix::Column &c) {
  return c * v;
}
mymtx::vector operator*(const mymtx::matrix::Column &c, double &coef) {
  mymtx::vector v(c.size);
  for (size_t i = 0; i < c.size; i++) {
    v[i] = c[i] * coef;
  }
  return v;
}
mymtx::vector operator*(double &coef, const mymtx::matrix::Column &c) {
  return c * coef;
}
mymtx::vector operator+(const mymtx::matrix::Column &c,
                        const mymtx::vector &v) {
  mymtx::vector vr(c.size);
  for (size_t i = 0; i < c.size; i++)
    vr[i] = c[i] + v[i];
  return vr;
}
mymtx::vector operator+(const mymtx::vector &v,
                        const mymtx::matrix::Column &c) {
  return c + v;
}
mymtx::vector operator-(const mymtx::matrix::Column &c,
                        const mymtx::vector &v) {
  mymtx::vector vr(c.size);
  for (size_t i = 0; i < c.size; i++)
    vr[i] = c[i] - v[i];
  return vr;
}
mymtx::vector operator-(const mymtx::vector &v,
                        const mymtx::matrix::Column &c) {
  mymtx::vector vr(c.size);
  for (size_t i = 0; i < c.size; i++)
    vr[i] = v[i] + c[i];
  return vr;
}
