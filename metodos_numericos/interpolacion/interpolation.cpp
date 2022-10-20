#include "interpolation.hpp"
#include "matrix_like/funcion_matriz.hpp"
#include "matrix_like/matrix.hpp"
#include <algorithm>
#include <cstdio>
#include <ratio>
PolyFunction interpolate_line(const mymtx::vector &X, const mymtx::vector &Y) {
  return interpolate_poly(X, Y, 1);
}
PolyFunction interpolate_poly(const mymtx::vector &X, const mymtx::vector &Y,
                              unsigned int grade) {
  mymtx::vector ab(grade + 1);
  mymtx::vector yx(grade + 1);
  mymtx::matrix A(grade + 1, grade + 1);
  double xpow;

  for (size_t i = 0; i <= grade; i++) {
    for (size_t j = 0; j < X.size; j++)
      yx[i] += Y[j] * pow(X[j], i);
  }
  for (size_t i = 0; i <= grade * 2; i++) {
    for (size_t j = 0; j < X.size; j++) {
      xpow += std::pow(X[j], i);
    }
    size_t current = i, qty = current + 1;
    size_t start = 0;
    if (current >= A.shape_x) {
      start = current - A.shape_x;
      current = A.shape_x - 1;
    }
    for (size_t j = start; j < qty && j < A.shape_y; ++j, --current) {
      A(current, j) = xpow;
    }
  }
  // gauss(A,ab,yx);
  crout(A, A, A);
  solucion_crout(A, ab, yx);
  return PolyFunction(ab);
}
MultiFunctionWrapper interpolate_funcs(const mymtx::vector &X,
                                       const mymtx::vector &Y,
                                       std::vector<FunctionWrapper *> fns) {
  mymtx::vector cs(fns.size());
  mymtx::vector ys(cs.size);
  mymtx::matrix A(cs.size, cs.size);
  for (size_t i = 0; i < fns.size(); i++) {
    for (size_t j = 0; j < X.size; j++) {
      A(i, 0) += fns[i]->eval(X[j]) * fns[0]->eval(X[j]);
      ys[i] += fns[i]->eval(X[j]);
    }
    for (size_t j = 1; j < fns.size(); j++) {
      for (size_t k = 0; k < X.size; k++) {
        A(i, j) += fns[i]->eval(X[k]) * fns[j]->eval(X[k]);
      }
    }
  }
  crout(A, A, A);
  solucion_crout(A, cs, ys);
  return MultiFunctionWrapper(fns, cs);
}
PolyFunction interpolate_poly_2(const mymtx::vector &X,
                                const mymtx::vector &Y) {
  mymtx::matrix A(X.size, X.size);
  mymtx::vector ys = Y;
  mymtx::vector as(Y.size);
  for (size_t i = 0; i < X.size; i++) {
    for (size_t j = 0; j < X.size; j++) {
      A(i, j) = std::pow(X[i], j);
    }
  }
  mymtx::matrix Q(A.shape_y, A.shape_x);
  mymtx::matrix R = Q;
  qr_decomposition(A, Q, R);
  solve_qr(Q, R, as, ys);
  return PolyFunction(as);
}
LagramFunction interpolate_lagram(const mymtx::vector &X,
                                  const mymtx::vector &Y) {
  return LagramFunction(X, Y);
}
NewtonPolyFunction interpolate_newton(const mymtx::vector &X,
                                      const mymtx::vector &Y) {
  return NewtonPolyFunction(X, Y);
}
double PolyFunction::eval(const double &x) {
  return const_cast<PolyFunction *>(this)->eval(x);
}
double PolyFunction::eval(const double &x) const {
  double sum = 0;
  for (auto a = A.begin(); a != A.end(); ++a)
    sum += *a * std::pow(x, a.get_col());

  return sum;
}
PolyFunction::PolyFunction(const mymtx::vector &v) : A(mymtx::vector(v)) {}
double MultiFunctionWrapper::eval(const double &x) {
  return const_cast<MultiFunctionWrapper *>(this)->eval(x);
}
double MultiFunctionWrapper::eval(const double &x) const {
  double sum = 0;
  auto c = coef.begin();
  for (auto i = fns.begin(); i != fns.end(); ++i, ++c) {
    sum += (*i)->eval(x) * (*c);
  }
  return sum;
}
MultiFunctionWrapper::MultiFunctionWrapper(std::vector<FunctionWrapper *> fns,
                                           mymtx::vector coef)
    : fns(fns), coef(coef) {}
LagramFunction::LagramFunction(const mymtx::vector &X, const mymtx::vector &Y)
    : X(mymtx::vector(X)), Y(mymtx::vector(Y)) {}
double LagramFunction::eval(const double &x) {
  return const_cast<LagramFunction *>(this)->eval(x);
}
double LagramFunction::eval(const double &x) const {
  auto Li = [&](const double x, const unsigned int i) {
    double prod = 1;
    for (size_t j = 0; j < X.size; j++) {
      if (i == j)
        continue;
      prod *= (x - X[j]) / (X[i] - X[j]);
    }
    return prod;
  };
  double sum = 0;
  for (size_t i = 0; i < X.size; i++) {
    sum += Y[i] * Li(x, i);
  }
  return sum;
}
double NewtonPolyFunction::NewtonPolyFunction::eval(const double &x) {
  return const_cast<NewtonPolyFunction *>(this)->eval(x);
}
double NewtonPolyFunction::NewtonPolyFunction::eval(const double &x) const {
  double acumulative_prod = 1;
  double sum = as[0];
  for (size_t i = 1; i < X.size; i++) {
    acumulative_prod *= (x - X[i - 1]);
    sum += acumulative_prod * as[i];
  }
  return sum;
}
double divided_difference(mymtx::matrix &computed, const mymtx::vector &X,
                          const mymtx::vector &Y, unsigned int from,
                          unsigned int to) {
  if (from == to)
    computed(from, to) = Y[from];
  if (!computed(from, to) && from != to) {
    computed(from, to) = (divided_difference(computed, X, Y, from + 1, to) -
                          divided_difference(computed, X, Y, from, to - 1)) /
                         (X[to] - X[from]);
  }
  return computed(from, to);
};
NewtonPolyFunction::NewtonPolyFunction(const mymtx::vector &X,
                                       const mymtx::vector &Y)
    : X(mymtx::vector(X)), Y(mymtx::vector(Y)), as(mymtx::vector(X.size)) {
  mymtx::matrix computed(X.size, X.size);
  divided_difference(computed, X, Y, 0, X.size - 1);
  as = computed[0];
}

unsigned int binarySearch(const std::vector<point> &arr, double x);
LineSpline::LineSpline(const std::vector<point> &points) : points(points) {
  std::sort(this->points.begin(), this->points.end(),
            [](const point &a, const point &b) { return a.x < b.x; });
}
double LineSpline::eval(const double &x) {
  return static_cast<const LineSpline *>(this)->eval(x);
}
double LineSpline::eval(const double &x) const {
  unsigned int n = binarySearch(points, x);
  const double m =
      (points[n + 1].y - points[n].y) / (points[n + 1].x - points[n].x);
  return points[n].y + m * (x - points[n].x);
}

CuadraticSpline::CuadraticSpline(const std::vector<point> &points)
    : points(points), Si(mymtx::vector(points.size())) {
  std::sort(this->points.begin(), this->points.end(),
            [](const point &a, const point &b) { return a.x < b.x; });
  mymtx::vector ti(Si.size - 1);
  mymtx::matrix A(ti.size, ti.size);
  for (size_t i = 0; i < ti.size; i++) {
    ti[i] = points[i + 1].y - points[i].y;
    auto Ai = (points[i + 2].x - points[i + 1].x) / 2;
    if (i > 0)
      A(i, i - 1) = Ai;
    A(i, i) = Ai;
  }
  auto subSi = Si.subvector(1, Si.size - 1);
  crout(A, A, A);
  solucion_crout(A, subSi, ti);
}
double CuadraticSpline::eval(const double &x) {
  return static_cast<const CuadraticSpline *>(this)->eval(x);
}
double CuadraticSpline::eval(const double &x) const {
  unsigned int n = binarySearch(points, x);
  const double dx = x - points[n].x;
  const double h = points[n + 1].x - points[n].x;
  const double a = (Si[n + 1] - Si[n]) / (2 * h);
  const double b = Si[n];
  const double c = points[n].y;
  return a * dx * dx + b * dx + c;
}

CubicSpline::CubicSpline(const std::vector<point> &points)
    : points(points), Si(mymtx::vector(points.size())) {
  std::sort(this->points.begin(), this->points.end(),
            [](const point &a, const point &b) { return a.x < b.x; });
  mymtx::vector ti(points.size() - 2);
  mymtx::matrix A(ti.size, ti.size);
  double hinext, hinow, tinext, tinow;
  hinow = points[1].x - points[0].x;
  tinow = points[0].y - points[1].y;
  hinext = points[2].x - points[1].x;
  tinext = points[2].y - points[1].y;
  for (size_t i = 0; i < ti.size; i++) {
    hinow = points[i + 2].x - points[i + 1].x;
    tinow = points[i + 2].y - points[i + 1].y;

    ti[i] = tinext / hinext - tinow / hinow;

    if (i > 0)
      A(i, i - 1) = hinow / 6;
    A(i, i) = (hinow + hinext) / 3;
    if (i < ti.size - 1)
      A(i, i + 1) = hinext / 6;

    hinow = hinext;
    tinow = tinext;
  }
  auto subSi = Si.subvector(1, Si.size - 2);
  crout(A, A, A);
  // conjugate_gradient(A, subSi, ti);
  solucion_crout(A, subSi, ti);
}
double CubicSpline::eval(const double &x) {
  return static_cast<const CubicSpline *>(this)->eval(x);
}
double CubicSpline::eval(const double &x) const {
  unsigned int n = binarySearch(points, x);
  const double h = points[n + 1].x - points[n].x;
  const double t = points[n + 1].y - points[n].y;
  const double dx = x - points[n].x;
  const double a = points[n].y;
  const double b = t / h - h / 3 * ((2 * Si[n] + Si[n + 1]) / 2);
  const double c = Si[n] / 2;
  const double d = (Si[n + 1] - Si[n]) / (3 * h);
  return a + b * dx + c * dx * dx + d * dx * dx * dx;
}

unsigned int binarySearch(const std::vector<point> &arr, double x) {
  unsigned int p = 0;
  unsigned int r = arr.size() - 1;
  if (x <= arr[p].x)
    return p;
  if (x >= arr[r].x)
    return r - 1;
  while (p < r) {
    const unsigned int mid = (p + r) / 2;
    if (arr[mid].x == x)
      return mid;
    if (arr[mid].x > x) {
      r = mid - 1;
      continue;
    }
    if (arr[mid].x < x)
      p = mid + 1;
  }
  if (x < arr[p].x)
    return p - 1;
  return p;
}

namespace mymtx {
vector map(const vector &v, FunctionWrapper *fn) {
  vector v_new(v.size);
  auto i = v.begin();
  for (auto j = v_new.begin(); j != v_new.end(); ++i, ++j)
    *j = fn->eval(*i);
  return v_new;
}
vector map(const vector &v, const FunctionWrapper &fn) {
  vector v_new(v.size);
  auto i = v.begin();
  for (auto j = v_new.begin(); j != v_new.end(); ++i, ++j)
    *j = fn.eval(*i);
  return v_new;
}
} // namespace mymtx
