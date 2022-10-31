#include "interpolation.hpp"
#include "function_wrapper/function_wrapper.hpp"
#include "matrix_like/matrix.hpp"
#include <cmath>
#include <cstddef>
#include <functional>
#include <tuple>
#include <vector>

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
  // TODO: programacion dinamica para m
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
  // TODO: programacion dinamica para a
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
    hinext = points[i + 2].x - points[i + 1].x;
    tinext = points[i + 2].y - points[i + 1].y;

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
  // TODO: programacion dinamica para b y d
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

FiniteElement::FiniteElement(const std::vector<point> &p_points,
                             const unsigned int nodes, const double lambda)
    : points(p_points), phi(mymtx::vector(nodes)) {
  std::sort(points.begin(), points.end(),
            [](const point &a, const point &b) { return a.x < b.x; });
  increment = ((points[points.size() - 1].x - points[0].x) / nodes);
  mymtx::matrix N_mtx(nodes, nodes);
  mymtx::vector yi(nodes);
  // compute Ni and ys
  size_t current_point_i = 0;
  double current_x = points[0].x;
  for (size_t diagonal = 0; diagonal < nodes - 1; ++diagonal) {
    N_mtx(diagonal, diagonal) += lambda / increment;
    N_mtx(diagonal + 1, diagonal + 1) += lambda / increment;
    N_mtx(diagonal + 1, diagonal) -= lambda / increment;
    N_mtx(diagonal, diagonal + 1) -= lambda / increment;
    while (current_point_i < points.size() &&
           points[current_point_i].x < current_x + increment) {
      double Ni = 1 - (points[current_point_i].x - current_x) / increment;
      double Ni_p_1 = (points[current_point_i].x - current_x) / increment;
      N_mtx(diagonal, diagonal) += Ni * Ni;
      N_mtx(diagonal + 1, diagonal + 1) += Ni_p_1 * Ni_p_1;
      N_mtx(diagonal + 1, diagonal) += Ni * Ni_p_1;
      N_mtx(diagonal, diagonal + 1) += Ni * Ni_p_1;
      yi[diagonal] += Ni * points[current_point_i].y;
      yi[diagonal + 1] += Ni_p_1 * points[current_point_i].y;
      ++current_point_i;
    }
    current_x += increment;
  }
  factor_cholesky(N_mtx, N_mtx);
  solve_cholesky(N_mtx, phi, yi);
}
double FiniteElement::eval(const double &x) const {
  double first = points[0].x;
  size_t index{0};
  while (first + increment <= x) {
    first += increment;
    ++index;
  }
  index = std::min(index, points.size() - 2);
  return phi[index] + ((x - first) / increment) * (phi[index + 1] - phi[index]);
}
double FiniteElement::eval(const double &x) {
  return static_cast<const FiniteElement *>(this)->eval(x);
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
double _area_montecarlo_2d(FunctionWrapper &fn, const double x0,
                           const double x1, double max_eval, sign s);
double area_montecarlo(FunctionWrapper &fn, double x0, double x1) {
  double area = 0;
  double x_min;
  double val_min = 0;
  double x_max;
  double val_max = 0;
  const double delta = 0.05;
  if (x0 > x1)
    std::swap(x0, x1);
  // find min and max
  for (double x = x0; x <= x1 + EP; x += delta) {
    {
      const double cur_eval = fn(x);
      if (cur_eval > val_max) {
        val_max = cur_eval;
        x_max = x;
      }
      if (cur_eval < val_min) {
        val_min = cur_eval;
        x_min = x;
      }
    }
  }
  if (val_min < 0)
    area += _area_montecarlo_2d(fn, x0, x1, val_min, sign::NEGATIVE);
  if (val_max > 0)
    area += _area_montecarlo_2d(fn, x0, x1, val_max, sign::POSITIVE);
  return area;
}
double _area_montecarlo_2d(FunctionWrapper &fn, const double x0,
                           const double x1, double max_eval, sign s) {
  const unsigned long points_increment = 3000;
  unsigned long total_points = 0;
  unsigned long coutn_in = 0;
  double current_area = 0;
  double prev_area;
  const double sign = (max_eval < 0) ? -1 : 1;
  const double height = std::abs(max_eval * 2);
  const double total_area = (x1 - x0) * height * sign;
  auto &rng = randgen::get_randgen();
  std::vector<double> point(2);
  do {
    auto xs = rng.sample(points_increment, x0, x1);
    auto ys = rng.sample(points_increment, 0, height);
    auto xi = xs.begin();
    auto yi = ys.begin();
    for (; xi != xs.end(); ++xi, ++yi) {
      double vali = fn(*xi);
      if (*yi < sign * vali)
        coutn_in++;
    }
    total_points += points_increment;
    prev_area = current_area;
    current_area = ((double)(coutn_in * total_area)) / ((double)(total_points));
  } while (std::abs(current_area - prev_area) < EP);

  return current_area;
}
double _vol_montecarlo(MultiVarFunctionWrapper &fn, point p0, point p1,
                       double height);
double volum_montecarlo(MultiVarFunctionWrapper &fn, point p0, point p1) {
  double area = 0;
  std::vector<double> point(2);
  std::vector<double> point_min(2);
  double val_min = 0;
  std::vector<double> point_max(2);
  double val_max = 0;
  const double delta = 0.05;
  if (p0.x > p1.x)
    std::swap(p0.x, p1.x);
  if (p0.y > p1.y)
    std::swap(p0.y, p1.y);
  // find min and max
  for (double x = p0.x; x <= p1.x + EP; x += delta) {
    for (double y = p0.y; y <= p1.y + EP; y += delta) {
      point[0] = x;
      point[1] = y;
      const double cur_eval = fn(point)[0];
      if (cur_eval > val_max) {
        val_max = cur_eval;
        point_max = point;
      }
      if (cur_eval < val_min) {
        val_min = cur_eval;
        point_min = point;
      }
    }
  }
  val_max = fn(point_max)[0];
  val_min = fn(point_min)[0];
  if (val_max > 0)
    area += _vol_montecarlo(fn, p0, p1, val_max);
  if (val_min < 0)
    area += _vol_montecarlo(fn, p0, p1, val_min);

  return area;
}
double _vol_montecarlo(MultiVarFunctionWrapper &fn, point p0, point p1,
                       double max_eval) {
  const unsigned long points_increment = 10000;
  unsigned long total_points = 0;
  unsigned long coutn_in = 0;
  double current_area = 0;
  double prev_area;
  const double sign = (max_eval < 0) ? -1 : 1;
  const double height = std::abs(max_eval * 2);
  const double total_volum = (p1.y - p0.y) * (p1.x - p0.x) * height * sign;
  auto &rng = randgen::get_randgen();
  std::vector<double> point(2);
  do {
    auto xs = rng.sample(points_increment, p0.x, p1.x);
    auto ys = rng.sample(points_increment, p0.y, p1.y);
    auto zs = rng.sample(points_increment, 0, height);
    auto xi = xs.begin();
    auto yi = ys.begin();
    auto zi = zs.begin();
    for (; xi != xs.end(); ++xi, ++yi, ++zi) {
      point[0] = *xi;
      point[1] = *yi;
      double vali = fn(point)[0];
      if (*zi < sign * vali)
        coutn_in++;
    }
    total_points += points_increment;
    prev_area = current_area;
    current_area =
        ((double)(coutn_in * total_volum)) / ((double)(total_points));
  } while (std::abs(current_area - prev_area) < EP);

  return current_area;
}
double bisection(FunctionWrapper &funcion, double x_inferior,
                 double x_superior) {
  auto calc_err = [](double &anterior, double &actual) {
    return std::abs((actual - anterior)) / std::abs(actual);
  };
  int inner_iter{0};
  int *real_iter;
  if (x_inferior == x_superior)
    throw - 1;
  if (x_inferior > x_superior) {
    std::swap(x_inferior, x_superior);
  }
  double y_inferior = funcion(x_inferior);
  double y_superior = funcion(x_superior);
  if (std::abs(y_inferior) <= EP)
    return x_inferior;
  if (std::abs(y_superior) <= EP)
    return x_superior;
  double x_medio{}, y_medio{};
  real_iter = &inner_iter;

  while ((*real_iter) < 1000) {
    ++(*real_iter);
    if (y_inferior * y_superior > 0)
      throw - 1;
    x_medio = (x_inferior + x_superior) / 2;
    y_medio = funcion(x_medio);
    if (std::abs(y_medio) <= EP)
      return x_medio;
    if (y_medio * y_inferior > 0) {
      x_inferior = x_medio;
    } else {
      x_superior = x_medio;
    }
  }
  return x_medio;
}

std::vector<double> NewtonMultivar(MultiVarFunctionWrapper &fn,
                                   const std::vector<double> &start_guess) {
  Gradient dfn(&fn);
  auto current_guess = start_guess;
  for (size_t i = 0; i < 100; i++) {
    double v = fn(start_guess)[0];
    if (std::abs(v) < EP)
      return current_guess;
    auto dv = dfn(start_guess);
    double dv_norm = norm(dv);
    current_guess = current_guess - (dv * v / (dv_norm * dv_norm));
  }
  if (std::any_of(current_guess.begin(), current_guess.end(),
                  [](double v) { return std::isnan(v) || std::isinf(v); }))
    current_guess = start_guess;
  return current_guess;
}

double norm(std::vector<double> &vec) {
  double sum{0};
  for (auto &&v : vec) {
    sum += v * v;
  }
  return std::sqrt(sum);
}

double integral_newton_cotes(FunctionWrapper &fn, const double from,
                             const double to, const unsigned grade) {
  /* clang-format off */
  static const std::vector<std::vector<double>> cotes_constans{
		{1/2.0   , 1/2.0},
		{1/3.0   , 4/3.0    , 1/3.0},
		{3/8.0   , 9/8.0    , 9/8.0   , 3/8.0},
		{14/45.0  , 64/45.0  , 24/45.0 , 64/45.0  , 14/45.0},
		{5.0/288 * 19, 5.0/288 * 75, 5.0/288 * 50, 5.0/288 * 50, 5.0/288 * 75, 5.0/288 * 19 },
		{1.0/140 * 41, 1.0/140 * 216, 1.0/140 * 27, 1.0/140 * 272, 1.0/140 * 27, 1.0/140 * 216, 1.0/140 * 41},
	};
  /* clang-format on */
  double integral{0};
  const double increment = (to - from) / (grade);
  double x = from;
  for (auto &&c : cotes_constans[grade - 1]) {
    integral += c * fn(x);
    x += increment;
  }
  return integral * increment;
}
double romberg_method(FunctionWrapper &fn, const double from, const double to,
                      const int maxRows, const double tolerance) {
  return richardson_extrapolation(fn, integral_newton_cotes, from, to, maxRows,
                                  tolerance, 1);
}
double richardson_extrapolation(FunctionWrapper &f, aproximation_function aprox,
                                double from, double to, const int maxRows,
                                const double tolerance,
                                const unsigned int grade) {
  // only require current row and previous
  mymtx::vector buff1(maxRows);
  mymtx::vector buff2(maxRows);
  mymtx::vector *_row = &buff1;
  mymtx::vector *_row_prev = &buff2;
  double h = (to - from);
  (*_row_prev)[0] = aprox(f, from, from + h, grade);
  for (size_t i = 1; i < maxRows; i++) {
    // get references for cleaner code
    mymtx::vector &row = *_row;
    mymtx::vector &row_prev = *_row_prev;
    h /= 2;
    // initilize row (since is not 0 due to previous iterations)
    row[0] = 0;
    for (size_t j = 0; j < (2 << (i - 1)); j++) {
      row[0] += aprox(f, from + j * h, from + (j + 1) * h, grade);
    }
    row[0] /= 2;
    row[0] += row_prev[0] / 2;
    for (size_t j = 1; j <= i; j++) {
      int n_k = 2 << j;
      row[j] = row[j - 1] + (row[j - 1] - row_prev[j - 1]) / (n_k - 1);
    }
    if (std::abs(row[i] - row_prev[i - 1]) < tolerance) {
      return row[i];
    }
    // make current row previous and recycle previous as new current
    std::swap(_row, _row_prev);
  }
  // because pointers were swaped last computed row is _row_prev
  return (*_row_prev)[maxRows - 1];
}
struct gaussian_cuadrature_pair {
  double weight;
  double point;
};
double gaussian_cuadrature(FunctionWrapper &fn, const double from,
                           const double to, const unsigned grade) {
  typedef gaussian_cuadrature_pair gp;
  static const std::vector<std::vector<gp>> gauss_constans{
      {gp{2, 0}},
      {gp{1, -1 / std::sqrt(3)}, gp{1, 1 / std::sqrt(3)}},
      {
          gp{5 / 9.0, -3 / std::sqrt(5)},
          gp{8 / 9.0, 0},
          gp{5 / 9.0, 3 / std::sqrt(5)},
      },
      {
          gp{(18 - std::sqrt(30)) / 36,
             -std::sqrt(1 / 7.0 * (3 + 2 * std::sqrt(6 / 5.0)))},
          gp{(18 + std::sqrt(30)) / 36,
             -std::sqrt(1 / 7.0 * (3 - 2 * std::sqrt(6 / 5.0)))},
          gp{(18 + std::sqrt(30)) / 36,
             std::sqrt(1 / 7.0 * (3 - 2 * std::sqrt(6 / 5.0)))},
          gp{(18 - std::sqrt(30)) / 36,
             std::sqrt(1 / 7.0 * (3 + 2 * std::sqrt(6 / 5.0)))},
      }};
  const int points = (grade + 1) / 2;
  double aprox = 0.0;
  for (auto &&pair : gauss_constans[points]) {
    aprox += pair.weight * fn((to - from) / 2 * pair.point + (from + to / 2));
  }
  return aprox * (to - from) / 2;
}
