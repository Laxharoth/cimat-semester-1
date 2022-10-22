#include "function_wrapper/function_wrapper.hpp"
#include "interpolation.hpp"
#include "macros.hpp"
#include <cmath>
#include <vector>
class sine : public FunctionWrapper {
  double eval(const double &x) const { return std::sin(x); }
  double eval(const double &x) { return std::sin(x); }
};
class fffn : public MultiVarFunctionWrapper {
  std::vector<double> eval(const std::vector<double> &x) const {
    return std::vector<double>{x[0] + x[1]};
  }
  std::vector<double> eval(const std::vector<double> &x) {
    return std::vector<double>{x[0] + x[1]};
  }
};
int main(int argc, const char **argv) {
  sine fn;

  auto normal_density = [](const double x) {
    return 1 / (10 * M_2_SQRTPI) * std::exp(-.5 * (x * x) / 100);
  };
  auto testfn = [](const double x) { return x + x * std::sin(x / 2) / 3; };

  auto mfn = &testfn;

  mymtx::vector xs(150);
  randomize(xs);
  const mymtx::vector ys = mymtx::map(xs, *mfn);
  std::vector<point> points(xs.size);
  Count count(0, 50, 100);
  const mymtx::vector Xs = mymtx::map(mymtx::vector(101), &count);
  mymtx::vector::fwrite("vec_out/finite/xs.vec", Xs);
  mymtx::vector::fwrite("vec_out/finite/ys.vec", mymtx::map(Xs, *mfn));
  for (size_t i = 0; i < xs.size; i++) {
    points[i].x = xs[i];
    points[i].y = ys[i];
  }
  // line splines
  FiniteElement line(points, 100, 0.2);
  {
    mymtx::vector xs(30);
    randomize(xs);
    mymtx::vector Ys = mymtx::map(Xs, &line);
    double error = (mymtx::map(xs, testfn) - mymtx::map(xs, &line)).distance();
    mymtx::vector errorvec(100);
    randomize(errorvec);
    strm_out("Error (Diferencias Finitas):" << error);
    mymtx::vector::fwrite("vec_out/finite/finite.vec", Ys);
  }
  strm_out("area (sin(x)) from 0 to pi:" << area_montecarlo(fn, 0, M_PI));
  strm_out("area (sin(x)) from 0 to 2pi:" << area_montecarlo(fn, 0, 2 * M_PI));

  point p0{0, 0};
  point p1{10, 10};
  fffn mvfn;
  strm_out("volumen (x+y) from 0  to 10:" << volum_montecarlo(mvfn, p0, p1));

  return 0;
}