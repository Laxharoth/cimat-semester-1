#include "function_wrapper/function_wrapper.hpp"
#include "interpolation.hpp"
#include "macros.hpp"
#include "matrix_like/funcion_matriz.hpp"
#include "matrix_like/matrix.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

int main(int argc, const char **argv) {

  auto normal_density = [](const double x) {
    return 1 / (10 * M_2_SQRTPI) * std::exp(-.5 * (x * x) / 100);
  };
  auto testfn = [](const double x) { return x + x * std::sin(x / 2) / 3; };

  mymtx::vector xs(30);
normal : {
  randomize(xs);
  const mymtx::vector ys = mymtx::map(xs, normal_density);
  std::vector<point> points(xs.size);
  Count count(0, 50, 100);
  const mymtx::vector Xs = mymtx::map(mymtx::vector(101), &count);
  mymtx::vector::fwrite("vec_out/normal/xs.vec", Xs);
  mymtx::vector::fwrite("vec_out/normal/ys.vec",
                        mymtx::map(Xs, normal_density));
  for (size_t i = 0; i < xs.size; i++) {
    points[i].x = xs[i];
    points[i].y = ys[i];
  }
  // line splines
  LineSpline line(points);
  {
    mymtx::vector Ys = mymtx::map(Xs, &line);
    double error = (ys - mymtx::map(xs, &line)).distance();
    mymtx::vector errorvec(30);
    randomize(errorvec);
    strm_out(
        "Error (line normal density):"
        << (mymtx::map(errorvec, line) - mymtx::map(errorvec, normal_density))
               .distance());
    mymtx::vector::fwrite("vec_out/normal/line_spline.vec", Ys);
  }
  // cuad splines
  CuadraticSpline para(points);
  {
    mymtx::vector Ys = mymtx::map(Xs, &para);
    mymtx::vector errorvec(30);
    randomize(errorvec);
    strm_out(
        "Error (cuadratic normal density):"
        << (mymtx::map(errorvec, para) - mymtx::map(errorvec, normal_density))
               .distance());
    mymtx::vector::fwrite("vec_out/normal/cuad_spline.vec", Ys);
  }
  // cub splines
  CubicSpline cub(points);
  {
    mymtx::vector Ys = mymtx::map(Xs, &cub);
    mymtx::vector errorvec(30);
    randomize(errorvec);
    strm_out(
        "Error (cubic normal density):"
        << (mymtx::map(errorvec, cub) - mymtx::map(errorvec, normal_density))
               .distance());
    mymtx::vector::fwrite("vec_out/normal/cubic_spline.vec", Ys);
  }
}
fn : {
  randomize(xs);
  const mymtx::vector ys = mymtx::map(xs, testfn);
  std::vector<point> points(xs.size);
  Count count(0, 50, 100);
  const mymtx::vector Xs = mymtx::map(mymtx::vector(101), &count);
  mymtx::vector::fwrite("vec_out/interp2/xs.vec", Xs);
  mymtx::vector::fwrite("vec_out/interp2/ys.vec", mymtx::map(Xs, testfn));
  for (size_t i = 0; i < xs.size; i++) {
    points[i].x = xs[i];
    points[i].y = ys[i];
  }
  // line splines
  LineSpline line(points);
  {
    mymtx::vector Ys = mymtx::map(Xs, &line);
    mymtx::vector errorvec(30);
    randomize(errorvec);
    strm_out("Error (cuadratic x+sin(x/2)/3):"
             << (mymtx::map(errorvec, line) - mymtx::map(errorvec, testfn))
                    .distance());
    mymtx::vector::fwrite("vec_out/interp2/line_spline.vec", Ys);
  }
  // cuad splines
  CuadraticSpline para(points);
  {
    mymtx::vector Ys = mymtx::map(Xs, &para);
    mymtx::vector errorvec(30);
    randomize(errorvec);
    strm_out("Error (cuadratic x+sin(x/2)/3):"
             << (mymtx::map(errorvec, para) - mymtx::map(errorvec, testfn))
                    .distance());
    mymtx::vector::fwrite("vec_out/interp2/cuad_spline.vec", Ys);
  }
  // cub splines
  CubicSpline cub(points);
  {
    mymtx::vector Ys = mymtx::map(Xs, &cub);
    mymtx::vector errorvec(30);
    randomize(errorvec);
    strm_out("Error (cuadratic x+sin(x/2)/3):"
             << (mymtx::map(errorvec, cub) - mymtx::map(errorvec, testfn))
                    .distance());
    mymtx::vector::fwrite("vec_out/interp2/cubic_spline.vec", Ys);
  }
}
  return 0;
}