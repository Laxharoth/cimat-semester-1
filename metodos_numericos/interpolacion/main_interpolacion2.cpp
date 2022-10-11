#include "function_wrapper/function_wrapper.hpp"
#include "interpolation.hpp"
#include "matrix_like/funcion_matriz.hpp"
#include "matrix_like/matrix.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

int main(int argc, const char **argv) {

  auto testfn = [](const double x) { return x + x * std::sin(x / 2) / 3; };

  mymtx::vector xs(30);
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
    mymtx::vector::fwrite("vec_out/interp2/line_spline.vec", Ys);
  }
  // cuad splines
  CuadraticSpline para(points);
  {
    mymtx::vector Ys = mymtx::map(Xs, &para);
    mymtx::vector::fwrite("vec_out/interp2/cuad_spline.vec", Ys);
  }

  // cub splines
  CubicSpline cub(points);
  {
    mymtx::vector Ys = mymtx::map(Xs, &cub);
    mymtx::vector::fwrite("vec_out/interp2/cubic_spline.vec", Ys);
  }

  return 0;
}