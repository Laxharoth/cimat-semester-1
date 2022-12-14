#include "function_wrapper/function_wrapper.hpp"
#include "interpolation.hpp"
#include "matrix_like/funcion_matriz.hpp"
#include "matrix_like/matrix.hpp"
#include "matrix_like/print.cpp"

#include <cmath>
#include <cstdio>

class cosine : public FunctionWrapper {
  double eval(const double &x) const { return std::cos(x); }
  double eval(const double &x) { return std::cos(x); }
};
class sine : public FunctionWrapper {
  double eval(const double &x) const { return std::sin(x); }
  double eval(const double &x) { return std::sin(x); }
};
class ln : public FunctionWrapper {
  double eval(const double &x) const { return std::log(x); }
  double eval(const double &x) { return std::log(x); }
};
class tow : public FunctionWrapper {
  double eval(const double &x) const { return std::pow(x, x); }
  double eval(const double &x) { return std::pow(x, x); }
};
class id : public FunctionWrapper {
  double eval(const double &x) const { return std::exp(x); }
  double eval(const double &x) { return std::exp(x); }
};
class ex : public FunctionWrapper {
  double eval(const double &x) const { return (x); }
  double eval(const double &x) { return (x); }
};
class func_1 : public FunctionWrapper {
  double eval(const double &x) const override { return 1.0 / (1 + 25 * x * x); }
  double eval(const double &x) override {
    return const_cast<const func_1 *>(this)->eval(x);
  }
};
class func_2 : public FunctionWrapper {
  double eval(const double &x) const override {
    return std::abs(x) - x / 2 - x * x;
  }
  double eval(const double &x) override {
    return const_cast<const func_2 *>(this)->eval(x);
  }
};

#define FN1 "fn1"
#define FN2 "fn2"
#define SIZE_5 "5"
#define SIZE_100 "100"
#define SIZE_1000 "1000"

int main(int argc, const char **argv) {
  Count count(-1, 1, 1000);
  mymtx::vector X_graph = mymtx::map(mymtx::vector(1001), &count);
  auto rys = mymtx::map(X_graph, func_1());
  /* Sample 5*/ {
    mymtx::vector X_sample(5);
    randomize(X_sample);
    mymtx::vector Y_sample = mymtx::map(X_sample, func_1());
    mymtx::vector::fwrite("vec_out/" FN1 "xs_point" SIZE_5 ".vec", X_sample);
    mymtx::vector::fwrite("vec_out/" FN1 "ys_point" SIZE_5 ".vec", Y_sample);
    /* Sample END */
    mymtx::vector::fwrite("vec_out/ys.vec", rys);
    mymtx::vector::fwrite("vec_out/"
                          "xs.vec",
                          X_graph);
    /* fn_1 n = 5*/ {
      {
        const FunctionWrapper &fn = interpolate_line(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 5 line: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "line_ys" SIZE_5 ".vec", Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_poly(X_sample, Y_sample, 4);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 5 polinomio: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "poly1_ys" SIZE_5 ".vec", Y_graph);
      }
      {
        std::vector<FunctionWrapper *> fns;
        cosine fn1;
        sine fn2;
        ln fn3;
        tow fn4;
        id fn5;
        fns.push_back(&fn1);
        fns.push_back(&fn2);
        // fns.push_back(&fn3);fns.push_back(&fn4);
        fns.push_back(&fn5);
        const FunctionWrapper &fn = interpolate_funcs(X_sample, Y_sample, fns);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 5 fns: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "fns_ys" SIZE_5 ".vec", Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_poly_2(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 5 polinomio n: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "poli2_ys" SIZE_5 ".vec", Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_lagram(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 5 lagram: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "lagr_ys" SIZE_5 ".vec", Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_newton(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 5 newton: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "newt_ys" SIZE_5 ".vec", Y_graph);
      }
    }
  }
  /* Sample 100 */ {
    mymtx::vector X_sample(100);
    randomize(X_sample);
    normalize(X_sample);
    mymtx::vector Y_sample = mymtx::map(X_sample, func_1());
    mymtx::vector::fwrite("vec_out/" FN1 "xs_point" SIZE_100 ".vec", X_sample);
    mymtx::vector::fwrite("vec_out/" FN1 "ys_point" SIZE_100 ".vec", Y_sample);
    /* Sample END */
    mymtx::vector::fwrite("vec_out/"
                          "xs.vec",
                          X_graph);
    /* fn_1 n = 5*/ {
      {
        const FunctionWrapper &fn = interpolate_line(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 100 line: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "line_ys" SIZE_100 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_poly(X_sample, Y_sample, 4);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 100 polinomio: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "poly1_ys" SIZE_100 ".vec",
                              Y_graph);
      }
      {
        std::vector<FunctionWrapper *> fns;
        cosine fn1;
        sine fn2;
        ln fn3;
        tow fn4;
        id fn5;
        fns.push_back(&fn1);
        fns.push_back(&fn2);
        // fns.push_back(&fn3);fns.push_back(&fn4);
        fns.push_back(&fn5);
        const FunctionWrapper &fn = interpolate_funcs(X_sample, Y_sample, fns);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 100 fns: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "fns_ys" SIZE_100 ".vec", Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_poly_2(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 100 polinomio n: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "poli2_ys" SIZE_100 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_lagram(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 100 lagram: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "lagr_ys" SIZE_100 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_newton(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 100 newton: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "newt_ys" SIZE_100 ".vec",
                              Y_graph);
      }
    }
  }
  /* Sample 1000 */ {
    mymtx::vector X_sample(1000);
    randomize(X_sample);
    normalize(X_sample);
    mymtx::vector Y_sample = mymtx::map(X_sample, func_1());
    mymtx::vector::fwrite("vec_out/" FN1 "xs_point" SIZE_1000 ".vec", X_sample);
    mymtx::vector::fwrite("vec_out/" FN1 "ys_point" SIZE_1000 ".vec", Y_sample);
    /* Sample END */
    mymtx::vector::fwrite("vec_out/"
                          "xs.vec",
                          X_graph);
    /* fn_1 n = 5*/ {
      {
        const FunctionWrapper &fn = interpolate_line(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 1000 line: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "line_ys" SIZE_1000 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_poly(X_sample, Y_sample, 4);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 1000 polinomio: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "poly1_ys" SIZE_1000 ".vec",
                              Y_graph);
      }
      {
        std::vector<FunctionWrapper *> fns;
        cosine fn1;
        sine fn2;
        ln fn3;
        tow fn4;
        id fn5;
        fns.push_back(&fn1);
        fns.push_back(&fn2);
        // fns.push_back(&fn3);fns.push_back(&fn4);
        fns.push_back(&fn5);
        const FunctionWrapper &fn = interpolate_funcs(X_sample, Y_sample, fns);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 1000 fns: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "fns_ys" SIZE_1000 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_poly_2(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 1000 polinomio n: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "poli2_ys" SIZE_1000 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_lagram(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 1000 lagram: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "lagr_ys" SIZE_1000 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_newton(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn1-Error 1000 newton: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN1 "newt_ys" SIZE_1000 ".vec",
                              Y_graph);
      }
    }
  }
  rys = mymtx::map(X_graph, func_2());
  /* Sample 5*/ {
    mymtx::vector X_sample(5);
    randomize(X_sample);
    normalize(X_sample);
    mymtx::vector Y_sample = mymtx::map(X_sample, func_2());
    mymtx::vector::fwrite("vec_out/" FN2 "xs_point" SIZE_5 ".vec", X_sample);
    mymtx::vector::fwrite("vec_out/" FN2 "ys_point" SIZE_5 ".vec", Y_sample);
    /* Sample END */
    mymtx::vector::fwrite("vec_out/ys2.vec", mymtx::map(X_graph, func_2()));
    mymtx::vector::fwrite("vec_out/"
                          "xs.vec",
                          X_graph);
    /* fn_1 n = 5*/ {
      {
        const FunctionWrapper &fn = interpolate_line(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 5 line: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "line_ys" SIZE_5 ".vec", Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_poly(X_sample, Y_sample, 4);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 5 polinomio: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "poly1_ys" SIZE_5 ".vec", Y_graph);
      }
      {
        std::vector<FunctionWrapper *> fns;
        cosine fn1;
        sine fn2;
        ln fn3;
        tow fn4;
        id fn5;
        fns.push_back(&fn1);
        fns.push_back(&fn2);
        // fns.push_back(&fn3);fns.push_back(&fn4);
        fns.push_back(&fn5);
        const FunctionWrapper &fn = interpolate_funcs(X_sample, Y_sample, fns);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 5 fns: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "fns_ys" SIZE_5 ".vec", Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_poly_2(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 5 polinomio n: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "poli2_ys" SIZE_5 ".vec", Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_lagram(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 5 lagram: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "lagr_ys" SIZE_5 ".vec", Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_newton(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 5 newton: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "newt_ys" SIZE_5 ".vec", Y_graph);
      }
    }
  }
  /* Sample 100 */ {
    mymtx::vector X_sample(100);
    randomize(X_sample);
    normalize(X_sample);
    mymtx::vector Y_sample = mymtx::map(X_sample, func_2());
    mymtx::vector::fwrite("vec_out/" FN2 "xs_point" SIZE_100 ".vec", X_sample);
    mymtx::vector::fwrite("vec_out/" FN2 "ys_point" SIZE_100 ".vec", Y_sample);
    /* Sample END */
    mymtx::vector::fwrite("vec_out/"
                          "xs.vec",
                          X_graph);
    /* fn_1 n = 5*/ {
      {
        const FunctionWrapper &fn = interpolate_line(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 100 line: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "line_ys" SIZE_100 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_poly(X_sample, Y_sample, 4);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 100 polinomio: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "poly1_ys" SIZE_100 ".vec",
                              Y_graph);
      }
      {
        std::vector<FunctionWrapper *> fns;
        cosine fn1;
        sine fn2;
        ln fn3;
        tow fn4;
        id fn5;
        fns.push_back(&fn1);
        fns.push_back(&fn2);
        // fns.push_back(&fn3);fns.push_back(&fn4);
        fns.push_back(&fn5);
        const FunctionWrapper &fn = interpolate_funcs(X_sample, Y_sample, fns);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 100 fns: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "fns_ys" SIZE_100 ".vec", Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_poly_2(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 100 polinomio n: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "poli2_ys" SIZE_100 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_lagram(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 100 lagram: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "lagr_ys" SIZE_100 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_newton(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 100 newton: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "newt_ys" SIZE_100 ".vec",
                              Y_graph);
      }
    }
  }
  /* Sample 1000 */ {
    mymtx::vector X_sample(1000);
    randomize(X_sample);
    normalize(X_sample);
    mymtx::vector Y_sample = mymtx::map(X_sample, func_2());
    mymtx::vector::fwrite("vec_out/" FN2 "xs_point" SIZE_1000 ".vec", X_sample);
    mymtx::vector::fwrite("vec_out/" FN2 "ys_point" SIZE_1000 ".vec", Y_sample);
    /* Sample END */
    mymtx::vector::fwrite("vec_out/"
                          "xs.vec",
                          X_graph);
    /* fn_1 n = 5*/ {
      {
        const FunctionWrapper &fn = interpolate_line(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 1000 line: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "line_ys" SIZE_1000 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_poly(X_sample, Y_sample, 4);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 1000 polinomio: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "poly1_ys" SIZE_1000 ".vec",
                              Y_graph);
      }
      {
        std::vector<FunctionWrapper *> fns;
        cosine fn1;
        sine fn2;
        ln fn3;
        tow fn4;
        id fn5;
        fns.push_back(&fn1);
        fns.push_back(&fn2);
        // fns.push_back(&fn3);fns.push_back(&fn4);
        fns.push_back(&fn5);
        const FunctionWrapper &fn = interpolate_funcs(X_sample, Y_sample, fns);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 1000 fns: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "fns_ys" SIZE_1000 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_poly_2(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 1000 polinomio n: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "poli2_ys" SIZE_1000 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_lagram(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 1000 lagram: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "lagr_ys" SIZE_1000 ".vec",
                              Y_graph);
      }
      {
        const FunctionWrapper &fn = interpolate_newton(X_sample, Y_sample);
        mymtx::vector Y_graph = mymtx::map(X_graph, fn);
        printf("fn2-Error 1000 newton: %e\n", (Y_graph - rys).distance());
        mymtx::vector::fwrite("vec_out/" FN2 "newt_ys" SIZE_1000 ".vec",
                              Y_graph);
      }
    }
  }
  return 0;
}
