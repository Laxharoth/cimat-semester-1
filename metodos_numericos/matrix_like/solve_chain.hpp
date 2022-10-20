#include "funcion_matriz.hpp"
#include "matrix.hpp"
#include <cstring>
#include <map>
#include <vector>

namespace mymtx {
class solve_handler;
class cant_solve_exception : public std::exception {
  const char *what() const noexcept;
};
class solve_chain {
  std::vector<solve_handler *> handlers;
  std::map<const matrix *, signed int> memo_solved;
  solve_chain(std::vector<solve_handler *>);

public:
  static const char *GAUSS;
  static const char *CROUT;
  static const char *QR;
  static const char *CHOLESKY;
  void solve(const matrix &A, vector &x, const vector &b);
  void solve(const matrix &A, vector &x, const vector &b, bool is_symetric);
  static solve_chain build_chain(const std::vector<const char *> &);
  ~solve_chain();
};
class solve_handler {
protected:
  std::map<const matrix *, matrix> memo;

public:
  virtual void solve(const matrix &A, vector &x, const vector &b) = 0;
  virtual void set_symetric(bool);
  virtual ~solve_handler();
  void clear_memo();
};
class gauss_handler : public solve_handler {
  void solve(const matrix &A, vector &x, const vector &b);
};
class crout_handler : public solve_handler {
  void solve(const matrix &A, vector &x, const vector &b);
};
class qr_handler : public solve_handler {
  std::map<const matrix *, matrix> R_memo;
  void solve(const matrix &A, vector &x, const vector &b);
};
class cholesky_handler : public solve_handler {
  bool is_symetric = false;
  void solve(const matrix &A, vector &x, const vector &b) override;
  void set_symetric(bool) override;
};
} // namespace mymtx