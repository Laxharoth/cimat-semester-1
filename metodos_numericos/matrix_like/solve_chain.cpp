#include "solve_chain.hpp"
#include "funcion_matriz.hpp"
#include "matrix.hpp"

using namespace mymtx;

solve_chain::solve_chain(std::vector<solve_handler *> l) {
  for (auto &&i : l) {
    this->handlers.push_back(i);
  }
}
solve_chain::~solve_chain() {
  for (auto &&i : handlers) {
    delete i;
  }
}
solve_handler ::~solve_handler() {}
void solve_chain::solve(const matrix &A, vector &x, const vector &b) {
  return solve(A, x, b, false);
}
void solve_chain::solve(const matrix &A, vector &x, const vector &b,
                        bool is_symetric) {
  if (memo_solved.find(&A) != memo_solved.end()) {
    return handlers[memo_solved.find(&A)->second]->solve(A, x, b);
  }
  for (int i = 0; i < handlers.size(); ++i) {
    try {
      handlers[i]->set_symetric(is_symetric);
      handlers[i]->solve(A, x, b);
      memo_solved.emplace(&A, i);
      return;
    } catch (const std::exception &e) {
    }
  }
  throw cant_solve_exception();
}
void solve_handler::clear_memo() { memo.clear(); }
void solve_handler::set_symetric(bool) {}
void gauss_handler::solve(const matrix &A, vector &x, const vector &b) {
  auto cpyA = A;
  auto cpyb = b;
  gauss(cpyA, x, cpyb);
}
void crout_handler::solve(const matrix &A, vector &x, const vector &b) {
  if (memo.find(&A) == memo.end()) {
    matrix mem(A.shape_y, A.shape_y);
    crout(A, mem, mem);
    memo.insert(std::make_pair(&A, mem));
  }
  solucion_crout(memo.find(&A)->second, x, b);
}

void qr_handler::solve(const matrix &A, vector &x, const vector &b) {
  if (memo.find(&A) == memo.end()) {
    matrix Q(A.shape_y, A.shape_y);
    matrix R(A.shape_y, A.shape_y);
    qr_decomposition(A, Q, R);
    memo.insert(std::make_pair(&A, Q));
    R_memo.insert(std::make_pair(&A, R));
  }
  solve_qr(memo.at(&A), R_memo.at(&A), x, b);
}
void cholesky_handler::set_symetric(bool b) { is_symetric = b; }
void cholesky_handler::solve(const matrix &A, vector &x, const vector &b) {
  if (memo.find(&A) == memo.end()) {
    if (!is_symetric) {
      throw cant_factor_exception();
    }
    matrix mem(A.shape_y, A.shape_y);
    factor_cholesky(A, mem);
    memo.insert(std::make_pair(&A, mem));
  }
  solve_cholesky(memo.at(&A), x, b);
}
solve_chain
solve_chain::build_chain(const std::vector<const char *> &link_names) {
  std::vector<solve_handler *> handlers;
  for (auto &&link : link_names) {
    if (strcmp(link, solve_chain::GAUSS) == 0) {
      handlers.push_back(new gauss_handler());
    }
    if (strcmp(link, solve_chain::QR) == 0) {
      handlers.push_back(new qr_handler());
    }
    if (strcmp(link, solve_chain::CHOLESKY) == 0) {
      handlers.push_back(new cholesky_handler());
    }
    if (strcmp(link, solve_chain::CROUT) == 0) {
      handlers.push_back(new crout_handler());
    }
  }
  return solve_chain(handlers);
}
const char *cant_solve_exception::what() const noexcept {
  return "añsdjfsdlfa";
}
const char *solve_chain::GAUSS = "382c421e-55d0-4604-82aa-b32f193c5e2e";
const char *solve_chain::CROUT = "07746c48-f5ed-4d33-81a4-247073b645cc";
const char *solve_chain::QR = "72b55ac3-edee-4471-8e4f-2807e678c4dd";
const char *solve_chain::CHOLESKY = "b350d165-7295-42f2-8513-0b24539adcc7";
