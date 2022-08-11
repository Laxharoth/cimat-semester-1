#ifndef FUNCTION_PARSER_HPP
#define FUNCTION_PARSER_HPP

#include "my_string_functions.hpp"

#include <string>
#include <vector>
#include <cmath>

using std::string;

namespace function_parser {
static const string BAD_OPERATOR_ERROR = "84c8f7f6-4719-4b00-b416-aa7e1d86940d";
static const std::vector<std::string> unitary_operators{ "sin", "cos", "tan", "csc", "sec", "coth", "sinh", "cosh", "tanh", "csch", "sech", "coth", "ln", "log", "sqrt", "cbrt" };

class function_operator_strategy{
public:
    virtual double calc() const = 0;
};
class function_operator_strategy_identity : public function_operator_strategy{
    const double arg;
    const double *ptr_arg;
public:
    function_operator_strategy_identity(const double arg);
    function_operator_strategy_identity(const double *arg);
    double calc() const;
};
class function_operator_strategy_uniary : public function_operator_strategy{
protected:
    function_operator_strategy *arg;
public:
    explicit function_operator_strategy_uniary(function_operator_strategy *arg);
    ~function_operator_strategy_uniary();
    virtual double calc() const = 0;
};
class function_operator_strategy_binary  : public function_operator_strategy{
protected:
    function_operator_strategy *arg1;
    function_operator_strategy *arg2;
public:
    explicit function_operator_strategy_binary(function_operator_strategy *arg1, function_operator_strategy *arg2);
    ~function_operator_strategy_binary();
    virtual double calc() const = 0;
};
class FunctionParser{
    double *variable{nullptr};
    function_operator_strategy *func;
    public:
    FunctionParser();
    FunctionParser(const std::string &func_rep);
    ~FunctionParser();
    double calc(double x);
    void parse(const std::string &func_rep);
};
function_operator_strategy * parse_commands(std::vector<std::string>::iterator *current, std::vector<std::string>::iterator end, const double *variable);

function_operator_strategy * uniary_strategy_factory(function_operator_strategy *arg, const std::string &operator_keyword);
function_operator_strategy * binary_strategy_factory(function_operator_strategy *arg1, function_operator_strategy *arg2, const std::string &operator_keyword);
}


#endif /* FUNCTION_PARSER_HPP */
