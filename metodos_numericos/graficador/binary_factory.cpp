#include "function_parser.hpp"

namespace function_parser{
function_operator_strategy_binary::function_operator_strategy_binary(
function_operator_strategy *arg1, function_operator_strategy *arg2):arg1(arg1),arg2(arg2){}
function_operator_strategy_binary::~function_operator_strategy_binary(){ delete arg1; delete arg2; }
class plus_strategy : public function_operator_strategy_binary{
    using function_operator_strategy_binary::function_operator_strategy_binary;
    public: double calc() const { return arg1->calc() + arg2->calc(); }
};
class minus_strategy : public function_operator_strategy_binary{
    using function_operator_strategy_binary::function_operator_strategy_binary;
    public: double calc() const { return arg1->calc() - arg2->calc(); }
};
class multliplication_strategy : public function_operator_strategy_binary{
    using function_operator_strategy_binary::function_operator_strategy_binary;
    public: double calc() const { return arg1->calc() * arg2->calc(); }
};
class division_strategy : public function_operator_strategy_binary{
    using function_operator_strategy_binary::function_operator_strategy_binary;
    public: double calc() const { return arg1->calc() / arg2->calc(); }
};
class power_strategy : public function_operator_strategy_binary{
    using function_operator_strategy_binary::function_operator_strategy_binary;
    public: double calc() const { return std::pow(arg1->calc(), arg2->calc()); }
};

function_operator_strategy * binary_strategy_factory(function_operator_strategy *arg1, function_operator_strategy *arg2, const std::string &operator_keyword) {
    if(operator_keyword == "+") return new plus_strategy(arg1, arg2);
    if(operator_keyword == "-") return new minus_strategy(arg1, arg2);
    if(operator_keyword == "*") return new multliplication_strategy(arg1, arg2);
    if(operator_keyword == "/") return new division_strategy(arg1, arg2);
    if(operator_keyword == "^") return new power_strategy(arg1, arg2);
    throw operator_keyword[0];
    throw BAD_OPERATOR_ERROR;
}
}

