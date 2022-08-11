#include "function_parser.hpp"
using std::string;
using std::vector;

namespace function_parser {
function_operator_strategy_identity::function_operator_strategy_identity(const double arg):arg(arg),ptr_arg(nullptr){};
function_operator_strategy_identity::function_operator_strategy_identity(const double *ptr_arg):arg(0),ptr_arg(ptr_arg){};
double function_operator_strategy_identity::calc() const{
    if(ptr_arg == nullptr)return arg;
    return *ptr_arg;
}
function_operator_strategy_uniary::function_operator_strategy_uniary(
    function_operator_strategy *arg, const std::string &operator_keyword):arg(arg),operator_keyword(operator_keyword){}
function_operator_strategy_uniary::~function_operator_strategy_uniary(){
    delete arg;
}
double function_operator_strategy_uniary::calc() const{
    if(operator_keyword == "sin") return std::sin(arg->calc());
    if(operator_keyword == "cos") return std::cos(arg->calc());
    if(operator_keyword == "tan") return std::tan(arg->calc());
    if(operator_keyword == "csc") return 1 / std::sin(arg->calc());
    if(operator_keyword == "sec") return 1 / std::cos(arg->calc());
    if(operator_keyword == "coth") return 1 / std::tan(arg->calc());
    if(operator_keyword == "sinh") return std::sinh(arg->calc());
    if(operator_keyword == "cosh") return std::cosh(arg->calc());
    if(operator_keyword == "tanh") return std::tanh(arg->calc());
    if(operator_keyword == "csch") return 1 / std::sinh(arg->calc());
    if(operator_keyword == "sech") return 1 / std::cosh(arg->calc());
    if(operator_keyword == "coth") return 1 / std::tanh(arg->calc());
    if(operator_keyword == "ln") return std::log10(arg->calc());
    if(operator_keyword == "log") return std::log(arg->calc());
    if(operator_keyword == "sqrt") return std::sqrt(arg->calc());
    if(operator_keyword == "cbrt") return std::cbrt(arg->calc());
    throw BAD_OPERATOR_ERROR;
}
function_operator_strategy_binary::function_operator_strategy_binary(
    function_operator_strategy *arg1, function_operator_strategy *arg2,const std::string &operator_keyword)
    :arg1(arg1),arg2(arg2),operator_keyword(operator_keyword){}
function_operator_strategy_binary::~function_operator_strategy_binary(){
    delete arg1; delete arg2;
}
double function_operator_strategy_binary::calc() const {
    if(operator_keyword == "+") return arg1->calc() + arg2->calc();
    if(operator_keyword == "-") return arg1->calc() - arg2->calc();
    if(operator_keyword == "*") return arg1->calc() * arg2->calc();
    if(operator_keyword == "/") return arg1->calc() / arg2->calc();
    if(operator_keyword == "^") return  std::pow(arg1->calc(), arg2->calc());
    throw BAD_OPERATOR_ERROR;
}
FunctionParser::FunctionParser(const std::string &func_rep){
    this->variable = new double;
    auto v = my_string_functions::split(func_rep, " ");
    auto begin = v.begin();
    func = parse_commands(&begin, v.end(), this->variable);
}
FunctionParser::~FunctionParser(){
    delete variable; delete func;
}
double FunctionParser::calc(double x){
    *(this->variable) = x;
    this->func->calc();
}
function_operator_strategy* parse_commands(vector<string>::iterator *ptr_current, vector<string>::iterator end, const double *variable){
    auto &current = *ptr_current;
    std::string operator_keyword = "";
    bool is_unitary = false;
    function_operator_strategy *arg1;
    if(*current == "("){
        ++current;
        arg1 = parse_commands(ptr_current, end, variable);
        ++current;
        if(current == end || *current ==")") return arg1;
        std::string operator_keyword_binary = "";
        operator_keyword_binary = *current;
        ++current;
        return new function_operator_strategy_binary(
            arg1, 
            parse_commands(ptr_current, end, variable),
            operator_keyword_binary
        );
    }
    if( my_string_functions::contains( unitary_operators, *current ) ){
        is_unitary = true;
        operator_keyword = *current;
        ++current;
    }
    if(is_unitary)arg1 = new function_operator_strategy_uniary{ parse_commands(ptr_current, end, variable), operator_keyword };
    else if(*current == "x")arg1 = new function_operator_strategy_identity(variable);
    else arg1= new function_operator_strategy_identity(std::stod(*current));
    ++current;
    if(current == end || *current ==")") return arg1;
    std::string operator_keyword_binary = "";
    operator_keyword_binary = *current;
    ++current;
    auto t =  new function_operator_strategy_binary(
        arg1, 
        parse_commands(ptr_current, end, variable),
        operator_keyword_binary
    );
    return t;
}
}
