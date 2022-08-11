#include "function_parser.hpp"
using std::string;
using std::vector;

namespace function_parser {
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
