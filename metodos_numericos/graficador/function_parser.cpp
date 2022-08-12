#include "function_parser.hpp"
#include <iostream>
using namespace std;

using std::string;
using std::vector;
//identity, uniary  = 0
// +, - = 1
// *, / = 2
// ^ = 3
// () = 4
namespace function_parser {
function_operator_strategy* get_single_argument(vector<string>::iterator *ptr_current, vector<string>::iterator end, const double *variable);

function_operator_strategy* parse_commands(vector<string>::iterator *ptr_current, vector<string>::iterator end, const double *variable){
    auto &current = *ptr_current;
    char function_priority = 0;
    function_operator_strategy* parsed_function = nullptr;
    function_operator_strategy *arg1 = get_single_argument(ptr_current, end, variable);
    if(current == end || *current ==")"){ return arg1;}
    std::string operator_keyword_binary = *current;
    ++current;
    do{
        if(operator_keyword_binary == "^"){
            arg1 = binary_strategy_factory(arg1, get_single_argument(ptr_current, end, variable), "^");
            if(current == end || *current ==")")return arg1;
        }
        if(operator_keyword_binary == "*" || operator_keyword_binary == "/"){
            function_operator_strategy* arg2 = get_single_argument(ptr_current, end, variable);
            if(current == end || *current ==")") return binary_strategy_factory(arg1, arg2, "*");
            operator_keyword_binary = *current;
            ++current;
            if(operator_keyword_binary == "^"){
                arg2 = binary_strategy_factory(arg2, get_single_argument(ptr_current, end, variable), "^");
                operator_keyword_binary = *current;
            }
            arg1 = binary_strategy_factory(arg1, arg2, "*");
        }
    }while(current != end && *current !=")" && operator_keyword_binary != "+" && operator_keyword_binary != "-");
    if(current == end || *current ==")"){ return arg1;}
    return  binary_strategy_factory(arg1, parse_commands(ptr_current, end, variable), operator_keyword_binary);
}
function_operator_strategy* get_single_argument(vector<string>::iterator *ptr_current, vector<string>::iterator end, const double *variable){
    auto &current = *ptr_current;
    bool is_uniary = false;
    std::string operator_keyword = "";
    function_operator_strategy* arg1 = nullptr;
    if(*current == "("){
        ++current;
        arg1 = parse_commands(ptr_current, end, variable);
        ++current;
        return arg1;
    }
    if( my_string_functions::contains( unitary_operators, *current ) ){
        is_uniary = true;
        operator_keyword = *current;
        ++current;
    }
    if(is_uniary) arg1 = uniary_strategy_factory( get_single_argument(ptr_current, end, variable), operator_keyword );
    else if(*current == "x")arg1 = new function_operator_strategy_identity(variable);
    else arg1= new function_operator_strategy_identity(std::stod(*current));
    if(*current != ")")++current;
    return arg1;
}

FunctionParser::FunctionParser():func(nullptr){ this->variable = new double; }
FunctionParser::FunctionParser(const std::string &func_rep):func(nullptr){
    this->variable = new double;
    this->parse(func_rep);
}
void FunctionParser::parse(const std::string &func_rep){
    if(this->func != nullptr) delete this->func;
    auto v = my_string_functions::split(func_rep, " ");
    auto begin = v.begin();
    func = parse_commands(&begin, v.end(), this->variable);
}
FunctionParser::~FunctionParser(){
    delete variable; if(func != nullptr) delete func;
}
double FunctionParser::calc(double x){
    *(this->variable) = x;
    this->func->calc();
}
}
