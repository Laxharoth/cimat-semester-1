#ifndef PARSER_WRAPPER_CPP
#define PARSER_WRAPPER_CPP
#include "parser_wrapper.hpp"
#include "function_wrapper.cpp"
FunctionParserAdapter_to_FWrapper::FunctionParserAdapter_to_FWrapper(FunctionParser *parser):parser(parser) {}
double FunctionParserAdapter_to_FWrapper::eval(const double &x){ return parser->Eval(&x); }

#endif /* PARSER_WRAPPER_CPP */
