#ifndef PARSER_WRAPPER_HPP
#define PARSER_WRAPPER_HPP

#include "function_wrapper.hpp"
#include "../fparser/fparser.hh"

class FunctionParserAdapter_to_FWrapper : public FunctionWrapper{
    FunctionParser *parser;
    public:
    FunctionParserAdapter_to_FWrapper(FunctionParser *parser);
    double eval(const double &x) const;
    double eval(const double &x);
};

#endif /* PARSER_WRAPPER_CPP */
