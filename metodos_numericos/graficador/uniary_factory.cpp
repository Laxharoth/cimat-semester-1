#include "function_parser.hpp"

namespace function_parser{
function_operator_strategy_identity::function_operator_strategy_identity(const double arg):arg(arg),ptr_arg(nullptr){};
function_operator_strategy_identity::function_operator_strategy_identity(const double *ptr_arg):arg(0),ptr_arg(ptr_arg){};
double function_operator_strategy_identity::calc() const{
    if(ptr_arg == nullptr)return arg;
    return *ptr_arg;
}
function_operator_strategy_uniary::function_operator_strategy_uniary(function_operator_strategy *arg):arg(arg){}
function_operator_strategy_uniary::~function_operator_strategy_uniary(){
    delete arg;
}

class sin_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return std::sin(arg->calc());}
};
class cos_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return std::cos(arg->calc());}
};
class tan_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return std::tan(arg->calc());}
};
class csc_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return 1/std::sin(arg->calc());}
};
class sec_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return 1 / std::cos(arg->calc());}
};
class cot_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return 1 / std::tan(arg->calc());}
};
class sinh_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return std::sinh(arg->calc());}
};
class cosh_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return std::cosh(arg->calc());}
};
class tanh_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return std::tanh(arg->calc());}
};
class csch_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return 1/ std::sinh(arg->calc());}
};
class sech_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return 1/ std::cosh(arg->calc());}
};
class coth_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return 1/ std::tanh(arg->calc());}
};
class ln_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return std::log(arg->calc());}
};
class log_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return std::log10(arg->calc());}
};
class sqrt_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return std::sqrt(arg->calc());}
};
class cbrt_strategy : public function_operator_strategy_uniary{
    using function_operator_strategy_uniary::function_operator_strategy_uniary;    
    public:
    double calc() const { return std::cbrt(arg->calc());}
};

function_operator_strategy *uniary_strategy_factory(function_operator_strategy *arg, const std::string &operator_keyword){
    if(operator_keyword == "sin") return new sin_strategy(arg);
    if(operator_keyword == "cos") return new cos_strategy(arg);
    if(operator_keyword == "tan") return new tan_strategy(arg);
    if(operator_keyword == "csc") return new csc_strategy(arg);
    if(operator_keyword == "sec") return new sec_strategy(arg);
    if(operator_keyword == "coth") return new cot_strategy(arg);
    if(operator_keyword == "sinh") return new sinh_strategy(arg);
    if(operator_keyword == "cosh") return new cosh_strategy(arg);
    if(operator_keyword == "tanh") return new tanh_strategy(arg);
    if(operator_keyword == "csch") return new csch_strategy(arg);
    if(operator_keyword == "sech") return new sech_strategy(arg);
    if(operator_keyword == "coth") return new coth_strategy(arg);
    if(operator_keyword == "ln") return new ln_strategy(arg);
    if(operator_keyword == "log") return new log_strategy(arg);
    if(operator_keyword == "sqrt") return new sqrt_strategy(arg);
    if(operator_keyword == "cbrt") return new cbrt_strategy(arg);
    throw BAD_OPERATOR_ERROR;
}
}
