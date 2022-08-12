#include "my_string_functions.hpp"

namespace my_string_functions{
std::string trim(std::string str){
    const char* typeOfWhitespaces = " \t\n\r\f\v";
    str.erase(str.find_last_not_of(typeOfWhitespaces) + 1);
    str.erase(0,str.find_first_not_of(typeOfWhitespaces));
    return str;
}
std::vector<std::string> split(std::string str, std::string substr){
    int index=0;
    int end_index{};
    std::vector<std::string> tokens;
    while( ( end_index = str.find(substr, index) ) >= 0 ){
        tokens.push_back(str.substr(index, end_index - index));
        index = end_index + substr.length();
    }
    tokens.push_back(str.substr(index));
    return tokens;
}

bool contains(const std::vector<std::string>& tokens, const std::string& str){
    for(auto i = tokens.cbegin(); i != tokens.cend(); ++i){
        if(str == *i)return true;
    }
    return false;
}
std::string rm_double_space(std::string str){
    size_t double_space = str.find("  ");
    while( double_space != std::string::npos ){
        str.erase(double_space, 1);
        double_space = str.find("  ");
    }
    return str;
}
}
