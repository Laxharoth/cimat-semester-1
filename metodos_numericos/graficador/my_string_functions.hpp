#ifndef MY_STRING_FUNCTIONS_HPP
#define MY_STRING_FUNCTIONS_HPP

#include <vector>
#include <string>
namespace my_string_functions{
std::string trim(std::string str);
std::vector<std::string> split(std::string str, std::string substr);
bool contains(const std::vector<std::string>& tokens, const std::string& str);
}
#endif /* MY_STRING_FUNCTIONS_HPP */
