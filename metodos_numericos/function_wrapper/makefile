function_wrapper.o: function_wrapper.hpp function_wrapper.cpp parser_wrapper.cpp parser_wrapper.hpp
	g++ -O2 -c parser_wrapper.cpp -o function_wrapper.o

test: test.run
	./test.run

test.run: test.cpp function_wrapper.hpp function_wrapper.cpp
	g++ -g -o test.run test.cpp function_wrapper.cpp
