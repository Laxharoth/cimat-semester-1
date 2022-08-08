programa: main.o simple_function_analizer.o
	g++ -o programa main.o simple_function_analizer.o

simple_function_analizer.o: simple_function_analizer.cpp point.hpp simple_function_analizer.hpp
	g++ -c simple_function_analizer.cpp

main.o: main.cpp point.hpp simple_function_analizer.hpp
	g++ -c main.cpp