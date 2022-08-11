gkt_include = `pkg-config --cflags --libs gtk+-3.0`
programa: simple_function_analizer.o main.o canvas_facade.o function_parser.o my_string_function.o
	g++ -o programa.run main.o simple_function_analizer.o canvas_facade.o $(gkt_include)
	./programa.run

simple_function_analizer.o: simple_function_analizer.cpp point.hpp simple_function_analizer.hpp
	g++ -c simple_function_analizer.cpp

main.o: main.cpp point.hpp simple_function_analizer.hpp
	g++ $(gkt_include) -c main.cpp

canvas_facade.o: canvas_cairo_facade.cpp canvas_facade.hpp point.hpp
	g++ $(gkt_include) -c canvas_cairo_facade.cpp -o canvas_facade.o
function_parser.o: function_parser.cpp function_parser.hpp my_string_function.hpp
	g++ -c function_parser.cpp
my_string_function.o: my_string_function.cpp my_string_function.hpp
	g++ -c my_string_function.cpp
one_go: simple_function_analizer.o main.o canvas_facade.o
	g++ $(gkt_include) -o programa main.cpp simple_function_analizer.cpp canvas_facade.cpp $(gkt_include)

clean:
	rm *.o