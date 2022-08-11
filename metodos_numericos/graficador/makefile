gkt_include = `pkg-config --cflags --libs gtk+-3.0`
object = `simple_function_analizer.o main.o canvas_facade.o function_parser.o my_string_function.o binary_factory.o uniary_factory.o` 
programa: $(object)
	g++ -o programa.run $(object) $(gkt_include)
	./programa.run
simple_function_analizer.o: simple_function_analizer.cpp point.hpp simple_function_analizer.hpp
	g++ -c simple_function_analizer.cpp
main.o: main.cpp point.hpp simple_function_analizer.hpp
	g++ $(gkt_include) -c main.cpp
canvas_facade.o: canvas_cairo_facade.cpp canvas_facade.hpp point.hpp
	g++ $(gkt_include) -c canvas_cairo_facade.cpp -o canvas_facade.o
function_parser.o: function_parser.cpp function_parser.hpp my_string_function.hpp uniary_factory.cpp binary_factory.cpp
	g++ -c function_parser.cpp uniary_factory.cpp binary_factory.cpp
binary_factory.o : function_parser.hpp binary_factory.cpp:
	g++ -c binary_factory.cpp
uniary_factory.o : function_parser.hpp uniary_factory.cpp:
	g++ -c uniary_factory.cpp
my_string_function.o: my_string_function.cpp my_string_function.hpp
	g++ -c my_string_function.cpp
one_go: simple_function_analizer.o main.o canvas_facade.o
	g++ $(gkt_include) -o programa main.cpp simple_function_analizer.cpp canvas_facade.cpp $(gkt_include)

clean:
	rm *.o