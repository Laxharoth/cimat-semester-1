gkt_include = `pkg-config --cflags --libs gtk+-3.0`
programa: simple_function_analizer.o main.o canvas_facade.o
	g++ -o programa.run main.o simple_function_analizer.o canvas_facade.o $(gkt_include)
	./programa.run

simple_function_analizer.o: simple_function_analizer.cpp point.hpp simple_function_analizer.hpp
	g++ -c simple_function_analizer.cpp

main.o: main.cpp point.hpp simple_function_analizer.hpp
	g++ $(gkt_include) -c main.cpp

canvas_facade.o: canvas_facade.cpp canvas_facade.hpp point.hpp
	g++ $(gkt_include) -c canvas_facade.cpp

one_go: simple_function_analizer.o main.o canvas_facade.o
	g++ $(gkt_include) -o programa main.cpp simple_function_analizer.cpp canvas_facade.cpp $(gkt_include)

clean:
	rm *.o