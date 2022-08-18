#Para ejecutar el programa se compila el con el comando make
#Al terminar la compilación se ejecutara el programa
#para ejecutar el programa de forma subsecuente se puede ejecutar como ./programa.run
#Para compilar se requiere que la libreria gtk3 esté instalada
gkt_include = `pkg-config --cflags --libs gtk+-3.0`
object_files = simple_function_analizer.o main.o canvas_facade.o function_wrapper.o fparser.o
exec: programa.run
	@echo Los archivos estan actualizados. Ejecutando programa
	./programa.run
programa.run: $(object_files)
	g++ -o programa.run $(object_files) $(gkt_include)
simple_function_analizer.o: simple_function_analizer.cpp point.hpp simple_function_analizer.hpp
	g++ -c simple_function_analizer.cpp
main.o: main.cpp point.hpp simple_function_analizer.hpp
	g++ $(gkt_include) -c main.cpp
canvas_facade.o: canvas_cairo_facade.cpp canvas_facade.hpp canvas_facade.cpp point.hpp
	g++ $(gkt_include) -c canvas_cairo_facade.cpp -o canvas_facade.o
my_string_functions.o: my_string_functions.hpp my_string_functions.cpp
	g++ -c my_string_functions.cpp
function_wrapper.o: function_wrapper/function_wrapper.hpp function_wrapper/function_wrapper.cpp
	$(MAKE) -C function_wrapper
	mv function_wrapper/function_wrapper.o .
fparser.o: fparser/fparser.hh fparser/fpconfig.hh fparser/fparser.cc
	$(MAKE) -C fparser
	mv fparser/fparser.o .
clean:
	@rm *.run
	@rm *.o
	@rm *.glade~