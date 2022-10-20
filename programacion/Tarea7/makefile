#Para obtener la tabla de frecuencias agregar la flag -DPRINT_TABLE a la compilacion
huffman.run: main_huffman.c huffman.h huffman.c quicksort.c quicksort.h tree.c tree.h queue.h queue.c
	gcc -g main_huffman.c huffman.c tree.c queue.c quicksort.c -o huffman.run -lm
img_multiprocess.run: img.multiprocess_main.c quicksort.c median.c pgm1.c
	gcc -g -DWINDOW_SIZE=7 -DFORK_SIZE=4 img.multiprocess_main.c quicksort.c -o img_multiprocess.run
img_multiprocess_time.run: img.multiprocess_main.c quicksort.c median.c pgm1.c get_time.h
	gcc -g -DMEASURE_TIME img.multiprocess_main.c quicksort.c -o img_multiprocess_time.run -lm
mtx.run: matrix_mult.c get_time.h
	gcc matrix_mult.c -o mtx.run -fopenmp
clean:
	rm -f *.run
	rm -f *.out
