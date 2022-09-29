#include "huffman.h"

#include <stdlib.h>
#include <stdio.h>

long fileSize(const char *file_name);
int main(int argc, char const *argv[]){
    {
        /*  compress */{
            const char *filename = "fractal_tree.ascii.pgm";
            long file_size = fileSize(filename);
            byte *buffer = malloc(file_size);
            FILE *f = fopen(filename,"rb");
            fread(buffer,1,file_size,f);
            fclose(f);
            compress(buffer,file_size,"fractal_tree.ascii.pgm" ".compresed");
            free(buffer);
        }
        /*  compress */
        {
            const char *filename = "fractal_tree.ascii.pgm" ".compresed";
            long file_size = fileSize(filename);
            byte *buffer = malloc(file_size);
            FILE *f = fopen(filename,"rb");
            fread(buffer,1,file_size,f);
            fclose(f);
            decompress(buffer,"decompresed." "fractal_tree.ascii.pgm");
            free(buffer);
        }
    }
    return 0;
}
long fileSize(const char *file_name){
    FILE* fp = fopen(file_name, "r");
    if (fp == NULL) {
        printf("File Not Found!\n");
        exit(-1);
    }
    fseek(fp, 0L, SEEK_END);
    long res = ftell(fp);
    fclose(fp);
    return res;
}
