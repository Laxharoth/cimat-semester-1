#include <stdio.h>
#include <stdlib.h>

#include "quicksort.c"

#define ABC_SIZE 26

struct letter_position_cout{
    int count;
    char position;
};
typedef struct letter_position_cout letter_position;

int cmp_letter_count( void* a, void* b);
void swap_letter_count( void* a, void* b );
void handle_my_file( const char* filename,unsigned int *filesize, unsigned int *n_lines , unsigned int *char_length, int *abcedary_count, int *abcedary_precedence_count);
void print_letter_percent( int *abcedary_count);
void print_precedende_letter_percent( int *abcedary_precedence_count, letter_position *position);

int cmp_letter_count( void* a, void* b){
    letter_position* aa = (letter_position*)a;
    letter_position* bb = (letter_position*)b;
    return aa->count - bb->count;
}
void swap_letter_count( void* a, void* b ){
    letter_position c;
    letter_position* aa = (letter_position*)a;
    letter_position* bb = (letter_position*)b;
    c.count = aa->count;
    c.position = aa->position;
    aa->count = bb->count;
    aa->position = bb->position;
    bb->count = c.count;
    bb->position = c.position;
}

int main(int argc, char const *argv[]){
    const char* filename = "Don-Quijote-Ingles.txt";
    unsigned int filesize = 0;
    unsigned int char_length = 0;
    unsigned int n_lines;
    int abcedary_count[ABC_SIZE] = {0};
    int abcedary_precedence_count[ABC_SIZE*ABC_SIZE] = {0};
    handle_my_file( filename, &filesize, &char_length, &n_lines, abcedary_count, abcedary_precedence_count);
    printf("Filename: %s\n", filename);
    printf("filesize: %d\n", filesize);
    printf("char_length: %d\n", char_length);
    printf("lines: %d\n", n_lines);
    print_letter_percent(abcedary_count);
    letter_position positions[ABC_SIZE];
    for (unsigned int i = 0; i < ABC_SIZE; i++){
        positions[i].count = abcedary_count[i];
        positions[i].position = i;
    }
    quicksort( positions, 0, ABC_SIZE, sizeof(letter_position), cmp_letter_count, swap_letter_count );
    print_precedende_letter_percent(abcedary_precedence_count, positions);
    return 0;
}

void handle_my_file( const char* filename,unsigned int *filesize, unsigned int *n_lines, 
        unsigned int *char_length, int *abcedary_count, int *abcedary_precedence_count){
    char ch;
    unsigned int count_position;
    unsigned int prev_count_position;
    FILE *fp = fopen(filename,"r");
    *n_lines = 1;
    if (fp == NULL){
        perror("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }
    while((ch = fgetc(fp)) != EOF){
        (*char_length)++;
        (*filesize)+=sizeof(char);
        if( ch == '\n' ) (*n_lines)++;
        if( ch >= 'a' && ch <= 'z' ){
            prev_count_position = count_position;
            count_position = ch - 'a';
        }
        else if( ch >= 'A' && ch <= 'Z' ){
            prev_count_position = count_position;
            count_position = ch - 'A';
        }
        else{ continue; }
        abcedary_count[count_position]++;
        // printf("%d\n", abcedary_count[count_position]);
        //    # start of the bidimensional array                                            #
        (* ( (get_address(abcedary_precedence_count, count_position, ABC_SIZE)) +
        //  # offset       # 
            prev_count_position ) )++;
    }   
    fclose(fp);
}

void print_letter_percent( int *abcedary_count){
    printf("letter count:\n");
    double sum = 0;
    for (unsigned int i = 0; i < ABC_SIZE; i++){
        sum+=abcedary_count[i];
    }
    for (unsigned int i = 0; i < ABC_SIZE; i++){
        printf("\t-%c: %lf \n", 'a'+i, (double)(abcedary_count[i])/sum );
    }
}
void print_precedende_letter_percent( int *abcedary_precedence_count, letter_position *position){
    printf("most used letters by precedent letter\n");
    double sum_per_letter[ABC_SIZE];
    for (unsigned int i = ABC_SIZE-1; i >= ABC_SIZE - 10; --i){
        printf( "\t%c:\n", position[i].position + 'a' );
        for(unsigned int j = 0; j< ABC_SIZE; ++j){
            printf("\t\t-%c: %lf \n", 'a'+j, (double)(abcedary_precedence_count[position[i].position * ABC_SIZE +j])/(double)(position[i].count) );
        }
    }
}