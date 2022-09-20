#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#if defined(_WIN32) || defined(WIN32)
#define CLEAR system("cls");
#else
#define CLEAR system("clear");
#endif

#define handle_statement(what) what;\
if (errno == EDOM) {\
    printf("ERROR: %s\n",strerror(errno));\
    exit(errno);\
}

struct id2481f504625843a1bf7b14bfc94a4025 {
char calif_promedio[2]; //Formato USA: {A+ > A > A- > B+ > B > B- > ...> F }
char grupo; // A, B, C, D
char turno; // M o V
};
typedef struct id2481f504625843a1bf7b14bfc94a4025 Info;
struct ide97e430423144c33aeb2b09da89864fb {
char *nombre;
int edad;
Info *E;
};
typedef struct ide97e430423144c33aeb2b09da89864fb Estudiante;

Estudiante create_student(char *nombre, int edad, char *calif_promedio, char grupo, char turno);
void delete_estudiante(Estudiante *student);
void print_student(Estudiante *student);
void print_students(Estudiante *student,int count);
void print_by_cmp(Estudiante *student,int count, int (*cmp)(void *a, void *b));
void print_count_group(Estudiante *student,int count);
void print_count_turn( Estudiante *student,int count);
void swap_stnt(void *a, void *b);
int grade_to_int(char grade[2]);
int cmp_name(void *a, void *b);
int cmp_age(void *a, void *b);
int cmp_grade(void *a, void *b);
void save_by_name(Estudiante *student,int count);
void remove_estudiante(Estudiante **student,int *count);
void add_estudiante(Estudiante **student,int *count, int *max_size);
void flush_buffer();