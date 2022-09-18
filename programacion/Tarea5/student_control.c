#include "student_struct.h"
#include "quicksort.c"

int main(int argc, char const *argv[]){
    char c;
    char name_buffer[100],grade_buffer[3], group_buffer,turn_buffer;
    int age_buffer, inserted = 0, max_size = 0;
    if(argc < 2){
        printf("ERROR: cannot open file\n\n");
        return -1;
    }
    FILE *file = fopen(argv[1], "r");
    if (file == NULL) {
        printf("ERROR: cannot open file\n\n");
        return -1;
    }
    Estudiante *stds=NULL;
    int eof = 0;
    while(1){
        handle_statement(eof=fscanf(file,"%s %d %s %c %c", name_buffer, &age_buffer, grade_buffer, &group_buffer, &turn_buffer));
        if(eof==EOF) break;
        inserted++;  max_size++;
        handle_statement( stds = realloc( stds, inserted*sizeof(Estudiante)) );
        stds[inserted-1] = create_student( name_buffer, age_buffer, grade_buffer, group_buffer, turn_buffer );
    }
    fclose(file);
    int choice = 0;
    while(choice != '9'){
        CLEAR
        printf("1) Desplegar archivo\n");
        printf("2) Ordenar por nombre:\n");
        printf("3) Ordenar por Edad;\n");
        printf("4) Ordenar por Promedio\n");
        printf("5) Número de estudiantes por grupo\n");
        printf("6) Número de estudiantes por turno\n");
        printf("7) Baja estudiante (por nombre)\n");
        printf("8) Alta estudiante (por nombre)\n");
        printf("9) Salir\n");
        printf("Seleccion: \n");
        handle_statement(choice = getchar());
        flush_buffer();
        switch (choice){
        case '1': print_students(stds,inserted);break;
        case '2': print_by_cmp(stds,inserted,cmp_name);break;
        case '3': print_by_cmp(stds,inserted,cmp_age);break;
        case '4': print_by_cmp(stds,inserted,cmp_grade);break;
        case '5': print_count_group(stds,inserted);break;
        case '6': print_count_turn(stds,inserted);break;
        case '7': remove_estudiante(&stds,&inserted);break;
        case '8': add_estudiante(&stds,&inserted,&max_size);break;
        case '9': save_by_name(stds, inserted); break;
        default: printf("Selection %c (%d) incorrect",choice,choice); getchar(); break;
        }
    }
    save_by_name(stds,inserted);
    for (int i = 0; i < inserted; i++){
        delete_estudiante(&stds[i]);
    }
    handle_statement(free(stds));
    return 0;
}
Estudiante create_student(char *nombre, int edad, char *calif_promedio, char grupo, char turno){
    handle_statement(unsigned name_size = strlen(nombre));
    handle_statement(char *name = malloc(name_size+1));
    handle_statement(strcpy(name,nombre));
    handle_statement(Info *E = malloc(sizeof(Info)));
    handle_statement(memcpy(E->calif_promedio,calif_promedio,2));
    E->grupo = grupo;
    E->turno = turno;
    Estudiante student;
    student.nombre = name;
    student.edad = edad;
    student.E = E;
    return student;
}
void delete_estudiante(Estudiante *student){
    free(student->E);
    free(student->nombre);
}
void print_student(Estudiante *student){
    printf("Nombre:%s\n",student->nombre);
    printf("Edad:%d\n",student->edad);
    printf("Promedio:%c%c\n",student->E->calif_promedio[0],student->E->calif_promedio[1]);
    printf("Grupo:%c\n",student->E->grupo);
    printf("Turno:%c\n",student->E->turno);
}
void print_students(Estudiante *student,int count){
    for (int i = 0; i < count; i++) print_student(student+i);
    printf("Enter:\n");
    handle_statement(getchar());
}
void print_by_cmp(Estudiante *student,int count, int (*cmp)(void *a, void *b)){
    Estudiante studentcpy[count];
    for (size_t i = 0; i < count; i++) studentcpy[i] = student[i];
    quicksort(studentcpy, 0, count, sizeof(Estudiante),cmp,swap_stnt);
    print_students(studentcpy,count);
}
void swap_stnt(void *a, void *b){
    Estudiante *aa = (Estudiante *)a;
    Estudiante *bb = (Estudiante *)b;
    Estudiante c = *aa;
    *aa = *bb;
    *bb = c;
}
int grade_to_int(char grade[2]){
    if(grade[0]=='A'&&grade[1]=='+')return 17;
    if(grade[0]=='A'&&grade[1]==0)return 16;
    if(grade[0]=='A'&&grade[1]=='-')return 15;
    if(grade[0]=='B'&&grade[1]=='+')return 14;
    if(grade[0]=='B'&&grade[1]==0)return 13;
    if(grade[0]=='B'&&grade[1]=='-')return 12;
    if(grade[0]=='C'&&grade[1]=='+')return 11;
    if(grade[0]=='C'&&grade[1]==0)return 10;
    if(grade[0]=='C'&&grade[1]=='-')return  9;
    if(grade[0]=='D'&&grade[1]=='+')return  8;
    if(grade[0]=='D'&&grade[1]==0)return  7;
    if(grade[0]=='D'&&grade[1]=='-')return  6;
    if(grade[0]=='E'&&grade[1]=='+')return  5;
    if(grade[0]=='E'&&grade[1]==0)return  4;
    if(grade[0]=='E'&&grade[1]=='-')return  3;
    if(grade[0]=='F'&&grade[1]=='+')return  2;
    if(grade[0]=='F'&&grade[1]==0)return  1;
    if(grade[0]=='F'&&grade[1]=='-')return  0;
    return -1;
}
int cmp_name(void *a, void *b){
    handle_statement(int res = strcmp(((Estudiante*)a)->nombre , ((Estudiante*)b)->nombre));
    return res;
}
int cmp_age(void *a, void *b){
    return ((Estudiante*)a)->edad - ((Estudiante*)b)->edad;
}
int cmp_grade(void *a, void *b){
    return  grade_to_int(((Estudiante*)a)->E->calif_promedio) - 
            grade_to_int(((Estudiante*)b)->E->calif_promedio) ;
}
void save_by_name(Estudiante *student,int count){
    quicksort(student, 0, count, sizeof(Estudiante),cmp_name,swap_stnt);
    FILE *fp = fopen("alumnos.2.txt", "w");
    for (size_t i = 0; i < count; i++){
        Estudiante *st = &student[i];
        if(st->E->calif_promedio[1]==0)
        fprintf(fp,"%s %d %c %c %c\n",
            st->nombre, st->edad, st->E->calif_promedio[0], st->E->grupo, st->E->turno);
        else
        fprintf(fp,"%s %d %c%c %c %c\n",
            st->nombre, st->edad, st->E->calif_promedio[0],st->E->calif_promedio[1], st->E->grupo, st->E->turno);
    }
    fclose(fp);
}
void print_count_group(Estudiante *student,int count){
    int groups[255]={0};
    printf("Estudiantes por grupo\n");
    for(int i=0;i<count;i++){
        groups[student[i].E->grupo]++;
    }
    for (unsigned char i = 0; i < 255; i++){
        if(groups[i]==0)continue;
        printf("%c:%d\n",i,groups[i]);
    }
    printf("Enter:\n");
    handle_statement(getchar());
}
void print_count_turn( Estudiante *student,int count){
    int groups[255]={0};
    printf("Estudiantes por turno\n");
    for(int i=0;i<count;i++){
        groups[student[i].E->turno]++;
    }
    for (unsigned char i = 0; i < 255; i++){
        if(groups[i]==0)continue;
        printf("%c:%d\n",i,groups[i]);
    }
    printf("Enter:\n");
    handle_statement(getchar());
}
void remove_estudiante(Estudiante **student,int *count){
    char name[100];
    int can_remove;
    printf("Nombre:");scanf("%s",name);
    flush_buffer();
    unsigned int i;
    for (i = 0; i < *count; i++){
        handle_statement(can_remove = strcmp(name,(*student)[i].nombre));
        if(can_remove == 0) break;
    }
    if(can_remove != 0) return;
    --(*count);
    swap_stnt((*student)+i,(*student)+(*count));
    delete_estudiante((*student)+(*count));
}
void add_estudiante(Estudiante **student,int *count, int *max_size){
    printf("Ingresa estudiante (%%s %%d %%s %%c %%c):");
    char name_buffer[100],grade_buffer[3], group_buffer,turn_buffer;
    int age_buffer;
    handle_statement(scanf("%s %d %s %c %c", name_buffer, &age_buffer, grade_buffer, &group_buffer, &turn_buffer));
    flush_buffer();
    (*count)++;
    if(*count > *max_size){
        handle_statement( *student = realloc( *student, (*count)*sizeof(Estudiante)) );
        *max_size = *count;
    }
    (*student)[(*count)-1] = create_student( name_buffer, age_buffer, grade_buffer, group_buffer, turn_buffer );
}
void flush_buffer(){
    char c  = '1';
    while (c != EOF && c != '\n') handle_statement(c = getchar()) ;
}
