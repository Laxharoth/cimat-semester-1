#include <sys/types.h>
#include <sys/wait.h>
#include <sys/shm.h>
#include <unistd.h> 
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include "get_time.h"

#ifndef Cores
#define CORES 4
#endif // CORES
#define print_error() printf("ERROR: %s",strerror(errno));
#define END_OK 0
#define ERROR_NOT_SUPPORT   -1
#define ERROR_SHMGET       -11
#define ERROR_SHMAT        -12
#define ERROR_SHMCTL       -13
#define ERROR_FORK         -21
#define ERROR_WAIT         -22
#define ERROR_NOT_IN_ARR   -23

#define REPETITIONS 20

struct _shape{
    unsigned shape_x;unsigned shape_y;
}; typedef struct _shape shape;
double *omp_mutl(double *A, double *B, shape shape_A, shape shape_B, shape *shape_C);
double *fork_mutl(double *A, double *B, shape shape_A, shape shape_B, shape *shape_C);

int main(){
    double *A=malloc(sizeof(double)*1000*1000);
    shape shape_A={1000,1000};
    double *B=malloc(sizeof(double)*1000*1000);
    shape shape_B={1000,1000};
    double *C; shape shape_C;
    // omp multiplication
    long time=0;
    for (size_t i = 0; i < REPETITIONS; i++)
    {
        time += GetTimeStamp(
            C = omp_mutl(A,B,shape_A,shape_B,&shape_C);
        );
    }
    printf("Open MP time: %ld microseconds\n",time/REPETITIONS);
    free(C);
    time=0;
    for (size_t i = 0; i < REPETITIONS; i++)
    {
        time += GetTimeStamp(
            C = fork_mutl(A,B,shape_A,shape_B,&shape_C);
        );
    }
    printf("Open MP time: %ld microseconds\n",time/REPETITIONS);
    
    free(C);
}

double *omp_mutl(double *A, double *B, shape shape_A, shape shape_B, shape *shape_C){
    shape_C->shape_y = shape_A.shape_y;
    shape_C->shape_x = shape_B.shape_x;
    double *C = calloc(shape_C->shape_y*shape_C->shape_x, sizeof(double));
    #pragma omp parallel for
    for (size_t i = 0; i < shape_A.shape_y; i++){
        for( size_t k = 0; k < shape_B.shape_y; ++k){
            for(size_t j = 0; j < shape_B.shape_x; ++j){
                C[i*shape_C->shape_x + j] += A[i*shape_A.shape_x+k] * 
                                             B[k*shape_B.shape_x+j];
            }
        }
    }
    return C;
}
double *fork_mutl(double *A, double *B, shape shape_A, shape shape_B, shape *shape_C){
    #define RELEASE() if( shmctl(shmid,IPC_RMID,0) == -1){ print_error(); exit(ERROR_SHMCTL); }
    shape_C->shape_y = shape_A.shape_y;
    shape_C->shape_x = shape_B.shape_x;
    key_t key = 5678;
    int shmid, pid=0;
    double *shtC;
    if ((shmid = shmget(key, 
            sizeof(double)*shape_C->shape_y*shape_C->shape_x,
            IPC_CREAT | 0666)) < 0) {
        RELEASE();
        print_error();
        exit(ERROR_SHMGET);
    }
    if ((shtC = shmat(shmid, NULL, 0)) == (void*)-1) {
        RELEASE();
        print_error();
        exit(ERROR_SHMAT);
    }
    memset(shtC, 0, sizeof(double)*shape_C->shape_y*shape_C->shape_x);
    //create children
    unsigned int child_n;
    for (child_n = 0; child_n < CORES; child_n++){
        if((pid=fork()) == -1){
            RELEASE();
            print_error();
            exit(ERROR_FORK);
        }
        //prevent child of forking
        if(pid == 0){ usleep(100);break; }
    }
    if (pid == 0){
        int start = (int)(shape_A.shape_y * (((double)child_n)/ CORES));
        int end   = (int)(shape_A.shape_y * (((double)(child_n+1))/ CORES));
        for (size_t i = start; i < end; i++){
        for( size_t k = 0; k < shape_B.shape_y; ++k){
            for(size_t j = 0; j < shape_B.shape_x; ++j){
                shtC[i*shape_C->shape_x + j] += A[i*shape_A.shape_x+k] * 
                                                B[k*shape_B.shape_x+j];
            }
        }}
        exit(END_OK);
    }
    int status; pid_t wpid;
    for (size_t i = 0; i < CORES; i++){
        if((wpid = wait(&status)) == (pid_t)-1){
            RELEASE();
            print_error();
            exit(ERROR_WAIT);
        }else if(WEXITSTATUS(status)){
            printf("Exited status %d\n", WEXITSTATUS(status));
        }
    }
    double *C = malloc(sizeof(double)*shape_C->shape_y*shape_C->shape_x);
    memcpy(C, shtC, sizeof(double)*shape_C->shape_y*shape_C->shape_x);
    RELEASE();
    return C;
}
