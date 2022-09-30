#include <sys/types.h>
#include <sys/wait.h>
#include <sys/shm.h>
#include <unistd.h> 
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "median.c"
#include "pgm1.c"

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#define print_error() printf("ERROR: %s",strerror(errno));

#define END_OK 0
#define ERROR_NOT_SUPPORT   -1
#define ERROR_SHMGET       -11
#define ERROR_SHMAT        -12
#define ERROR_SHMCTL       -13
#define ERROR_FORK         -21
#define ERROR_WAIT         -22
#define ERROR_NOT_IN_ARR   -23

// pgm1 se modifico para usar apuntador simple
void fork_replace(char *mtx,int dim_y, int dim_x, filter fnc, int window_size, size_t forks);
int main(int argc, char const *argv[]){
    int cols, rows;
    char *img = (char *)pgmRead("fractal_tree.ascii.pgm",&rows, &cols);
    // replace( img, rows,cols,median, 9 );
    fork_replace( img, rows,cols,median, 9, 2);
    pgmWrite("fractal_tree.2.ascii.pgm",rows,cols,(unsigned char*)img,"");
    free(img);
    return 0;
}

void fork_replace(char *mtx,int dim_y, int dim_x, filter fnc, int window_size, size_t forks){
    #define RELEASE() if( shmctl(shmid,IPC_RMID,0) == -1){ print_error(); exit(ERROR_SHMCTL); }
    if(forks == 1) return replace(mtx,dim_y,dim_x,fnc,window_size);
    int start_x, end_x;
    int start_y, end_y;
    key_t key = 5678;
    pid_t *internal_forks_id;
    char *shrmtx;
    int shmid, pid=0;
    //declare shared memory
    if ((shmid = shmget(key, 
            forks * sizeof(pid_t) + sizeof(char) *dim_x*dim_y, 
            IPC_CREAT | 0666)) < 0) {
        RELEASE();
        print_error();
        exit(ERROR_SHMGET);
    }
    if ((internal_forks_id = shmat(shmid, NULL, 0)) == (void*)-1) {
        RELEASE();
        print_error();
        exit(ERROR_SHMAT);
    }
    //copy image in shared memory
    shrmtx = (char*)(internal_forks_id + forks);
    memcpy(shrmtx,mtx,sizeof(char) *dim_x*dim_y);
    for (size_t i = 0; i < forks; i++){
        if((pid=fork()) == -1){
            RELEASE();
            print_error();
            exit(ERROR_FORK);
        }
        //prevent child of forking
        if(pid == 0){ usleep(100);break; }
        internal_forks_id[i] = pid;
    }
    //child process
    if(pid == 0){
        int child_num = 0;
        //get child number
        while(child_num<forks && 
              getpid()!=internal_forks_id[child_num])++child_num;
        if(child_num>=forks) exit(ERROR_NOT_IN_ARR);
        switch (forks){
            case 2:
                start_x = 0; end_x = dim_x;
                start_y = child_num*(dim_y / 2); end_y = (child_num+1)*(dim_y / 2);
                break;
            case 3:
                start_x = 0; end_x = dim_x;
                start_y = child_num*(dim_y / 3); end_y = (child_num+1)*(dim_y / 3);
                break;
            case 4:
                start_x = (child_num/2)*(dim_x / 2);end_x = (child_num/2 +1)*(dim_x / 2);
                start_y = (child_num/2)*(dim_y / 2);end_y = (child_num/2 +1)*(dim_y / 2);
                break;
            case 9:
                start_x = (child_num/3)*(dim_x / 3);end_x = (child_num/3 +1)*(dim_x / 3);
                start_y = (child_num/3)*(dim_y / 3);end_y = (child_num/3 +1)*(dim_y / 3);
                break;
            default:
                printf("case not supported");
                exit(ERROR_SHMCTL);
                break;
        }
        for (size_t i = start_y*child_num; i < end_y; ++i){
            for (size_t j = start_x; j < end_x; ++j){
                shrmtx[i*dim_x+j] = fnc(shrmtx,i,j,window_size,dim_y,dim_x);
            }
        }
        exit(END_OK);
    }    
    int status;
    pid_t wpid;
    for (size_t i = 0; i < forks; i++){
        if((wpid = wait(&status)) == (pid_t)-1){
            // in this case we decide to not finish the program even if 
            // there is an error
            print_error();
        }
    }
    //re-copy from shrmtx to mtx
    memcpy(mtx,shrmtx,sizeof(char) *dim_x*dim_y);
    //Release shared memory
    RELEASE();
}

