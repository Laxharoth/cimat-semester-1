#include "get_time.h"
#include "median.c"
#include "pgm1.c"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#ifdef _WIN32
#include <windows.h>
#else
#endif

#define print_error() printf("ERROR: %s", strerror(errno));

#define END_OK 0
#define ERROR_NOT_SUPPORT -1
#define ERROR_SHMGET -11
#define ERROR_SHMAT -12
#define ERROR_SHMCTL -13
#define ERROR_FORK -21
#define ERROR_WAIT -22
#define ERROR_NOT_IN_ARR -23

#ifndef WINDOW_SIZE
#define WINDOW_SIZE 3
#endif // !WINDOW_SIZE

#ifndef FORK_SIZE
#define FORK_SIZE 2
#endif // !FORK_SIZE

// pgm1 se modifico para usar apuntador simple
void fork_replace(char *mtx, int dim_y, int dim_x, filter fnc, int window_size,
                  size_t forks);
int main(int argc, char const *argv[]) {
  int cols, rows;
  char *img = (char *)pgmRead("fractal_tree.ascii.pgm", &rows, &cols);
// replace( img, rows,cols,median, 9 );
#ifndef MEASURE_TIME
  fork_replace(img, rows, cols, median, WINDOW_SIZE, FORK_SIZE);
  pgmWrite("fractal_tree.2.ascii.pgm", rows, cols, (unsigned char *)img, "");
#else
  int window_size[] = {3, 5, 7, 9};
  int forks[] = {1, 2, 3, 4, 9};
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 5; j++) {
      long time = (long)GetTimeStamp(
          fork_replace(img, rows, cols, median, window_size[i], forks[j]););
      printf("window:%d, forks:%d, time:%lf\n", window_size[i], forks[j],
             (double)time * (1e-6));
    }
  }
#endif // DEBUG
  free(img);
  return 0;
}

void fork_replace(char *mtx, int dim_y, int dim_x, filter fnc, int window_size,
                  size_t forks) {
#define RELEASE()                                                              \
  if (shmctl(shmid, IPC_RMID, 0) == -1) {                                      \
    print_error();                                                             \
    exit(ERROR_SHMCTL);                                                        \
  }
  if (forks == 1)
    return replace(mtx, dim_y, dim_x, fnc, window_size);
  int start_x, end_x;
  int start_y, end_y;
  key_t key = 5678;
  char *shrmtx;
  int shmid, pid = 0;
  // declare shared memory
  if ((shmid = shmget(key, sizeof(char) * dim_x * dim_y, IPC_CREAT | 0666)) <
      0) {
    RELEASE();
    print_error();
    exit(ERROR_SHMGET);
  }
  if ((shrmtx = shmat(shmid, NULL, 0)) == (void *)-1) {
    RELEASE();
    print_error();
    exit(ERROR_SHMAT);
  }
  memcpy(shrmtx, mtx, sizeof(char) * dim_x * dim_y);
  unsigned int child_n;
  // create children
  for (child_n = 0; child_n < forks; child_n++) {
    if ((pid = fork()) == -1) {
      RELEASE();
      print_error();
      exit(ERROR_FORK);
    }
    // prevent child of forking
    if (pid == 0) {
      usleep(100);
      break;
    }
  }
  // child process
  if (pid == 0) {
    const int child_id = child_n;
    switch (forks) {
    case 2:
      start_x = 0;
      end_x = dim_x;
      start_y = child_id * (dim_y / 2);
      end_y = (child_id + 1) * (dim_y / 2);
      break;
    case 3:
      start_x = 0;
      end_x = dim_x;
      start_y = child_id * (dim_y / 3);
      end_y = (child_id + 1) * (dim_y / 3);
      break;
    case 4:
      start_x = (child_id % 2) * (dim_x / 2);
      end_x = (child_id % 2 + 1) * (dim_x / 2);
      start_y = (child_id / 2) * (dim_y / 2);
      end_y = (child_id / 2 + 1) * (dim_y / 2);
      break;
    case 9:
      start_x = (child_id % 3) * (dim_x / 3);
      end_x = (child_id % 3 + 1) * (dim_x / 3);
      start_y = (child_id / 3) * (dim_y / 3);
      end_y = (child_id / 3 + 1) * (dim_y / 3);
      break;
    default:
      printf("case not supported");
      exit(ERROR_SHMCTL);
      break;
    }
    for (size_t i = start_y; i < end_y; ++i) {
      for (size_t j = start_x; j < end_x; ++j) {
        shrmtx[i * dim_x + j] = fnc(shrmtx, i, j, window_size, dim_y, dim_x);
      }
    }
    exit(END_OK);
  }
  int status;
  pid_t wpid;
  for (size_t i = 0; i < forks; i++) {
    if ((wpid = wait(&status)) == (pid_t)-1) {
      // in this case we decide to not finish the program even if
      // there is an error
      print_error();
    }
  }
  // re-copy from shrmtx to mtx
  memcpy(mtx, shrmtx, sizeof(char) * dim_x * dim_y);
  // Release shared memory
  RELEASE();
}
