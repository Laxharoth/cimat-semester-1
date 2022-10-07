#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define bool int
#define true 1
#define false 0
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

struct _grid_meta {
  unsigned int width;
  unsigned int height;
};
struct _position {
  unsigned int y;
  unsigned int x;
};

typedef struct _grid_meta grid_meta;
typedef struct _position position;
typedef char byte;

/**
 * @brief Get the char in the grid, if out of bounds return 0
 *
 * @param holder where the char is stored
 * @param grid the grid
 * @param grid_dim the dimensions of the grid
 * @param y row
 * @param x column
 */
char get_at(char *grid, grid_meta grid_dim, unsigned int y, unsigned int x);
void set_at(char new_char, char *grid, grid_meta grid_dim, unsigned int y,
            unsigned int x);
bool is_free_at(char *grid, grid_meta grid_dim, unsigned int y, unsigned int x);
bool can_move_from(char *grid, grid_meta grid_dim, unsigned int y,
                   unsigned int x);
void print_grid(char *grid, grid_meta grid_dim);

byte generate_step();
int main(int argc, char const *argv[]) {
  time_t t;
  srand((unsigned)time(&t));
  char grid[100];
  for (int i = 0; i < 100; i++)
    grid[i] = '.';
  grid_meta meta;
  meta.height = 10;
  meta.width = 10;
  char current_char = '.';
  position current_pos;
  position new_pos;
  current_pos.y = 0;
  current_pos.x = 0;
  // initial position
  grid[0] = ++current_char;
  while (1) {
    if (!can_move_from(grid, meta, current_pos.y, current_pos.x))
      break;
    do {
      new_pos = current_pos;
      byte move = generate_step();
      switch (move) {
      case UP:
        new_pos.y--;
        break;
      case DOWN:
        new_pos.y++;
        break;
      case LEFT:
        new_pos.x--;
        break;
      case RIGHT:
        new_pos.x++;
        break;
      }
    } while (!is_free_at(grid, meta, new_pos.y, new_pos.x));
    current_pos = new_pos;
    set_at(++current_char, grid, meta, new_pos.y, new_pos.x);
  }
  print_grid(grid, meta);
  return 0;
}

char get_at(char *grid, grid_meta grid_dim, unsigned int y, unsigned int x) {
  if (x < 0 || x >= grid_dim.width || y < 0 || y >= grid_dim.height) {
    return 0;
  }
  return grid[y * grid_dim.width + x];
}
void set_at(char new_char, char *grid, grid_meta grid_dim, unsigned int y,
            unsigned int x) {
  if (x < 0 || x >= grid_dim.width || y < 0 || y >= grid_dim.height) {
    return;
  }
  grid[y * grid_dim.width + x] = new_char;
}

byte generate_step() { return rand() % 4; }

bool is_free_at(char *grid, grid_meta grid_dim, unsigned int y,
                unsigned int x) {
  if (get_at(grid, grid_dim, y, x) == '.')
    return true;
  return false;
}
bool can_move_from(char *grid, grid_meta grid_dim, unsigned int y,
                   unsigned int x) {
  return is_free_at(grid, grid_dim, y - 1, x) ||
         is_free_at(grid, grid_dim, y + 1, x) ||
         is_free_at(grid, grid_dim, y, x - 1) ||
         is_free_at(grid, grid_dim, y, x + 1);
}
void print_grid(char *grid, grid_meta grid_dim) {
  for (unsigned int i = 0; i < grid_dim.height; i++) {
    for (unsigned int j = 0; j < grid_dim.width; j++) {
      printf("%c", get_at(grid, grid_dim, i, j));
    }
    printf("\n");
  }
}