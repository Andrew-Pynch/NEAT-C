#ifndef RENDERER_H
#define RENDERER_H

#include "common.h"
#include "game.h"
#include <GLFW/glfw3.h>
#include <stdbool.h>

void init_renderer();
GLFWwindow *create_window(int width, int height, const char *title);
void draw_grid(Population *pop, CellType grid[WIDTH][HEIGHT]);
void render_frame(GLFWwindow *window, Population *pop,
                  CellType grid[WIDTH][HEIGHT]);
void cleanup_renderer(GLFWwindow *window);
bool should_close_window(GLFWwindow *window);
void process_input(GLFWwindow *window);
bool check_collision(int x1, int y1, int size1, int x2, int y2, int size2);

#endif // RENDERER_H
