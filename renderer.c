#include "renderer.h"
#include "common.h"
#include <stdio.h>

void error_callback(int error, const char *description) {
    fprintf(stderr, "Error: %s\n", description);
}

void key_callback(GLFWwindow *window, int key, int scancode, int action,
                  int mods) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

void init_renderer() {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glfwSetErrorCallback(error_callback);
    if (!glfwInit()) {
        fprintf(stderr, "Failed to initialize GLFW\n");
        exit(EXIT_FAILURE);
    }
}

GLFWwindow *create_window(int width, int height, const char *title) {
    GLFWmonitor *primary = glfwGetPrimaryMonitor();
    const GLFWvidmode *mode = glfwGetVideoMode(primary);

    glfwWindowHint(GLFW_RED_BITS, mode->redBits);
    glfwWindowHint(GLFW_GREEN_BITS, mode->greenBits);
    glfwWindowHint(GLFW_BLUE_BITS, mode->blueBits);
    glfwWindowHint(GLFW_REFRESH_RATE, mode->refreshRate);

    GLFWwindow *window = glfwCreateWindow(width, height, title, NULL, NULL);
    if (!window) {
        fprintf(stderr, "Failed to create GLFW window\n");
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    // Position the window on the primary monitor
    int monitorX, monitorY;
    glfwGetMonitorPos(primary, &monitorX, &monitorY);
    glfwSetWindowPos(window, monitorX + (mode->width - width) / 2,
                     monitorY + (mode->height - height) / 2);

    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, key_callback);
    return window;
}

void draw_cell(int x, int y, int cell_width, float r, float g, float b) {
    float half_width = (float)cell_width / 2;
    float cell_x = x * CELL_SIZE + (CELL_SIZE / 2);
    float cell_y = y * CELL_SIZE + (CELL_SIZE / 2);

    glColor3f(r, g, b);
    glBegin(GL_QUADS);
    glVertex2f(cell_x - half_width, cell_y - half_width);
    glVertex2f(cell_x + half_width, cell_y - half_width);
    glVertex2f(cell_x + half_width, cell_y + half_width);
    glVertex2f(cell_x - half_width, cell_y + half_width);
    glEnd();
}

void draw_grid(Population *pop, CellType grid[WIDTH][HEIGHT]) {
    glClear(GL_COLOR_BUFFER_BIT);

    // Draw empty cells and food
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            switch (grid[x][y]) {
            case EMPTY:
                draw_cell(x, y, DEFAULT_CELL_SIZE, 0.0f, 0.0f, 0.0f); // Black
                break;
            case FOOD:
                draw_cell(x, y, DEFAULT_CELL_SIZE, 0.0f, 1.0f, 0.0f); // Green
                break;
            case PREDATOR:
                // We'll draw predators separately
                break;
            }
        }
    }

    // Draw predators
    for (int i = 0; i < pop->num_predators; i++) {
        Predator *pred = &pop->predators[i];
        if (pred->health > 0) {
            draw_cell(pred->x, pred->y, pred->cell_size, 1.0f, 0.0f,
                      0.0f); // Red
        }
    }
}

void render_frame(GLFWwindow *window, Population *pop,
                  CellType grid[WIDTH][HEIGHT]) {
    glViewport(0, 0, WIDTH * CELL_SIZE, HEIGHT * CELL_SIZE);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, WIDTH * CELL_SIZE, 0, HEIGHT * CELL_SIZE, -1, 1);
    glMatrixMode(GL_MODELVIEW);

    draw_grid(pop, grid);

    glfwSwapBuffers(window);
    glfwPollEvents();
}

void cleanup_renderer(GLFWwindow *window) {
    glfwDestroyWindow(window);
    glfwTerminate();
}

bool should_close_window(GLFWwindow *window) {
    return glfwWindowShouldClose(window);
}

void process_input(GLFWwindow *window) { glfwPollEvents(); }

bool check_collision(int x1, int y1, int size1, int x2, int y2, int size2) {
    // Calculate the boundaries of each object
    int left1 = x1 - size1 / 2;
    int right1 = x1 + size1 / 2;
    int top1 = y1 + size1 / 2;
    int bottom1 = y1 - size1 / 2;

    int left2 = x2 - size2 / 2;
    int right2 = x2 + size2 / 2;
    int top2 = y2 + size2 / 2;
    int bottom2 = y2 - size2 / 2;

    // Check for overlap
    if (left1 > right2 || right1 < left2) {
        return false; // No horizontal overlap
    }

    if (bottom1 > top2 || top1 < bottom2) {
        return false; // No vertical overlap
    }

    return true;
}
