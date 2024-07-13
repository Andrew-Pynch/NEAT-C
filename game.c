/*
 * NEAT-C is a C implementation of the NEAT algorithm, a genetic algorithm for
 * evolving NN toppologies.
 *
 * This includes a custom environment for the NEAT algorithm, where predators
 * compete for natural resources in an evolving environment.
 */
#include "game.h"
#include "common.h"
#include "neat.h"
#include "renderer.h"

#include <GLFW/glfw3.h>
#include <assert.h>
#include <math.h>
#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

CellType grid[WIDTH][HEIGHT];
int current_generation = 0;
int current_game_step = 0;

void init_grid() {
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            grid[x][y] = EMPTY;
        }
    }
}

void init_food() {
    int food_placed = 0;
    while (food_placed < NUM_FOOD) {
        int food_x = rand() % WIDTH;
        int food_y = rand() % HEIGHT;
        if (grid[food_x][food_y] == EMPTY) {
            set_cell(food_x, food_y, FOOD);
            food_placed++;
        }
    }
}

void set_cell(int x, int y, CellType type) { grid[x][y] = type; }
int distance(int x1, int y1, int x2, int y2) {
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

int find_nearest_food(int x, int y) {
    int min_distance = MAX_DISTANCE;
    for (int food_x = 0; food_x < WIDTH; food_x++) {
        for (int food_y = 0; food_y < HEIGHT; food_y++) {
            if (grid[food_x][food_y] == FOOD) {
                int dist = distance(x, y, food_x, food_y);
                if (dist < min_distance) {
                    min_distance = dist;
                }
            }
        }
    }

    return min_distance;
}

void set_adjacent_tiles(int x, int y, CellType adjacent_tiles[8]) {
    /*
     * =  =  =
     * = x,y =
     * =  =  =
     * */
    // up left
    adjacent_tiles[0] = (y < HEIGHT - 1) ? grid[x][y + 1] : EMPTY;
    // up
    adjacent_tiles[1] = (y < HEIGHT - 1) ? grid[x][y + 1] : EMPTY;
    // up right
    adjacent_tiles[2] = (y < HEIGHT - 1) ? grid[x][y + 1] : EMPTY;
    // left
    adjacent_tiles[3] = (x > 0) ? grid[x - 1][y] : EMPTY;
    // right
    adjacent_tiles[4] = (x < WIDTH - 1) ? grid[x + 1][y] : EMPTY;
    // down left
    adjacent_tiles[5] = (y > 0) ? grid[x][y - 1] : EMPTY;
    // down
    adjacent_tiles[6] = (y > 0) ? grid[x][y - 1] : EMPTY;
    // down right
    adjacent_tiles[7] = (y > 0) ? grid[x][y - 1] : EMPTY;
}

void consume_food(CellType grid[WIDTH][HEIGHT], int x, int y, Predator *pred) {
    bool is_pred_colliding =
        check_collision(pred->x, pred->y, pred->cell_size, x, y, CELL_SIZE);
}

void update_grid(Population *pop) {
    int food_eaten_this_step = 0;
    // test

    // Clear the grid of predators (but keep food in place)
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            if (grid[x][y] == PREDATOR) {
                grid[x][y] = EMPTY;
            }
        }
    }

    for (int i = 0; i < pop->num_predators; i++) {
        Predator *pred = &pop->predators[i];
        int old_x = pred->x;
        int old_y = pred->y;

        Action action = get_action(pop, pred);

        switch (action) {
        case MOVE_UP_LEFT:
            pred->x = (pred->x - 1 + WIDTH) % WIDTH;
            pred->y = (pred->y - 1 + HEIGHT) % HEIGHT;
            break;
        case MOVE_UP:
            pred->y = (pred->y + 1) % HEIGHT;
            break;
        case MOVE_UP_RIGHT:
            pred->x = (pred->x + 1) % WIDTH;
            pred->y = (pred->y - 1 + HEIGHT) % HEIGHT;
            break;
        case MOVE_LEFT:
            pred->x = (pred->x - 1 + WIDTH) % WIDTH;
            break;
        case MOVE_RIGHT:
            pred->x = (pred->x + 1) % WIDTH;
            break;
        case MOVE_DOWN_LEFT:
            pred->x = (pred->x - 1 + WIDTH) % WIDTH;
            pred->y = (pred->y + 1) % HEIGHT;
            break;
        case MOVE_DOWN:
            pred->y = (pred->y - 1 + HEIGHT) % HEIGHT;
            break;
        case MOVE_DOWN_RIGHT:
            pred->x = (pred->x + 1) % WIDTH;
            pred->y = (pred->y + 1) % HEIGHT;
            break;
        case STAY:
            break;
        }

        // Check for food collisions within MAX_CELL_SIZE range
        int range = pred->cell_size / 2;
        for (int dx = -range; dx <= range; dx++) {
            for (int dy = -range; dy <= range; dy++) {
                int check_x = (pred->x + dx + WIDTH) % WIDTH;
                int check_y = (pred->y + dy + HEIGHT) % HEIGHT;
                if (grid[check_x][check_y] == FOOD) {
                    if (check_collision(pred->x, pred->y, pred->cell_size,
                                        check_x, check_y, DEFAULT_CELL_SIZE)) {
                        pred->food_eaten++;
                        pred->health =
                            fmin(HEALTH_POINTS,
                                 pred->health + HEALTH_RESTORED_PER_FOOD);
                        grid[check_x][check_y] = EMPTY;
                        food_eaten_this_step++;
                    }
                }
            }
        }

        // Handle collisions with other predators
        for (int j = 0; j < pop->num_predators; j++) {
            if (i != j &&
                check_collision(pred->x, pred->y, pred->cell_size,
                                pop->predators[j].x, pop->predators[j].y,
                                pop->predators[j].cell_size)) {
                // both predators lose 1 health point from fighting
                pred->health--;
                pop->predators[j].health--;

                // If one predator dies, the other gains a kill
                if (pred->health <= 0 && pop->predators[j].health > 0) {
                    pop->predators[j].enemies_slain++;
                } else if (pop->predators[j].health <= 0 && pred->health > 0) {
                    pred->enemies_slain++;
                }
            }
        }

        // Update cell size based on food eaten
        float growth = log2f(pred->food_eaten * GROWTH_FACTOR + 1);
        pred->cell_size =
            MIN_CELL_SIZE + (MAX_CELL_SIZE - MIN_CELL_SIZE) *
                                (growth / log2f(GROWTH_FACTOR * NUM_FOOD + 1));

        // Ensure cell size stays within bounds
        pred->cell_size =
            fmaxf(MIN_CELL_SIZE, fminf(MAX_CELL_SIZE, pred->cell_size));

        // Decrease health over time (incentivize exploration)
        pred->health = fmax(0.0f, pred->health - 0.0001f);

        // Set the grid cell
        if (pred->health > 0) {
            set_cell(pred->x, pred->y, PREDATOR);
        }
    }

    // Replenish only the food that was eaten this step
    if (food_eaten_this_step > 0) {
        int current_food_count = get_grid_current_food_count();
        int food_to_place =
            fmin(food_eaten_this_step, NUM_FOOD - current_food_count);
        int attempts = 0;
        while (food_to_place > 0 && attempts < 1000) {
            int food_x = rand() % WIDTH;
            int food_y = rand() % HEIGHT;

            if (grid[food_x][food_y] == EMPTY) {
                set_cell(food_x, food_y, FOOD);
                food_to_place--;
            }
            attempts++;
        }
    }
}

int get_grid_current_food_count() {
    int food_count = 0;
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            if (grid[x][y] == FOOD) {
                food_count++;
            }
        }
    }
    return food_count;
}

void simulation(Population *pop, bool render, bool force_render) {
    log_string("Starting simulation");

    GLFWwindow *window = NULL;
    if (render || force_render) {
        log_string("Initializing renderer");
        init_renderer();
        window =
            create_window(WIDTH * CELL_SIZE, HEIGHT * CELL_SIZE, "NEAT Game");
        log_string("Renderer initialized");
    }

    init_grid();
    init_food();

    current_generation = 0;
    while (current_generation < MAX_GENERATIONS &&
           (!render || !should_close_window(window))) {
        // Print population statistics at the start of each generation
        print_population_stats(pop);

        // Run 100 game steps for this generation
        for (pop->current_game_step = 0; pop->current_game_step < 200;
             pop->current_game_step++) {
            log_string("Updating grid for step");
            log_int(current_game_step);
            update_grid(pop);

            if (render) {
                render_frame(window, pop, grid);
            }

            usleep(50000); // Sleep for 50ms between steps (adjust as needed)
        }

        log_string("Evaluating fitness");
        evaluate_fitness(pop);

        // Pause every 10 generations to let the simulation play out
        if (current_generation % VISUALIZE_EVERY_GENERATION == 0 ||
            current_generation == 0 && current_generation > 0) {
            for (int i = 0; i < 200; i++) {
                // Run for 200 steps without evolving
                update_grid(pop);
                if (render || force_render) {
                    render_frame(window, pop, grid);
                }
                usleep(50000);
            }
        }

        log_string("Speciating population");
        speciate(pop);

        log_string("Reproducing population");
        reproduce(pop);

        current_generation++;
        printf("Generation %d completed\n", current_generation);
    }

    if (render || force_render) {
        cleanup_renderer(window);
    }
}

void print_population_stats(Population *pop) {
    float total_fitness = 0.0f;
    float min_fitness = INFINITY;
    float max_fitness = -INFINITY;
    int total_connections = 0;
    int total_nodes = 0;

    for (int i = 0; i < pop->num_predators; i++) {
        Genome *genome = pop->predators[i].genome;
        float fitness = genome->fitness;

        total_fitness += fitness;
        min_fitness = fminf(min_fitness, fitness);
        max_fitness = fmaxf(max_fitness, fitness);

        total_connections += genome->num_connections;
        total_nodes += genome->num_nodes;
    }

    float avg_fitness = total_fitness / pop->num_predators;
    float avg_connections = (float)total_connections / pop->num_predators;
    float avg_nodes = (float)total_nodes / pop->num_predators;

    printf("Generation %d Statistics:\n", current_generation);
    printf("  Min Fitness: %.2f\n", min_fitness);
    printf("  Max Fitness: %.2f\n", max_fitness);
    printf("  Avg Fitness: %.2f\n", avg_fitness);
    printf("  Avg Connections: %.2f\n", avg_connections);
    printf("  Avg Nodes: %.2f\n", avg_nodes);
    printf("  Number of Species: %d\n", pop->num_species);
}

void segfault_handler(int signum) {
    printf("Segfault handler called with signal %d\n", signum);
    exit(EXIT_FAILURE);
}

int main(void) {
    signal(SIGSEGV, segfault_handler);
    log_string("Starting program");

    Population pop;
    log_string("Initializing population");
    initialize_population(&pop, INITIAL_POPULATION_SIZE);

    log_string("Starting simulation");
    simulation(&pop, DEBUG, true); // Set force_render to true

    log_string("Cleaning up population");
    cleanup_population(&pop);

    log_string("Program completed successfully");
    exit(EXIT_SUCCESS);
}
