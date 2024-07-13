#ifndef NEAT_GAME_H
#define NEAT_GAME_H

#include "common.h"
#include "neat.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define DEBUG false
// game params
#define MAX_DISTANCE (WIDTH + HEIGHT)

// NEAT params

// x, y coordinates, health, food_eaten, enemies_slain, adjacent tile
// info
#define EPSILON 0.1
#define NUM_INPUTS 13
#define NUM_OUTPUTS                                                            \
    9 // Up, Down, Left, Right, Up-Left, Up-Right, Down-Left, Down-Right, Stay
#define MAX_GENERATIONS 11

#define C1 1.0 // excess genes coef
#define C2 1.0 // disjoint genes coef
#define C3 0.4 // weight difference coef
#define COMPATIBILITY_THRESHOLD 5.0
#define MAX_CONNECTION_ATTEMPTS 50

#define WEIGHT_MUTATION_RATE 0.8
#define WEIGHT_PERTURBATION_CHANCE 0.9
#define ADD_CONNECTION_RATE 0.5
#define ADD_NODE_RATE 0.5

#define VISUALIZE_EVERY_GENERATION 5
#define GRID_MEMORY_SIZE 10 // number of last positions to remember

void init_grid();
void init_food();
void set_cell(int x, int y, CellType type);
int distance(int x1, int y1, int x2, int y2);
int find_nearest_food(int x, int y);
void set_adjacent_tiles(int x, int y, CellType adjacent_tiles[8]);
void update_grid(Population *pop);
void simulation(Population *pop, bool render, bool force_render);
void print_population_stats(Population *pop);
int get_grid_current_food_count();
void consume_food(CellType grid[WIDTH][HEIGHT], int x, int y, Predator *pred);

#endif // NEAT_GAME_H
