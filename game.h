#ifndef NEAT_GAME_H
#define NEAT_GAME_H

#include <GLFW/glfw3.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define DEBUG false
// game params
#define WIDTH 128
#define HEIGHT 128
#define CELL_SIZE 5
#define MAX_DISTANCE (WIDTH + HEIGHT)
#define NUM_FOOD 100
#define INITIAL_POPULATION_SIZE 20
#define HEALTH_POINTS 100.0f
#define HEALTH_RESTORED_PER_FOOD 10.0f

// NEAT params

// x, y coordinates, health, food_eaten, enemies_slain, adjacent tile
// info
#define EPSILON 0.1
#define NUM_INPUTS 13
#define NUM_OUTPUTS                                                            \
    9 // Up, Down, Left, Right, Up-Left, Up-Right, Down-Left, Down-Right, Stay
#define MAX_GENERATIONS 1

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

typedef enum { EMPTY, PREDATOR, FOOD } CellType;
typedef enum {
    MOVE_UP,
    MOVE_DOWN,
    MOVE_LEFT,
    MOVE_RIGHT,
    MOVE_UP_LEFT,
    MOVE_UP_RIGHT,
    MOVE_DOWN_LEFT,
    MOVE_DOWN_RIGHT,
    STAY
} Action;

typedef struct {
    int from_node;
    int to_node;
    float weight;
    bool enabled;
    int innovation_number;
} Connection;

typedef struct {
    int id;
    enum { INPUT, HIDDEN, OUTPUT } type;
} Node;

typedef struct {
    Node *nodes;
    Connection *connections;
    int num_nodes;
    int num_connections;
    float fitness;
} Genome;

typedef struct {
    Genome *genomes;
    int population_size;
    int offspring_count;
    float adjusted_fitness;
} Species;

typedef struct {
    int x;
    int y;
    float health;
    int food_eaten;
    int enemies_slain;
    Genome *genome;
} Predator;

typedef struct {
    Species *species;
    int num_species;
    int generation;
    Predator *predators;
    int num_predators;
} Population;

void error_callback(int error, const char *description);
void key_callback(GLFWwindow *window, int key, int scancode, int action,
                  int mods);
void init_grid();
void init_food();
void draw_grid();
void set_cell(int x, int y, CellType type);
int distance(int x1, int y1, int x2, int y2);
int find_nearest_food(int x, int y);
void log_string(char *message);
void log_int(int message);
void set_adjacent_tiles(int x, int y, CellType adjacent_tiles[8]);
Action get_action(Predator *predator);
void initialize_population(Population *pop, int pop_size);
void evaluate_fitness(Population *pop);
float calculate_average_weight_diff(Genome *genome1, Genome *genome2);
float calculate_compatibility(Genome *genome1, Genome *genome2);
void clear_species(Population *pop);
void speciate(Population *pop);
void reproduce(Population *pop);
int get_max_innovation(Genome *genome);
Genome *crossover(Genome *parent1, Genome *parent2);
void mutate(Genome *genome);
void add_node_mutation(Genome *genome);
bool is_node_connected(Genome *genome, int from_node, int to_node);
void add_connection_mutation(Genome *genome);
void update_grid(Population *pop);
void simulation(Population *pop, bool render, bool force_render);
void render_frame(GLFWwindow *window);
Connection *find_connection_by_innovation(Genome *genome,
                                          int innovation_number);
void cleanup_population(Population *pop);
void cleanup_opengl(GLFWwindow *window);
void add_connection_to_child(Genome *child, Connection *connection);
float activate(float x);
void evaluate_network(Genome *genome, float *inputs, float *outputs);
Genome *clone_genome(Genome *original);
bool is_valid_node_id(Genome *genome, int node_id);
void print_population_stats(Population *pop);
int get_grid_current_food_count();

#endif // NEAT_GAME_H
