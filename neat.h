#ifndef NEAT_H
#define NEAT_H

#include "common.h"
#include <stdbool.h>
#include <stdlib.h>

// NEAT-specific structures and constants
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
    int cell_size;
} Predator;

typedef struct {
    Species *species;
    int num_species;
    int generation;
    Predator *predators;
    int num_predators;
    int current_game_step;
    int current_generation;
} Population;

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

void initialize_population(Population *pop, int pop_size);
void initialize_genome(Genome *genome, int num_inputs, int num_outputs);
Genome *crossover(Genome *parent1, Genome *parent2);
void mutate(Genome *genome);
void add_node_mutation(Genome *genome);
void add_connection_mutation(Genome *genome);
float calculate_compatibility(Genome *genome1, Genome *genome2);
float calculate_average_weight_diff(Genome *genome1, Genome *genome2);
void clear_species(Population *pop);
void speciate(Population *pop);
void reproduce(Population *pop);
int get_max_innovation(Genome *genome);
bool is_node_connected(Genome *genome, int from_node, int to_node);
void evaluate_fitness(Population *pop);
Action get_action(Population *pop, Predator *predator);
Connection *find_connection_by_innovation(Genome *genome,
                                          int innovation_number);
void cleanup_population(Population *pop);
void add_connection_to_child(Genome *child, Connection *connection);
void evaluate_network(Population *pop, Genome *genome, float *inputs,
                      float *outputs);
Genome *clone_genome(Genome *original);
bool is_valid_node_id(Genome *genome, int node_id);
float activate(float x);

#endif // NEAT_H
