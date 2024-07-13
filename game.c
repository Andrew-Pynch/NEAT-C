/*
 * NEAT-C is a C implementation of the NEAT algorithm, a genetic algorithm for
 * evolving NN toppologies.
 *
 * This includes a custom environment for the NEAT algorithm, where predators
 * compete for natural resources in an evolving environment.
 */
#include "game.h"
#include "common.h"
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

void safe_free(void **pp) {
    if (pp != NULL && *pp != NULL) {
        // printf("Freeing memory at address: %p\n", *pp);
        free(*pp);
        *pp = NULL;
    }
}

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

void log_string(char *message) {
    if (DEBUG) {
        printf("%s\n", message);
    }
}

void log_int(int message) {
    if (DEBUG) {
        printf("%d\n", message);
    }
}

void log_float(float message) {
    if (DEBUG) {
        printf("%f\n", message);
    }
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

Action get_action(Predator *pred) {
    float inputs[NUM_INPUTS];
    float outputs[NUM_OUTPUTS];

    // Get information about adjacent tiles
    CellType adjacent_tiles[8];
    set_adjacent_tiles(pred->x, pred->y, adjacent_tiles);

    // Set inputs based on adjacent tiles
    for (int i = 0; i < 8; i++) {
        switch (adjacent_tiles[i]) {
        case EMPTY:
            inputs[i] = 0.0f;
            break;
        case FOOD:
            inputs[i] = 0.5f;
            break;
        case PREDATOR:
            inputs[i] = 1.0f;
            break;
        }
    }

    // Add predator's current health as an input
    inputs[8] = (float)pred->health / HEALTH_POINTS;

    // Add predator's current position
    inputs[9] = (float)pred->x / WIDTH;
    inputs[10] = (float)pred->y / HEIGHT;

    // Add food eaten and enemies slain
    inputs[11] =
        (float)pred->food_eaten / 10.0f; // Normalize assuming max 10 food eaten
    inputs[12] = (float)pred->enemies_slain /
                 5.0f; // Normalize assuming max 5 enemies slain

    evaluate_network(pred->genome, inputs, outputs);

    // Epsilon-greedy exploration
    if ((float)rand() / RAND_MAX < EPSILON) {
        return rand() % NUM_OUTPUTS;
    }

    // Choose action based on highest output
    int max_index = 0;
    for (int i = 1; i < NUM_OUTPUTS; i++) {
        if (outputs[i] > outputs[max_index]) {
            max_index = i;
        }
    }

    return (Action)max_index;
}

void initialize_population(Population *pop, int pop_size) {
    // log pop size
    log_string("Initializing population with size:");
    log_int(pop_size);

    pop->species = NULL;
    pop->num_species = 0;
    pop->num_predators = pop_size;
    pop->predators = malloc(sizeof(Predator) * pop_size);
    pop->generation = 0;

    log_string("Allocated memory for predators");

    for (int i = 0; i < pop_size; i++) {
        Predator *pred = &pop->predators[i];
        pred->x = rand() % WIDTH;
        pred->y = rand() % HEIGHT;
        pred->health = HEALTH_POINTS;
        pred->food_eaten = 0;
        pred->enemies_slain = 0;
        pred->genome = malloc(sizeof(Genome));
        pred->cell_size = DEFAULT_CELL_SIZE;

        log_string("Initializing genome for predator");
        log_int(i);

        Genome *genome = pred->genome;

        // init nodes
        pred->genome->num_nodes = NUM_INPUTS + NUM_OUTPUTS;
        pred->genome->nodes = malloc(sizeof(Node) * pred->genome->num_nodes);

        // setup the input for the nodes
        for (int j = 0; j < NUM_INPUTS; j++) {
            pred->genome->nodes[j].id = j;
            pred->genome->nodes[j].type = INPUT;
        }

        // setup for output nodes
        for (int j = NUM_INPUTS; j < pred->genome->num_nodes; j++) {
            pred->genome->nodes[j].id = j; // TODO: Convert this to UUID?
            pred->genome->nodes[j].type = OUTPUT;
        }

        // init the connections for the pred->genome
        pred->genome->num_connections = NUM_INPUTS * NUM_OUTPUTS;
        pred->genome->connections =
            malloc(sizeof(Connection) * pred->genome->num_connections);

        int cxn_idx = 0;
        for (int in = 0; in < NUM_INPUTS; in++) {
            for (int out = NUM_INPUTS; out < pred->genome->num_nodes; out++) {
                pred->genome->connections[cxn_idx].from_node = in;
                pred->genome->connections[cxn_idx].to_node = out;
                // Random weight between -1 and 1
                pred->genome->connections[cxn_idx].weight =
                    ((float)rand() / RAND_MAX) * 4 - 2;
                pred->genome->connections[cxn_idx].enabled = true;
                // Simple initial innovation numbering
                pred->genome->connections[cxn_idx].innovation_number = cxn_idx;
                cxn_idx++;
            }
        }

        // init fitness
        pred->genome->fitness = 0.0f;

        log_string("Initialized genome for predator");
    }

    log_string("Finished initializing population");
}

void evaluate_fitness(Population *pop) {
    for (int i = 0; i < pop->num_predators; i++) {
        Predator *pred = &pop->predators[i];

        float food_multiplier = 1.0f + (pred->food_eaten * 0.5f);
        float food_score = pred->food_eaten * 100.0f * food_multiplier;
        float health_score = pred->health;
        float combat_score = pred->enemies_slain * 2.0f;

        // New scoring components
        float exploration_score = 0.0f;
        float diverse_action_score = 0.0f;
        float time_alive_score =
            current_game_step * 0.1f; // Reward for surviving longer

        // Calculate exploration score
        static int recent_positions[GRID_MEMORY_SIZE][2] = {0};
        static int position_index = 0;
        bool new_position = true;
        for (int j = 0; j < GRID_MEMORY_SIZE; j++) {
            if (recent_positions[j][0] == pred->x &&
                recent_positions[j][1] == pred->y) {
                new_position = false;
                break;
            }
        }
        if (new_position) {
            exploration_score += 5.0f;
            recent_positions[position_index][0] = pred->x;
            recent_positions[position_index][1] = pred->y;
            position_index = (position_index + 1) % GRID_MEMORY_SIZE;
        }

        // Calculate diverse action score
        static int action_counts[5] = {0}; // Assuming 5 possible actions
        Action last_action = get_action(pred);
        action_counts[last_action]++;
        float action_entropy = 0.0f;
        for (int j = 0; j < 5; j++) {
            if (action_counts[j] > 0) {
                float p = (float)action_counts[j] / current_game_step;
                action_entropy -= p * log2f(p);
            }
        }
        diverse_action_score = action_entropy * 10.0f;

        // Proximity-based food score
        int nearest_food_distance = find_nearest_food(pred->x, pred->y);
        float proximity_score = 100.0f / (nearest_food_distance + 1) * 1.5f;

        // Penalize constant upward movement
        float up_movement_penalty = 0.0f;
        if (last_action == MOVE_UP) {
            up_movement_penalty = -5.0f;
        }

        // Calculate final fitness
        pred->genome->fitness = food_score + health_score + combat_score +
                                exploration_score + diverse_action_score +
                                time_alive_score + proximity_score +
                                up_movement_penalty;
    }
}

float calculate_average_weight_diff(Genome *genome1, Genome *genome2) {
    float total_weight_diff = 0.0f;
    int matching_genes = 0;

    for (int i = 0; i < genome1->num_connections; i++) {
        for (int j = 0; j < genome2->num_connections; j++) {
            if (genome1->connections[i].innovation_number ==
                genome2->connections[j].innovation_number) {
                total_weight_diff += fabs(genome1->connections[i].weight -
                                          genome2->connections[j].weight);
                matching_genes++;
                break;
            }
        }
    }

    return (matching_genes > 0) ? (total_weight_diff / matching_genes) : 0.0f;
}

float calculate_compatibility(Genome *genome1, Genome *genome2) {
    // 1. Count excess and disjoint genes
    int max_innovation1 = 0;
    int max_innovation2 = 0;
    int i = 0;
    int j = 0;
    int excess = 0;
    int disjoint = 0;

    // max innovation number for each genome
    for (int k = 0; k < genome1->num_connections; k++) {
        if (genome1->connections[k].innovation_number > max_innovation1) {
            max_innovation1 = genome1->connections[k].innovation_number;
        }
    }

    for (int k = 0; k < genome2->num_connections; k++) {
        if (genome2->connections[k].innovation_number > max_innovation2) {
            max_innovation2 = genome2->connections[k].innovation_number;
        }
    }

    // compare connections
    while (i < genome1->num_connections && j < genome2->num_connections) {
        int genome1_innovation_number =
            genome1->connections[i].innovation_number;
        int genome2_innovation_number =
            genome2->connections[j].innovation_number;
        if (genome1_innovation_number == genome2_innovation_number) {
            // matching genes
            i++;
            j++;
        } else if (genome1_innovation_number < genome2_innovation_number) {
            // gene in genome1 but not in genome2
            if (genome1_innovation_number >
                fmin(max_innovation1, max_innovation2)) {
                excess++;
            } else {
                disjoint++;
            }
            i++;
        } else {
            // gene in genome2 but not in genome1
            if (genome2_innovation_number >
                fmin(max_innovation1, max_innovation2)) {
                excess++;
            } else {
                disjoint++;
            }
            j++;
        }
    }

    // Count remaining genes in genome1
    while (i < genome1->num_connections) {
        if (genome1->connections[i].innovation_number >
            fmin(max_innovation1, max_innovation2)) {
            excess++;
        } else {
            disjoint++;
        }
        i++;
    }

    // Count remaining genes in genome2
    while (j < genome2->num_connections) {
        if (genome2->connections[j].innovation_number >
            fmin(max_innovation1, max_innovation2)) {
            excess++;
        } else {
            disjoint++;
        }
        j++;
    }
    // 2. Calculate average weight difference of matching genes
    float average_weight_diff = calculate_average_weight_diff(genome1, genome2);
    // 3. Apply the compatibility formula
    int n = fmax(genome1->num_connections, genome2->num_connections);
    if (n < 20) {
        n = 1; // normalize for really tiny genomes
    }

    float compatibility =
        (C1 * excess / n) + (C2 * disjoint / n) + (C3 * average_weight_diff);

    return compatibility;
}

void clear_species(Population *pop) {
    for (int s = 0; s < pop->num_species; s++) {
        Species *species = &pop->species[s];
        if (species->population_size > 0) {
            Genome representative = species->genomes[0];
            safe_free((void **)&species->genomes);
            species->genomes = malloc(sizeof(Genome));
            species->genomes[0] = representative;
            species->population_size = 1;
            species->adjusted_fitness = 0.0f;
            species->offspring_count = 0;
        }
    }
}

void speciate(Population *pop) {
    // Free existing species before reallocating
    for (int s = 0; s < pop->num_species; s++) {
        safe_free((void **)&pop->species[s].genomes);
    }
    safe_free((void **)&pop->species);
    pop->num_species = 0;
    pop->species = NULL;

    // Temporary array to hold all genomes
    Genome *all_genomes = malloc(sizeof(Genome) * INITIAL_POPULATION_SIZE);
    int total_genomes = 0;
    for (int s = 0; s < pop->num_species; s++) {
        for (int g = 0; g < pop->species[s].population_size; g++) {
            all_genomes[total_genomes++] = pop->species[s].genomes[g];
        }
    }

    // speciate each genome
    for (int i = 0; i < total_genomes; i++) {
        Genome *genome = &all_genomes[i];
        bool found_species = false;

        for (int s = 0; s < pop->num_species; s++) {
            Species *species = &pop->species[s];
            // could make this a random gene from the species genome
            Genome *representative = &species->genomes[0];

            float compatibility =
                calculate_compatibility(genome, representative);

            if (compatibility > COMPATIBILITY_THRESHOLD) {
                // Add to this species
                species->population_size++;
                species->genomes =
                    realloc(species->genomes,
                            sizeof(Genome) * species->population_size);
                species->genomes[species->population_size - 1] =
                    *genome; // Shallow copy
                found_species = true;
                break;
            }
        }

        if (!found_species) {
            // create a new one!
            pop->num_species++;
            pop->species =
                realloc(pop->species, sizeof(Species) * pop->num_species);
            Species *new_species = &pop->species[pop->num_species - 1];
            new_species->genomes = malloc(sizeof(Genome));
            new_species->genomes[0] = *genome; // Shallow copy
            new_species->population_size = 1;
        }
    }

    safe_free((void **)&all_genomes);
}

Genome *crossover(Genome *parent1, Genome *parent2) {
    log_string("Starting crossover");
    Genome *child = malloc(sizeof(Genome));
    if (child == NULL) {
        log_string("Failed to allocate memory for child genome");
        return NULL;
    }
    memset(child, 0, sizeof(Genome));

    // print both parens
    log_string("\n\nParent 1:");
    log_int(parent1->fitness);

    log_string("Parent 2:");
    log_int(parent2->fitness);
    Genome *fitter_parent;
    Genome *less_fit_parent;
    if (parent1->fitness > parent2->fitness) {
        fitter_parent = parent1;
        less_fit_parent = parent2;
    } else {
        fitter_parent = parent2;
        less_fit_parent = parent1;
    }

    log_string("\n\nFitter parent:");
    log_int(fitter_parent->fitness);

    // Handle nodes
    child->num_nodes = fmax(parent1->num_nodes, parent2->num_nodes);
    child->nodes = malloc(sizeof(Node) * child->num_nodes);

    log_string("Copying nodes from fitter parent to child");
    // Copy nodes from fitter parent
    memcpy(child->nodes, fitter_parent->nodes,
           sizeof(Node) * fitter_parent->num_nodes);

    log_string("Adding nodes from less fit parent to child");
    // Add any additional nodes from the less fit parent
    for (int i = fitter_parent->num_nodes; i < child->num_nodes; i++) {
        child->nodes[i] = less_fit_parent->nodes[i];
    }

    log_string("Copying connections from fitter parent to child");
    // Handle connections
    int max_innovation =
        fmax(get_max_innovation(parent1), get_max_innovation(parent2));
    child->connections =
        malloc(sizeof(Connection) * (parent1->num_connections +
                                     parent2->num_connections)); // Worst case
    child->num_connections = 0;

    for (int i = 0; i <= max_innovation; i++) {
        Connection *gene1 = find_connection_by_innovation(parent1, i);
        Connection *gene2 = find_connection_by_innovation(parent2, i);

        if (gene1 && gene2) {
            // Both parents have this gene
            Connection *selected_gene = (rand() % 2 == 0) ? gene1 : gene2;
            add_connection_to_child(child, selected_gene);
        } else if (gene1 && fitter_parent == parent1) {
            // Only fitter parent (parent1) has this gene
            add_connection_to_child(child, gene1);
        } else if (gene2 && fitter_parent == parent2) {
            // Only fitter parent (parent2) has this gene
            add_connection_to_child(child, gene2);
        }
        // Genes unique to less fit parent are not inherited
    }

    log_string("Enabling/disabling connections from fitter parent to child");
    // Handle enabling/disabling of connections
    for (int i = 0; i < child->num_connections; i++) {
        Connection *child_gene = &child->connections[i];
        Connection *parent1_gene = find_connection_by_innovation(
            parent1, child_gene->innovation_number);
        Connection *parent2_gene = find_connection_by_innovation(
            parent2, child_gene->innovation_number);

        if ((!parent1_gene || !parent1_gene->enabled) &&
            (!parent2_gene || !parent2_gene->enabled)) {
            child_gene->enabled = false;
        } else if ((parent1_gene && !parent1_gene->enabled) ||
                   (parent2_gene && !parent2_gene->enabled)) {
            child_gene->enabled =
                (rand() % 100 < 75); // 75% chance of being enabled
        } else {
            child_gene->enabled = true;
        }
    }

    log_string("Checking if child has any valid connections");
    if (child->num_connections == 0) {
        // No valid connections, cleanup and return NULL
        safe_free((void **)&child->nodes);
        safe_free((void **)&child->connections);
        safe_free((void **)&child);
        return NULL;
    }

    log_string("Finished crossover successfully");

    return child;
}

void reproduce(Population *pop) {
    log_string("Reproducing population");

    float total_adjusted_fitness = 0.0f;
    int total_offspring = 0;

    log_string("Calculating total adjusted fitness");
    // Calculate total adjusted fitness
    for (int s = 0; s < pop->num_species; s++) {
        Species *species = &pop->species[s];
        float species_avg_fitness = 0.0f;
        for (int g = 0; g < species->population_size; g++) {
            species_avg_fitness += species->genomes[g].fitness;
        }
        species_avg_fitness /= species->population_size;

        species->adjusted_fitness =
            species_avg_fitness / species->population_size;
        total_adjusted_fitness += species->adjusted_fitness;
    }
    log_float(total_adjusted_fitness);

    log_string("Determining offspring count for each species");
    // Determine offspring count for each species
    if (total_adjusted_fitness > 0 && pop->num_species > 0) {
        for (int s = 0; s < pop->num_species; s++) {
            Species *species = &pop->species[s];
            species->offspring_count =
                (int)((species->adjusted_fitness / total_adjusted_fitness) *
                      INITIAL_POPULATION_SIZE);
            total_offspring += species->offspring_count;
        }
        log_int(total_offspring);
    } else {
        // Handle the case where total_adjusted_fitness is zero or there are no
        // species
        log_string(
            "Warning: Total adjusted fitness is zero or no species exist");
        // Create a single species with all individuals if no species exist
        if (pop->num_species == 0) {
            pop->num_species = 1;
            pop->species = malloc(sizeof(Species));
            pop->species[0].population_size = pop->num_predators;
            pop->species[0].genomes =
                malloc(sizeof(Genome) * pop->num_predators);
            for (int i = 0; i < pop->num_predators; i++) {
                pop->species[0].genomes[i] = *(pop->predators[i].genome);
            }
        }

        // Distribute offspring equally among species
        int offspring_per_species = INITIAL_POPULATION_SIZE / pop->num_species;
        for (int s = 0; s < pop->num_species; s++) {
            pop->species[s].offspring_count = offspring_per_species;
            total_offspring += offspring_per_species;
        }
        log_int(total_offspring);
    }

    log_string("Adjusting offspring counts to match desired population size");
    // Adjust offspring counts to match desired population size
    while (total_offspring < INITIAL_POPULATION_SIZE) {
        pop->species[rand() % pop->num_species].offspring_count++;
        total_offspring++;
    }
    log_int(total_offspring);

    log_string("Removing offspring from species if necessary");
    while (total_offspring > INITIAL_POPULATION_SIZE) {
        int s = rand() % pop->num_species;
        if (pop->species[s].offspring_count > 0) {
            pop->species[s].offspring_count--;
            total_offspring--;
        }
    }
    log_int(total_offspring);

    // Create new population
    log_string("Creating new population");
    Predator *new_population =
        malloc(sizeof(Predator) * INITIAL_POPULATION_SIZE);
    int new_population_size = 0;

    for (int s = 0; s < pop->num_species; s++) {
        Species *species = &pop->species[s];
        for (int i = 0; i < species->offspring_count; i++) {
            Predator *new_pred = &new_population[new_population_size];
            new_pred->x = rand() % WIDTH;
            new_pred->y = rand() % HEIGHT;
            new_pred->health = 100.0f;
            new_pred->food_eaten = 0;
            new_pred->enemies_slain = 0;
            new_pred->cell_size = DEFAULT_CELL_SIZE;

            if (species->population_size == 1 || rand() % 100 < 25) {
                // Asexual reproduction (clone and mutate)
                new_pred->genome = clone_genome(
                    &species->genomes[rand() % species->population_size]);
                log_string("Mutating new predator genome");
                mutate(new_pred->genome);
            } else {
                // Sexual reproduction (crossover and mutate)
                log_string("Crossing over genomes");
                Genome *parent1 =
                    &species->genomes[rand() % species->population_size];
                Genome *parent2 =
                    &species->genomes[rand() % species->population_size];
                new_pred->genome = crossover(parent1, parent2);
                if (new_pred->genome != NULL) {
                    log_string("Mutating crossd predator genome from parent1 "
                               "and parent2");
                    mutate(new_pred->genome);
                } else {
                    // If crossover failed, clone a parent instead
                    log_string("Cloning genome since crossover failed");
                    new_pred->genome = clone_genome(parent1);
                    mutate(new_pred->genome);
                }
            }
            new_population_size++;
        }
    }

    log_string("Replacing old population with new population");
    // Replace old population with new population
    for (int i = 0; i < pop->num_predators; i++) {
        safe_free((void **)&pop->predators[i].genome->nodes);
        safe_free((void **)&pop->predators[i].genome->connections);
        safe_free((void **)&pop->predators[i].genome);
    }
    safe_free((void **)&pop->predators);

    pop->predators = new_population;
    pop->num_predators = INITIAL_POPULATION_SIZE;

    // Increment generation counter
    pop->generation++;

    log_string("Finished reproduction");
}

int get_max_innovation(Genome *genome) {
    int max_innovation = 0;
    for (int i = 0; i < genome->num_connections; i++) {
        if (genome->connections[i].innovation_number > max_innovation) {
            max_innovation = genome->connections[i].innovation_number;
        }
    }
    return max_innovation;
}

Connection *find_connection_by_innovation(Genome *genome,
                                          int innovation_number) {
    for (int i = 0; i < genome->num_connections; i++) {
        if (genome->connections[i].innovation_number == innovation_number) {
            return &genome->connections[i];
        }
    }
    return NULL;
}

void cleanup_population(Population *pop) {
    for (int i = 0; i < pop->num_predators; i++) {
        safe_free((void **)&pop->predators[i].genome->nodes);
        safe_free((void **)&pop->predators[i].genome->connections);
        safe_free((void **)&pop->predators[i].genome);
    }
    safe_free((void **)&pop->predators);

    for (int s = 0; s < pop->num_species; s++) {
        // We don't need to free individual genomes in species, as they are
        // shallow copies
        safe_free((void **)&pop->species[s].genomes);
    }
    safe_free((void **)&pop->species);
}

void add_connection_to_child(Genome *child, Connection *connection) {
    child->connections[child->num_connections] = *connection;
    child->num_connections++;
}

void add_node_mutation(Genome *genome) {
    log_string("Starting add node mutation");
    log_string("Number of connections:");
    log_int(genome->num_connections);

    if (genome->num_connections == 0) {
        log_string("No connections to split, exiting add_node_mutation");
        return;
    }

    // choose a random connection to split
    int connection_to_split = rand() % genome->num_connections;
    Connection *old_connection = &genome->connections[connection_to_split];
    log_string("Connection to split:");
    log_int(connection_to_split);

    // disable old connection
    old_connection->enabled = false;
    log_string("Old connection disabled");

    log_string("Adding new node");
    // add a new node
    genome->num_nodes++;
    Node *new_nodes = realloc(genome->nodes, sizeof(Node) * genome->num_nodes);
    if (new_nodes == NULL) {
        log_string("Failed to allocate memory for new node");
        genome->num_nodes--;
        return;
    }
    genome->nodes = new_nodes;

    Node *new_node = &genome->nodes[genome->num_nodes - 1];
    new_node->id = genome->num_nodes - 1;
    new_node->type = HIDDEN;
    log_string("New node added");

    log_string("Adding 2 new connections");
    // add 2 new connections
    size_t new_connections_size =
        sizeof(Connection) * (genome->num_connections + 2);
    Connection *new_connections =
        realloc(genome->connections, new_connections_size);
    if (new_connections == NULL) {
        log_string("Failed to allocate memory for new connections");
        genome->num_nodes--;
        return;
    }
    genome->connections = new_connections;
    genome->num_connections += 2;
    log_string("2 new connections added");

    log_string("Creating new connections");
    // connection from source to new node
    Connection *new_connection1 =
        &genome->connections[genome->num_connections - 2];
    new_connection1->from_node = old_connection->from_node;
    new_connection1->to_node = new_node->id;
    new_connection1->weight = 1.0;
    new_connection1->enabled = true;
    new_connection1->innovation_number = genome->num_connections - 2;

    log_string("Creating new connection");
    // Connection from new node to destination
    Connection *new_connection2 =
        &genome->connections[genome->num_connections - 1];
    new_connection2->from_node = new_node->id;
    new_connection2->to_node = old_connection->to_node;
    new_connection2->weight = old_connection->weight; // Keep the old weight
    new_connection2->enabled = true;
    new_connection2->innovation_number = genome->num_connections - 1;

    log_string("Finished add node mutation");
}

bool is_node_connected(Genome *genome, int from_node, int to_node) {
    if (genome->nodes[from_node].type == OUTPUT ||
        genome->nodes[to_node].type == INPUT) {
        return false;
    }

    for (int i = 0; i < genome->num_connections; i++) {
        if (genome->connections[i].from_node == from_node &&
            genome->connections[i].to_node == to_node) {
            return true;
        }
    }
    return false;
}

void add_connection_mutation(Genome *genome) {
    log_string("Starting add_connection_mutation");
    int from_node, to_node;
    bool connection_exists;
    int attempts = 0;

    do {
        // Choose random nodes
        from_node = rand() % genome->num_nodes;
        to_node = rand() % genome->num_nodes;

        attempts++;
    } while (is_node_connected(genome, from_node, to_node) &&
             attempts < MAX_CONNECTION_ATTEMPTS);

    if (attempts == MAX_CONNECTION_ATTEMPTS) {
        log_string("Failed to add new connection after max attempts");
        return;
    }

    // create a new connection
    genome->num_connections++;
    Connection *new_connections = realloc(
        genome->connections, sizeof(Connection) * genome->num_connections);
    if (new_connections == NULL) {
        log_string("Failed to allocate memory for new connection");
        genome->num_connections--;
        return;
    }
    genome->connections = new_connections;

    // get new alloc connection (empty)
    Connection *new_connection =
        &genome->connections[genome->num_connections - 1];

    new_connection->from_node = from_node;
    new_connection->to_node = to_node;
    new_connection->weight = ((float)rand() / RAND_MAX) * 2 - 1;
    new_connection->enabled = true;
    // innovation number is cxn count
    new_connection->innovation_number = genome->num_connections - 1;

    log_string("Successfully added new connection");
}

void mutate(Genome *genome) {
    log_string("Starting mutation");
    log_string("Number of nodes:");
    log_int(genome->num_nodes);
    log_string("Number of connections:");
    log_int(genome->num_connections);

    // mutate the weights
    for (int i = 0; i < genome->num_connections; i++) {
        if ((float)rand() / RAND_MAX < WEIGHT_PERTURBATION_CHANCE) {
            // Perturb weight
            genome->connections[i].weight +=
                ((float)rand() / RAND_MAX) * 0.4 - 0.2;
        } else {
            // Assign new random weight
            genome->connections[i].weight = ((float)rand() / RAND_MAX) * 4 - 2;
        }
    }

    // add new connection mutation
    if ((float)rand() / RAND_MAX < ADD_CONNECTION_RATE) {
        log_string("Attempting to add new connection");
        add_connection_mutation(genome);
    }

    // add the node mutation
    if ((float)rand() / RAND_MAX < ADD_NODE_RATE) {
        log_string("Attempting to add new node");
        add_node_mutation(genome);
    }

    log_string("Finished mutation");
    log_string("Final number of nodes:");
    log_int(genome->num_nodes);
    log_string("Final number of connections:");
    log_int(genome->num_connections);
}

void consume_food(CellType grid[WIDTH][HEIGHT], int x, int y, Predator *pred) {
    bool is_pred_colliding =
        check_collision(pred->x, pred->y, pred->cell_size, x, y, CELL_SIZE);
}

void update_grid(Population *pop) {
    int food_eaten_this_step = 0;

    // Clear the grid of predators (but keep food in place)
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            if (grid[x][y] == PREDATOR) {
                // test
                grid[x][y] = EMPTY;
            }
        }
    }

    for (int i = 0; i < pop->num_predators; i++) {
        Predator *pred = &pop->predators[i];
        int old_x = pred->x;
        int old_y = pred->y;
        // test

        Action action = get_action(pred);

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
        int range = MAX_CELL_SIZE / 2;
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

Genome *clone_genome(Genome *original) {
    Genome *clone = malloc(sizeof(Genome));
    clone->num_nodes = original->num_nodes;
    clone->num_connections = original->num_connections;
    clone->nodes = malloc(sizeof(Node) * clone->num_nodes);
    clone->connections = malloc(sizeof(Connection) * clone->num_connections);

    memcpy(clone->nodes, original->nodes, sizeof(Node) * clone->num_nodes);
    for (int i = 0; i < clone->num_connections; i++) {
        clone->connections[i] = original->connections[i];
    }

    clone->fitness = original->fitness;

    return clone;
}

bool is_valid_node_id(Genome *genome, int node_id) {
    return node_id >= 0 && node_id < genome->num_nodes;
}

float activate(float x) {
    // Avoid division by zero
    if (x < -700.0f)
        return 0.0f;
    if (x > 700.0f)
        return 1.0f;
    return 1.0f / (1.0f + expf(-x));
}

void evaluate_network(Genome *genome, float *inputs, float *outputs) {
    log_string("Evaluating network");
    log_string("Current generation:");
    log_int(current_generation);
    log_string("Current game step:");
    log_int(current_game_step);

    log_string("Genome structure:");
    log_string("Number of nodes:");
    log_int(genome->num_nodes);
    log_string("Number of connections:");
    log_int(genome->num_connections);

    if (genome == NULL || genome->nodes == NULL ||
        genome->connections == NULL) {
        log_string("Error: Invalid genome structure");
        return;
    }

    float *node_values = calloc(genome->num_nodes, sizeof(float));
    if (node_values == NULL) {
        log_string("Error: Failed to allocate memory for node_values");
        return;
    }
    assert(genome->num_nodes >= NUM_INPUTS + NUM_OUTPUTS);
    assert(genome->num_connections > 0);
    // set the input values
    for (int i = 0; i < NUM_INPUTS; i++) {
        node_values[i] = inputs[i];
    }

    // evaluate hidden -> output layers
    for (int i = NUM_INPUTS; i < genome->num_nodes; i++) {
        float sum = 0.0f;
        for (int j = 0; j < genome->num_connections; j++) {
            Connection *conn = &genome->connections[j];
            if (conn->enabled) {
                if (conn->from_node < 0 ||
                    conn->from_node >= genome->num_nodes) {
                    log_string("Error: Invalid from_node in connection");
                    log_int(conn->from_node);
                    continue;
                }
                if (conn->to_node != i) {
                    continue;
                }
                sum += node_values[conn->from_node] * conn->weight;
            }
        }
        node_values[i] = activate(sum);
    }

    // set outputs
    for (int i = 0; i < NUM_OUTPUTS; i++) {
        outputs[i] = node_values[genome->num_nodes - NUM_OUTPUTS + i];
        log_string("Network output:");
        log_int(i);
        log_float(outputs[i]);
    }

    safe_free((void **)&node_values);
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
        for (current_game_step = 0; current_game_step < 100;
             current_game_step++) {
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
