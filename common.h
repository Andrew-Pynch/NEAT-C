// common.h
#ifndef COMMON_H
#define COMMON_H

#define WIDTH 32
#define HEIGHT 32

#define CELL_SIZE 5
#define DEFAULT_CELL_SIZE 4
#define MIN_CELL_SIZE 1
#define MAX_CELL_SIZE 10
#define GROWTH_FACTOR 0.5f

#define NUM_FOOD 10
#define HEALTH_RESTORED_PER_FOOD 1000.0f

#define INITIAL_POPULATION_SIZE 10
#define HEALTH_POINTS 100.0f

#define MAX_DISTANCE (WIDTH + HEIGHT)

typedef enum { EMPTY, PREDATOR, FOOD } CellType;

#endif // COMMON_H