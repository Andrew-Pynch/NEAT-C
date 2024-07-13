#!/bin/bash

# Compile with debugging symbols
gcc -g -c game.c -o game.o
gcc -g -o game game.o -lglfw -lGL -lm

# Run with Valgrind
# valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./game

./game
# If you want to save the Valgrind output to a file, use this command instead:
# valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./game 2> valgrind_output.txt
