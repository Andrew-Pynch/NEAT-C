#!/bin/bash

# Function to compile the project
compile_project() {
    echo "Compiling project..."
    gcc -g -c game.c -o game.o
    gcc -g -c renderer.c -o renderer.o
    gcc -g -c neat.c -o neat.o
    gcc -g -c common.c -o common.o
    gcc -g -o game game.o renderer.o neat.o common.o -lglfw -lGL -lm
    if [ $? -eq 0 ]; then
        echo "Compilation successful"
    else
        echo "Compilation failed"
        return 1
    fi
}

# Function to run the game
run_game() {
    ./game
}

# Check if nodemon is installed
if ! command -v nodemon &> /dev/null; then
    echo "nodemon is not installed. Please install it using: npm install -g nodemon"
    exit 1
fi

# Initial compilation
compile_project && run_game

# Use nodemon to watch for changes and recompile
# nodemon --watch '*.c' --watch '*.h' --ext 'c,h' --exec 'bash -c "compile_project && run_game"'

# Uncomment the following lines if you want to use Valgrind
# valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./game
# If you want to save the Valgrind output to a file, use this command instead:
# valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./game 2> valgrind_output.txt
