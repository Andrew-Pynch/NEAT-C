#!/bin/bash

# Function to compile the project
compile_project() {
    echo "Compiling project..."
    gcc -g -c game.c -o game.o
    gcc -g -c renderer.c -o renderer.o
    gcc -g -o game game.o renderer.o -lglfw -lGL -lm
    if [ $? -eq 0 ]; then
        echo "Compilation successful"
    else
        echo "Compilation failed"
        return 1
    fi
}

run_game() {
    echo "Running game..."
    ./game
}


# Check if nodemon is installed
if ! command -v nodemon &> /dev/null; then
    echo "nodemon is not installed. Please install it using: npm install -g nodemon"
    exit 1
else
    echo "nodemon is installed at: $(which nodemon)"
fi

# Initial compilation
echo "Performing initial compilation..."
if compile_project; then
    echo "Initial compilation successful. Running game..."
    run_game
else
    echo "Initial compilation failed. Exiting."
    exit 1
fi

echo "Setting up nodemon watch..."
# Use nodemon to watch for changes and recompile
nodemon --watch '*.c' --watch '*.h' --ext 'c,h' --exec 'bash -c "echo Recompiling...; compile_project && run_game"'

echo "Script reached end. This line should not be printed if nodemon is running correctly."
