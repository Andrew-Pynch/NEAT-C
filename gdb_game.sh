#!/bin/bash

# Compile with debugging symbols
gcc -g -c game.c -o game.o
gcc -g -o game game.o -lglfw -lGL -lm

# Check if compilation was successful
if [ $? -ne 0 ]; then
    echo "Compilation failed. Exiting."
    exit 1
fi

# Create a GDB command file
echo "run
bt full
quit" > gdb_commands.txt

# Run GDB with the command file
gdb -x gdb_commands.txt --args ./game

# Clean up the temporary command file
rm gdb_commands.txt
