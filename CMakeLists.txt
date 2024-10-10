# This file aids in integrating the project with various development tools.
# It's not strictly necessary for creating the final executable.

# manual DEVELOPMENT    compilation use    "gcc -Wall -pedantic -O0 -std=c99 main.c -o factor"
# manual PRODUCTION     compilation use    "gcc -Wall -pedantic -O3 -std=c99 main.c -o factor"

cmake_minimum_required(VERSION 3.20)
project(qs C)
add_executable(qs main.c)

target_compile_options(qs PRIVATE -Wall -pedantic -g -O3 -std=c99)
