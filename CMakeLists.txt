cmake_minimum_required(VERSION 3.16)
project(MyProgram)

# Set the C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)

# Add include directory
include_directories(include)

# Specify the path to the Eigen library
set(EIGEN_INCLUDE_DIR "./eigen")

# Collect all source files in the "src" directory
file(GLOB CPP_SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
file(GLOB C_SOURCES "${PROJECT_SOURCE_DIR}/src/*.c")

# Set the output directory for the executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

# Build the main executable
add_executable(program ${CPP_SOURCES} ${C_SOURCES} main.cpp)

# Include the Eigen library directory
target_include_directories(program PRIVATE ${EIGEN_INCLUDE_DIR})

