CC = gcc
CXX = g++
CFLAGS = -c
CXXFLAGS = -c -std=c++20 -I eigen-3.4.0 -I googletest
LDFLAGS =

OBJ_DIR = obj
INCLUDE_DIR = include
SRC_DIR = src

SRCS = $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(filter %.cpp,$(SRCS))) $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(filter %.c,$(SRCS)))

all: program

program: $(OBJS) main.o
	$(CXX) $(LDFLAGS) $(OBJS) main.o -o program

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INCLUDE_DIR)/%.h
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) $< -o $@

main.o: main.cpp $(wildcard $(INCLUDE_DIR)/*.h)
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) main.cpp -o main.o

clean:
	rm -f $(OBJ_DIR)/*.o program
