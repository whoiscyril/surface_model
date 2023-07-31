#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <string>
#include "Struct_Atom.h" // Include the Atom struct definition here
#include <vector>

// int echo_input(std::string filename);
std::vector<Atom> get_input_coordinates(std::string filename);
std::vector<Specie> get_input_species(std::string filename);
std::vector<Buckingham> get_input_buckingham(std::string filename);

#endif