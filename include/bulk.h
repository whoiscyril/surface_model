#ifndef BULK_H
#define BULK_H
#include <vector> // Include the header for std::vector
#include "Struct_Atom.h"


double bulk_energy(std::vector<Atom> atom_positions, std::vector<Specie> species, std::vector<Buckingham> potentials);
#endif