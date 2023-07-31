//Calculate single point energy for bulk structure
#include "Struct_Atom.h"
#include "Input_Parser.h"
#include<iostream>

double bulk_energy(std::vector<Atom> atom_positions, std::vector<Specie> species, std::vector<Buckingham> potentials)
{
	double energy = 0.0;
	std::string filename= "input.in";
	//find the largest potential cut-off

	double max_pot_cutoff = 0.0;
	std::vector<Buckingham> buckingham_potentials = get_input_buckingham(filename);

	for (const auto& pot : buckingham_potentials)
	{
        max_pot_cutoff = std::max(max_pot_cutoff, pot.cut_off2);
	}

	//Now get lattice constants/vectors

	return energy;
}
