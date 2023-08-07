//Calculate single point energy for bulk structure
#include "Struct_Atom.h"
#include "Input_Parser.h"
#include <iostream>
#include <cmath>
#include <math.h>
#include <Simple_Math.h>
#include <Energy.h>

double bulk_energy(std::vector<Atom> atom_positions, std::vector<Specie> species, std::vector<Buckingham> potentials)
{
	double energy = 0.0;
	std::string filename= "input.in";
	std::vector<double> lattice_constants(6, 0);
	//find the largest potential cut-off

	double max_pot_cutoff = 0.0;
	std::vector<Buckingham> buckingham_potentials = get_input_buckingham(filename);

	for (const auto& pot : buckingham_potentials)
	{
        max_pot_cutoff = std::max(max_pot_cutoff, pot.cut_off2);
	}

	//Now get lattice constants/vectors
	lattice_constants = get_lattice_constants(filename);

	// for (const auto& val : lattice_constants)
	// {
	// 	std::cout << val << std::endl;
	// }
	double a, b, c, alpha, beta, gamma;

	a = lattice_constants[0];
	b = lattice_constants[1];
	c = lattice_constants[2];
	alpha = lattice_constants[3];
	beta = lattice_constants[4];
	gamma = lattice_constants[5];
		

	//Find the number of repetition in x, y, z directions according to potential cut-offs

	int xrep, yrep, zrep;

	xrep = ceil(max_pot_cutoff/a);
	yrep = ceil(max_pot_cutoff/b);
	zrep = ceil(max_pot_cutoff/c);


	//Now create Periodic boundary conditions


	std::vector<Atom> bulk_supercell;
	for (auto& elem : atom_positions)
	{
		for (int i = -xrep; i <= xrep; ++i)
		{
			for (int j = -yrep; j<= yrep; ++j)
			{
				for (int k = -zrep; k<= zrep; ++k)
				{
					if (i == 0 && j == 0 && k == 0)
					{
						continue;
					}
					Atom atom;
					atom.index = elem.index;
					atom.label = elem.label;
					atom.type = elem.type;
					atom.x = elem.x + i;
					atom.y = elem.y + j;
					atom.z = elem.z + k;

					bulk_supercell.push_back(atom);
				}
			}
		}
	}

	//Now Convert fractional to cartesian
	std::vector<Atom> cell = atom_positions;
	for (auto& elem : cell)
	{
		elem.x *= a;
		elem.y *= b;
		elem.z *= c;
	}

	for (auto& elem : bulk_supercell)
	{
		elem.x *= a;
		elem.y *= b;
		elem.z *= c;
	}

	// std::cout << bulk_supercell.size() << std::endl;
	// for (auto& elem : bulk_supercell)
	// {
	// 	std::cout << elem.x << " " << elem.y << " " << elem.z << std::endl;
	// }


	//Calculate the short-range potential energy

	double short_range_energy = 0.0;
	short_range_energy = calc_short_range_buckingham_potential(cell, bulk_supercell, buckingham_potentials);

	/*
	CALCULATE ELECTROSTATIC ENERGY
	*/

	//Short-range real space contribution

	double eta = 0.0;
	double volume = 0.0;
	volume = a * b * c;
	// std::cout << volume << " " << eta << std::endl;

	//Add charge to input coord

	for (auto& elem1 : cell)
	{
		for (const auto& elem2 : species)
		{
			if (elem1.label == elem2.label && elem1.type == elem2.type)
			{
				elem1.q = elem2.charge;
			} 
		}
	}
	for (auto& elem1 : bulk_supercell)
	{
		for (const auto& elem2 : species)
		{
			if (elem1.label == elem2.label && elem1.type == elem2.type)
			{
				elem1.q = elem2.charge;
			} 
		}
	}


	// for (const auto& elem : cell)
	// {
	// 	std::cout << elem.label << " " << elem.type << " " << elem.x << " " << elem.y << " " << elem.z << " " << elem.q << std::endl;
	// }

}
