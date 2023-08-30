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
    double total_energy = 0.;
    double short_range_energy = 0.0;
    double electrostatic_energy = 0.0;

    //Calculate the short-range potential energy
    UnitCell unitcell("input.in");
    short_range_energy = calc_short_range_buckingham_potential(unitcell);
    electrostatic_energy = calc_electrostatics_3D(unitcell);
    total_energy = short_range_energy + electrostatic_energy;
    // std::cout << short_range_energy << " " << electrostatic_energy <<" " << total_energy << std::endl;
    return total_energy;


}
