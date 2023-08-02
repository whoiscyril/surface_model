#include <cmath>
#include <math.h>
#include <vector>
#include <Simple_Math.h>
#include <Struct_Atom.h>
#include <iostream>



double calc_short_range_buckingham_potential (const std::vector<Atom>& cell,
                                              const std::vector<Atom>& bulk_supercell,
                                              const std::vector<Buckingham>& buckingham_potentials)
{
	double short_range_energy = 0.0;

	std::vector<Atom> core_atoms;
	std::vector<Atom> shell_atoms;
	std::vector<Atom> supercell_core_atoms;
	std::vector<Atom> supercell_shell_atoms;

	//Split core and shell

	for (const auto& elem : cell)
	{
		if (elem.type == "core")
			{
				core_atoms.push_back(elem);
			}
		else if (elem.type == "shel")
		{
			shell_atoms.push_back(elem);
		}
	}

	for (const auto& elem : bulk_supercell)
	{
		if (elem.type == "core")
			{
				supercell_core_atoms.push_back(elem);
			}
		else if (elem.type == "shel")
		{
			supercell_shell_atoms.push_back(elem);
		}
	}


	//Now calculate interactions between unit cell and images
for (const auto& elem1 : cell) {
    for (const auto& elem2 : bulk_supercell) {
        double distance = get_distance(elem1, elem2);
        for (const auto& pot : buckingham_potentials) {
            if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                distance <= pot.cut_off2) {
                double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
                short_range_energy += energy;
            } else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                       ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                       distance <= pot.cut_off2) {
                double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
                short_range_energy += energy;
            }
        }
    }
}

// Calculate interactions within the unit cell
for (const auto& elem1 : cell) {
    for (const auto& elem2 : cell) {
        if (&elem1 != &elem2) { // Ensure we don't calculate self-interaction
            double distance = get_distance(elem1, elem2);
            for (const auto& pot : buckingham_potentials) {
                if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                    ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                    distance <= pot.cut_off2 && distance > 0.1) {
                    double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
                    short_range_energy += energy;
                } else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                           ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                           distance <= pot.cut_off2 && distance > 0.1) {
                    double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
                    short_range_energy += energy;
                }
            }
        }
    }
}
	short_range_energy *= 0.5;
    return short_range_energy;
}
