#ifndef ENERGY_H
#define ENERGY_H
#include <vector>
double calc_pairwise_potential (const Atom& atom1, const Atom& atom2, const Buckingham& potential);

double calc_short_range_buckingham_potential (const std::vector<Atom>& cell,
        const std::vector<Atom>& bulk_supercell,
        const std::vector<Buckingham>& buckingham_potentials);

double calc_electrostatics_3D(const std::vector<Atom>& cell);

#endif