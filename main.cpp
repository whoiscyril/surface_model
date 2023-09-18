#include "bulk.h"
#include "Energy.h"
#include "Input_Parser.h"
#include "Struct_Atom.h"
#include "distributor.h"
#include "derivatives.h"
#include <Eigen/Dense>
#include <iostream>
int main(int argc, char const *argv[])
{
    std::string filename = "input.in";
    std::vector<double> lattice_constants;
    UnitCell unitcell(filename);

    double energy;
    energy = calc_electrostatics_3D(unitcell) + calc_short_range_buckingham_potential(unitcell);
    std::cout << energy << std::endl;


    // calc_electrostatics_3D(unitcell);
    // calc_forces(unitcell);
    // calc_strain_deriv(unitcell);
    // calc_lattice_deriv(unitcell);


    // double energy = calc_short_range_buckingham_potential(unitcell);
    // std::cout << energy << std::endl;
    // //Creating the bulk model
    // bulk_energy(get_input_coordinates(filename), get_input_species(filename), get_input_buckingham(filename));

    // lattice_constants=get_lattice_constants(filename);

    // UnitCell unitcell;

    // // populate_lattice_const(unitcell, lattice_constants);

    // // // // std::cout << unitcell.lattice_constants << std::endl;

    // // // // for (const auto& val : lattice_constants)
    // // // // {
    // // // //     std::cout << val << std::endl;

    // // // // }
    // // // // UnitCell unitcell;

    // // // std::cout << unitcell.lattice_angles << std::endl;

    // std::vector<Atom> atoms;
    // atoms = get_input_coordinates(filename);
    // populate_coordinates(unitcell, atoms);

    // for (const auto& elem : unitcell.coordinates_frac)
    // {
    //     std::cout << elem.label << " " << elem.type << " " << elem.x << " " << elem.y <<" " <<elem.z << std::endl;
    // }
    // // std::cout << unitcell.coordinates_frac.type << std::endl;

    // // The members of unitcell should now be populated based on the data from the parser
    // UnitCell unitcell(filename);
    // for (const auto& elem : unitcell.buckingham_potentials)
    // std::cout << elem.atom1_label <<" " << elem.atom1_type <<" " <<elem.atom2_label <<" " << elem.atom2_type<<" " <<
    // elem.A << " " << elem.rho <<" " << elem.C << " " <<
    //             elem.cut_off1 << " " << elem.cut_off2 << std::endl;
    // // std::cout << unitcell.lattice_angles << std::endl;

    // for (const auto& elem : unitcell.coordinates_cart)
    // {
    //     std::cout << elem.label << " " << elem.type << " " << elem.x << " " << elem.y <<" " <<elem.z << " " << elem.q<< std::endl;
    // }
    // std::cout << unitcell.lattice_vectors << std::endl;
    // std::cout << unitcell.reciprocal_vectors << std::endl;
    // std:: cout << unitcell.volume << std::endl;
    return 0;
}
