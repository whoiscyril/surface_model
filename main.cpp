#include "bulk.h"
#include "Input_Parser.h"
#include "Struct_Atom.h"
#include "distributor.h"
#include <Eigen/Dense>
#include <iostream>
int main(int argc, char const *argv[])
{
    std::string filename = "input.in";
    std::vector<double> lattice_constants;
    // //Creating the bulk model
    // // bulk_energy(get_input_coordinates(filename), get_input_species(filename), get_input_buckingham(filename));
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
    UnitCell unitcell(filename);

    // std::cout << unitcell.lattice_constants << std::endl;
    // std::cout << unitcell.lattice_angles << std::endl;

    // for (const auto& elem : unitcell.coordinates_cart)
    // {
    //     std::cout << elem.label << " " << elem.type << " " << elem.x << " " << elem.y <<" " <<elem.z << " " << elem.q<< std::endl;
    // }
    std::cout << unitcell.lattice_vectors << std::endl;
    std::cout << unitcell.reciprocal_vectors << std::endl;
    std:: cout << unitcell.volume << std::endl;
    return 0;
}
