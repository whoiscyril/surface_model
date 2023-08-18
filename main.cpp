#include "bulk.h"
#include "Input_Parser.h"
#include <Eigen/Dense>
#include <iostream>
int main(int argc, char const *argv[])
{
    std::string filename = "input.in";
    Eigen::Vector3d lattice_constants;
    //Creating the bulk model
    // bulk_energy(get_input_coordinates(filename), get_input_species(filename), get_input_buckingham(filename));
    // lattice_constants=get_lattice_constants(filename);

    // for (const auto& val : lattice_constants)
    // {
    //     std::cout << val << std::endl;

    // }

    return 0;
}
