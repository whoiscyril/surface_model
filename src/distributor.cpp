#include "Input_Parser.h"
#include "Simple_Math.h"
#include <Eigen/Dense>
#include <Struct_Atom.h>
#include <iostream>


//Methods to distribute all input parameters into corresponding structs;

//Populate lattice constants;

Eigen::Vector3d populate_lattice_const(UnitCell& unitcell, std::vector<double>& input)
{
    Eigen::Vector3d constants;
    for (size_t i = 0; i < 3; ++i)
    {
        std::cout << input[i] << std::endl;
        constants[i] = input[i];
    }
    //     for (size_t i = 3; i < input.size(); ++i)
    // {
    //     unitcell.lattice_angles[i-3] = input[i];
    // }
    // std::cout << unitcell.lattice_angles << std::endl;
    return constants;
}

//Now populate coordinates and convert them into cartesian

// void populate_coordinates(UnitCell& unitcell, std::vector<Atom>& input)
// {
// 	for (const auto& elem : input)
// 	{
//         unitcell.coordinates_frac.push_back(elem);
// 	}


// 	//Populate cartesian coordinates by first setting lattice vectors, and then transform it into cartesian

// 	Eigen::Matrix3d lattice_vectors;
// 	//function to populate lattice vectors;
// 	lattice_vectors << unitcell.lattice_constants[0], unitcell.lattice_constants[1],unitcell.lattice_constants[2], 0, 0, 0, 0, 0, 0;
// 	std::cout << lattice_vectors << std::endl;
// }