#include "Struct_Atom.h"
#include <fstream>
#include <iostream>
#include <sstream> // Include this header for std::istringstream
#include <stdexcept>
#include <array>

UnitCell::UnitCell(const std::string& filename)
{
    //First populate lattice constants and angles

    std::string line;
    std::ifstream in(filename, std::ios_base::in);

    if (!in)
    {
        throw std::runtime_error("error opening file: " + filename);

    }

    while(std::getline(in,line))
    {
        if (line.compare(0,4,"cell") == 0 && !line.empty())
        {
            while (std::getline(in, line) && !line.empty())
            {
                std::istringstream iss(line);
                for (size_t i = 0; i < 3; ++i)
                {
                    double val;
                    if (iss >> val)
                    {
                        lattice_constants[i] = val;
                    }
                }
                for (size_t i = 0; i < 3; ++i)
                {
                    double val;
                    if (iss >> val)
                    {
                        lattice_angles[i] = val;
                    }
                }
            }
        }
    }


    in.clear();
    in.seekg(0);
    //Now populate fractional coordinates

    Atom atom;
    while (std::getline(in, line))
    {
        // std::cout << line << std::endl;
        if (line.compare(0, 4, "frac") == 0 && !line.empty())
        {
            while (std::getline(in, line) && !line.empty())
            {
                std::istringstream iss(line);
                Atom atom;

                if (iss >> atom.label >> atom.type >> atom.x >> atom.y >> atom.z >> atom.q && !line.empty())
                {
                    coordinates_frac.push_back(atom);
                }
                else if (line.empty())
                {
                    break;
                }
                else
                {
                    throw std::runtime_error("error parsing atom data: " + line);
                }
            }
        }
    }


    //Now calculate lattice vector so can convert from frac to cart
//Obtuse angle not accounted for;

    const double& a = lattice_constants[0];
    const double& b = lattice_constants[1];
    const double& c = lattice_constants[2];

    const double alpha_rad = lattice_angles[0] * M_PI / 180.0;
    const double beta_rad = lattice_angles[1] * M_PI / 180.0;
    const double gamma_rad = lattice_angles[2] * M_PI / 180.0;

    lattice_vectors << a, 0, 0,
                    b * std::cos(gamma_rad), b * std::sin(gamma_rad), 0,
                    c * std::cos(beta_rad),
                    c * (std::cos(alpha_rad) - std::cos(beta_rad) * std::cos(gamma_rad)) /
                    std::sin(gamma_rad),
                    c * std::sqrt(1.0 - std::cos(alpha_rad) * std::cos(alpha_rad) -
                                  std::cos(beta_rad) * std::cos(beta_rad) -
                                  std::cos(gamma_rad) * std::cos(gamma_rad) +
                                  2.0 * std::cos(alpha_rad) * std::cos(beta_rad) * std::cos(gamma_rad)) /
                    std::sin(gamma_rad);

// std::cout << "Lattice Vectors:" << std::endl;
// std::cout << lattice_vectors << std::endl;

    for (auto& elem : coordinates_frac)
    {
        Eigen::Vector3d frac;
        frac << elem.x, elem.y, elem.z; // Use commas to separate values

        // Convert fractional to Cartesian for the current atom
        Eigen::Vector3d cart = lattice_vectors * frac;

        // Create a new atom for the Cartesian coordinates
        Atom cartesian_atom = elem; // Copy the atom data
        cartesian_atom.x = cart[0];
        cartesian_atom.y = cart[1];
        cartesian_atom.z = cart[2];

        // Add the Cartesian atom to the coordinates_cart vector
        coordinates_cart.push_back(cartesian_atom);
    }


//Calculate volume;

    volume = lattice_vectors.row(0).dot(lattice_vectors.row(1).cross(lattice_vectors.row(2)));


//Populate reciprocal lattice vectors;
    reciprocal_vectors.row(0) = 2. * M_PI * lattice_vectors.row(1).cross(lattice_vectors.row(2)) / volume;
    reciprocal_vectors.row(1) = 2. * M_PI * lattice_vectors.row(2).cross(lattice_vectors.row(0)) / volume;
    reciprocal_vectors.row(2) = 2. * M_PI * lattice_vectors.row(0).cross(lattice_vectors.row(1)) / volume;


//Populate short-range potential information

    in.clear();
    in.seekg(0);

    line = ' ';

    while (std::getline(in, line))
    {
        if (line.compare(0, 10, "buckingham") == 0 && !line.empty())
        {
            while(std::getline(in, line))
            {
                std::istringstream iss(line);
                Buckingham potential;
                if (iss >> potential.atom1_label >> potential.atom1_type >> potential.atom2_label >> potential.atom2_type >> potential.A >> potential.rho >> potential.C >> potential.cut_off1 >> potential.cut_off2)
                {
                    buckingham_potentials.push_back(potential);
                }
                else if (line.empty())
                {
                    break;
                }
                else
                {
                    throw std::runtime_error("Error parsing atom potentials");
                }
            }
        }


    }

}