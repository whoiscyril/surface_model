#include "Input_Parser.h"
#include "Struct_Atom.h"
#include "distributor.h"
#include <fstream>
#include <iostream>
#include <sstream> // Include this header for std::istringstream
#include <stdexcept>
#include <array>

int echo_input(std::string filename)
{
    std::string line;
    std::ifstream in(filename, std::ios_base::in);

    if (!in)
    {
        std::cerr << "Error opening file " << filename << std::endl;
        return 1; // Indicate an error occurred
    }

    while(in.peek() != EOF)
    {
        std::getline(in, line); // Read a line from the input file and store it in 'line'
        // std::cout << line << std::endl; // Output the read line to the console
    }
    in.close();
    return 0;
}

std::vector<Atom> get_input_coordinates(std::string filename)
{
    std::vector<Atom> input_coordinates;
    Atom atom;

    std::string line;
    std::ifstream in(filename, std::ios_base::in);

    if (!in)
    {
        throw std::runtime_error("Error opening file: " + filename);

    }

    while (std::getline(in, line))
    {
        // Check if the line starts with "frac"
        if (line.compare(0, 4, "frac") == 0 && !line.empty())
        {
            while (std::getline(in, line))
            {
                std::istringstream iss(line);
                Atom atom;

                // Parse the atom data and store it in the Atom struct
                if (iss >> atom.label >> atom.type >> atom.x >> atom.y >> atom.z)
                {
                    input_coordinates.push_back(atom);
                }
                else if (line.empty())
                {
                    break;
                }
                else
                {
                    throw std::runtime_error("Error parsing atom data: " + line);
                }
            }
        }
    }

// Add index to input_coordinates
    int ctn = 0;
    for (auto& atom : input_coordinates)
    {
        ctn += 1;
        atom.index = ctn;
        // std :: cout << atom.index << " " << atom.label << " " << atom.type << " " << atom.x << " " << atom.y << " " << atom.z << std::endl;
    }

    in.close();
    return input_coordinates;

}

std::vector<Specie> get_input_species(std::string filename)
{
    std::vector<Specie> input_species;
    std::string line;
    std::ifstream in(filename, std::ios_base::in);

    if (!in)
    {
        throw std::runtime_error("Error opening file: " + filename);

    }

    while (std::getline(in, line))
    {
        // Check if the line starts with "frac"
        if (line.compare(0, 7, "species") == 0 && !line.empty())
        {
            while (std::getline(in, line))
            {
                std::istringstream iss(line);
                Specie specie;

                // Parse the atom data and store it in the Atom struct
                if (iss >> specie.label >> specie.type >> specie.charge)
                {
                    input_species.push_back(specie);
                }
                else if (line.empty())
                {
                    break;
                }
                else
                {
                    throw std::runtime_error("Error parsing atom species: " + line);
                }
            }
        }
    }
    // for (const auto& specie : input_species)
    // {
    // 	std::cout << specie.label << " " << specie.type << " " << specie.charge << std::endl;
    // }
    in.close();
    return input_species;

}

std::vector<Buckingham> get_input_buckingham(std::string filename)
{
    std::vector<Buckingham> input_buckingham_potentials;
    std::string line;
    std::ifstream in(filename, std::ios_base::in);
    if (!in)
    {
        throw std::runtime_error("Error opening file: " + filename);

    }
    while (std::getline(in, line))
    {
        // Check if the line starts with "frac"
        if (line.compare(0, 10, "buckingham") == 0 && !line.empty())
        {
            while (std::getline(in, line))
            {
                std::istringstream iss(line);
                Buckingham potential;

                // Parse the atom data and store it in the Atom struct
                if (iss >> potential.atom1_label >> potential.atom1_type >> potential.atom2_label >> potential.atom2_type >> potential.A >> potential.rho >> potential.C >> potential.cut_off1 >> potential.cut_off2)
                {
                    input_buckingham_potentials.push_back(potential);
                }
                else if (line.empty())
                {
                    break;
                }
                else
                {
                    throw std::runtime_error("Error parsing atom potentials: " + line);
                }
            }
        }
    }

// for (const auto& potential : input_buckingham_potentials)
// {
// 			std:: cout<<" " <<potential.atom1_label <<" " <<potential.atom1_type << " " <<potential.atom2_label <<" " <<potential.atom2_type <<" " <<potential.A << " " <<potential.rho <<" " <<potential.C <<" " <<potential.cut_off1 <<" " <<potential.cut_off2 << std::endl;
// }

    in.close();
    return input_buckingham_potentials;

}

std::vector<double> get_lattice_constants(std::string filename)
{
    std::string line;
    std::vector<double> lattice_constants;
    std::ifstream in(filename, std::ios_base::in);

    if (!in)
    {
        throw std::runtime_error("Error opening file: " + filename);
    }

    int ctn = 0;

    while (std::getline(in, line))
    {
        if (line.compare(0, 4, "cell") == 0 && !line.empty())
        {
            while (std::getline(in, line) && !line.empty())
            {
                ctn += 1;
            }
            break; // Exit the outer loop once the required lines are counted
        }
    }

    in.clear();  // Clear the EOF flag and error flags
    in.seekg(0); // Move the file pointer to the beginning of the file

    while (std::getline(in, line))
    {
        if (line.compare(0, 4, "cell") == 0 && !line.empty())
        {
            if (ctn == 1)
            {
                while (std::getline(in, line) && !line.empty())
                {
                    std::istringstream iss(line);
                    double val;
                    while (iss >> val)
                    {
                        lattice_constants.push_back(val);
                    }
                }
                in.close();

                return lattice_constants;
            }
            else if (ctn == 3)
            {
                while (std::getline(in, line) && !line.empty())
                {
                    std::istringstream iss(line);
                    double val;
                    while (iss >> val)
                    {
                        lattice_constants.push_back(val);
                    }
                }
                in.close();
                //need fixing
                return lattice_constants;
            }
            else
            {
                in.close();
                throw std::runtime_error("Invalid cell structure input");
            }
        }
    }

    in.close();
    throw std::runtime_error("Cell structure not found in the file");
}
