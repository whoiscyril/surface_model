#ifndef STRUCT_ATOM_H
#define STRUCT_ATOM_H
#include <string>
#include <tuple>
#include <Eigen/Dense>


struct Atom
{
    int value;
    int index;
    std::string label;
    std::string type;
    double x, y, z, q, o;
    double rad;
    double dp_x, dp_y, dp_z;

    bool operator<(const Atom& other) const
    {
        return std::tie(label, x, y, z) < std::tie(other.label, other.x, other.y, other.z);
    }
};
struct Buckingham
{
    std::string atom1_label;
    std::string atom1_type;
    std::string atom2_label;
    std::string atom2_type;
    double A;
    double rho;
    double C;
    double cut_off1;
    double cut_off2;
};

struct UnitCell
{
    Eigen::Vector3d lattice_constants;
    Eigen::Vector3d lattice_angles;
    std::vector<Atom> coordinates_cart;
    std::vector<Atom> coordinates_frac;
    Eigen::Matrix3d lattice_vectors;
    Eigen::Matrix3d reciprocal_vectors;

    std::vector<Buckingham> buckingham_potentials;
    double volume;

    UnitCell(const std::string& filename);
};

struct Specie
{
    std::string label;
    std::string type;
    double charge;
};


#endif // 