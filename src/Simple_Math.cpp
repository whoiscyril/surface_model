#include "Simple_Math.h"
#include <math.h>
#include <cmath>
#include<string>
#include <Eigen/Dense>
#include <vector>
#include <iostream>


std::tuple<double, double, double> rad2deg(const std::tuple<double, double, double>& angles)
{
    double roll_r, pitch_r, yaw_r;
    std::tie(roll_r, pitch_r, yaw_r) = angles;
    double roll_d = roll_r * 180.0 / M_PI;
    double pitch_d = pitch_r * 180.0 / M_PI;
    double yaw_d = yaw_r * 180.0 / M_PI;
    return std::make_tuple(roll_d, pitch_d, yaw_d);
}


double get_distance(const Atom& atom1, const Atom& atom2)
{
    double dx = atom1.x - atom2.x;
    double dy = atom1.y - atom2.y;
    double dz = atom1.z - atom2.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

double absolute_sum(const std::vector<double>& vec)
{
    double sum = 0;
    for (const double& val : vec)
    {
        sum += std::abs(val);
    }
    return sum;
}

double absolute_sum(const Eigen::Vector3d& vec)
{
    return std::abs(vec[0]) + std::abs(vec[1]) + std::abs(vec[2]);
}

double vmag(const std::vector<double> v)
{
    double sum=0.0;
    for (const auto& val : v)
    {
        sum += pow(val,2);
    }
    return sqrt(sum);
}
double vmag(const Eigen::Vector3d v)
{
    return v.norm();
}
double len(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}
double dp(double x, double y, double z, double h, double k, double l)
{
    return x*h+y*k+z*l;
}

Eigen::Matrix3d frac2cart(const std::vector<Atom>& cell, Eigen::Matrix3d& lvecs)
{
    //first calculate angles
    Eigen::MatrixXd transformation_matrix(3,3);
    Eigen::VectorXd temp(9);
    Eigen::Vector3d a1, a2, a3;
    double alpha, beta, gamma;
    double reci_alpha;
    a1.setZero();
    a2.setZero();
    a3.setZero();

    a1 = lvecs.row(0);
    a2 = lvecs.row(1);
    a3 = lvecs.row(2);

    alpha = acos(a2.normalized().dot(a3.normalized()));
    beta = acos(a1.normalized().dot(a3.normalized()));
    gamma = acos(a1.normalized().dot(a2.normalized()));

    reci_alpha = acos( (cos(beta) * cos(gamma) - cos(alpha)) / (sin(beta) * sin(gamma)));

    transformation_matrix <<
                          a1.norm(), a2.norm() * std::cos(gamma), a3.norm() * std::cos(beta),
                                  0.0, a2.norm() * std::sin(gamma), -a3.norm() * std::sin(beta) * std::cos(reci_alpha),
                                  0.0, 0.0, a3.norm() * std::sin(beta) * std::sin(reci_alpha);





    // std::cout << a1.norm() << " " << a2.norm() << " " << a3.norm()<< std::endl;
    // std::cout << alpha * 180. / M_PI<< " " << beta * 180. / M_PI<< " " << gamma * 180./M_PI<< std::endl;

    // std::cout << transformation_matrix << std::endl;
    return transformation_matrix;

}
//function that takes in lattice constants and angles and return lattice vectors;

Eigen::Matrix3d lattice_vectors(Eigen::Vector3d& lattice_constants, Eigen::Vector3d& lattice_angles)
{
    Eigen::Matrix3d lattice_vectors;

    double cos_alpha = std::cos(lattice_angles[0] * M_PI / 180.0);
    double cos_beta = std::cos(lattice_angles[1] * M_PI / 180.0);
    double cos_gamma = std::cos(lattice_angles[2] * M_PI / 180.0);
    double sin_gamma = std::sin(lattice_angles[2] * M_PI / 180.0);

    lattice_vectors << lattice_constants[0], 0, 0,
                    lattice_constants[1] * cos_gamma, lattice_constants[1] * sin_gamma, 0,
                    lattice_constants[2] * cos_beta, lattice_constants[2] * (cos_alpha - cos_beta * cos_gamma) / sin_gamma,
                    lattice_constants[2] * std::sqrt(1.0 - cos_alpha * cos_alpha - cos_beta * cos_beta - cos_gamma * cos_gamma +
                            2.0 * cos_alpha * cos_beta * cos_gamma) / sin_gamma;


    return lattice_vectors;
}
