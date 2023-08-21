#ifndef SIMPLE_MATH_H
#define SIMPLE_MATH_H
#include <string>
#include <math.h>
#include "Struct_Atom.h"
#include <Eigen/Dense>



std::tuple<double, double, double> rad2deg(const std::tuple<double, double, double>& angles);
double absolute_sum(const std::vector<double>& vec);
double absolute_sum(const Eigen::Vector3d& vec);
double get_distance(const Atom& atom1, const Atom& atom2);
double vmag(const std::vector<double> v);
double vmag(const Eigen::Vector3d v);
double rad2deg(double radian);
double len(double x, double y, double z);
double dp(double x, double y, double z, double h, double k, double l);
Eigen::Matrix3d frac2cart(const std::vector<Atom>& cell, Eigen::Matrix3d& lvecs);
Eigen::Matrix3d lattice_vectors(Eigen::Vector3d& lattice_constants, Eigen::Vector3d& lattice_angles);
#endif // MATRIX_FUNCTIONS_H	