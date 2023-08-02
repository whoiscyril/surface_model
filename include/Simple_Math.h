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

#endif // MATRIX_FUNCTIONS_H	