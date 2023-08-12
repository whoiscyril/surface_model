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