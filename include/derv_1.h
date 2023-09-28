#ifndef DERV_1_H
#define DERV_1_H
#include <cmath>
#include <Struct_Atom.h>
#include <Eigen/Dense>
Eigen::Vector3d buck_derv1(Atom atom1, Atom atom2, Eigen::Vector3d n, Buckingham pot);
double buck_derv1_scalar (Atom atom1, Atom atom2, UnitCell unitcell_init);
double buck_derv2_scalar (Atom atom1, Atom atom2, UnitCell unitcell_init);
double derv1_scalar(Atom atom1, Atom atom2, Eigen::Vector3d n, Buckingham pot);
double derv2_scalar(Atom atom1, Atom atom2, Eigen::Vector3d n, Buckingham pot);

#endif