#include <derv_2.h>
#include <cmath>
#include <Struct_Atom.h>
#include <Eigen/Dense>

double buck_derv2(Atom atom1, Atom atom2, Eigen::Vector3d n, Buckingham pot);
{
	
	Eigen::Vector3d rij, rijn;
	rij.setZero();
	rijn = rij + n;
	double intact[5];
	intact[0] = -1./rijn.norm() + 
}
