#include <derv_1.h>
#include <cmath>
#include <Struct_Atom.h>
#include <Eigen/Dense>


Eigen::Vector3d buck_derv1(Atom atom1, Atom atom2, Eigen::Vector3d n, Buckingham pot)
{
	Eigen::Vector3d derv1;
	derv1.setZero();
	Eigen::Vector3d rij;
	Eigen::Vector3d rijn;
	double intact[2];
	rij.setZero();
	rijn.setZero();


	rij << atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z;
	rijn = rij + n;
	intact[0] = (-pot.A/pot.rho) * exp(-rijn.norm()/pot.rho);
	intact[1] = 6. * pot.C / pow(rijn.norm(), 7.);

	if ((pot.atom1_type == atom1.type && pot.atom1_label == atom1.label) &&
		(pot.atom2_type == atom2.type && pot.atom2_label == atom2.label) &&
		rijn.norm() < pot.cut_off2)
	{
		derv1 += (rijn/rijn.norm()) * (intact[0] + intact[1]);
	}
	else if ((pot.atom2_type == atom1.type && pot.atom2_label == atom1.label) &&
			 (pot.atom1_type == atom2.type && pot.atom1_label == atom2.label) &&
			 rijn.norm() < pot.cut_off2)
	{
		derv1 += (rijn/rijn.norm()) * (intact[0] + intact[1]);

	}
return derv1;

}

double buck_derv1_scalar (Atom atom1, Atom atom2, Eigen::Vector3d n, Buckingham pot)
{
	double derv1;
	Eigen::Vector3d rij;
	Eigen::Vector3d rijn;
	double intact[2];
	rij.setZero();
	rijn.setZero();


	rij << atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z;
	rijn = rij + n;
	intact[0] = (-pot.A/pot.rho) * exp(-rijn.norm()/pot.rho);
	intact[1] = 6. * pot.C / pow(rijn.norm(), 7.);

	if ((pot.atom1_type == atom1.type && pot.atom1_label == atom1.label) &&
		(pot.atom2_type == atom2.type && pot.atom2_label == atom2.label) &&
		rijn.norm() < pot.cut_off2)
	{
		derv1 +=  (intact[0] + intact[1]);
	}
	else if ((pot.atom2_type == atom1.type && pot.atom2_label == atom1.label) &&
			 (pot.atom1_type == atom2.type && pot.atom1_label == atom2.label) &&
			 rijn.norm() < pot.cut_off2)
	{
		derv1 +=  (intact[0] + intact[1]);

	}
return derv1;

}

double buck_derv2_scalar (Atom atom1, Atom atom2, Eigen::Vector3d n, Buckingham pot)
{
	double derv2;
	Eigen::Vector3d rij;
	Eigen::Vector3d rijn;
	double intact[2];
	rij.setZero();
	rijn.setZero();


	rij << atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z;
	rijn = rij + n;
	intact[0] = (pot.A/pot.rho/pot.rho) * exp(-rijn.norm()/pot.rho);
	intact[1] = -42. * pot.C / pow(rijn.norm(), 8.);

	if ((pot.atom1_type == atom1.type && pot.atom1_label == atom1.label) &&
		(pot.atom2_type == atom2.type && pot.atom2_label == atom2.label) &&
		rijn.norm() < pot.cut_off2)
	{
		derv2 +=  (intact[0] + intact[1]);
	}
	else if ((pot.atom2_type == atom1.type && pot.atom2_label == atom1.label) &&
			 (pot.atom1_type == atom2.type && pot.atom1_label == atom2.label) &&
			 rijn.norm() < pot.cut_off2)
	{
		derv2 +=  (intact[0] + intact[1]);

	}
return derv2;

}