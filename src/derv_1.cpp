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

// double buck_derv1_scalar (Atom atom1, Atom atom2, Eigen::Vector3d n, Buckingham pot)
// {
// 	double derv1;
// 	Eigen::Vector3d rij;
// 	Eigen::Vector3d rijn;
// 	double intact[2];
// 	rij.setZero();
// 	rijn.setZero();

// 	rij << atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z;
// 	rijn = rij + n;
// 	intact[0] = (-pot.A/pot.rho) * exp(-rijn.norm()/pot.rho);
// 	intact[1] = 6. * pot.C / pow(rijn.norm(), 7.);

// 	if ((pot.atom1_type == atom1.type && pot.atom1_label == atom1.label) &&
// 		(pot.atom2_type == atom2.type && pot.atom2_label == atom2.label) &&
// 		rijn.norm() < pot.cut_off2)
// 	{
// 		derv1 +=  (intact[0] + intact[1]);
// 	}
// 	else if ((pot.atom2_type == atom1.type && pot.atom2_label == atom1.label) &&
// 			 (pot.atom1_type == atom2.type && pot.atom1_label == atom2.label) &&
// 			 rijn.norm() < pot.cut_off2)
// 	{
// 		derv1 +=  (intact[0] + intact[1]);

// 	}
// return derv1;

// }

double buck_derv1_scalar (Atom atom1, Atom atom2, UnitCell unitcell_init)
{
    double derv1;
    std::vector<Buckingham> buck = unitcell_init.buckingham_potentials;
    std::vector<Atom> atoms = unitcell_init.coordinates_cart;
    Eigen::Matrix3d lattice_vectors = unitcell_init.lattice_vectors;
    Eigen::Vector3d n;
    Eigen::Vector3d rij;
    Eigen::Vector3d rijn;
    double intact[2];
    rij.setZero();
    rijn.setZero();
    n.setZero();
    Eigen::Vector3d v1 = lattice_vectors.row(0);
    Eigen::Vector3d v2 = lattice_vectors.row(1);
    Eigen::Vector3d v3 = lattice_vectors.row(2);

    double cutoff;
    for (const auto& elem : buck)
    {
        if (elem.cut_off2 > cutoff)
        {
            cutoff = elem.cut_off2;
        }
        else
        {
            continue;
        }
    }

    int maxx = ceil(cutoff/v1.norm());
    int maxy = ceil(cutoff/v2.norm());
    int maxz = ceil(cutoff/v3.norm());

    rij << atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z;

    for (int i = -maxx; i <= maxx; ++i)
    {
        for (int j = -maxy; j <= maxy; ++j)
        {
            for (int k = -maxz; k <= maxz; ++k)
            {
                n = i * v1 + j * v2 + k * v3;
                rijn = rij + n;

                if (i ==0 && j == 0 && k == 0)
                {
                    for (const auto pot : buck)
                    {
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
                    }
                }
                else
                {
                    for (const auto pot : buck)
                    {
                        {
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
                        }
                    }
                }
            }
        }
    }

    return derv1;
}

double buck_derv2_scalar (Atom atom1, Atom atom2, UnitCell unitcell_init)
{
    double derv2;
    std::vector<Buckingham> buck = unitcell_init.buckingham_potentials;
    std::vector<Atom> atoms = unitcell_init.coordinates_cart;
    Eigen::Matrix3d lattice_vectors = unitcell_init.lattice_vectors;
    Eigen::Vector3d n;
    Eigen::Vector3d rij;
    Eigen::Vector3d rijn;
    double intact[2];
    rij.setZero();
    rijn.setZero();
    n.setZero();
    Eigen::Vector3d v1 = lattice_vectors.row(0);
    Eigen::Vector3d v2 = lattice_vectors.row(1);
    Eigen::Vector3d v3 = lattice_vectors.row(2);

    double cutoff;
    for (const auto& elem : buck)
    {
        if (elem.cut_off2 > cutoff)
        {
            cutoff = elem.cut_off2;
        }
        else
        {
            continue;
        }
    }

    int maxx = ceil(cutoff/v1.norm());
    int maxy = ceil(cutoff/v2.norm());
    int maxz = ceil(cutoff/v3.norm());

    rij << atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z;

    for (int i = -maxx; i <= maxx; ++i)
    {
        for (int j = -maxy; j <= maxy; ++j)
        {
            for (int k = -maxz; k <= maxz; ++k)
            {
                n = i * v1 + j * v2 + k * v3;
                rijn = rij + n;

                if (i ==0 && j == 0 && k == 0)
                {
                    for (const auto pot : buck)
                    {
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
                    }
                }
                else
                {
                    for (const auto pot : buck)
                    {
                        {
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
                        }
                    }
                }
            }
        }
    }

    return derv2;
}

double derv1_scalar(Atom atom1, Atom atom2, Eigen::Vector3d n, Buckingham pot)
{
	double result;
	Eigen::Vector3d rij, rijn;
	rij.setZero();
	rijn.setZero();

	rij << atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z;
	rijn = rij - n;

	result = -pot.A/pot.rho * exp(-rijn.norm() / pot.rho) + (6.* pot.C / pow(rijn.norm(), 7.));

	return result;

}
double derv2_scalar(Atom atom1, Atom atom2, Eigen::Vector3d n, Buckingham pot)
{
	double result;
	Eigen::Vector3d rij, rijn;
	rij.setZero();
	rijn.setZero();
	
	rij << atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z;
	rijn = rij - n;

	result = (pot.A/pot.rho / pot.rho) * exp(-rijn.norm() / pot.rho) - (42.* pot.C / pow(rijn.norm(), 8.));

	return result;

}