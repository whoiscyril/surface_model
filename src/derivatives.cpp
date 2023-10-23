#include <cmath>
#include <Energy.h>
#include <iostream>
#include <Struct_Atom.h>
#include <derv_1.h>

void internal_derv2(UnitCell& unitcell_init)
{

    int natom = unitcell_init.coordinates_cart.size();
    std::vector<Atom> atoms = unitcell_init.coordinates_cart;
    double cutoff;
    //Max buck cut off
    std::vector<Buckingham> buck = unitcell_init.buckingham_potentials;
    Eigen::Matrix3d lattice_vectors = unitcell_init.lattice_vectors;
    Eigen::Vector3d v1, v2, v3;
    Eigen::Vector3d n;
    Eigen::Vector3d total_buck_derv1;
    Eigen::MatrixXd derv2(3*natom, 3*natom);
    derv2.setZero();
    n.setZero();
    v1 = lattice_vectors.row(0);
    v2 = lattice_vectors.row(1);
    v3 = lattice_vectors.row(2);
    Eigen::Vector3d a1 = unitcell_init.lattice_vectors.row(0);
    Eigen::Vector3d a2 = unitcell_init.lattice_vectors.row(1);
    Eigen::Vector3d a3 = unitcell_init.lattice_vectors.row(2);

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

    //Max images in 3 directions;

    double xmax = ceil(cutoff/v1.norm());
    double ymax = ceil(cutoff/v2.norm());
    double zmax = ceil(cutoff/v3.norm());

    Eigen::Vector3d rij, rijn;

    //First put coordinates into row and column vectors
    std::vector<double> row_coordinates;
    std::vector<double> col_coordinates;

    for (const auto & elem : atoms)
    {
        row_coordinates.push_back(elem.x);
        row_coordinates.push_back(elem.y);
        row_coordinates.push_back(elem.z);
    }

    //Now compute distance and put it in a matrix
    Eigen::MatrixXd dist_mat(3*natom, 3*natom);
    dist_mat.setZero();

    for (int i = 0; i < 3*natom; ++i)
    {
        for (int j = 0; j < 3*natom ; ++j)
        {
            dist_mat(j, i) = row_coordinates[i] - row_coordinates[j];
        }
    }
    // std::cout << dist_mat << std::endl;

//Rewrite second derivatives;
    Eigen::Matrix3d temp, temp_diag;
    temp.setZero();
    temp_diag.setZero();
    for (int i = 0; i <atoms.size(); ++i)
    {
        Atom elem1 = atoms[i];
        temp.setZero();
        temp_diag.setZero();
        for (int j = 0; j < atoms.size(); ++j)
        {
            Atom elem2 =atoms[j];
            Eigen::Vector3d rij;
            rij << elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z;
            double r_norm = rij.norm();
            double r_sqr = r_norm * r_norm;
            double temp_derv_1 = 0.;
            double temp_derv_2 = 0.;

            for (int ii = -xmax; ii <= xmax; ++ii)
            {
                for (int jj = -ymax; jj <= ymax; ++jj)
                {
                    for (int kk = -zmax; kk <= zmax; ++kk)
                    {
                        n = ii * a1 + jj * a2 + kk * a3;
                        rijn = rij + n;
                        r_norm = rijn.norm();
                        r_sqr = r_norm * r_norm;
                        if (ii == 0 && jj == 0 && kk == 0)
                        {

                            if (i != j)
                            {
                                for (const auto pot : buck)
                                {
                                    temp_derv_1 = derv1_scalar(elem1,elem2,n,pot);
                                    temp_derv_2 = derv2_scalar(elem1,elem2,n,pot);
                                }
                                for (int ia = 0; ia < 3; ++ia)
                                {
                                    for (int jb = 0; jb < 3; ++jb)
                                    {
                                        if (ia == jb)
                                        {
                                            temp(ia,jb) += (-1. * (temp_derv_1 / r_norm) - (((rijn[ia] * rijn[jb]) / r_sqr) * (temp_derv_2 - temp_derv_1 / r_norm))) ;
                                        }
                                        else
                                        {
                                            temp(ia,jb) += (-0. * (temp_derv_1 / r_norm) - (((rijn[ia] * rijn[jb]) / r_sqr) * (temp_derv_2 - temp_derv_1 / r_norm))) ;
                                        }
                                    }
                                }
                            }
                            else if (i == j)
                            {

                                for (int k = 0; k < atoms.size(); ++k)
                                {
                                    Atom elem3 = atoms[k];
                                    Eigen::Vector3d rik, rikn;
                                    rik << elem1.x - elem3.x, elem1.y - elem3.y, elem1.z - elem3.z;
                                    rikn = rik + n;
                                    double rk_norm = rikn.norm();
                                    double rk_sqr = rk_norm * rk_norm;
                                    for (const auto pot : buck)
                                    {
                                        temp_derv_1 = derv1_scalar(elem1,elem3,n,pot);
                                        temp_derv_2 = derv2_scalar(elem1,elem3,n,pot);
                                    }
                                    if (i != k)
                                    {
                                        for (int ia = 0; ia < 3; ++ia)
                                        {
                                            for (int jb = 0; jb < 3; ++jb)
                                            {
                                                if (ia == jb)
                                                {
                                                    temp_diag(ia,jb) += (1. * temp_derv_1 / rk_norm) + ((rikn[ia] * rikn[jb]) / rk_sqr) * (temp_derv_2 - temp_derv_1 / rk_norm);

                                                }
                                                else
                                                {
                                                    temp_diag(ia,jb) += (0. * temp_derv_1 / rk_norm) + ((rikn[ia] * rikn[jb]) / rk_sqr) * (temp_derv_2 - temp_derv_1 / rk_norm);

                                                }
                                            }
                                        }
                                    }
                                }

                            }

                        }
                        else
                        {

                            if (i != j)
                            {
                                for (const auto pot : buck)
                                {
                                    temp_derv_1 =   derv1_scalar(elem1, elem2, n, pot);
                                    temp_derv_2 =   derv2_scalar(elem1, elem2, n, pot);
                                }
                                for (int ia = 0; ia < 3; ++ia)
                                {
                                    for (int jb = 0; jb < 3; ++jb)
                                    {
                                        if (ia == jb)
                                        {
                                            temp(ia,jb) +=  (-1. * (temp_derv_1 / r_norm) - (((rijn[ia] * rijn[jb]) / r_sqr) * (temp_derv_2 - temp_derv_1 / r_norm))) ;
                                        }
                                        else
                                        {
                                            temp(ia,jb) +=  (-0. * (temp_derv_1 / r_norm) - (((rijn[ia] * rijn[jb]) / r_sqr) * (temp_derv_2 - temp_derv_1 / r_norm))) ;
                                        }
                                    }
                                }
                            }
                            else if ( i == j )
                            {

                                for (int k = 0; k < atoms.size(); ++k)
                                {
                                    Atom elem3 = atoms[k];
                                    Eigen::Vector3d rik, rikn;
                                    rik << elem1.x - elem3.x, elem1.y - elem3.y, elem1.z - elem3.z;
                                    rikn = rik + n;
                                    double rk_norm = rikn.norm();
                                    double rk_sqr = rk_norm * rk_norm;
                                    for (const auto pot : buck)
                                    {
                                        temp_derv_1 =   derv1_scalar(elem1, elem3, n, pot);
                                        temp_derv_2 =   derv2_scalar(elem1, elem3, n, pot);
                                    }
                                    if (i != k)
                                    {
                                        for (int ia = 0 ; ia < 3; ++ia)
                                        {
                                            for (int jb = 0; jb < 3; ++jb)
                                            {
                                                if (ia == jb)
                                                {
                                                    temp_diag(ia,jb) += (1. * temp_derv_1 / rk_norm) + ((rikn[ia] * rikn[jb]) / rk_sqr) * (temp_derv_2 - temp_derv_1 / rk_norm);
                                                }
                                                else
                                                {
                                                    temp_diag(ia,jb) += (0. * temp_derv_1 / rk_norm) + ((rikn[ia] * rikn[jb]) / rk_sqr) * (temp_derv_2 - temp_derv_1 / rk_norm);

                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // temp_derv_1 = buck_derv1_scalar(elem1, elem2, unitcell_init);
            // temp_derv_2 = buck_derv2_scalar(elem1, elem2, unitcell_init);
            // // std::cout << temp_derv_1 << std::endl;

            // if(j != i)
            // {
            //     for (int ii = 0; ii < 3; ++ii)
            //     {
            //         for (int jj = 0; jj < 3; ++jj)
            //         {
            //             if (ii == jj)
            //             {
            //                 temp(ii,jj) = -1. * (temp_derv_1 / r_norm) - (((rij[ii] * rij[jj]) / r_sqr) * (temp_derv_2 - temp_derv_1 / r_norm));
            //             }
            //             else
            //             {
            //                 temp(ii,jj) = 0. * (temp_derv_1 / r_norm) - (((rij[ii] * rij[jj]) / r_sqr) * (temp_derv_2 - temp_derv_1 / r_norm));
            //             }
            //             // std::cout << rij[ii] << " " << rij[jj] << std::endl;
            //         }
            //     }
            // }
        }
        // std::cout << temp_diag << std::endl;
    }

//Second derivatives;
    // double temp_derv_1;
    // double temp_derv_2;
    // double result[atoms.size()][3][atoms.size()][3];
    // Eigen::MatrixXd temp(3, 3);
    // temp.setZero();
    // derv2.setZero();
    // for (int i = 0; i <atoms.size(); ++i)
    // {
    //     Atom elem1 = atoms[i];
    //     Eigen::Vector3d el1_coord;
    //     el1_coord << elem1.x, elem1.y, elem1.z;
    //     for (int j = 0; j < atoms.size(); ++j)
    //     {
    //         Atom elem2 = atoms[j];
    //         Eigen::Vector3d el2_coord;
    //         el2_coord << elem2.x, elem2.y, elem2.z;
    //         rij << elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z;
    //         for (int ii = -xmax; ii <= xmax; ++ii)
    //         {
    //             for (int jj = -ymax; jj <= xmax; ++jj)
    //             {
    //                 for (int kk = -zmax; kk <= zmax; ++kk)
    //                 {
    //                     n = ii * a1 + jj * a2 + kk * a3;
    //                     rijn = rij - n;
    //                     double r_norm = rijn.norm();
    //                     double r_sqr = r_norm * r_norm;
    //                     //Central image
    //                     if (ii == 0 && jj == 0 && kk == 0)
    //                     {
    //                         //Calculate psi'
    //                         for (const auto& pot : buck)
    //                         {
    //                             if ((pot.atom1_type == elem1.type && pot.atom1_label == elem1.label) &&
    //                                     (pot.atom2_type == elem2.type && pot.atom2_label == elem2.label) &&
    //                                     rijn.norm() < pot.cut_off2)
    //                             {
    //                                 temp_derv_1 += derv1_scalar(elem1, elem2, n, pot);
    //                                 temp_derv_2 += derv2_scalar(elem1, elem2, n, pot);
    //                             }
    //                             else if ((pot.atom2_type == elem1.type && pot.atom2_label == elem1.label) &&
    //                                      (pot.atom2_type == elem2.type && pot.atom1_label == elem2.label) &&
    //                                      rijn.norm() < pot.cut_off2)
    //                             {
    //                                 temp_derv_1 += derv1_scalar(elem1, elem2, n, pot);
    //                                 temp_derv_2 += derv2_scalar(elem1, elem2, n, pot);
    //                             }

    //                         }
    //                     }

    //                     //other images;
    //                     else
    //                     {
    //                         for (const auto& pot : buck)
    //                         {
    //                             if ((pot.atom1_type == elem1.type && pot.atom1_label == elem1.label) &&
    //                                     (pot.atom2_type == elem2.type && pot.atom2_label == elem2.label) &&
    //                                     rijn.norm() < pot.cut_off2)
    //                             {
    //                                 temp_derv_1 += derv1_scalar(elem1, elem2, n, pot);
    //                                 temp_derv_2 += derv2_scalar(elem1, elem2, n, pot);
    //                             }
    //                             else if ((pot.atom2_type == elem1.type && pot.atom2_label == elem1.label) &&
    //                                      (pot.atom2_type == elem2.type && pot.atom1_label == elem2.label) &&
    //                                      rijn.norm() < pot.cut_off2)
    //                             {
    //                                 temp_derv_1 += derv1_scalar(elem1, elem2, n, pot);
    //                                 temp_derv_2 += derv2_scalar(elem1, elem2, n, pot);
    //                             }

    //                         }
    //                     }

    //                     if (j != i)
    //                     {
    //                         for (int ia = 0; ia < 3; ++ia)
    //                         {
    //                             for (int jb = 0; jb < 3; ++jb)
    //                             {
    //                                 if (ia == jb)
    //                                 {
    //                                     // temp(ia,jb) = -1. * (temp_derv_1 / r_norm) - (((rijn[ia] * rijn[jb]) / r_sqr) * (temp_derv_2 - temp_derv_1 / r_norm));
    //                                 }
    //                                 else
    //                                 {
    //                                     temp(ia,jb) = -(((rijn[ia] * rijn[jb]) / r_sqr) * (temp_derv_2 - temp_derv_1 / r_norm));
    //                                                                             std::cout << temp << std::endl;
    //                                 }
    //                             }
    //                         }
    //                     }
    //                     // else
    //                     // {
    //                     //     for (int k = 0; k < atoms.size(); ++k)
    //                     //     {
    //                     //         Atom elem3 = atoms[k];
    //                     //         Eigen::Vector3d rik;
    //                     //         Eigen::Vector3d rikn;
    //                     //         rik << elem1.x - elem3.x, elem1.y - elem3.y, elem1.z - elem3.z;
    //                     //         for (int kx = -xmax; kx <= xmax; ++kx)
    //                     //         {
    //                     //             for (int ky = -ymax; ky <= ymax; ++ky)
    //                     //             {
    //                     //                 for (int kz = -zmax; kz <= zmax; ++kz)
    //                     //                 {
    //                     //                     Eigen::Vector3d nk;
    //                     //                     nk = kx * a1 + ky * a2 + kz * a3;
    //                     //                     rikn = rik - nk;
    //                     //                     double rk_norm = rikn.norm();
    //                     //                     double rk_sqr = rk_norm * rk_norm;
    //                     //                     if (kx == 0 && ky == 0 && kz == 0)
    //                     //                     {

    //                     //                     }
    //                     //                     else
    //                     //                     {
    //                     //                     if (k != i)
    //                     //                     {
    //                     //                         for (int ia = 0; ia < 3; ++ia)
    //                     //                         {
    //                     //                             for (int jb = 0; jb < 3; ++jb)
    //                     //                             {
    //                     //                                 if (ia == jb)
    //                     //                                 {
    //                     //                                     temp(ia,jb) += -1. * (temp_derv_1 / rk_norm) - (((rikn[ia] * rikn[jb]) / rk_sqr) * (temp_derv_2 - temp_derv_1 / rk_norm));
    //                     //                                 }
    //                     //                                 else
    //                     //                                 {
    //                     //                                     temp(ia,jb) += 0. * (temp_derv_1 / rk_norm) - (((rikn[ia] * rikn[jb]) / rk_sqr) * (temp_derv_2 - temp_derv_1 / rk_norm));
    //                     //                                 }
    //                     //                             }
    //                     //                         }
    //                     //                     }
    //                     //                 }

    //                     //                 }
    //                     //             }
    //                     //         }

    //                     //     }
    //                     // }
    //                 }
    //             }
    //         }
    //     }
    // }
    std::cout << temp_diag << std::endl;
//First derivative
    for (int i = 0; i <atoms.size(); ++i)
    {
        Atom elem1 = atoms[i];
        total_buck_derv1.setZero();
        for (int j = 0; j<atoms.size(); ++j)
        {
            Atom elem2 = atoms[j];
            for (int ii = -xmax; ii <= xmax; ++ii)
            {
                for (int jj = -ymax; jj <= ymax; ++jj)
                {
                    for (int kk = -zmax; kk <= zmax ; ++kk)
                    {
                        n = ii * a1 + jj * a2 + kk * a3;

                        if (ii == 0 && jj == 0 && kk == 0)
                        {
                            for (const auto pot : buck)
                            {
                                total_buck_derv1 += buck_derv1(elem1, elem2, n, pot);
                            }
                        }
                        else
                        {
                            for (const auto pot : buck)
                            {
                                total_buck_derv1 += buck_derv1(elem1, elem2, n, pot);
                            }
                        }
                    }
                }
            }
        }
        total_buck_derv1 = lattice_vectors * total_buck_derv1;
    }
}

void calc_lattice_deriv(UnitCell& unitcell_init)
{
    Eigen::Matrix3d strain_deriv = unitcell_init.strain_deriv;
    Eigen::Matrix3d lattice_deriv;
    lattice_deriv.setZero();
    Eigen::Matrix3d lattice_vectors = unitcell_init.lattice_vectors;
    Eigen::Matrix3d lattice_vectors_inverse = lattice_vectors.inverse();

    for (int i=0; i < 3; ++i)
    {

        lattice_deriv.col(0)(0) += strain_deriv.col(0)(i) * lattice_vectors_inverse.col(0)(i);
        lattice_deriv.col(0)(1) += strain_deriv.col(1)(i) * lattice_vectors_inverse.col(0)(i);
        lattice_deriv.col(0)(2) += strain_deriv.col(2)(i) * lattice_vectors_inverse.col(0)(i);

        lattice_deriv.col(1)(0) += strain_deriv.col(0)(i) * lattice_vectors_inverse.col(1)(i);
        lattice_deriv.col(1)(1) += strain_deriv.col(1)(i) * lattice_vectors_inverse.col(1)(i);
        lattice_deriv.col(1)(2) += strain_deriv.col(2)(i) * lattice_vectors_inverse.col(1)(i);

        lattice_deriv.col(2)(0) += strain_deriv.col(0)(i) * lattice_vectors_inverse.col(2)(i);
        lattice_deriv.col(2)(1) += strain_deriv.col(1)(i) * lattice_vectors_inverse.col(2)(i);
        lattice_deriv.col(2)(2) += strain_deriv.col(2)(i) * lattice_vectors_inverse.col(2)(i);

    }
    unitcell_init.lattice_deriv = lattice_deriv;

}

void calc_strain_deriv(UnitCell& unitcell_init)
{
    UnitCell unitcell_wstrain;
    Eigen::Matrix3d deriv;
    deriv.setZero();

    // //Now populate strain derivatives due to buckingham potentials;

    std::vector<Atom> atoms = unitcell_init.coordinates_cart;
    Eigen::Vector3d lattice_constants = unitcell_init.lattice_constants;
    Eigen::Matrix3d lattice_vectors = unitcell_init.lattice_vectors;
    Eigen::Vector3d v1, v2, v3;
    v1 = lattice_vectors.row(0);
    v2 = lattice_vectors.row(1);
    v3 = lattice_vectors.row(2);
    Eigen::Vector3d a1 = unitcell_init.lattice_vectors.row(0);
    Eigen::Vector3d a2 = unitcell_init.lattice_vectors.row(1);
    Eigen::Vector3d a3 = unitcell_init.lattice_vectors.row(2);
    double cutoff = 0;
    std::vector<Buckingham> buck = unitcell_init.buckingham_potentials;

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

    Eigen::Vector3d rij;
    Eigen::Vector3d n;
    Eigen::Vector3d rijn;

    for (auto& elem1 : atoms)
    {
        for (auto& elem2 : atoms)
        {
            rij << elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z;
            // std::cout << rij.norm() << std::endl;
            for (int i = -maxx; i <= maxx; ++i)
            {
                for (int j = -maxy; j <= maxy; ++j)
                {
                    for (int k = -maxz; k <= maxz; ++k)
                    {
                        n = i * a1 + j * a2+ k * a3;
                        rijn = rij + n;

                        if (i == 0 && j ==0 && k == 0)
                        {
                            for (const auto& pot : unitcell_init.buckingham_potentials)
                            {
                                double expterm = (-pot.A / pot.rho) * exp(-rij.norm()/pot.rho) + 6. * pot.C / pow(rij.norm(), 7.);
                                // std::cout << expterm << std::endl;
                                if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                                        ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                                        rij.norm() <= pot.cut_off2 && rij.norm() > 0.1)
                                {
                                    deriv.col(0)(0) += rij[0] * rij[0] * expterm / rij.norm();
                                    deriv.col(1)(0) += rij[0] * rij[1] * expterm / rij.norm();
                                    deriv.col(1)(1) += rij[1] * rij[1] * expterm / rij.norm();
                                    deriv.col(2)(0) += rij[0] * rij[2] * expterm / rij.norm();
                                    deriv.col(2)(1) += rij[1] * rij[2] * expterm / rij.norm();
                                    deriv.col(2)(2) += rij[2] * rij[2] * expterm / rij.norm();

                                }
                                else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                                         ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                                         rij.norm() <= pot.cut_off2 && rij.norm() > 0.1)
                                {
                                    deriv.col(0)(0) += rij[0] * rij[0] * expterm / rij.norm();
                                    deriv.col(1)(0) += rij[0] * rij[1] * expterm / rij.norm();
                                    deriv.col(1)(1) += rij[1] * rij[1] * expterm / rij.norm();
                                    deriv.col(2)(0) += rij[0] * rij[2] * expterm / rij.norm();
                                    deriv.col(2)(1) += rij[1] * rij[2] * expterm / rij.norm();
                                    deriv.col(2)(2) += rij[2] * rij[2] * expterm / rij.norm();
                                }
                            }
                        }
                        else
                        {
                            for (const auto pot : unitcell_init.buckingham_potentials)
                            {
                                double expterm = (-pot.A / pot.rho) * exp(-rijn.norm()/pot.rho) + 6. * pot.C / pow(rijn.norm(), 7.);
                                if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                                        ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                                        rijn.norm() <= pot.cut_off2)
                                {
                                    deriv.col(0)(0) += rijn[0] * rijn[0] * expterm / rijn.norm();
                                    deriv.col(1)(0) += rijn[0] * rijn[1] * expterm / rijn.norm();
                                    deriv.col(1)(1) += rijn[1] * rijn[1] * expterm / rijn.norm();
                                    deriv.col(2)(0) += rijn[0] * rijn[2] * expterm / rijn.norm();
                                    deriv.col(2)(1) += rijn[1] * rijn[2] * expterm / rijn.norm();
                                    deriv.col(2)(2) += rijn[2] * rijn[2] * expterm / rijn.norm();
                                }
                                else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                                         ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                                         rijn.norm() <= pot.cut_off2)
                                {
                                    deriv.col(0)(0) += rijn[0] * rijn[0] * expterm / rijn.norm();
                                    deriv.col(1)(0) += rijn[0] * rijn[1] * expterm / rijn.norm();
                                    deriv.col(1)(1) += rijn[1] * rijn[1] * expterm / rijn.norm();
                                    deriv.col(2)(0) += rijn[0] * rijn[2] * expterm / rijn.norm();
                                    deriv.col(2)(1) += rijn[1] * rijn[2] * expterm / rijn.norm();
                                    deriv.col(2)(2) += rijn[2] * rijn[2] * expterm / rijn.norm();
                                }
                            }
                        }
                    }

                }
            }
        }
    }
    deriv = deriv * 0.5;
    // std::cout << deriv << std::endl;

    //Real contribution
    const double toeV = 14.39964390675221758120;

    int natom = unitcell_init.coordinates_frac.size();

    Eigen::Vector3d g1 = unitcell_init.reciprocal_vectors.row(0);
    Eigen::Vector3d g2 = unitcell_init.reciprocal_vectors.row(1);
    Eigen::Vector3d g3 = unitcell_init.reciprocal_vectors.row(2);

    double V = unitcell_init.volume;

    long double kappa = pow(((natom * 1. * M_PI*M_PI*M_PI)/ V / V), 1./6.);
    double rcut = pow( -log(10E-17)/ (kappa * kappa),1./2.);
    double kcut = 2.*kappa*sqrt((-log(10E-17)));

    int nmax_x = ceil(rcut/a1.norm());
    int nmax_y = ceil(rcut/a2.norm());
    int nmax_z = ceil(rcut/a3.norm());
    int kmax_x = ceil(kcut/g1.norm());
    int kmax_y = ceil(kcut/g2.norm());
    int kmax_z = ceil(kcut/g3.norm());

    // double kappa = 1./2.39958;
    // double rcut = 15.013;
    // double kcut = 5.21468;

    // int nmax_x = 3;
    // int nmax_y = 3;
    // int nmax_z = 2;
    // int kmax_x = 2;
    // int kmax_y = 4;
    // int kmax_z = 4;
    rij.setZero();
    n.setZero();
    rijn.setZero();

    for ( auto& elem1 : atoms)
    {
        for ( auto& elem2 : atoms)
        {
            rij << elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z;
            for (int i = -nmax_x; i<= nmax_x; ++i)
            {
                for (int j = -nmax_y; j <= nmax_y; ++j)
                {
                    for (int k = -nmax_z; k <= nmax_z; ++k)
                    {
                        n = i * a1 + j * a2 + k * a3;
                        rijn = rij + n;
                        double r_norm = rijn.norm();
                        double r_sqr = r_norm * r_norm;
                        if (rijn.norm() <= rcut)
                        {
                            double intact = toeV*(0.5*elem1.q*elem2.q)*((-2.*kappa/sqrt(M_PI))*(exp(-r_sqr*kappa*kappa)/r_norm)-(erfc(r_norm*kappa)/r_sqr))/r_norm;
                            if (i == 0 && j == 0 && k ==0)
                            {
                                if (rijn.norm() > 0.1)
                                {
                                    deriv.col(0)(0) += rijn[0] * rijn[0] * intact ;
                                    deriv.col(1)(0) += rijn[0] * rijn[1] * intact ;
                                    deriv.col(1)(1) += rijn[1] * rijn[1] * intact ;
                                    deriv.col(2)(0) += rijn[0] * rijn[2] * intact ;
                                    deriv.col(2)(1) += rijn[1] * rijn[2] * intact ;
                                    deriv.col(2)(2) += rijn[2] * rijn[2] * intact ;
                                }
                            }
                            else
                            {
                                deriv.col(0)(0) += rijn[0] * rijn[0] * intact ;
                                deriv.col(1)(0) += rijn[0] * rijn[1] * intact ;
                                deriv.col(1)(1) += rijn[1] * rijn[1] * intact ;
                                deriv.col(2)(0) += rijn[0] * rijn[2] * intact ;
                                deriv.col(2)(1) += rijn[1] * rijn[2] * intact ;
                                deriv.col(2)(2) += rijn[2] * rijn[2] * intact ;
                            }
                        }
                    }
                }
            }
        }
    }
    // std::cout << deriv << std::endl;

    // std::cout <<  real_deriv[0]  <<" " <<  real_deriv[1] << " " <<  real_deriv[2] << std::endl;

    //Reciprocal Contribution - 3 terms; volume, r vector, k vector;
    rij.setZero();
    n.setZero();
    rijn.setZero();
    Eigen::Vector3d kvecs;
    kvecs.setZero();
    double intact[5];

    for ( auto& elem1 : atoms)
    {
        for ( auto& elem2 : atoms)
        {
            rij << elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z;
            for (int i = -kmax_x; i<= kmax_x; ++i)
            {
                for (int j = -kmax_y; j <= kmax_y; ++j)
                {
                    for (int k = -kmax_z; k <= kmax_z; ++k)
                    {
                        kvecs.setZero();
                        kvecs = (g1 * i) + (g2 * j) + (g3 * k);
                        double kk = kvecs.norm() * kvecs.norm();
                        intact[0] = toeV * (2. * M_PI / V) * elem1.q * elem2.q;
                        intact[1] = -2. * exp(-kk/4./kappa/kappa)/kk/kk * cos(kvecs.dot(rij));
                        intact[2] = -1./2. * exp(-kk/4./kappa/kappa)/(kk * kappa * kappa) * cos(kvecs.dot(rij));
                        intact[3] = -exp(-kk/4./kappa/kappa)/kk * sin(kvecs.dot(rij));
                        intact[4] = -exp(-kk/4./kappa/kappa)/kk * sin(kvecs.dot(rij));
                        if (kvecs.norm() <= kcut)
                        {

                            if (i == 0 && j == 0 && k ==0)
                            {
                                continue;
                            }
                            else
                            {
                                // wrt volume;
                                deriv.col(0)(0) += toeV*(-2. * M_PI * elem1.q * elem2.q / V / kk) * exp(-kk / (4. * kappa * kappa)) * cos(kvecs.dot(rij));
                                deriv.col(1)(1) += toeV*(-2. * M_PI * elem1.q * elem2.q / V / kk) * exp(-kk / (4. * kappa * kappa)) * cos(kvecs.dot(rij));
                                deriv.col(2)(2) += toeV*(-2. * M_PI * elem1.q * elem2.q / V / kk) * exp(-kk / (4. * kappa * kappa)) * cos(kvecs.dot(rij));

                                // wrt rvec;
                                deriv.col(0)(0) += intact[0] * intact[4] * kvecs[0] * rij[0];
                                deriv.col(1)(0) += intact[0] * intact[4] * kvecs[0] * rij[1];
                                deriv.col(1)(1) += intact[0] * intact[4] * kvecs[1] * rij[1];
                                deriv.col(2)(0) += intact[0] * intact[4] * kvecs[0] * rij[2];
                                deriv.col(2)(1) += intact[0] * intact[4] * kvecs[1] * rij[2];
                                deriv.col(2)(2) += intact[0] * intact[4] * kvecs[2] * rij[2];

                                //wrt kvec;
                                deriv.col(0)(0) += -kvecs[0] * intact[0] * (kvecs[0] * (intact[1] + intact[2]) + rij[0] * intact[3]);
                                deriv.col(1)(0) += -kvecs[0] * intact[0] * (kvecs[1] * (intact[1] + intact[2]) + rij[1] * intact[3]);
                                deriv.col(1)(1) += -kvecs[1] * intact[0] * (kvecs[1] * (intact[1] + intact[2]) + rij[1] * intact[3]);
                                deriv.col(2)(0) += -kvecs[0] * intact[0] * (kvecs[2] * (intact[1] + intact[2]) + rij[2] * intact[3]);
                                deriv.col(2)(1) += -kvecs[1] * intact[0] * (kvecs[2] * (intact[1] + intact[2]) + rij[2] * intact[3]);
                                deriv.col(2)(2) += -kvecs[2] * intact[0] * (kvecs[2] * (intact[1] + intact[2]) + rij[2] * intact[3]);

                            }
                        }
                    }
                }
            }
        }
    }
    deriv(1,0) = deriv(0,1);
    deriv(2,0) = deriv(0,2);
    deriv(2,1) = deriv(1,2);

    // std::cout << deriv << std::endl;
    unitcell_init.strain_deriv = deriv;
}

//Function that takes in energy and coordinates to compute the forces on each ions - rigid ion model only;

UnitCell calc_forces(UnitCell unitcell_init)
{
    UnitCell unitcell_wforces = unitcell_init;

    double total_energy = calc_short_range_buckingham_potential(unitcell_init) + calc_electrostatics_3D(unitcell_init);
    // std::cout << total_energy << std::endl;
    std::vector<Atom> atoms = unitcell_init.coordinates_cart;
    Eigen::Vector3d lattice_constants = unitcell_init.lattice_constants;
    Eigen::Matrix3d lattice_vectors = unitcell_init.lattice_vectors;
    Eigen::Vector3d v1, v2, v3;
    v1 = lattice_vectors.row(0);
    v2 = lattice_vectors.row(1);
    v3 = lattice_vectors.row(2);
    Eigen::Vector3d a1 = unitcell_init.lattice_vectors.row(0);
    Eigen::Vector3d a2 = unitcell_init.lattice_vectors.row(1);
    Eigen::Vector3d a3 = unitcell_init.lattice_vectors.row(2);

    double a = lattice_constants[0];
    double b = lattice_constants[1];
    double c = lattice_constants[2];

    //Find max force cut-off, should be same as max buck cut off
    double cutoff = 0;
    std::vector<Buckingham> buck = unitcell_init.buckingham_potentials;

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

    Eigen::Vector3d rij;
    Eigen::Vector3d n;
    Eigen::Vector3d rijn;

    Eigen::Vector3d pot_deriv;
    pot_deriv.setZero();

    for (auto& elem1 : atoms)
    {
        pot_deriv.setZero();

        for (auto& elem2 : atoms)
        {
            rij << elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z;
            // std::cout << rij.norm() << std::endl;
            for (int i = -maxx; i <= maxx; ++i)
            {
                for (int j = -maxy; j <= maxy; ++j)
                {
                    for (int k = -maxz; k <= maxz; ++k)
                    {
                        n = i * a1 + j * a2+ k * a3;
                        rijn = rij + n;

                        if (i == 0 && j ==0 && k == 0)
                        {
                            for (const auto& pot : unitcell_init.buckingham_potentials)
                            {
                                double expterm = (-pot.A / pot.rho) * exp(-rij.norm()/pot.rho) + 6. * pot.C / pow(rij.norm(), 7.);
                                if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                                        ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                                        rij.norm() <= pot.cut_off2 && rij.norm() > 0.1)
                                {
                                    pot_deriv += rij/rij.norm() * expterm;

                                }
                                else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                                         ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                                         rij.norm() <= pot.cut_off2 && rij.norm() > 0.1)
                                {
                                    pot_deriv += rij/rij.norm() * expterm;

                                }
                            }
                        }
                        else
                        {
                            for (const auto pot : unitcell_init.buckingham_potentials)
                            {
                                double expterm = (-pot.A / pot.rho) * exp(-rijn.norm()/pot.rho) + 6. * pot.C / pow(rijn.norm(), 7.);
                                if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                                        ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                                        rijn.norm() <= pot.cut_off2)
                                {
                                    pot_deriv += (rijn/rijn.norm()) * expterm;

                                }
                                else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                                         ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                                         rijn.norm() <= pot.cut_off2)
                                {
                                    pot_deriv += (rijn/rijn.norm()) * expterm;
                                }
                            }
                        }
                    }

                }
            }
        }
        // pot_deriv *= 0.5;
        pot_deriv = lattice_vectors * pot_deriv;

        elem1.fx = pot_deriv[0];
        elem1.fy = pot_deriv[1];
        elem1.fz = pot_deriv[2];
        // std::cout << pot_deriv[0]   <<" " << pot_deriv[1]  << " " << pot_deriv[2]  << std::endl;

    }

    //Now calculate electrostatic contributions to force
    const double toeV = 14.39964390675221758120;
    int natom = unitcell_init.coordinates_frac.size();

    Eigen::Vector3d g1 = unitcell_init.reciprocal_vectors.row(0);
    Eigen::Vector3d g2 = unitcell_init.reciprocal_vectors.row(1);
    Eigen::Vector3d g3 = unitcell_init.reciprocal_vectors.row(2);

    double V = unitcell_init.volume;

    long double kappa = pow(((natom * 1. * M_PI*M_PI*M_PI)/ V / V), 1./6.);
    double rcut = pow( -log(10E-17)/ (kappa * kappa),1./2.);
    double kcut = 2.*kappa*sqrt((-log(10E-17)));

    // double kappa = sqrt(0.159532);
    // double rcut = 13.1605;
    // double kcut = 4.199059;

    int nmax_x = ceil(rcut/a1.norm());
    int nmax_y = ceil(rcut/a2.norm());
    int nmax_z = ceil(rcut/a3.norm());
    int kmax_x = ceil(kcut/g1.norm());
    int kmax_y = ceil(kcut/g2.norm());
    int kmax_z = ceil(kcut/g3.norm());

    rij.setZero();
    n.setZero();
    rijn.setZero();

    //Real part contribution Note calculated derivatives.
    Eigen::Vector3d real_deriv;
    for ( auto& elem1 : atoms)
    {
        real_deriv.setZero();
        for ( auto& elem2 : atoms)
        {
            rij << elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z;
            for (int i = -nmax_x; i<= nmax_x; ++i)
            {
                for (int j = -nmax_y; j <= nmax_y; ++j)
                {
                    for (int k = -nmax_z; k <= nmax_z; ++k)
                    {
                        n = i * a1 + j * a2 + k * a3;
                        rijn = rij + n;
                        double r_norm = rijn.norm();
                        double r_sqr = r_norm * r_norm;
                        // Eigen::Vector3d prefactor = rijn/pow(rijn.norm(), 3.);
                        // double erfcterm = erfc(k * rijn.norm()) + (2. * kappa / sqrt(M_PI) ) * rijn.norm() * exp(-kappa*kappa*rijn.norm() * rijn.norm());
                        if (rijn.norm() <= rcut)
                        {
                            if (i == 0 && j == 0 && k ==0)
                            {
                                if (rijn.norm() > 0.1)
                                {
                                    real_deriv += rijn * toeV * ( elem1.q * elem2.q) *
                                                  ((-2. * kappa / sqrt(M_PI)) *
                                                   (exp(-r_sqr * kappa * kappa) / r_norm) -
                                                   (erfc(r_norm *kappa) / r_sqr)) / r_norm;
                                    // real_deriv += rijn* toeV*(0.5*elem1.q*elem2.q)*((-2. * kappa/sqrt(M_PI))*(exp(-r_sqr*kappa * kappa)/rnorm)-(erfc(rnorm * kappa )/r_sqr))/rnorm;
                                    // real_deriv +=  (-rijn / rijn.norm()) * (-0.5 * elem1.q * elem2.q / (rijn.norm() * rijn.norm()))
                                    // * ( (2. * kappa * rijn.norm() * exp(-kappa * kappa * rijn.norm() * rijn.norm())) / sqrt(M_PI) + erfc(kappa * rijn.norm()));
                                    // real_deriv += toeV * elem1.q * elem2.q * prefactor * erfcterm;
                                }
                            }
                            else
                            {
                                real_deriv += rijn * toeV * ( elem1.q * elem2.q) *
                                              ((-2. * kappa / sqrt(M_PI)) *
                                               (exp(-r_sqr * kappa * kappa) / r_norm) -
                                               (erfc(r_norm *kappa) / r_sqr)) / r_norm;
                                // real_deriv += rijn* toeV*(0.5*elem1.q*elem2.q)*((-2. * kappa/sqrt(M_PI))*(exp(-r_sqr*kappa * kappa)/rnorm)-(erfc(rnorm * kappa )/r_sqr))/rnorm;
                                // real_deriv += (-rijn / rijn.norm()) * (-0.5 * elem1.q * elem2.q / (rijn.norm() * rijn.norm()))
                                // * ( (2. * kappa * rijn.norm() * exp(-kappa * kappa * rijn.norm() * rijn.norm())) / sqrt(M_PI) + erfc(kappa * rijn.norm()));
                                // real_deriv += toeV * elem1.q * elem2.q * prefactor * erfcterm;

                            }
                        }
                    }
                }
            }
        }
        // std::cout <<  real_deriv[0]  <<" " <<  real_deriv[1] << " " <<  real_deriv[2] << std::endl;
        real_deriv = lattice_vectors * real_deriv;
        elem1.fx += real_deriv[0];
        elem1.fy += real_deriv[1];
        elem1.fz += real_deriv[2];
    }
    rijn.setZero();

    //Reciprocal contribution, again derivatives, not forces;
    Eigen::Vector3d reci_deriv;
    Eigen::Vector3d kvecs;
    kvecs.setZero();
    for ( auto& elem1 : atoms)
    {
        reci_deriv.setZero();
        for ( auto& elem2 : atoms)
        {
            rijn << (elem1.x - elem2.x), (elem1.y - elem2.y), (elem1.z - elem2.z);
            // std::cout << rijn[0] << " " << rijn[1] << " " << rijn[2] << std::endl;
            for (int i = -kmax_x; i <= kmax_x; ++i)
            {
                for (int j = -kmax_y; j <= kmax_y; ++j)
                {
                    for (int k = -kmax_z; k <= kmax_z; ++k)
                    {
                        kvecs = i * g1 + j * g2 + k * g3;
                        // std::cout << i << " " << j << " " << k << " " <<  kvecs[0] << " " << kvecs[1] << " " << kvecs[2] << std::endl;
                        double kk = kvecs.norm() * kvecs.norm();

                        if (kvecs.norm() <= kcut)
                        {
                            if (i == 0 && j == 0 && k == 0)
                            {
                                continue;
                            }
                            else
                            {
                                reci_deriv += kvecs*toeV * ((4. * M_PI)/ V) * elem1.q * elem2.q * exp(-kk/(4 * kappa * kappa))/kk * -sin(kvecs.dot(rijn));
                                // reci_deriv += toeV *  kvecs * toeV *((2.*M_PI)/V)*elem1.q*elem2.q*exp(-0.25 * kk/kappa/kappa)/kk * -sin(kvecs.dot(rijn));
                                // reci_deriv += -((2. * M_PI * elem1.q * elem2.q) / (V * kk)) * kvecs
                                // * exp(-kk/4. * kappa * kappa) * -sin(kvecs.dot(rijn));
                                // std::cout << elem1.label<< " "<< elem2.label << " " << rijn.norm() << std::endl;
                            }
                        }
                    }
                }
            }

        }
        // std::cout << reci_deriv[0] <<" " << reci_deriv[1]  << " " << reci_deriv[2]  << std::endl;
        reci_deriv = lattice_vectors * reci_deriv;
        elem1.fx += reci_deriv[0];
        elem1.fy += reci_deriv[1];
        elem1.fz += reci_deriv[2];
    }
    // for ( auto& elem : atoms)
    // {
    //     std::cout << elem.fx <<" " << elem.fy << " " << elem.fz << std::endl;
    // }

    return unitcell_wforces;
}
