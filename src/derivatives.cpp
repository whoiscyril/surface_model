#include <cmath>
#include <Energy.h>
#include <iostream>
#include <Struct_Atom.h>

//Function that takes in energy and coordinates to compute the forces on each ions - rigid ion model only;

UnitCell calc_forces(UnitCell unitcell_init)
{
    UnitCell unitcell_wforces = unitcell_init;


    double total_energy = calc_short_range_buckingham_potential(unitcell_init) + calc_electrostatics_3D(unitcell_init);
    // std::cout << total_energy << std::endl;
    std::vector<Atom> atoms = unitcell_init.coordinates_cart;
    Eigen::Vector3d lattice_constants = unitcell_init.lattice_constants;

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


    int maxx = ceil(cutoff/a);
    int maxy = ceil(cutoff/b);
    int maxz = ceil(cutoff/c);


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

            for (int i = -maxx; i <= maxx; ++i)
            {
                for (int j = -maxy; j <= maxy; ++j)
                {
                    for (int k = -maxz; k <= maxz; ++k)
                    {
                        rij << elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z;
                        n << i * a, j * b, k * c;
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
                                    pot_deriv += -rij/rij.norm() * expterm;


                                }
                                else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                                         ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                                         rij.norm() <= pot.cut_off2 && rij.norm() > 0.1)
                                {
                                    pot_deriv += -rij/rij.norm() * expterm;


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
                                    pot_deriv += (-rijn/rijn.norm()) * expterm;

                                }
                                else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                                         ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                                         rijn.norm() <= pot.cut_off2)
                                {
                                    pot_deriv += (-rijn/rijn.norm()) * expterm;
                                }
                            }
                        }
                    }    

                }
            }
        }
        // std::cout << pot_deriv[0] * a <<" " << pot_deriv[1] * b << " " << pot_deriv[2] * c << std::endl;

    }

    //both of these contributions below were missing a factor of 2;

    //Now calculate electrostatic contributions to force
    const double toeV = 14.39964390675221758120;
    int natom = unitcell_init.coordinates_frac.size();
    Eigen::Vector3d a1 = unitcell_init.lattice_vectors.row(0);
    Eigen::Vector3d a2 = unitcell_init.lattice_vectors.row(1);
    Eigen::Vector3d a3 = unitcell_init.lattice_vectors.row(2);

    Eigen::Vector3d g1 = unitcell_init.reciprocal_vectors.row(0);
    Eigen::Vector3d g2 = unitcell_init.reciprocal_vectors.row(1);
    Eigen::Vector3d g3 = unitcell_init.reciprocal_vectors.row(2);

    double V = unitcell_init.volume;

    // long double kappa = pow(((natom * 1. * M_PI*M_PI*M_PI)/ V / V), 1./6.);
    // double rcut = pow( -log(10E-17)/ (kappa * kappa),1./2.);
    // double kcut = 2.*kappa*sqrt((-log(10E-17)));

    double kappa = 1./2.04154;
    double rcut = 12.7729;
    double kcut = 6.12922;

    int nmax_x = ceil(rcut/a1.norm());
    int nmax_y = ceil(rcut/a2.norm());
    int nmax_z = ceil(rcut/a3.norm());
    int kmax_x = ceil(kcut/g1.norm());
    int kmax_y = ceil(kcut/g2.norm());
    int kmax_z = ceil(kcut/g3.norm());

    rij.setZero();
    n.setZero();
    rijn.setZero();

        //Real part contribution Note calculated forces, not derivatves, hence negative sign in front of rijn.
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
                        n << i * a, j * b, k * c;
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
                                    real_deriv += rijn * toeV * (0.5 * elem1.q * elem2.q) *
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
                                    real_deriv += rijn * toeV * (0.5 * elem1.q * elem2.q) *
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
        std::cout <<  real_deriv[0]  <<" " <<  real_deriv[1] << " " <<  real_deriv[2] << std::endl;

    }
    rijn.setZero();

    //Reciprocal contribution, again forces, not derivatives;
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
                                reci_deriv += kvecs*toeV * ((2. * M_PI)/ V) * elem1.q * elem2.q * exp(-kk/(4 * kappa * kappa))/kk * -sin(kvecs.dot(rijn));
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
        std::cout << reci_deriv[0] <<" " << reci_deriv[1]  << " " << reci_deriv[2]  << std::endl;

    }


    return unitcell_wforces;
}