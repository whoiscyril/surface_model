#include <cmath>
#include <math.h>
#include <vector>
#include <Simple_Math.h>
#include <Struct_Atom.h>
#include <iostream>
#include <iomanip>

//Calculate buckingham potential contributions to energy

double calc_short_range_buckingham_potential (const UnitCell unitcell)
{
    //First find out maximum cut-off in the potentials;
    double max_pot_cutoff = 0.;
    for (const auto elem : unitcell.buckingham_potentials)
    {
        if (elem.cut_off2 > max_pot_cutoff)
        {
            max_pot_cutoff = elem.cut_off2;
        }
        else
        {
            continue;
        }
    }

    //setting global parameters
    double short_range_energy = 0.;
    Eigen::Vector3d a1, a2, a3;
    a1 = unitcell.lattice_vectors.row(0);
    a2 = unitcell.lattice_vectors.row(1);
    a3 = unitcell.lattice_vectors.row(2);

    int nmaxx = ceil(max_pot_cutoff/unitcell.lattice_constants[0]);
    int nmaxy = ceil(max_pot_cutoff/unitcell.lattice_constants[1]);
    int nmaxz = ceil(max_pot_cutoff/unitcell.lattice_constants[2]);
    // std::cout << nmaxx << " " << nmaxy << " " << nmaxz << std::endl;

    Eigen::Vector3d rij;
    Eigen::Vector3d n;
    Eigen::Vector3d rijn;
    for (const auto elem1 : unitcell.coordinates_cart)
    {
        for (const auto elem2 : unitcell.coordinates_cart)
        {
            rij << elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z;
            for (int i = -nmaxx; i <= nmaxx+1; ++i)
            {
                for (int j = -nmaxy; j <= nmaxy+1; ++j)
                {
                    for (int k = -nmaxz; k<= nmaxz+1; ++k)
                    {
                        if (i == 0 && j == 0 && k == 0)
                        {
                            // std::cout << i << " " << j << " " << k << std::endl;
                            //Calculate interactions within the unit cell;
                            for (const auto pot : unitcell.buckingham_potentials)
                            {
                                if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                                        ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                                        rij.norm() <= pot.cut_off2 && rij.norm() > 0.1)
                                {
                                    short_range_energy += pot.A * exp(-1. / pot.rho * rij.norm()) - pot.C / pow(rij.norm(), 6.);
                                    // std::cout << short_range_energy << std::endl;

                                }
                                else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                                         ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                                         rij.norm() <= pot.cut_off2 && rij.norm() > 0.1)
                                {
                                    short_range_energy += pot.A * exp(-1. / pot.rho * rij.norm()) - pot.C / pow(rij.norm(), 6.);
                                    // std::cout << short_range_energy << std::endl;

                                }
                            }
                        }
                        else
                        {
                            n = i * a1 + j * a2+ k * a3;
                            rijn = rij + n;
                            for (const auto pot : unitcell.buckingham_potentials)
                            {
                                if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                                        ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                                        rijn.norm() <= pot.cut_off2)
                                {
                                    short_range_energy += pot.A * exp(-1. / pot.rho * rijn.norm()) - pot.C / pow(rijn.norm(), 6.);
                                    // std::cout << short_range_energy << std::endl;

                                }
                                else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                                         ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                                         rijn.norm() <= pot.cut_off2)
                                {
                                    short_range_energy += pot.A * exp(-1. / pot.rho * rijn.norm()) - pot.C / pow(rijn.norm(), 6.);
                                }
                            }

                        }
                    }
                }
            }
        }

    }
    // std::cout << 0.5 * short_range_energy << std::endl;
    return 0.5 * short_range_energy;
}
// double calc_short_range_buckingham_potential (const std::vector<Atom>& cell,
//         const std::vector<Atom>& bulk_supercell,
//         const std::vector<Buckingham>& buckingham_potentials)
// {
//     double short_range_energy = 0.0;

//     std::vector<Atom> core_atoms;
//     std::vector<Atom> shell_atoms;
//     std::vector<Atom> supercell_core_atoms;
//     std::vector<Atom> supercell_shell_atoms;

//     //Split core and shell

//     for (const auto& elem : cell)
//     {
//         if (elem.type == "core")
//         {
//             core_atoms.push_back(elem);
//         }
//         else if (elem.type == "shel")
//         {
//             shell_atoms.push_back(elem);
//         }
//     }

//     for (const auto& elem : bulk_supercell)
//     {
//         if (elem.type == "core")
//         {
//             supercell_core_atoms.push_back(elem);
//         }
//         else if (elem.type == "shel")
//         {
//             supercell_shell_atoms.push_back(elem);
//         }
//     }

//     //Now calculate interactions between unit cell and images
//     for (const auto& elem1 : cell)
//     {
//         for (const auto& elem2 : bulk_supercell)
//         {
//             double distance = get_distance(elem1, elem2);
//             for (const auto& pot : buckingham_potentials)
//             {
//                 if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
//                         ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
//                         distance <= pot.cut_off2)
//                 {
//                     double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
//                     short_range_energy += energy;
//                 }
//                 else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
//                          ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
//                          distance <= pot.cut_off2)
//                 {
//                     double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
//                     short_range_energy += energy;
//                 }
//             }
//         }
//     }

// // Calculate interactions within the unit cell
//     for (const auto& elem1 : cell)
//     {
//         for (const auto& elem2 : cell)
//         {
//             if (&elem1 != &elem2)   // Ensure we don't calculate self-interaction
//             {
//                 double distance = get_distance(elem1, elem2);
//                 for (const auto& pot : buckingham_potentials)
//                 {
//                     if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
//                             ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
//                             distance <= pot.cut_off2 && distance > 0.1)
//                     {
//                         double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
//                         short_range_energy += energy;
//                     }
//                     else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
//                              ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
//                              distance <= pot.cut_off2 && distance > 0.1)
//                     {
//                         double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
//                         short_range_energy += energy;
//                     }
//                 }
//             }
//         }
//     }
//     short_range_energy *= 0.5;
//     return short_range_energy;
// }

double calc_electrostatics_3D(UnitCell unitcell)
{
    //Setting global parameters
    double total_electrostatic_energy_3D = 0.;
    const double toeV = 14.39964390675221758120;
    using std::endl;
    using std::cout;
    int natom = unitcell.coordinates_frac.size();

    // cout << unitcell.lattice_vectors << endl;
    // Setting real space parameters
    Eigen::Vector3d a1 = unitcell.lattice_vectors.row(0);
    Eigen::Vector3d a2 = unitcell.lattice_vectors.row(1);
    Eigen::Vector3d a3 = unitcell.lattice_vectors.row(2);
    Eigen::Vector3d g1 = unitcell.reciprocal_vectors.row(0);
    Eigen::Vector3d g2 = unitcell.reciprocal_vectors.row(1);
    Eigen::Vector3d g3 = unitcell.reciprocal_vectors.row(2);
    double V = unitcell.volume;
    double c1 = unitcell.lattice_constants[0];
    double c2 = unitcell.lattice_constants[1];
    double c3 = unitcell.lattice_constants[2];

    long double kappa = pow(((natom * 1. * M_PI*M_PI*M_PI)/ V / V), 1./6.);
    double rcut = pow( -log(10E-17)/ (kappa * kappa),1./2.);
    double kcut = 2.*kappa*sqrt((-log(10E-17)));
    // long double kappa = 1/2.24278;
    // double rcut = 14.032;
    // double kcut = 5.57926;

    int nmax_x = ceil(rcut/a1.norm());
    int nmax_y = ceil(rcut/a2.norm());
    int nmax_z = ceil(rcut/a3.norm());
    int kmax_x = ceil(kcut/g1.norm());
    int kmax_y = ceil(kcut/g2.norm());
    int kmax_z = ceil(kcut/g3.norm());

    // int nmax_x = 3;
    // int nmax_y = 3;
    // int nmax_z = 2;
    // int kmax_x = 2;
    // int kmax_y = 4;
    // int kmax_z = 4;
    //real space contribution

    double real_energy = 0.;

    Eigen::Vector3d rij;
    Eigen::Vector3d n;
    Eigen::Vector3d rijn;
    for (const auto& elem1 : unitcell.coordinates_cart)
    {
        for (const auto& elem2 : unitcell.coordinates_cart)
        {
            rij << elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z;
            for (int i = -nmax_x; i <= nmax_x; i++)
            {
                for (int j = -nmax_y; j <= nmax_y; j++)
                {
                    for (int k = -nmax_z; k <= nmax_z; k++)
                    {
                        n.setZero();
                        n = i * a1 + j * a2 + k * a3;
                        rijn = rij - n;
                        if (n.norm() < rcut)
                        {
                            // std::cout << i << " " << j << " " << k << " "<< "Trans "<< n[0] << " " << n[1] << " " << n[2];
                            if (i ==0 && j==0 && k ==0 )
                            {
                                if (rijn.norm() > 0.1)
                                {
                                    real_energy += 0.5 * elem1.q * elem2.q * erfc(kappa * rijn.norm()) / rijn.norm();
                                    // std::cout << rijn.norm() << " " << erfc(rijn.norm() * kappa) << std::endl;
                                    // std::cout << " Rij " << rijn[0] <<" " <<rijn[1] << " " << rijn[2] << std::endl;

                                }

                            }
                            else
                            {
                                real_energy += 0.5 * elem1.q * elem2.q * erfc(kappa * rijn.norm()) / rijn.norm();
                                // std::cout << rijn.norm() << " " << erfc(rijn.norm() * kappa) << std::endl;
                                // std::cout << " Rij "<<rijn[0] <<" " <<rijn[1] << " " << rijn[2] << std::endl;

                            }
                        }

                    }
                }
            }

        }
    }

    //Reciprocal contribution

    Eigen::Vector3d distvecs;
    Eigen::Vector3d kvecs;

    long double reciprocal_energy = 0.;
    for (const auto& elem1 : unitcell.coordinates_cart)
    {
        for (const auto& elem2 : unitcell.coordinates_cart)
        {
            rijn << (elem1.x - elem2.x), (elem1.y - elem2.y), (elem1.z - elem2.z);
            for (int i = -kmax_x; i <= kmax_x; i++)
            {
                for (int j =-kmax_y; j <= kmax_y; j++)
                {
                    for (int k=-kmax_z; k <= kmax_z; k++)
                    {

                        kvecs << 0, 0, 0;
                        kvecs = (g1 * i) + (g2 * j) + (g3 * k);

                        long double kk = kvecs.norm() * kvecs.norm();

                        if (kvecs.norm() <= kcut)
                        {

                            if (i == 0 && j == 0 && k == 0)
                            {
                                continue;
                            }
                            else
                            {

                                reciprocal_energy +=  (2. * M_PI / V) * elem1.q * elem2.q * std::exp(-kk/(4.*kappa * kappa)) * std::cos(kvecs.dot(rijn))/ kk;

                            }
                        }

                    }
                }
            }
        }
    }

    //Self energy contribution
    double self_energy = 0.;
    for (const auto& elem : unitcell.coordinates_cart)
    {
        self_energy += (-kappa / sqrt(M_PI)) * elem.q * elem.q;
    }
    // cout <<"self " <<self_energy * 14.39964390675221758120 << endl;
    // cout << "reci " << reciprocal_energy * 14.39964390675221758120  << endl;
//
    // cout <<"real " <<real_energy * 14.39964390675221758120 << endl;
    // cout << "reci + self " << (reciprocal_energy + self_energy) * 14.39964390675221758120 <<endl;

    // for (const auto& elem : unitcell.coordinates_cart)
    // {
    //     std::cout << elem.x << " " << elem.y << " " << elem.z << std::endl;
    // }

    // cout << "total "<<std::setprecision(10) <<(real_energy+reciprocal_energy + self_energy) * 14.39964390675221758120 <<endl;
    total_electrostatic_energy_3D = self_energy + reciprocal_energy + real_energy;
    // std::cout << total_electrostatic_energy_3D * toeV << std::endl;
    return total_electrostatic_energy_3D * toeV;

}