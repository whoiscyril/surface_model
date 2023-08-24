#include <cmath>
#include <Energy.h>
#include <iostream>
#include <Struct_Atom.h>

//Function that takes in energy and coordinates to compute the forces on each ions - rigid ion model only;

UnitCell calc_forces(UnitCell unitcell_init)
{
    UnitCell unitcell_wforces = unitcell_init;
    //         for (auto& elem : unitcell_wforces.coordinates_cart)
    // {
    // 	std::cout << elem.fx << " " << elem.fy << " " << elem.fz << std::endl;
    // }

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

    double pot_derivx = 0.;
    double pot_derivy = 0.;
    double pot_derivz = 0.;
    // int maxx = ceil(cutoff/a);
    // int maxy = ceil(cutoff/b);
    // int maxz = ceil(cutoff/c);
    int maxx = 0;
    int maxy = 0;
    int maxz = 0;

    Eigen::Vector3d rij;
    Eigen::Vector3d n;
    Eigen::Vector3d rijn;

    for (auto& elem1 : atoms)
    {
        // std::cout << pot_derivx << " " << pot_derivy << " " << pot_derivz << std::endl;
        pot_derivx = 0;
        pot_derivy = 0;
        pot_derivz = 0;
        // std::cout << pot_derivx << " " << pot_derivy << " " << pot_derivz << std::endl;

        for (auto& elem2 : atoms)
        {
            rij << elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z;
            // std::cout << rij << std::endl;
            // std::cout << elem1.label << " " << elem2.label << " " << std::endl;
            for (int i = -maxx; i <= maxx; ++i)
            {
                for (int j = -maxy; j <= maxy; ++j)
                {
                    for (int k = -maxz; k <= maxz; ++k)
                    {

                        if (i == 0 && j ==0 & k == 0)
                        {
                            for (const auto& pot : unitcell_init.buckingham_potentials)
                            {
                                double expterm = (-pot.A / pot.rho) * exp(-rij.norm()/pot.rho) + 6. * pot.C / pow(rij.norm(), 7.);
                                if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                                        ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                                        rij.norm() <= pot.cut_off2 && rij.norm() > 0.1)
                                {
                                    pot_derivx += (- elem1.x / rij.norm()) * expterm;
                                    pot_derivy += (- elem1.y / rij.norm()) * expterm;
                                    pot_derivz += (- elem1.z / rij.norm()) * expterm;
                                    // std::cout << " Yes1" << " " << elem1.label << " " << elem2.label <<" " << pot_derivx << " " << pot_derivy << " " << pot_derivz << std::endl;
                                    // std::cout << pot_derivx << " " << pot_derivy << " " << pot_derivz << std::endl;
                                    // std::cout << pot.A <<" " << pot.rho << " " << pot.C << std::endl;
                                    std::cout << (- elem1.x / rij.norm()) << std::endl;

                                }
                                else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                                         ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                                         rij.norm() <= pot.cut_off2 && rij.norm() > 0.1)
                                {
                                    pot_derivx += (- elem1.x / rij.norm()) * expterm;
                                    pot_derivy += (- elem1.y / rij.norm()) * expterm;
                                    pot_derivz += (- elem1.z / rij.norm()) * expterm;
                                    // std::cout << " Yes2" << " " << elem1.label << " " << elem2.label <<" " << pot_derivx << " " << pot_derivy << " " << pot_derivz << std::endl;
                                    // std::cout << pot_derivx << " " << pot_derivy << " " << pot_derivz << std::endl;
                                                             std::cout << (- elem1.x / rij.norm()) << std::endl;

                                }
                                // std::cout <<elem1.label << " " << elem2.label <<" " << pot_derivx << " " << pot_derivy << " " << pot_derivz << std::endl;
                            }
                            // std::cout <<elem1.label << " " << elem2.label <<" " << pot_derivx << " " << pot_derivy << " " << pot_derivz << std::endl;
                        }
                        // else
                        // {
                        //     n << i * a, j * b, k * c;
                        //     rijn << rij + n;
                        //     for (const auto pot : unitcell_init.buckingham_potentials)
                        //     {
                        //         double expterm = (-pot.A / pot.rho) * exp(-rijn.norm()/pot.rho) + 6. * pot.C / pow(rijn.norm(), 7.);

                        //         if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                        //                 ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                        //                 rijn.norm() <= pot.cut_off2)
                        //         {
                        //             pot_derivx += (- elem1.x / rijn.norm()) * expterm;
                        //             pot_derivy += (- elem1.y / rijn.norm()) * expterm;
                        //             pot_derivz += (- elem1.z / rijn.norm()) * expterm;
                        //             // std::cout << pot_derivx << " " << pot_derivy << " " << pot_derivz << std::endl;
                        //             // std:: cout << rijn.norm() << std::endl;

                        //         }
                        //         else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                        //                  ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                        //                  rijn.norm() <= pot.cut_off2)
                        //         {
                        //             pot_derivx += (- elem1.x / rijn.norm()) * expterm;
                        //             pot_derivy += (- elem1.y / rijn.norm()) * expterm;
                        //             pot_derivz += (- elem1.z / rijn.norm()) * expterm;
                        //         }
                        //     }
                        // }
                    }

                }
            }
            // std::cout <<elem1.label << " " << elem2.label <<" " << pot_derivx << " " << pot_derivy << " " << pot_derivz << std::endl;
        }
        elem1.fx = pot_derivx;
        elem1.fy = pot_derivy;
        elem1.fz = pot_derivz;
   	    // std::cout <<elem1.label << " " << elem2.label <<" " << elem1.fx << " " << elem1.fy << " " << elem1.fz << std::endl;


    }

    // for (const auto& elem : unitcell_init.coordinates_cart)
    // {
    //     std::cout << elem.fx << " " << elem.fy << " " << elem.fz << std::endl;
    // }

    // std::cout << pot_derivx << " " << pot_derivy << " " << pot_derivz << std::endl;
    return unitcell_wforces;
}