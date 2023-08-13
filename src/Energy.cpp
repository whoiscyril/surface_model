#include <cmath>
#include <math.h>
#include <vector>
#include <Simple_Math.h>
#include <Struct_Atom.h>
#include <iostream>
#include <iomanip>



double calc_short_range_buckingham_potential (const std::vector<Atom>& cell,
        const std::vector<Atom>& bulk_supercell,
        const std::vector<Buckingham>& buckingham_potentials)
{
    double short_range_energy = 0.0;

    std::vector<Atom> core_atoms;
    std::vector<Atom> shell_atoms;
    std::vector<Atom> supercell_core_atoms;
    std::vector<Atom> supercell_shell_atoms;

    //Split core and shell

    for (const auto& elem : cell)
    {
        if (elem.type == "core")
        {
            core_atoms.push_back(elem);
        }
        else if (elem.type == "shel")
        {
            shell_atoms.push_back(elem);
        }
    }

    for (const auto& elem : bulk_supercell)
    {
        if (elem.type == "core")
        {
            supercell_core_atoms.push_back(elem);
        }
        else if (elem.type == "shel")
        {
            supercell_shell_atoms.push_back(elem);
        }
    }


    //Now calculate interactions between unit cell and images
    for (const auto& elem1 : cell)
    {
        for (const auto& elem2 : bulk_supercell)
        {
            double distance = get_distance(elem1, elem2);
            for (const auto& pot : buckingham_potentials)
            {
                if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                        ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                        distance <= pot.cut_off2)
                {
                    double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
                    short_range_energy += energy;
                }
                else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                         ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                         distance <= pot.cut_off2)
                {
                    double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
                    short_range_energy += energy;
                }
            }
        }
    }

// Calculate interactions within the unit cell
    for (const auto& elem1 : cell)
    {
        for (const auto& elem2 : cell)
        {
            if (&elem1 != &elem2)   // Ensure we don't calculate self-interaction
            {
                double distance = get_distance(elem1, elem2);
                for (const auto& pot : buckingham_potentials)
                {
                    if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                            ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                            distance <= pot.cut_off2 && distance > 0.1)
                    {
                        double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
                        short_range_energy += energy;
                    }
                    else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                             ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                             distance <= pot.cut_off2 && distance > 0.1)
                    {
                        double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
                        short_range_energy += energy;
                    }
                }
            }
        }
    }
    short_range_energy *= 0.5;
    return short_range_energy;
}

double calc_electrostatics_3D(const std::vector<Atom>& cell, std::vector<double> lattice_constants)
{

    //First convert lattice constants from vectors to Eigen matrix

    Eigen::Matrix3d lvecs;
    Eigen::Matrix3d trans;
    lvecs.setZero();
    lvecs = Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(lattice_constants.data());

    trans = frac2cart(cell, lvecs);

    for (const auto& elem : cell)
    {
        Eigen::Vector3d frac;
        Eigen::Vector3d cart;
        frac << elem.x, elem.y, elem.z;
        cart = trans * frac;
        std::cout << frac << std::endl;
        std::cout << cart << std::endl;

    }




    // std::cout << lvecs << std::endl;
    Eigen::Vector3d a1, a2, a3;
    Eigen::Vector3d g1, g2, g3;
    a1.setZero();
    a2.setZero();
    a3.setZero();

    a1 = lvecs.row(0);
    a2 = lvecs.row(1);
    a3 = lvecs.row(2);

    double V = a1.dot(a2.cross(a3));

    g1 = (2. * M_PI / V) * (a2.cross(a3));
    g2 = (2. * M_PI / V) * (a3.cross(a1));
    g3 = (2. * M_PI / V) * (a1.cross(a2));

    double total_energy = 0.0;
    double real_energy = 0.0;
    double reciprocal_energy =0.0;
    double self_energy = 0.0;
    double rmax = 0.;
    double gmax = 0.;
    double kappa = 0.;
    int real_max_image = 0;
    int reciprocal_max_image = 0;


    double accuracy = 10E-30;
    int N = cell.size();
    double w = 1.0;

    // std::cout << N << std::endl;



    kappa = pow((N * w * pow(M_PI, 3.0)) / pow(V, 2.0), -1.0/6.0);
    rmax = pow(-log(accuracy) / pow(kappa,2.0),1.0/2.0);
    gmax = 2.0 * kappa * pow((-log(accuracy)),1.0/2.0);

    int rmaxx = ceil(rmax/a1.norm());
    int rmaxy = ceil(rmax/a2.norm());
    int rmaxz = ceil(rmax/a3.norm());
    int recimaxx = ceil(gmax/g1.norm());
    int recimaxy = ceil(gmax/g2.norm());
    int recimaxz = ceil(gmax/g3.norm());

    real_max_image = ceil(rmax/4.);
    reciprocal_max_image = ceil(gmax/(2.0 * M_PI /4.));




    //Real space energy

    // std::cout << real_max_image << " " << gmax<< std::endl;

    for (const auto& elem1 : cell)
    {
        for (const auto& elem2 : cell)
        {
            for (int ri = -rmaxx; ri <= rmaxx ; ++ri)
            {
                for (int rj = -rmaxy; rj <= rmaxy; ++rj)
                {
                    for (int rk = -rmaxz; rk <= rmaxz; ++rk)
                    {
                        double distancex = elem1.x - elem2.x + ri * a1.norm();
                        double distancey = elem1.y - elem2.y + rj * a2.norm();
                        double distancez = elem1.z - elem2.z + rk * a3.norm();
                        if (len(distancex, distancey, distancez) <= rmax)
                        {
                            if(ri == 0 && rj ==0 && rk==0)
                            {
                                continue;
                            }
                            else
                            {
                                // std::cout << ri << " " << rj << " " << rk << std::endl;
                                // double distancex = elem1.x - elem2.x + ri * 5.0;
                                // double distancey = elem1.y - elem2.y + rj * 5.0;
                                // double distancez = elem1.z - elem2.z + rk * 5.0;
                                real_energy += 0.5*elem1.q * elem2.q * std::erfc(kappa * len(distancex, distancey, distancez))/len(distancex, distancey, distancez);
                            }
                        }

                    }
                }
            }
        }
    }

    // Reciprocal space energy
    double gx, gy, gz;
    for (const auto& elem1 : cell)
    {
        for(const auto& elem2 : cell)
        {
            for (int rrx = -recimaxx ; rrx <= recimaxx; ++rrx)
            {
                for (int rry = -recimaxy; rry <= recimaxy; ++rry)
                {
                    for (int rrz = -recimaxz; rrz <= recimaxz; ++rrz)
                    {
                        gx = 2. * M_PI * rrx/a1.norm();
                        gy = 2. * M_PI * rry/a2.norm();
                        gz = 2. * M_PI * rrz/a3.norm();
                        if(len(gx, gy, gz) <=gmax)
                        {
                            if(rrx == 0 && rry == 0 && rrz ==0)
                            {
                                continue;
                            }
                            else
                            {
                                // gx = 2. * M_PI * rrx/5.;
                                // gy = 2. * M_PI * rry/5.;
                                // gz = 2. * M_PI * rrz/5.;
                                double gg = dp(gx,gy,gz,gx,gy,gz);
                                reciprocal_energy += 0.5 / (M_PI * V) * elem1.q * elem2 .q * 4. * M_PI * M_PI / gg * exp(-gg/(4.*kappa*kappa)) * cos(dp(gx, gy, gz, elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z));
                            }
                        }
                    }
                }
            }
        }
    }

    for (const auto& elem : cell)
    {
        self_energy += -kappa / sqrt(M_PI) * elem.q * elem.q;
    }


    total_energy = (real_energy + self_energy + reciprocal_energy) * 14.39964390675221758120;
    // std::cout <<std::setprecision(10) << total_energy  <<std::endl;
    // std::cout << kappa <<std::endl;
    return total_energy;


}