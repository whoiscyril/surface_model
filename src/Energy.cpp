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


double calc_electrostatics_3D( std::vector<Atom>& cell, std::vector<double> lattice_constants)
{
    //Setting global parameters
    double total_electrostatic_energy_3D = 0.;
    const double toeV = 14.39964390675221758120;
    using std::endl;
    using std::cout;
    int natom = cell.size();
    //First convert lattice constants from vectors to Eigen matrix

    Eigen::Matrix3d lvecs;
    Eigen::Matrix3d trans;
    lvecs.setZero();
    lvecs = Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(lattice_constants.data());
    trans = frac2cart(cell, lvecs);


    //Setting real space parameters
    Eigen::Vector3d a1 = lvecs.row(0);
    Eigen::Vector3d a2 = lvecs.row(1);
    Eigen::Vector3d a3 = lvecs.row(2);


    //Setting reciprocal space parameters
    long double V = a1.dot(a2.cross(a3));



    Eigen::Vector3d g1 = (2 * M_PI / V) * (a2.cross(a3));
    Eigen::Vector3d g2 = (2 * M_PI / V) * (a3.cross(a1));
    Eigen::Vector3d g3 = (2 * M_PI / V) * (a1.cross(a2));

    long double kappa = pow(((natom * 1. * M_PI*M_PI*M_PI)/ V / V), 1./6.);
    double rcut = pow( -log(10E-17)/ (kappa * kappa),1./2.);
    double kcut = 2.*kappa*sqrt((-log(10E-17)));
 


    int nmax_x = ceil(rcut/a1.norm());
    int nmax_y = ceil(rcut/a2.norm());
    int nmax_z = ceil(rcut/a3.norm());
    int kmax_x = ceil(kcut/g1.norm());
    int kmax_y = ceil(kcut/g2.norm());
    int kmax_z = ceil(kcut/g3.norm());



    //real space contribution

    double real_energy = 0.;
    
    Eigen::Vector3d rij;
    Eigen::Vector3d n;
    Eigen::Vector3d rijn;
    for (const auto& elem1 : cell)
    {
        for (const auto& elem2 : cell)
        {

            for (int i = -nmax_x; i <= nmax_x; i++)
            {
                for (int j = -nmax_y; j <= nmax_y; j++)
                {
                    for (int k = -nmax_z; k <= nmax_z; k++)
                    {

                        rijn << (elem1.x - elem2.x + i) * a1 + (elem1.y - elem2.y + j) * a2 + (elem1.z - elem2.z + k) * a3;

                        if (rijn.norm() < rcut)
                        {
                            if (i ==0 && j==0 && k ==0 )
                            {
                                if (rijn.norm() > 0.1)
                                {
                                    real_energy += 0.5 * elem1.q * elem2.q * erfc(kappa * rijn.norm()) / rijn.norm();
   
                                }

                            }
                            else 
                            {
                                real_energy += 0.5 * elem1.q * elem2.q * erfc(kappa * rijn.norm()) / rijn.norm();

                            }
                        }

                    }
                }
            }

        }
    }
    Eigen::Vector3d distvecs;
    Eigen::Vector3d kvecs;
    // cout << kcut << endl;
    //Reciprocal contribution
    long double reciprocal_energy = 0.;
    for (const auto& elem1 : cell)
    {
        for (const auto& elem2 : cell)
        {
            rijn << (elem1.x - elem2.x) * a1 + (elem1.y - elem2.y) * a2 + (elem1.z - elem2.z) * a3;
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
    for (const auto& elem : cell)
    {
        self_energy += (-kappa / sqrt(M_PI)) * elem.q * elem.q;
    }
    // cout <<"self " <<self_energy * 14.39964390675221758120 << endl;
    // cout << "reci " << reciprocal_energy * 14.39964390675221758120  << endl;

    // cout <<"real " <<real_energy * 14.39964390675221758120 << endl;
    // cout << "reci + self " << (reciprocal_energy + self_energy) * 14.39964390675221758120 <<endl;
    // cout << "total "<<std::setprecision(10) <<(real_energy+reciprocal_energy + self_energy) * 14.39964390675221758120 <<endl;

    return total_electrostatic_energy_3D * toeV;






}