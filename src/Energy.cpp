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

// double calc_electrostatics_3D(const std::vector<Atom>& cell, std::vector<double> lattice_constants)
// {
//     using std::endl;
//     using std::cout;
//     //First convert lattice constants from vectors to Eigen matrix

//     Eigen::Matrix3d lvecs;
//     Eigen::Matrix3d trans;
//     lvecs.setZero();
//     lvecs = Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(lattice_constants.data());
//     trans = frac2cart(cell, lvecs);

//     // for (const auto& elem : cell)
//     // {
//     //     Eigen::Vector3d frac;
//     //     Eigen::Vector3d cart;
//     //     frac << elem.x, elem.y, elem.z;
//     //     cart = trans * frac;
//     //     std::cout << frac << std::endl;
//     //     std::cout << cart << std::endl;

//     // }
//    // std::cout << lvecs << std::endl;
//     Eigen::Vector3d a1, a2, a3;
//     Eigen::Vector3d g1, g2, g3;
//     a1.setZero();
//     a2.setZero();
//     a3.setZero();

//     a1 = lvecs.row(0);
//     a2 = lvecs.row(1);
//     a3 = lvecs.row(2);
//     // cout << a1 << endl;


//     double V = a1.dot(a2.cross(a3));
//     // cout << V << endl;
//     g1 = (2. * M_PI / V) * (a2.cross(a3));
//     g2 = (2. * M_PI / V) * (a3.cross(a1));
//     g3 = (2. * M_PI / V) * (a1.cross(a2));


//     double total_energy = 0.0;
//     double real_energy = 0.0;
//     double reciprocal_energy =0.0;
//     double self_energy = 0.0;
//     double rmax = 0.;
//     double gmax = 0.;
//     double kappa = 0.;
//     double accuracy = 10E-30;
//     int N = cell.size();
//     double w = 1.0;

//     // std::cout << N << std::endl;
//     // kappa = 0.536656;
//     kappa = pow((N * w * pow(M_PI, 3.0)) / pow(V, 2.0), 1.0/6.0);
//     rmax = pow(-log(accuracy) / pow(kappa,2.0),1.0/2.0);
//     gmax = 2.0 * kappa * pow((-log(accuracy)),1.0/2.0);

//     // cout << kappa << " " << rmax << " " << gmax << endl;

//     int rmaxx = ceil(rmax/a1.norm());
//     int rmaxy = ceil(rmax/a2.norm());
//     int rmaxz = ceil(rmax/a3.norm());
//     int recimaxx = ceil(gmax/g1.norm());
//     int recimaxy = ceil(gmax/g2.norm());
//     int recimaxz = ceil(gmax/g3.norm());

//     // cout << rmaxx << endl;


//     // std::cout << rmax << std::endl;


//     //Real space energy

//     // std::cout << real_max_image << " " << gmax<< std::endl;

//     for (const auto& elem1 : cell)
//     {
//         for (const auto& elem2 : cell)
//         {
//             for (int ri = -rmaxx; ri <= rmaxx ; ri++)
//             {
//                 for (int rj = -rmaxy; rj <= rmaxy; rj++)
//                 {
//                     for (int rk = -rmaxz; rk <= rmaxz; rk++)
//                     {
//                         Eigen::Vector3d rij;
//                         Eigen::Vector3d n;
//                         Eigen::Vector3d rijn;
//                         rij << elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z;
//                         n << ri, rj, rk;
//                         rijn = rij + n;
//                         rijn = trans * rijn;
//                         if (rijn.norm() <= rmax)
//                         {   

//                             if(ri == 0 && rj ==0 && rk==0)
//                             {
//                                 if (rij.norm() != 0.)
//                                 {   
//                                     real_energy += 0.5*elem1.q * elem2.q * std::erfc(kappa * rijn.norm())/rijn.norm();

//                                 }
//                             }
//                             else
//                             {

//                                 real_energy += 0.5*elem1.q * elem2.q * std::erfc(kappa * rijn.norm())/rijn.norm();
//                             }
//                         }

//                     }
//                 }
//             }
//         }
//     }       

//     // Reciprocal space energy

//     for (const auto& elem1 : cell)
//     {
//         for(const auto& elem2 : cell)
//         {
//             for (int rrx = -recimaxx ; rrx <= recimaxx; rrx++)
//             {
//                 for (int rry = -recimaxy; rry <= recimaxy; rry++)
//                 {
//                     for (int rrz = -recimaxz; rrz <= recimaxz; rrz++)
//                     {
//                         // cout << rrx << " " << rry << " " << rrz << endl;

//                         Eigen::Vector3d kvec =  rrx * g1 +  rry * g2 + rrz * g3;
//                         // cout << kvec << endl;
//                         Eigen::Vector3d dist_vec;
//                         // kvec = trans * kvec;
//                         dist_vec << elem1.x - elem2.x, elem1.y - elem2.y, elem1.z - elem2.z;
//                         dist_vec = trans * dist_vec;
//                         // cout << kvec.norm() << endl;
//                         if(kvec.norm() <=gmax)
//                         {
//                             // cout << kvec.norm() << endl;
//                             if(rrx == 0 && rry == 0 && rrz ==0)
//                             {
//                                 if (dist_vec.norm() < 0.1)
//                                 {
//                                     continue;
//                                 }
//                                 continue;
//                             }
//                             else
//                             {

//                                 reciprocal_energy += 0.5 / (M_PI * V) * elem1.q * elem2 .q * 4. * M_PI * M_PI / (kvec.norm()*kvec.norm()) * (exp(-kvec.norm() * kvec.norm())/(4.*kappa*kappa)) * cos(kvec.dot(dist_vec));
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     for (const auto& elem : cell)
//     {
//         self_energy += -kappa / sqrt(M_PI) * elem.q * elem.q;
//     }


//     total_energy = (real_energy + self_energy + reciprocal_energy) * 14.39964390675221758120;
//     std::cout <<std::setprecision(10) << total_energy  <<std::endl;
//     cout << real_energy * 14.39964390675221758120 << " " << reciprocal_energy*14.39964390675221758120  << " " <<self_energy * 14.39964390675221758120 << endl;
//     cout << kappa << endl;
//     // std::cout << kappa <<std::endl;
//     return total_energy;


// }

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
    // cout << V << endl;


    // cout << V << endl;
    Eigen::Vector3d g1 = (2 * M_PI / V) * (a2.cross(a3));
    Eigen::Vector3d g2 = (2 * M_PI / V) * (a3.cross(a1));
    Eigen::Vector3d g3 = (2 * M_PI / V) * (a1.cross(a2));

    double kappa = 1./2.04154;
    // double kappa = sqrt(0.220395);
    // long double kappa = pow(((natom * 0.8 * M_PI*M_PI*M_PI)/ V / V), 1./6.);
    // double rcut = pow( -log(10E-17)/ (kappa * kappa),1./2.);
    // double kcut = 2.*kappa*sqrt((-log(10E-17)));
    double rcut = 50.;
    double kcut = 50.;


    // cout << "rmax" << rcut << endl;
    // cout << "gmax" << kcut << endl;

    int nmax_x = ceil(rcut/a1.norm());
    int nmax_y = ceil(rcut/a2.norm());
    int nmax_z = ceil(rcut/a3.norm());
    // int kmax_x = ceil(kcut/g1.norm());
    // int kmax_y = ceil(kcut/g2.norm());
    // int kmax_z = ceil(kcut/g3.norm());
    int kmax_x = 6;
    int kmax_y = 6;
    int kmax_z = 6;
    // cout << g1 << endl;
    // cout << g2 << endl;
    // cout << g3 << endl;


    // cout << rcut << endl;
    // cout << kcut << endl;

    // cout << natom << endl;

    // cout << nmax_x << " " << nmax_y << " " <<nmax_z<<endl;
    // cout << kmax_x << " " <<kmax_y << " " <<kmax_z<<endl;



    //real space contribution

    double real_energy = 0.;
    
    Eigen::Vector3d rij;
    Eigen::Vector3d n;
    Eigen::Vector3d rijn;
    for (const auto& elem1 : cell)
    {
        for (const auto& elem2 : cell)
        {
            rij << (elem1.x - elem2.x), (elem1.y - elem2.y), (elem1.z - elem2.z);
            // cout << rij << endl;
            // cout << trans * rij << endl;
            for (int i = -nmax_x; i <= nmax_x; i++)
            {
                for (int j = -nmax_y; j <= nmax_y; j++)
                {
                    for (int k = -nmax_z; k <= nmax_z; k++)
                    {
                        // double distx = elem1.x - elem2.x + i * 5.;
                        // double disty = elem1.y - elem2.y + j * 5.;
                        // double distz = elem1.z - elem2.z + k * 5.;
                        n << i, j, k;
                        rijn = rij + n;
                        rijn = trans * rijn;
                        // rijn = rij + n;
                        // double rijn = len(distx, disty, distz);
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
    cout << kcut << endl;
    //Reciprocal contribution
    long double reciprocal_energy = 0.;
    for (const auto& elem1 : cell)
    {
        for (const auto& elem2 : cell)
        {
            rijn << (elem1.x - elem2.x) * a1 + (elem1.y - elem2.y) * a2 + (elem1.z - elem2.z) * a3;

            // cout << " --------after-------" << endl;
            // cout <<trans*rij << endl;
            for (int i = -kmax_x; i <= kmax_x; i++)
            {
                for (int j =-kmax_y; j <= kmax_y; j++)
                {
                    for (int k=-kmax_z; k <= kmax_z; k++)
                    {
                        // double kx = (2. * M_PI / 5.) * i;
                        // double ky = (2. * M_PI / 5.) * j;
                        // double kz = (2. * M_PI / 5.) * k;
                        // cout << distvecs << endl;
                        // kvecs << kx, ky, kz;
                        // cout << distvecs << endl;
                        kvecs << 0, 0, 0;
                        kvecs = (g1 * i) + (g2 * j) + (g3 * k);
                        // kvecs << g1[0] * i + g2[0] * i + g3[0]*i, g1[0] * j + g2[0] * j + g3[0]*j, g1[0] * k + g2[0] * k + g3[0]*k;
                        // kvecs[0] = i * g1[0] + j * g2[0] + k * g3[0];
                        // kvecs[1] = i * g1[1] + j * g2[1] + k * g3[1];
                        // kvecs[2] = i * g1[2] + j * g2[2] + k * g3[2];

                        // cout << "=====begin=======" << endl;
                        // cout << g1 << endl;
                        // cout << g1 * i << endl;
                        // cout << "=====finish=======" << endl;
                        // cout << kvecs << endl;
                        long double kk = kvecs[0] * kvecs[0] + kvecs[1] * kvecs[1] + kvecs[2] * kvecs[2];
                        // double kk = kvecs.norm() * kvecs.norm();
                        // cout << kvecs.norm() << " " <<kk << endl;
                        // double kijn = len(kx, ky, kz);
                        // cout << "=================" << endl;
                        // cout << kvecs << endl;
                        if (kvecs.norm() <= kcut)
                        {
                            
                            if (i == 0 && j == 0 && k == 0)
                            {
                                continue;
                            }
                            else
                            {
                                // cout << kvecs[0] << " " << kvecs[1] << " " << kvecs[2] << endl;
                                // double kk = dp(kx,ky,kz,kx,ky,kz);
                                // reciprocal_energy += 0.5 / (M_PI * V) * (elem1.q * elem2.q * 4. * M_PI * M_PI / kijn*kijn) * exp(-kijn*kijn/(4.*kappa*kappa)) * cos(kx*(elem1.x - elem2.x) + ky * (elem1.y - elem2.y) + kz * (elem1.z - elem2.z));
                                // cout << reciprocal_energy << endl;
                                // cout << rijn[0] << " " << rijn[1] << " " << rijn[2] << endl;
                                // cout << kvecs[0] << " " << kvecs[1] << " " << kvecs[2] << " " << rijn[0] << " " << rijn[1] << " " << rijn[2] << " " << kvecs.dot(rijn) << endl;
                                reciprocal_energy +=  (2. * M_PI / V) * elem1.q * elem2.q * std::exp(-kk/(4.*kappa * kappa)) * std::cos(kvecs.dot(rijn))/ kk;
                                // cout << -kk/(4.*kappa * kappa) << endl;
                                // cout << reciprocal_energy * toeV << " " <<  i <<" " << j << " " << k << endl;
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
    cout <<"self " <<self_energy * 14.39964390675221758120 << endl;
    cout << "reci " << reciprocal_energy * 14.39964390675221758120  << endl;

    cout <<"real " <<real_energy * 14.39964390675221758120 << endl;
    cout << "reci + self " << (reciprocal_energy + self_energy) * 14.39964390675221758120 <<endl;
    cout << "total "<<std::setprecision(10) <<(real_energy+reciprocal_energy + self_energy) * 14.39964390675221758120 <<endl;
    // cout << kappa << "  " << 1/kappa << endl;
    // total_electrostatic_energy_3D = (real_energy + reciprocal_energy + self_energy) * toeV;
    // cout<<total_electrostatic_energy_3D << endl;
    return total_electrostatic_energy_3D;






}