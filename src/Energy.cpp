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
for (const auto& elem1 : cell) {
    for (const auto& elem2 : bulk_supercell) {
        double distance = get_distance(elem1, elem2);
        for (const auto& pot : buckingham_potentials) {
            if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                distance <= pot.cut_off2) {
                double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
                short_range_energy += energy;
            } else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                       ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                       distance <= pot.cut_off2) {
                double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
                short_range_energy += energy;
            }
        }
    }
}

// Calculate interactions within the unit cell
for (const auto& elem1 : cell) {
    for (const auto& elem2 : cell) {
        if (&elem1 != &elem2) { // Ensure we don't calculate self-interaction
            double distance = get_distance(elem1, elem2);
            for (const auto& pot : buckingham_potentials) {
                if (((pot.atom1_type == elem1.type) && (pot.atom1_label == elem1.label)) &&
                    ((pot.atom2_type == elem2.type) && (pot.atom2_label == elem2.label)) &&
                    distance <= pot.cut_off2 && distance > 0.1) {
                    double energy = pot.A * exp(-1 / pot.rho * distance) - pot.C / pow(distance, 6);
                    short_range_energy += energy;
                } else if (((pot.atom2_type == elem1.type) && (pot.atom2_label == elem1.label)) &&
                           ((pot.atom1_type == elem2.type) && (pot.atom1_label == elem2.label)) &&
                           distance <= pot.cut_off2 && distance > 0.1) {
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

double calc_electrostatics_3D(const std::vector<Atom>& cell)
{
	double total_energy = 0.0;
	double real_energy = 0.0;
	double reciprocal_energy =0.0;
	double self_energy;

	double rmax;
	double gmax;
	double kappa;

	int real_max_image;
	int reciprocal_max_image;


	double accuracy = 10E-30;
	int N = 5;
	double w = 1.0;
	double V = 125.0;

	kappa = 1.4;
	// kappa = pow((N * w * pow(M_PI, 3.0)) / pow(V, 2.0), 1.0/6.0);
	// rmax = pow(-log(accuracy) / pow(kappa,2.0) ,1.0/2.0);
	// gmax = 2.0 * kappa * pow((-log(accuracy)),1.0/2.0);
	rmax = 100;
	gmax = 100;

	// real_max_image = ceil(rmax/5.0);
	// reciprocal_max_image = ceil(gmax/(2.0 * M_PI / 5.0));

	real_max_image = 100;
	reciprocal_max_image = 100;
	int rx, ry, rz;
	int rrx, rry, rrz;


	//Real space energy

	// std::cout << real_max_image << " " << gmax<< std::endl;

	for (const auto& elem1 : cell)
	{
		for (const auto& elem2 : cell)
		{
			for (int ri = -real_max_image; ri <= real_max_image ; ++ri)
			{
				for (int rj = -real_max_image; rj <= real_max_image; ++rj)
				{
					for (int rk = -real_max_image; rk <= real_max_image; ++rk)
					{							
							double distancex = elem1.x - elem2.x + ri * 5.0;
							double distancey = elem1.y - elem2.y + rj * 5.0;
							double distancez = elem1.z - elem2.z + rk * 5.0;
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
			for (int rrx = -reciprocal_max_image ; rrx <= reciprocal_max_image; ++rrx)
			{
				for (int rry = -reciprocal_max_image; rry <= reciprocal_max_image; ++rry)
				{
					for (int rrz = -reciprocal_max_image; rrz <= reciprocal_max_image; ++rrz)
					{
							gx = 2. * M_PI * rrx/5.;
							gy = 2. * M_PI * rry/5.;
							gz = 2. * M_PI * rrz/5.;
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
	std::cout <<std::setprecision(10) << total_energy  <<std::endl;
	// std::cout << kappa <<std::endl;
	return total_energy;


}