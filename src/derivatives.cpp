#include <cmath> 
#include <Energy.h>
#include <iostream>
#include <Struct_Atom.h>








//Function that takes in energy and coordinates to compute the forces on each ions - rigid ion model only;

UnitCell calc_forces(UnitCell unitcell_init)
{
	UnitCell unitcell_wforces;
	double total_energy = calc_short_range_buckingham_potential(unitcell_init) + calc_electrostatics_3D(unitcell_init);
	// std::cout << total_energy << std::endl;
	std::vector<Atom> atoms = unitcell_init.coordinates_cart;

	for (auto& elem : atoms)
	{

	}
	return unitcell_wforces;
}