#ifndef DERIVATIVES_H
#define DERIVATIVES_H
#include <cmath> 
#include <Energy.h>
#include <iostream>
#include <Struct_Atom.h>

UnitCell calc_forces(UnitCell unitcell_init);
void calc_strain_deriv(UnitCell& unitcell_init);
void calc_lattice_deriv(UnitCell& unitcell_init);
void internal_derv2(UnitCell& unitcell);

#endif