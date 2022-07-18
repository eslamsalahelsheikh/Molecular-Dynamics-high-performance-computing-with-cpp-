#ifndef MOLECULAR_DYNAMICS_ENERGY_H
#define MOLECULAR_DYNAMICS_ENERGY_H

# include "atoms.h"
# include "lj_direct_summation.h"

double KineticEnergy(Atoms &atoms, double mass);
double TotalEnergy(Atoms &atoms, double mass);

#endif //MOLECULAR_DYNAMICS_ENERGY_H
