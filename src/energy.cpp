#include "../include/Molecular-Dynamics/energy.h"
# include <iostream>


double KineticEnergy(Atoms &atoms, double mass) {
    double kinetic_energy = 0.0;
    for (int i = 0; i < atoms.nb_atoms(); ++i) {
        kinetic_energy += 0.5 * mass * pow(atoms.velocities.col(i).matrix().norm(), 2);
    }
    std::cout << "Kinetic energy: " << kinetic_energy << std::endl;
    return kinetic_energy;
}


double TotalEnergy(Atoms &atoms, double mass) {
    double potential_energy = lj_direct_summation(atoms);
    double kinetic_energy = KineticEnergy(atoms, mass);
    return potential_energy + kinetic_energy;
}