# include "../include/Molecular-Dynamics/verlet.h"

void verlet_step1(Atoms &atoms, double timestep, double mass) {
    atoms.velocities += 0.5 * atoms.forces * timestep / mass;   // update velocities
    atoms.positions += atoms.velocities * timestep;            // update positions
}

void verlet_step2(Atoms &atoms, double timestep, double mass) {
    atoms.velocities += 0.5 * atoms.forces * timestep / mass;   // update velocities
}

