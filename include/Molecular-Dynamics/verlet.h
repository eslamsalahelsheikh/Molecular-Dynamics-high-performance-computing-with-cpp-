#ifndef __VERLET_H
#define __VERLET_H

# include <iostream>
# include "atoms.h"

/**
 * @brief Step1 to update positions and velocities using Verlet algorithm
 * @param atoms Atoms class object
 * @param timestep Timestep
 * @param mass Mass of the atoms
 */
void verlet_step1(Atoms &atoms, double timestep, double mass);

/**
 * @brief Step2 to velocities using Verlet algorithm
 * @param atoms Atoms class object
 * @param timestep Timestep
 * @param mass Mass of the atoms
 */
void verlet_step2(Atoms &atoms, double timestep, double mass);

#endif // VERLET_H