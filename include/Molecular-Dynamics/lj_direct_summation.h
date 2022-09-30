#ifndef MOLECULAR_DYNAMICS_LJ_DIRECT_SUMMATION_H
#define MOLECULAR_DYNAMICS_LJ_DIRECT_SUMMATION_H

# include "atoms.h"
# include "neighbors.h"

/**
* @brief Compute the forces and energies of the Lennard-Jones potential
* @param atoms Atoms object
* @param epsilon Depth of the potential well
* @param sigma Distance at which the potential is a minimum
* @return The potential energy
*/
double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

/**
* @brief Compute the forces and energies of the Lennard-Jones potential using a neighbor list
* @param atoms Atoms object
* @param neighbor_list NeighborList object
* @param epsilon Depth of the potential well
* @param sigma Distance at which the potential is a minimum
* @return The potential energy
*/
double
lj_direct_summation_neighbors(Atoms &atoms, NeighborList &neighbor_list, double epsilon = 1.0, double sigma = 1.0);

#endif //MOLECULAR_DYNAMICS_LJ_DIRECT_SUMMATION_H
