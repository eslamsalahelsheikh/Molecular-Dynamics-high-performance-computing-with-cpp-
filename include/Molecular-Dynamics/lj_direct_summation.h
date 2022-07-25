#ifndef MOLECULAR_DYNAMICS_LJ_DIRECT_SUMMATION_H
#define MOLECULAR_DYNAMICS_LJ_DIRECT_SUMMATION_H

# include "atoms.h"
# include "neighbors.h"

double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0) ;
double lj_direct_summation_neighbors(Atoms &atoms, NeighborList &neighbor_list, double epsilon = 1.0, double sigma = 1.0);

#endif //MOLECULAR_DYNAMICS_LJ_DIRECT_SUMMATION_H
