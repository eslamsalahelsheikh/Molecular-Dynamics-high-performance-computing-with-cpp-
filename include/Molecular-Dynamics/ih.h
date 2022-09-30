#ifndef MOLECULARDYNAMICS_IH_H
#define MOLECULARDYNAMICS_IH_H

# include "vector.h"
#include "../../../include/Molecular-Dynamics/atoms.h"
#include <fstream>
#include <iostream>
#include <cstring>
#include <memory>

/**
 * @brief  generate a cluster of atoms
 * 
 * @param n layer number
 * @param d atomic distance between atoms
 * @return Positions_t the positions of atoms
 */
Positions_t generate_cluster(int n, double d);

#endif //MOLECULARDYNAMICS_IH_H
