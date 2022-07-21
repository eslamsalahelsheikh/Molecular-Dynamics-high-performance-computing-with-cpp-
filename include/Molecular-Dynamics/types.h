#ifndef MYPROJECT_TYPES_H
#define MYPROJECT_TYPES_H

#include "Eigen/Dense"

// Type aliases
using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
using Energies_t = Eigen::Array3Xd;
using Names_t = std::vector<std::string>;
using Masses_t = Eigen::ArrayXd;
#endif //MYPROJECT_TYPES_H
