#ifndef MOLECULAR_DYNAMICS_ATOMS_H
#define MOLECULAR_DYNAMICS_ATOMS_H
#include "types.h"

struct Atoms {
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Energies_t energies;
    Atoms(Positions_t &p)
            : positions{p},
              velocities{3, p.cols()},
              forces{3, p.cols()},
              energies{3, p.cols()} {
        velocities.setZero();
        forces.setZero();
        energies.setZero();
    }
    Atoms(const Positions_t &p, const Velocities_t &v) :
            positions{p}, velocities{v}, forces{3, p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
    }

    size_t nb_atoms() const {
        return positions.cols();
    }
};

#endif //MOLECULAR_DYNAMICS_ATOMS_H
