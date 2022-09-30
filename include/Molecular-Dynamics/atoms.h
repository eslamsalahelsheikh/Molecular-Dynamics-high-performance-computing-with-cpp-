#ifndef MOLECULAR_DYNAMICS_ATOMS_H
#define MOLECULAR_DYNAMICS_ATOMS_H

#include "types.h"
#include <iostream>

struct Atoms {
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Energies_t energies;
    Names_t names;
    Masses_t masses;

    // number of atoms as input
    Atoms() {}

    Atoms(const int &nb_atoms) :
            positions{3, nb_atoms}, velocities{3, nb_atoms}, forces{3, nb_atoms}, masses{nb_atoms}, names(nb_atoms) {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
        for (int i = 0; i < nb_atoms; i++) { names[i] = "H"; }
    }

    Atoms(Positions_t &p)
            : positions{p},
              velocities{3, p.cols()},
              forces{3, p.cols()},
              energies{3, p.cols()},
              masses{p.cols()} {
        velocities.setZero();
        forces.setZero();
        energies.setZero();
        masses.setOnes();
    }

    Atoms(const Positions_t &p, const Velocities_t &v) :
            positions{p}, velocities{v}, forces{3, p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
    }

    Atoms(const Names_t &names, const Positions_t &poses) :
            names{names}, positions{poses}, velocities{3, poses.cols()}, forces{3, poses.cols()}, masses{poses.cols()} {
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
    }

    int nb_atoms() const {
        return static_cast<int>(positions.cols());
    }

    void resize(const int local_length) {
        if (local_length <= 0) throw std::runtime_error("local_length must be positive");
        positions.conservativeResizeLike(Positions_t::Zero(3, local_length));
        velocities.conservativeResizeLike(Velocities_t::Zero(3, local_length));
        forces.conservativeResizeLike(Forces_t::Zero(3, local_length));
        energies.conservativeResizeLike(Energies_t::Zero(3, local_length));
        masses.conservativeResizeLike(Masses_t::Ones(local_length));
    }
};

inline Atoms lattice(int nx, int ny, int nz, double sigma) {
    Atoms atoms(nx * ny * nz);
    atoms.positions.setZero();
    atoms.velocities.setZero();
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int index = i + j * nx + k * nx * ny;
                atoms.positions.col(index) << i * sigma, j * sigma, k * sigma;
            }
        }
    }
    return atoms;
}

#endif //MOLECULAR_DYNAMICS_ATOMS_H
