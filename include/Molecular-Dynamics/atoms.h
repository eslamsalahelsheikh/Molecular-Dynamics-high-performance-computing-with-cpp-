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

    /**
     * @brief Construct a new Atoms object
     * 
     * @param nb_atoms number of atoms
     */
    Atoms(const int &nb_atoms) :
            positions{3, nb_atoms}, velocities{3, nb_atoms}, forces{3, nb_atoms}, masses{nb_atoms}, names(nb_atoms) {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
        for (int i = 0; i < nb_atoms; i++) { names[i] = "H"; }
    }
    /**
     * @brief Construct a new Atoms object
     * @param p Positions of atoms
     */
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
    /**
     * @brief Construct a new Atoms object
     * @param p Positions of atoms
     * @param v Velocities of atoms
     */
    Atoms(const Positions_t &p, const Velocities_t &v) :
            positions{p}, velocities{v}, forces{3, p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
    }
    /**
     * @brief Construct a new Atoms object
     * @param names Names of atoms
     * @param poses Positions of atoms
     */
    Atoms(const Names_t &names, const Positions_t &poses) :
            names{names}, positions{poses}, velocities{3, poses.cols()}, forces{3, poses.cols()}, masses{poses.cols()} {
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
    }
    /**
     * @return The number of atoms
     */
    int nb_atoms() const {
        return static_cast<int>(positions.cols());
    }

    /**
     * @brief Resize the atoms object to the given number of atoms
     * @param local_length New number of atoms
     */
    void resize(const int local_length) {
        if (local_length <= 0) throw std::runtime_error("local_length must be positive");
        positions.conservativeResizeLike(Positions_t::Zero(3, local_length)); // resize and keep the old values
        velocities.conservativeResizeLike(Velocities_t::Zero(3, local_length)); // resize and keep the old values
        forces.conservativeResizeLike(Forces_t::Zero(3, local_length)); // resize and keep the old values
        energies.conservativeResizeLike(Energies_t::Zero(3, local_length)); // resize and keep the old values
        masses.conservativeResizeLike(Masses_t::Ones(local_length));    // resize and keep the old values
    }
};

/**
 * @brief Create a cubic lattice of atoms
 * @param nx Lattic size in x direction
 * @param ny Lattic size in y direction
 * @param nz Lattic size in z direction
 * @param sigma Lattice constant
 * @return Atoms object with the given lattice
 */
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
