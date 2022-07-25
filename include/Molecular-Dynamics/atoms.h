#ifndef MOLECULAR_DYNAMICS_ATOMS_H
#define MOLECULAR_DYNAMICS_ATOMS_H
#include "types.h"

struct Atoms {
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Energies_t energies;
    Names_t names;
    Masses_t masses;
    // number of atoms as input
    Atoms(const int &nb_atoms) :
            positions{3,nb_atoms}, velocities{3,nb_atoms}, forces{3, nb_atoms}, masses{nb_atoms}, names(nb_atoms) {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
        for(int i=0; i<nb_atoms; i++){names[i] = "H";}
    }
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
    Atoms(const Names_t &names, const Positions_t &poses) :
            names{names} ,positions{poses}, velocities{3, poses.cols()}, forces{3, poses.cols()}, masses{poses.cols()}{
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
    }
    size_t nb_atoms() const {
        return positions.cols();
    }
};

inline Atoms lattice(int nx, int ny, int nz, double sigma){
    Atoms atoms(nx*ny*nz);
    atoms.positions.setZero();
    atoms.velocities.setZero();
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int index = i + j*nx + k*nx*ny;
                atoms.positions.col(index) << i*sigma, j*sigma, k*sigma;
            }
        }
    }
    return atoms;
}
#endif //MOLECULAR_DYNAMICS_ATOMS_H
