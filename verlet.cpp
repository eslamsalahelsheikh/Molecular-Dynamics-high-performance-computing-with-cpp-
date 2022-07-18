# include "verlet.h"
# include "lj_direct_summation.h"

void verlet_step1(Atoms &atoms, double timestep, double mass) {
    atoms.velocities += 0.5 * atoms.forces * timestep / mass;

    atoms.positions += atoms.velocities * timestep;
    std::cout << "new positions" << atoms.positions << std::endl;
}
void verlet_step2(Atoms &atoms, double timestep, double mass) {
    atoms.velocities += 0.5 * atoms.forces * timestep / mass;
}

void compute_force(Atoms &atoms, float sigma, float epsilon) {
    Eigen::Array3d r{atoms.positions.col(1) - atoms.positions.col(0)};
    Eigen::Array3Xd sig_over_r = sigma / r;
    lj_direct_summation(atoms, epsilon, sigma);

    Energies_t V_new = 4 * epsilon * (pow(sig_over_r, 12) - pow(sig_over_r, 6));  // potential energy between two atoms with distance r
    Energies_t F_new = 48 * epsilon * pow(sigma, 12) / pow(r, 13)
                     - 24 * epsilon * pow(sigma, 6) / pow(r, 7);  // force between two atoms with distance r

    std::cout << "new energy " << V_new << std::endl;
    atoms.forces.col(0) = F_new;
    atoms.forces.col(1) = -1 * F_new;
    std::cout << "forces " << atoms.forces << std::endl;
}

int main() {
    int nb_steps = 5;
//     double timestep = 1.0e-15; // 1 fs
    double mass = 1.0;
    float epsilon = 3.0;   // epsilon for Argon
    float sigma = 3.0; // sigma for Argon
    double timestep = 0.001 * sqrt(mass * pow(sigma, 2) / epsilon);
    int nb_atoms=2;
    Positions_t positions(3, nb_atoms);
    positions.col(0) << 0, 0, 0;
    positions.col(1) << 1, 1, 1;
    Velocities_t velocities(3, nb_atoms);
    velocities << 0.01, 0.01, 0.01, 0.01, 0.01, 0.01;
    Atoms atoms{positions,velocities};
    std::cout << "atom1 pos: " << atoms.positions.col(0) << std::endl;
    std::cout << "atom2 pos: " << atoms.positions.col(1) << std::endl;


//    atoms.forces << 0.01, 0.01, 0.01,
//                   0.01, 0.01, 0.01;

    for (int i = 0; i < nb_steps; ++i) {
        std::cout << "Step: " << i << std::endl;
        Forces_t forces_old = atoms.forces;

        verlet_step1(atoms, timestep, mass);
        compute_force(atoms,sigma,epsilon);
        verlet_step2(atoms, timestep, mass);
//        std::cout << atoms.positions << std::endl;
//        std::cout << atoms.forces - forces_old << std::endl;

    }

    std::cout << "Finished: Hoooooooooooooooray!!!!" << std::endl;

    return 0;
}