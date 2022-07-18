# include "../include/Molecular-Dynamics/verlet.h"


int main() {
    // Initialize atoms parameters
    double mass = 1.0;
    double sigma = 1.0;
    double epsilon = 1.0;
    int total_time = 100* sqrt(mass * pow(sigma,2) / epsilon);
    double time_step = total_time/100000;
    // Reading initial positions and velocities from xyz file
    auto [names, positions, velocities]{read_xyz_with_velocities("../data/lj54.xyz")};

    Atoms atoms{positions,velocities};
    // Main simulation loop
    for (int i = 0; i < total_time; ++i) {
        std::cout << "Step: " << i << std::endl;

        verlet_step1(atoms, time_step, mass);
        double total_energy = TotalEnergy(atoms, mass);
        std::cout << "Total energy: " << total_energy << std::endl;
        verlet_step2(atoms, time_step, mass);

    }

    std::cout << "Finished: Hoooooooooooooooray!!!!" << std::endl;

    return 0;
}