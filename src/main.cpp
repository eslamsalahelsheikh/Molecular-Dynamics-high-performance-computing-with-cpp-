# include "../include/Molecular-Dynamics/verlet.h"


int main() {
    // Initialize atoms parameters
    double mass = 1.0;
    double sigma = 1.0;
    double epsilon = 1.0;
    double time = sqrt(mass * pow(sigma,2) / epsilon);
    double total_time = 100* time;
    std::cout << "total_time: " << total_time << std::endl;
    double time_step = time/1000.0;
    std::cout << "time_step: " << time_step << std::endl;
    // Reading initial positions and velocities from xyz file
    auto [names, positions, velocities]{read_xyz_with_velocities("../data/lj54.xyz")};
    // Creating visualization file
    std::ofstream traj("../data/output/traj.xyz");

    // initializing atoms with poses and velocities
    Atoms atoms{positions,velocities};
    // Main simulation loop
    for (int i = 0; i < total_time; ++i) {
        if (i % 10 == 0 ){ write_xyz(traj, atoms);}
        verlet_step1(atoms, time_step, mass);
        double total_energy = TotalEnergy(atoms, mass);
        verlet_step2(atoms, time_step, mass);
    }
    traj.close();
    std::cout << "Finished: Hoooooooooooooooray!!!!" << std::endl;

    return 0;
}