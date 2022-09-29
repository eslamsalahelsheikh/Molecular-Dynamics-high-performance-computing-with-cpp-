# include "../../../include/Molecular-Dynamics/verlet.h"
# include "../../../include/Molecular-Dynamics/xyz.h"
# include "../../../include/Molecular-Dynamics/energy.h"

int main() {
    // Initialize atoms parameters
    double mass = 1.0;
    double sigma = 1.0;
    double epsilon = 1.0;
    double time = sqrt(mass * pow(sigma, 2) / epsilon);
    double total_time = 5000 * time;
    // Reading initial initial_positions and initial_velocities from xyz file
    auto [names, initial_positions, initial_velocities]{read_xyz_with_velocities("lj54.xyz")};
    // Creating visualization file milestones/04/output
    std::ofstream traj("/home/fr/fr_fr/fr_ee64/Molecular-Dynamics/output/milestone_04/traj.xyz");

    // Main simulation loop
    std::ofstream energy_file("/home/eslam/Desktop/Molecular-Dynamics/output/milestone_04/energies_vs_time_step.csv");
    energy_file << "time_step,total_energy,kinetic_energy,potential_energy" << std::endl;
    // Main simulation loop
    for (double time_step = 1e-4; time_step < 30e-3; time_step += 1e-3) {
        // initializing atoms with poses and initial_velocities
        Atoms atoms{initial_positions, initial_velocities};
        Energy energy(atoms, epsilon, sigma, mass); // initialize energy class
        std::cout << "time_step: " << time_step << std::endl;
        for (int i = 0; i < total_time; ++i) {
            if (i % 10 == 0) write_xyz(traj, atoms);    // writing to xyz file
            double old_total_energy{energy.get_total_energy()};    // storing the old total energy
            verlet_step1(atoms, time_step, mass);   // updating positions
            energy.energy_update(atoms, epsilon, sigma);    // updating energies
            verlet_step2(atoms, time_step, mass);   // updating velocities
            // check if the system reached equilibrium point to move on to the next time step
            if (abs(energy.get_total_energy() - old_total_energy) < 0.001 and i > 2) {
                std::cout << "quilibrium reached for this step , exiting now!!" << std::endl;
                energy_file << time_step << ", " << energy.get_total_energy() << ", " << energy.get_kinetic_energy()
                            << ", " << energy.get_potential_energy() << std::endl;
                break;
            }
        }
    }
    traj.close();
    std::cout << "Simulation finished successfully !!" << std::endl;

    return 0;
}