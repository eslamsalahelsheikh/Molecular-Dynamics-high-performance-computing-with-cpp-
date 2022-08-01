#include "../../../include/Molecular-Dynamics/simulation.h"


int main() {
    // Reading initial positions and velocities from xyz file
    auto [names, positions]{read_xyz("cluster_923.xyz")};
    // initializing atoms with poses
    Atoms atoms{positions};
    // Initialize simulation for given atoms
    Simulation simulation(atoms);
    // Run simulation
    simulation.initial_loop();
    // Deposit heat
    simulation.deposit_heat(atoms, 0.001);
    // Run relaxation loop
    double average_temp = simulation.relaxation_loop(5000);
    std::cout << "average_temp: " << average_temp << std::endl;

    return 0;
}