#include "../../../include/Molecular-Dynamics/simulation.h"
#include "../../../include/Molecular-Dynamics/ih.h"

int main() {
    // Initialize atoms class and getting simulation data
    Atoms atoms;
    SimulationData data;
    if (data.continue_old_experiment) {
        // Reading initial positions and velocities from old xyz file
        auto [names, positions, Velocities_t]{read_xyz_with_velocities(data.old_experiment_file)};
        atoms = Atoms{positions, Velocities_t};
    } else {
        // generating a cluster of atoms with layer numbers and atomic distance
        Positions_t positions = generate_cluster(data.layer_numbers, data.atomic_distance);
        atoms = Atoms{positions};
    }
    // Initialize simulation for given atoms
    Simulation simulation(atoms);

    // Run initial simulation until equilibrium is reached or continue old simulation
    if (data.continue_old_experiment)
        simulation.energy_update(atoms);
    else
        simulation.initial_loop();
    // Run heating and relaxation experiments
    for (int i = 0; i < data.expermint_num; i++) {
        // Deposit heat
        simulation.add_heat();
        // Run relaxation loop
        simulation.relaxation_loop(i);
    }
    return 0;
}