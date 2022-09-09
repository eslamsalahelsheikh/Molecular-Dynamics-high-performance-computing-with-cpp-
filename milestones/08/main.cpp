#include "../../../include/Molecular-Dynamics/simulation.h"
#include "../../../include/Molecular-Dynamics/ih.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    // Reading initial positions and velocities from xyz file

    Atoms atoms;
    SimulationData data;
    if (data.continue_old_experiment) {
        auto [names, positions, Velocities_t]{read_xyz_with_velocities(data.old_experiment_file)};
        atoms = Atoms{positions, Velocities_t};
    } else {
        Positions_t positions = generate_cluster(data.layer_numbers, data.atomic_distance);
        atoms = Atoms{positions};
    }
    // Initialize simulation for given atoms
    Simulation simulation(atoms);
    Domain domain(MPI_COMM_WORLD, {3.5, 4.5, 5.5}, {1, 1, 4}, {0, 0, 1});

    // Run initial simulation
    if (!data.continue_old_experiment) simulation.initial_loop(domain);
    else simulation.energy_update(atoms);

    // Run heating and relaxation experiments
    for (int i = 0; i < data.expermint_num; i++) {
        // Deposit heat
        simulation.add_heat();
        // Run relaxation loop
        simulation.relaxation_loop(i);
    }
    MPI_Finalize();
    return 0;
}